#include "MathUtils/Cartesian3D.h"
#include <cmath>
#include <utility>
#include <array>
#include <cassert>
#include "DetectorsBase/Constants.h"
#include "DetectorsBase/Utils.h"


/**********************************************************************
 *                                                                    *
 * Material data on the cell on cylindrical layer                     *
 * Cell is limited by 2 radial planes at phiMin,phiMax,               *
 * radii radMin,radMax and XY planes at zMin, zMax.                   *
 * This limits are defined in the cell container class                *
 *                                                                    *
 **********************************************************************/

struct MatCell {
  static constexpr int NParams = 2; // number of material parameters stored
  float mRhoAv = 0.f;    ///< 
  float mX2X0  = 0.f;    ///<
};


/**********************************************************************
 *                                                                    *
 * Cylindrical material layer                                         *
 *                                                                    *
 **********************************************************************/
class MatLayerCyl {

 public:
  
  MatLayerCyl(float rMin,float rMax,float zMin,float zMax,int nz,int nphi);
  ~MatLayerCyl() = default;
  
  // ----------------------- segmentation 
  void  initSegmentation(float rMin,float rMax,float zMin,float zMax,int nz,int nphi);
  float getRMin()        const {return mRMin;}
  float getRMax()        const {return mRMax;}
  float getZMin()        const {return mZMin;}
  float getZMax()        const {return mZMax;}
  int   getNZSlices()    const {return mNZSlices;}
  int   getNPhiSlices()  const {return mNPhiSlices;}

  float getRMin2()       const {return mRMin2;}
  float getRMax2()       const {return mRMax2;}

  
  // linearized cell ID
  int getCellID(int iphi,int iz)       const {return iphi*mNZSlices + iz;}
  
  // obtain material cell, cell ID must be valid
  const MatCell& getCell(int iphi, int iz)   const {return mCells[getCellID(iphi,iz)];}
  MatCell& getCell(int iphi, int iz)               {return mCells[getCellID(iphi,iz)];}
  
  // ---------------------- Z slice manipulation
  // convert Z to Zslice
  int getZSliceID(float z)             const {return z>mZMin && z<mZMax ? int((z-mZMin)*mDZInv) : -1;}

  // lower boundary of Z slice
  float getZSliceMin(int id)           const {return mZMin + id*mDZ;}

  // upper boundary of Z slice
  float getZSliceMax(int id)           const {return mZMin + (id+1)*mDZ;}

  // ---------------------- Phi slice manipulation (0:2pi convention, no check is done)
  // convert Phi to Zslice
  int getPhiSliceID(float phi)         const {
    assert(phi>=0.f); // temporary, to remove
    //    BringTo02Pi(phi);
    return phi*mDPhiInv;
  }
  
  // lower boundary of phi slice
  float getPhiSliceMin(int id)         const {return id*mDPhi;}

  // upper boundary of phi slice
  float getPhiSliceMax(int id)         const {return (id+1)*mDPhi;}
  
  
 protected:
  

  
  float mRMin = 0.f;       ///< min radius
  float mRMax = 0.f;       ///< max radius
  float mZMin = 0.f;       ///< min Z
  float mZMax = 0.f;       ///< max
  float mNZSlices = 0.f;   ///< number of Z slices
  float mNPhiSlices = 0.f; ///< number of phi slices
  //
  // auxiliary varaibles
  float mRMin2 = 0.f;      ///< squared min r
  float mRMax2 = 0.f;      ///< squared max r
  float mDZ = 0.f;         ///< Z slice thickness
  float mDZInv = 0.f;      ///< Z slice thickness inverse
  float mDPhi = 0.f;       ///< phi slice thickness
  float mDPhiInv = 0.f;    ///< phi slice thickness inverse

  std::vector<MatCell> mCells; ///< vector of mat.budget per cell
  
  ClassDefNV(MatLayerCyl,1);
};


//________________________________________________________________________________
MatLayerCyl::MatLayerCyl(float rMin,float rMax,float zMin,float zMax,int nz,int nphi)
{
  // main constructor
  initSegmentation(rMin,rMax,zMin,zMax,nz,nphi);
}

//________________________________________________________________________________
void MatLayerCyl::initSegmentation(float rMin,float rMax,float zMin,float zMax,int nz,int nphi)
{
  // recalculate aux parameters
  mRMin = rMin;
  mRMax = rMax;
  mZMin = zMin;
  mZMax = zMax;
  mNZSlices = nz;
  mNPhiSlices = nphi;
  
  assert(mRMin<mRMax);
  assert(mZMin<mZMax);
  assert(mNZSlices>0);
  assert(mNPhiSlices>0);

  mRMin2 = mRMin*mRMin;
  mRMax2 = mRMax*mRMax;
  
  mDZ = (mZMax-mZMin) / mNZSlices;
  mDZInv = 1.f/mDZInv;

  mDPhi = o2::Base::Constants::k2PI / mNPhiSlices;
  mDPhiInv = 1.f/mDPhiInv;
  //
  int nCells = mNZSlices*mNPhiSlices;
  mCells.reserve(nCells);
}


/**********************************************************************
 *                                                                    *
 * Ray parameterized via itst endpoints as                            *
 * Vi = Vi0 + t*(Vi1-Vi0), with Vi (i=0,1,2) for global X,Y,Z         *
 * and 0 < t < 1                                                      *
 *                                                                    *
 **********************************************************************/

class Ray {
  
 public:

  using CrossPar = std::pair<float,float>;
  
  // get material budget
  void getMatBudget(const MatLayerCyl& lr) const;

  
  static constexpr float InvalidT = -1e9;
  static constexpr float Tiny = 1e-9;

  Ray(const Point3D<float> point0,const Point3D<float> point1);

  bool  crossCircleR(float r2, CrossPar& cross) const;
  int   crossLayerR(const MatLayerCyl& lr, std::array<CrossPar,2>& cross) const;
  float crossRadial(float cs, float sn) const;
  float crossZ(float z) const;
  
  Point3D<float> getPos(float t) const {
    return Point3D<float>(mP0.X()+t*mDx,mP0.Y()+t*mDy,mP0.Z()+t*mDz);
  }
    
 private:
  
  Point3D<float> mP0;         ///< entrance point
  Point3D<float> mP1;         ///< exit point 
  float mDx = 0.f;            ///< X distance
  float mDy = 0.f;            ///< Y distance
  float mDz = 0.f;            ///< Z distance
  float mDist2 = 0.f;         ///< dist^2 between points in XY plane
  float mDist2i = 0.f;        ///< inverse dist^2 between points in XY plane
  float mXDxPlusYDy = 0.f;    ///< aux x0*DX+y0*DY
  float mXDxPlusYDyRed = 0.f; ///< aux (x0*DX+y0*DY)/mDist2
  float mXDxPlusYDy2 = 0.f;   ///< aux (x0*DX+y0*DY)^2  
  float mR02 = 0.f;           ///< radius^2 of mP0
 
 ClassDefNV(Ray,1);
};

//______________________________________________________
Ray::Ray(const Point3D<float> point0,const Point3D<float> point1)
:mP0(point0)
,mP1(point1)
  ,mDx(point1.X()-point0.X())
  ,mDy(point1.Y()-point0.Y())
  ,mDz(point1.Z()-point0.Z())
{
  mDist2 = mDx*mDx + mDy*mDy;
  mDist2i = mDist2>0 ? 1.f/mDist2 : 0.f;
  mXDxPlusYDy = point0.X()*mDx + point0.Y()*mDy;
  mXDxPlusYDyRed = -mXDxPlusYDy*mDist2i;
  mXDxPlusYDy2 = mXDxPlusYDy*mXDxPlusYDy;
  mR02 = point0.Perp2();
} 


//______________________________________________________
inline bool Ray::crossCircleR(float r2, CrossPar& cross) const
{
  // calculate parameters t of intersection with circle of radius r^2
  // calculated as solution of equation
  // t^2*mDist2 +- sqrt( mXDxPlusYDy^2 - mDist2*(mR02 - r^2) )
  // 
  float det = mXDxPlusYDy2 - mDist2*(mR02 - r2);
  if (det<0) return false;  // no intersection
  float detRed = std::sqrt(det)*mDist2i;
  cross.first  = mXDxPlusYDyRed + detRed; // (-mXDxPlusYDy + det)*mDist2i;
  cross.second = mXDxPlusYDyRed - detRed; // (-mXDxPlusYDy - det)*mDist2i;
  return true;
}

//______________________________________________________
inline int Ray::crossLayerR(const MatLayerCyl& lr, std::array<CrossPar,2>& cross) const
{
  // Calculate parameters t of intersection with cyl.layer of rmin^2 and rmax^2.
  // It is assumed rmin<rmax (not checked)
  // Calculated as solution of equation
  // t^2*mDist2 +- sqrt( mXDxPlusYDy^2 - mDist2*(mR02 - r^2) )
  // Region of valid t is 0:1.
  // Straigh line may have 2 crossings with cyl. layer 

  float detMax = mXDxPlusYDy2 - mDist2*(mR02 - lr.getRMax2());
  if (detMax<0) return 0;  // does not reach outer R, hence inner also
  float detMaxRed = std::sqrt(detMax)*mDist2i;  
  float tCross0Max = mXDxPlusYDyRed + detMaxRed; // largest possible t
  
  if (tCross0Max<0) { // max t is outside of the limiting point -> other t's also
    return 0;
  }

  int nCross = 0;
  float tCross0Min = mXDxPlusYDyRed - detMaxRed; // smallest possible t
  if (tCross0Min>1.f) { // min t is outside of the limiting point -> other t's also
    return nCross;
  }

  float detMin = mXDxPlusYDy2 - mDist2*(mR02 - lr.getRMin2());
  if (detMin<0) {  // does not reach inner R -> just 1 tangential crossing
    cross[nCross].first  = tCross0Min>0.f ? tCross0Min : 0.f;
    cross[nCross].second = tCross0Max<1.f ? tCross0Max : 1.f;
    return ++nCross;  
  }
  float detMinRed = std::sqrt(detMin)*mDist2i;
  float tCross1Max = mXDxPlusYDyRed + detMinRed;
  float tCross1Min = mXDxPlusYDyRed - detMinRed;

  if (tCross1Max<1.f) {
    cross[nCross].first  = tCross0Max<1.f ? tCross0Max:1.f;
    cross[nCross].second = tCross1Max>0.f ? tCross1Max:0.f;
    nCross++;
  }
  
  if (tCross1Min>-0.f) {
    cross[nCross].first  = tCross1Min<1.f ? tCross1Min:1.f;
    cross[nCross].second = tCross0Min>0.f ? tCross0Min:0.f;
    nCross++;
  }

  return nCross;
}

//______________________________________________________
inline float Ray::crossRadial(float cs, float sn) const
{
  // calculate t of crossing with radial line with inclination cosine and sine
  float den = mDx*sn - mDy*cs;
  if (std::abs(den)<Tiny) return InvalidT;
  return (mP0.Y()*cs - mP0.X()*sn)/den;
}

//______________________________________________________
inline float Ray::crossZ(float z) const
{
  // calculate t of crossing XY plane at Z
  return std::abs(mDz)>Tiny ? (z-mP0.Z())/mDz : InvalidT;
}


//______________________________________________________
void Ray::getMatBudget(const MatLayerCyl& lr) const
{
  
}
