// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Ray.cxx
/// \brief Call for the ray between start-end points for material budget estimate

#ifndef ALICEO2_RAY_H
#define ALICEO2_RAY_H

#include "Rtypes.h"
#include "MathUtils/Cartesian3D.h"
#include "MatLayerCyl.h"
#include "MathUtils/Utils.h"
#include <array>

/**********************************************************************
 *                                                                    *
 * Ray parameterized via its endpoints as                             *
 * Vi = Vi0 + t*(Vi1-Vi0), with Vi (i=0,1,2) for global X,Y,Z         *
 * and 0 < t < 1                                                      *
 *                                                                    *
 **********************************************************************/

class Ray {
  
 public:

  using CrossPar = std::pair<float,float>;
  
  
  static constexpr float InvalidT = -1e9;
  static constexpr float Tiny = 1e-9;

  Ray(const Point3D<float> point0,const Point3D<float> point1);
  
  int   crossLayer(const MatLayerCyl& lr);
  bool  crossCircleR(float r2, CrossPar& cross) const;

  float crossRadial(const MatLayerCyl& lr, int sliceID) const;
  float crossRadial(float cs, float sn) const;
  float crossZ(float z) const;

  const CrossPar& getCrossParams(int i) const { return mCrossParams[i];}
  
  Point3D<float> getPos(float t) const {
    return Point3D<float>(mP0.X()+t*mDx,mP0.Y()+t*mDy,mP0.Z()+t*mDz);
  }

  float getPhi(float t) const {
    float p = std::atan2(mP0.Y()+t*mDy,mP0.X()+t*mDx);
    o2::utils::BringTo02Pi(p);
    return p;
  }
  
  float getZ(float t) const { return mP0.Z()+t*mDz; }
  
  bool validateZRange(CrossPar& cpar, const MatLayerCyl& lr) const;
  
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
  std::array<CrossPar,2> mCrossParams; ///< parameters of crossing the layer

  ClassDefNV(Ray,1);
};


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
inline int Ray::crossLayer(const MatLayerCyl& lr)
{
  // Calculate parameters t of intersection with cyl.layer
  // Calculated as solution of equation for ray crossing with circles of r (rmin and rmax)
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

  float tCross0Min = mXDxPlusYDyRed - detMaxRed; // smallest possible t
  if (tCross0Min>1.f) { // min t is outside of the limiting point -> other t's also
    return 0;
  }
  float detMin = mXDxPlusYDy2 - mDist2*(mR02 - lr.getRMin2());
  if (detMin<0) {  // does not reach inner R -> just 1 tangential crossing
    mCrossParams[0].first  = tCross0Min>0.f ? tCross0Min : 0.f;
    mCrossParams[0].second = tCross0Max<1.f ? tCross0Max : 1.f;
    return validateZRange( mCrossParams[0], lr ) ;
  }
  int nCross = 0;
  float detMinRed = std::sqrt(detMin)*mDist2i;
  float tCross1Max = mXDxPlusYDyRed + detMinRed;
  float tCross1Min = mXDxPlusYDyRed - detMinRed;

  if (tCross1Max<1.f) {
    mCrossParams[0].first  = tCross0Max<1.f ? tCross0Max:1.f;
    mCrossParams[0].second = tCross1Max>0.f ? tCross1Max:0.f;
    if (validateZRange( mCrossParams[nCross], lr )) {
      nCross++;
    }
  }
  
  if (tCross1Min>-0.f) {
    mCrossParams[nCross].first  = tCross1Min<1.f ? tCross1Min:1.f;
    mCrossParams[nCross].second = tCross0Min>0.f ? tCross0Min:0.f;
    if (validateZRange( mCrossParams[nCross], lr )) {
      nCross++;
    }
  }

  return nCross;
}

inline float Ray::crossRadial(const MatLayerCyl& lr, int sliceID) const
{
  // calculate t of crossing with phimin of layer's slice sliceID
  const auto& cs = lr.getCosSinSlice(sliceID);
  return crossRadial(cs.first,cs.second);
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
inline bool Ray::validateZRange(CrossPar& cpar, const MatLayerCyl& lr) const
{
  // make sure that estimated crossing parameters are compatible
  // with Z coverage of the layer
  MatLayerCyl::RangeStatus zout0 = lr.isZOutside(getZ(cpar.first)), zout1 = lr.isZOutside(getZ(cpar.second));
  if (zout0==zout1) { // either both points outside w/o crossing or boht inside
    return zout0==MatLayerCyl::Within ? true : false;  
  }
  // at least 1 point is outside, but the there is a crossing
  if (zout0!=MatLayerCyl::Within) {
    cpar.first = crossZ(zout0==MatLayerCyl::Below ? lr.getZMin() : lr.getZMax());
  }
  if (zout1!=MatLayerCyl::Within) {
    cpar.second = crossZ(zout1==MatLayerCyl::Below ? lr.getZMin() : lr.getZMax());
  }
  return true;
}

#endif
