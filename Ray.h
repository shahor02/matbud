#include "MathUtils/Cartesian3D.h"
#include <cmath>
#include <utility>
#include <array>
#include <cassert>
#include "CommonConstants/MathConstants.h"
#include "MathUtils/Utils.h"
#include <TGeoManager.h>
#include <TFile.h>
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/TreeStreamRedirector.h"


#include <FairLogger.h>

class Ray;
class MatLayerCyl;
class MatLayerCylSet;

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
  float mRho = 0.f;      ///< mean density, g/cm^3
  float mX2X0  = 0.f;    ///< fraction of radiaton lenght
  MatCell operator+(const MatCell& rhs) {
    MatCell cell;
    cell.mRho = this->mRho + rhs.mRho;
    cell.mX2X0 = this->mX2X0 + rhs.mX2X0;
    return cell;
  }
  MatCell& operator+=(const MatCell& rhs) {
    mRho += rhs.mRho;
    mX2X0 += rhs.mX2X0;
    return *this;
  }
  void scale(float scale) {
    mRho *= scale;
    mX2X0 *= scale;
  }

  
  ClassDefNV(MatCell,1);
};


/**********************************************************************
 *                                                                    *
 * Cylindrical material layer                                         *
 *                                                                    *
 * Logical division is given by mNZBins and mNPhiBins                 *
 * but the actual number phi slices might be less if neighbouring     *
 * phi bins with similar properties are merged together. The actual   *
 * phi slice should be accessed via mPhiBin2Slice map.                *
 *                                                                    *
 **********************************************************************/
class MatLayerCyl {

 public:

  MatLayerCyl() = default;
  MatLayerCyl(float rMin,float rMax,float zMin,float zMax,short nz,short nphi);
  MatLayerCyl(float rMin,float rMax,float zHalfSpan, float dzMin,float drphiMin);
  ~MatLayerCyl() = default;
  
  // ----------------------- segmentation 
  void  initSegmentation(float rMin,float rMax,float zMin,float zMax,short nz,short nphi);
  float getRMin()        const {return mRMin;}
  float getRMax()        const {return mRMax;}
  float getZMin()        const {return mZMin;}
  float getZMax()        const {return mZMax;}
  short getNZBins()      const {return mNZBins;}
  short getNPhiBins()    const {return mNPhiBins;}
  // actual number of phi slices stored
  short getNPhiSlices()  const {return mNPhiSlices;}
  short getNPhiBinsInSlice(int i) const;
  
  float getRMin2()       const {return mRMin2;}
  float getRMax2()       const {return mRMax2;}

  void print(bool data=false) const;
  void populateFromTGeo(int ntrPerCell=10);
  void populateFromTGeo(short ip, short iz, int ntrPerCell);
  
  // linearized cell ID
  int getCellID(short iphi,short iz) const {
    return int(phiBin2Slice(iphi))*mNZBins + iz;
  }
  
  // obtain material cell, cell ID must be valid
  const MatCell& getCell(short iphi, short iz)  const {return mCells[getCellID(iphi,iz)];}
  MatCell& getCell(short iphi, short iz)              {return mCells[getCellID(iphi,iz)];}
  
  // ---------------------- Z slice manipulation
  // convert Z to Zslice
  int getZBinID(float z)        const {return z>mZMin && z<mZMax ? int((z-mZMin)*mDZInv) : -1;}

  // lower boundary of Z slice
  float getZBinMin(short id)  const {return mZMin + id*mDZ;}

  // upper boundary of Z slice
  float getZBinMax(short id)  const {return mZMin + (id+1)*mDZ;}

  // ---------------------- Phi slice manipulation (0:2pi convention, no check is done)
  // convert Phi (in 0:2pi convention) to PhiBinID
  short getPhiBinID(float phi)  const {
    assert(phi>=0.f); // temporary, to remove
    return short(phi*mDPhiInv);
  }
  short phiBin2Slice(short i)    const {return mPhiBin2Slice[i];}
  short getPhiSliceID(float phi)   const { return phiBin2Slice(getPhiBinID(phi));}
  short getNPhiBinsInSlice(short iSlice, short &binMin, short &binMax) const;
  
  // lower boundary of phi slice
  float getPhiBinMin(short id)     const {return id*mDPhi;}

  // upper boundary of phi slice
  float getPhiBinMax(short id)     const {return (id+1)*mDPhi;}
  
  const std::vector<MatCell> & getCells() const {return mCells;}

  void getMeanRMS(MatCell &mean, MatCell &rms) const;


  bool cellsDiffer(const MatCell& cellA, const MatCell& cellB, float maxRelDiff) const;
  bool canMergePhiSlices(int i,int j, float maxRelDiff=0.05, int maxDifferent=1) const;
  void optimizePhiSlices(float maxRelDiff=0.05);
  std::size_t getSize() const {
    return sizeof(*this) + getNPhiBins()*sizeof(short) + getNPhiSlices()*getNZBins()*sizeof(MatCell);
  }
  
 protected:
  
  float mRMin = 0.f;       ///< min radius
  float mRMax = 0.f;       ///< max radius
  float mZMin = 0.f;       ///< min Z
  float mZMax = 0.f;       ///< max
  short mNZBins = 0;       ///< number of Z bins
  short mNPhiBins = 0;     ///< number of phi bins (logical)
  short mNPhiSlices = 0;   ///< actual number of phi slices
  //
  // auxiliary varaibles
  float mRMin2 = 0.f;      ///< squared min r
  float mRMax2 = 0.f;      ///< squared max r
  float mDZ = 0.f;         ///< Z slice thickness
  float mDZInv = 0.f;      ///< Z slice thickness inverse
  float mDPhi = 0.f;       ///< phi slice thickness
  float mDPhiInv = 0.f;    ///< phi slice thickness inverse

  std::vector<short> mPhiBin2Slice; ///< mapping from analytical phi bin ID to real slice ID
  std::vector<MatCell> mCells; ///< vector of mat.budget per cell
  
  ClassDefNV(MatLayerCyl,1);
};

//==================================================================================
/**********************************************************************
 *                                                                    *
 * Set of cylindrical material layer                                  *
 *                                                                    *
 **********************************************************************/
class MatLayerCylSet {

 public:
  MatLayerCylSet() = default;
  ~MatLayerCylSet() = default;
  int getNLayers() const {return mLayers.size();}
  void addLayer(float rmin,float rmax, float zmax, float dz, float drphi);
  const MatLayerCyl& getLayer(int i) const {return mLayers[i];}
  const vector<MatLayerCyl>& getLayers() const {return mLayers;}

  float getRMin() const {return mRMin;}
  float getRMax() const {return mRMax;}
  float getZMax() const {return mZMax;}
  float getRMin2() const {return mRMin2;}
  float getRMax2() const {return mRMax2;}
  
  void print(bool data=false) const;
  void populateFromTGeo(int ntrPerCel=10);

  void dumpToTree(const std::string outName = "matBudTree.root") const;

  void writeToFile(std::string outFName = "matBudg.root", std::string name = "MatBud");
  static MatLayerCylSet* loadFromFile(std::string inpFName = "matBudg.root", std::string name = "MatBud");

  void optimizePhiSlices(float maxRelDiff=0.05);

  std::size_t getSize() const;

  MatCell getMatBudget(const Point3D<float> &point0,const Point3D<float> &point1) const;
  MatCell getMatBudget(float x0, float y0, float z0, float x1, float y1, float z1) const;
  
 protected:
  float mRMin = 0.f; ///< min radius
  float mRMax = 0.f; ///< max radius
  float mZMax  =0.f; ///< max Z span
  float mRMin2 = 0.f;
  float mRMax2 = 0.f; 
  vector<MatLayerCyl> mLayers; ///< set of cylinrical layers

  ClassDefNV(MatLayerCylSet,1);
};


//==================================================================================

//________________________________________________________________________________
MatLayerCyl::MatLayerCyl(float rMin,float rMax,float zMin,float zMax,short nz,short nphi)
{
  // main constructor
  initSegmentation(rMin,rMax,zMin,zMax,nz,nphi);
}

//________________________________________________________________________________
MatLayerCyl::MatLayerCyl(float rMin,float rMax,float zHalfSpan, float dzMin,float drphiMin)
{
  // main constructor
  if (dzMin<0.001f) dzMin = 0.001f;
  if (drphiMin<0.001f) drphiMin = 0.001f;
  float peri = (rMax+rMin)*o2::constants::math::PI;
  int nz = 2*zHalfSpan/dzMin, nphi = peri/drphiMin;
  initSegmentation(rMin,rMax,-zHalfSpan,zHalfSpan, nz>0 ? nz:1, nphi>0 ? nphi:1);
}

//________________________________________________________________________________
void MatLayerCyl::initSegmentation(float rMin,float rMax,float zMin,float zMax,short nz,short nphi)
{
  // recalculate aux parameters
  mRMin = rMin;
  mRMax = rMax;
  mZMin = zMin;
  mZMax = zMax;
  mNZBins = nz;
  mNPhiBins = nphi;
  
  assert(mRMin<mRMax);
  assert(mZMin<mZMax);
  assert(mNZBins>0);
  assert(mNPhiBins>0);

  mRMin2 = mRMin*mRMin;
  mRMax2 = mRMax*mRMax;
  
  mDZ = (mZMax-mZMin) / mNZBins;
  mDZInv = 1.f/mDZInv;

  mDPhi = o2::constants::math::TwoPI / mNPhiBins;
  mDPhiInv = 1.f/mDPhiInv;
  //
  mNPhiSlices = mNPhiBins;
  mPhiBin2Slice.resize(mNPhiBins);
  int nCells = int(mNZBins)*mNPhiSlices;
  mCells.resize(nCells);
  for (int i=mNPhiSlices;i--;) {
    mPhiBin2Slice[i] = i; // fill with trivial mapping
  }
}

//________________________________________________________________________________
inline short MatLayerCyl::getNPhiBinsInSlice(short iSlice, short &binMin, short &binMax) const
{
  // slow method to get number of phi bins for given phi slice
  int nb = 0;
  binMin = binMax = -1;
  for (int ib=mNPhiBins;ib--;) {
    if (mPhiBin2Slice[ib]==iSlice) {
      binMax<0 ? binMin = binMax=ib : binMin=ib;
      nb++;
      continue;
    }
    if (binMax>=0) {
      break; // no more bins since they are consecutive
    }
  }
  return nb;
}

//________________________________________________________________________________
void MatLayerCyl::populateFromTGeo(int ntrPerCell)
{
  /// populate layer with info extracted from TGeometry
  ntrPerCell = ntrPerCell>1 ? ntrPerCell : 1;
  for (int iz=mNZBins;iz--;) {
    for (int ip=mNPhiBins;ip--;) {
      populateFromTGeo(ip, iz ,ntrPerCell);
    }
  }
}

//________________________________________________________________________________
void MatLayerCyl::populateFromTGeo(short ip, short iz, int ntrPerCell)
{
  /// populate cell with info extracted from TGeometry
  
  float zmn = getZBinMin(iz), phmn = getPhiBinMin(ip), sn,cs;
  double meanRho = 0., meanX2X0 = 0., lgt = 0.;;
  for (int isz=ntrPerCell;isz--;) {
    float zs = zmn + (isz+0.5)*mDZ/ntrPerCell;
    for (int isp=ntrPerCell;isp--;) {
      o2::utils::sincosf(phmn + (isp+0.5)*mDPhi/ntrPerCell, sn,cs);
      auto bud = o2::Base::GeometryManager::MeanMaterialBudget(mRMin*cs,mRMin*sn,zs,mRMax*cs,mRMax*sn,zs);
      if (bud.length>0.) {
	meanRho += bud.length*bud.meanRho;
	meanX2X0 += bud.length*bud.meanX2X0;
	lgt += bud.length;
      }	  
    }
  }
  if (lgt>0.) {
    auto &cell = getCell(ip,iz);
    cell.mRho = meanRho/lgt;
    cell.mX2X0 = meanX2X0/lgt;
  }  
}

//________________________________________________________________________________
bool MatLayerCyl::canMergePhiSlices(int i,int j, float maxRelDiff, int maxDifferent) const
{
  if (std::abs(i-j)>1 || i==j || std::max(i,j)>=mNPhiSlices) {
    LOG(ERROR)<<"Only existing "<<mNPhiSlices<<" slices with difference of 1 can be merged, input is "
	      <<i<<" and "<<FairLogger::endl;
    return false;
  }
  int ndiff = 0; // number of different cells
  for (int iz=getNZBins();iz--;) {
    const auto &cellI = getCell(i,iz);
    const auto &cellJ = getCell(j,iz);
    if (cellsDiffer(cellI, cellJ, maxRelDiff)) {
      if (++ndiff>maxDifferent) {
	return false;
      }
    }
  }
  return true;
}

//________________________________________________________________________________
bool MatLayerCyl::cellsDiffer(const MatCell& cellA, const MatCell& cellB, float maxRelDiff) const
{
  /// check if the cells content is different within the tolerance
  float rav = 0.5*(cellA.mRho + cellB.mRho), xav = 0.5*(cellA.mX2X0 + cellB.mX2X0);
  float rdf = 0.5*(cellA.mRho - cellB.mRho), xdf = 0.5*(cellA.mX2X0 - cellB.mX2X0);
  if (rav>0 && std::abs(rdf/rav)>maxRelDiff) return true;
  if (xav>0 && std::abs(xdf/xav)>maxRelDiff) return true;
  return false;
}

//________________________________________________________________________________
void MatLayerCyl::optimizePhiSlices(float maxRelDiff)
{
  // merge compatible phi slices
  if (mNPhiSlices<mNPhiBins) {
    LOG(ERROR)<<mNPhiBins<<" phi bins were already merged to "<<mNPhiSlices<<" slices"<<FairLogger::endl;
    return;
  }
  short newSl = 0;
  for (int is=1;is<mNPhiSlices;is++) {
    if (!canMergePhiSlices(is-1,is,maxRelDiff)) {
      newSl++;
    }
    mPhiBin2Slice[is] = newSl;
  }
  LOG(INFO)<<newSl+1<<" slices out of "<<mNPhiBins<<FairLogger::endl;
  if (newSl+1==mNPhiSlices) {
    return;
  }
  newSl = 0;
  short slMin=0,slMax=0, is=0;
  while (is++<mNPhiSlices) {
    while (is<mNPhiSlices && mPhiBin2Slice[is]==newSl) { // select similar slices 
      slMax++;
      is++;
    }
    if (slMax>slMin || newSl!=slMin) {  // merge or shift slices
      mCells[newSl] = mCells[slMin];
      for (int ism=slMin+1;ism<=slMax;ism++) {
	mCells[newSl] += mCells[ism];
      }
      mCells[newSl].scale(1.f/(1.f+slMax-slMin));
      LOG(INFO)<<"mapping "<<slMin<<":"<<slMax<<" to new slice "<<newSl<<FairLogger::endl;
    }
    newSl++;
    slMin = slMax = is;
  }
  mNPhiSlices = newSl;
  LOG(INFO)<<"Updated Nslices = "<<mNPhiSlices<<FairLogger::endl;
  mCells.resize(mNPhiSlices);
}
  
//________________________________________________________________________________
void MatLayerCyl::print(bool data) const
{
  ///< print layer data
  float szkb = float(getSize())/1024;
  printf("Cyl.Layer %.3f<r<%.3f %+.3f<Z<%+.3f | Nphi: %5d (%d slices) Nz: %5d Size: %.3f KB\n",
	 mRMin,mRMax,mZMin,mZMax,mNPhiBins, getNPhiSlices() ,mNZBins, szkb);
  if (!data) {
    return;
  }
  for (int ip=0;ip<getNPhiSlices();ip++) {
    short ib0,ib1;
    short nb = getNPhiBinsInSlice(ip,ib0,ib1);
    printf("phi slice: %d (bins %d-%d %.4f:%.4f) ... [iz/<rho>/<x/x0>] \n",ip,ib0,ib1,mDPhi*ib0,mDPhi*(ib1+1));
    for (int iz=0;iz<mNZBins;iz++) {
      auto cell = getCell(ib0,iz);
      printf("%3d/%.2e/%.2e ",iz,cell.mRho,cell.mX2X0);
      if (((iz+1)%5)==0) {
	printf("\n");
      }
    }
    if (mNZBins%5) {
      printf("\n");
    }
  }
}

//________________________________________________________________________________
void MatLayerCyl::getMeanRMS(MatCell &mean, MatCell &rms) const
{
  // mean and RMS over layer
  mean.mRho = rms.mRho = 0.f;
  mean.mX2X0 = rms.mX2X0 = 0.f;
  for (int ip=mNPhiBins;ip--;) {
    for (int iz=mNZBins;iz--;) {
      const auto& cell = getCell(ip,iz);
      mean.mRho += cell.mRho;
      mean.mX2X0 += cell.mX2X0;
      rms.mRho += cell.mRho*cell.mRho;
      rms.mX2X0 += cell.mX2X0*cell.mX2X0;
    }
  }
  int nc = mNPhiBins*mNZBins;
  mean.mRho /= nc;
  mean.mX2X0 /= nc;  
  rms.mRho /= nc;
  rms.mX2X0 /= nc;  
  rms.mRho -= mean.mRho*mean.mRho;
  rms.mX2X0 -= mean.mX2X0*mean.mX2X0;
  rms.mRho = rms.mRho>0.f ? std::sqrt(rms.mRho) : 0.f;
  rms.mX2X0 = rms.mX2X0>0.f ? std::sqrt(rms.mX2X0) : 0.f;  
}


//==========================================================

//________________________________________________________________________________
void MatLayerCylSet::addLayer(float rmin,float rmax, float zmax, float dz, float drphi)
{
  // add new layer checking for overlaps
  assert(rmin<rmax && zmax>0 && dz>0 && drphi>0);
  for (int il=0;il<getNLayers();il++) {
    const auto &lr = getLayer(il);
    if (lr.getRMax()>rmin && rmax>lr.getRMin()) {
      LOG(FATAL)<<"new layer overlaps with layer "<<il<<FairLogger::endl;
    }
  }
  mLayers.emplace_back(rmin,rmax, zmax, dz, drphi);
  mRMin = mRMin>rmin ? rmin : mRMin;
  mRMax = mRMax<rmax ? rmax : mRMax;
  mZMax = mZMax<zmax ? zmax : mZMax;
  mRMin2 = mRMin*mRMin;
  mRMax2 = mRMax*mRMax;  
}

//________________________________________________________________________________
void MatLayerCylSet::print(bool data) const
{
  ///< print layer data
  for (int i=0;i<getNLayers();i++) {
    printf("#%3d | ",i);
    mLayers[i].print(data);
  }
  printf("%.2f < R < %.2f  %d layers with total size %.2f MB\n",mRMin,mRMax,getNLayers(),float(getSize())/1024/1024);
}

void MatLayerCylSet::optimizePhiSlices(float maxRelDiff)
{
  // merge similar phi slices
  for (int i=getNLayers();i--;) mLayers[i].optimizePhiSlices(maxRelDiff);
}

//________________________________________________________________________________
void MatLayerCylSet::populateFromTGeo(int ntrPerCel)
{
  ///< populate layers
  if (!gGeoManager || !gGeoManager->IsClosed()) {
    LOG(ERROR) << "No active geometry or geometry not yet closed!" << FairLogger::endl;
    return;
  }
  for (int i=0;i<getNLayers();i++) {
    int nz = mLayers[i].getNZBins(), np = mLayers[i].getNPhiBins();
    LOG(INFO)<<"populating layer " << i << " NZ: "<<nz<<" NPhi: "<<np<<FairLogger::endl;
    mLayers[i].populateFromTGeo(ntrPerCel);
  }
}

//________________________________________________________________________________
std::size_t MatLayerCylSet::getSize() const
{
  std::size_t sz = sizeof(*this);
  for (int i=0;i<getNLayers();i++) {
    sz += mLayers[i].getSize();
  }
  return sz;
}

//________________________________________________________________________________
void MatLayerCylSet::dumpToTree(const std::string outName) const
{
  /// dump per cell info to the tree
  
  o2::utils::TreeStreamRedirector dump(outName.data(),"recreate");
  for (int i=0;i<getNLayers();i++) {
    const auto & lr = getLayer(i);
    float r = 0.5*(lr.getRMin()+lr.getRMax());
    // per cell dump
    char merge = 0;
    int nphib = lr.getNPhiBins();
    for (int ip=0;ip<nphib;ip++) {
      float phi = 0.5*( lr.getPhiBinMin(ip)+lr.getPhiBinMax(ip) );
      float sn,cs;
      int ips = lr.phiBin2Slice(ip);
      merge = 0; // not mergeable
      if (ip+1<nphib) {
	int ips1 = lr.phiBin2Slice(ip+1);
	merge = ips==ips1 ? -1 : lr.canMergePhiSlices(ips,ips1); // -1 for already merged
      }
      else merge = -2; // last one
      o2::utils::sincosf(phi,sn,cs);
      float x = r*cs, y = r*sn;
      for (int iz=0;iz<lr.getNZBins();iz++) {
	float z = 0.5*( lr.getZBinMin(iz)+lr.getZBinMax(iz) );
	auto cell = lr.getCell(ip,iz);
	dump<<"cell"<<"ilr="<<i<<"r="<<r<<"phi="<<phi<<"x="<<x<<"y="<<y<<"z="<<z<<"ip="<<ip<<"ips="<<ips<<"iz="<<iz
	    <<"mrgnxt="<<merge<<"val="<<cell<<"\n";
      }
    }
    //
    // statistics per layer
    MatCell mean,rms;
    lr.getMeanRMS(mean,rms);
    dump<<"lay"<<"ilr="<<i<<"r="<<r<<"mean="<<mean<<"rms="<<rms<<"\n";
  }
}

//________________________________________________________________________________
void MatLayerCylSet::writeToFile(std::string outFName, std::string name)
{
  /// store to file
  
  TFile outf(outFName.data(), "recreate");
  if (outf.IsZombie()) {
    return;
  }
  if (name.empty()) {
    name = "matBud";
  }
  outf.WriteObjectAny(this, Class(), name.data());
  outf.Close();
}

//________________________________________________________________________________
MatLayerCylSet* MatLayerCylSet::loadFromFile(std::string inpFName, std::string name)
{
  if (name.empty()) {
    name = "MatBud";
  }
  TFile inpf(inpFName.data());
  if (inpf.IsZombie()) {
    LOG(ERROR) << "Failed to open input file " << inpFName << FairLogger::endl;
    return nullptr;
  }
  MatLayerCylSet* mb = reinterpret_cast<MatLayerCylSet*>(inpf.GetObjectChecked(name.data(),Class()));
  if (!mb) {
    LOG(ERROR) << "Failed to load " << name << " from " << inpFName << FairLogger::endl;
  }
  return mb;
}


//_________________________________________________________________________________________________
MatCell MatLayerCylSet::getMatBudget(float x0, float y0, float z0, float x1, float y1, float z1) const
{
  // get material budget traversed on the line between point0 and point1
  Point3D<float> point0(x0,y0,z0), point1(x1,y1,z1);
  return getMatBudget(point0,point1);
}

//_________________________________________________________________________________________________
MatCell MatLayerCylSet::getMatBudget(const Point3D<float> &point0,const Point3D<float> &point1) const
{
  // get material budget traversed on the line between point0 and point1
  MatCell rval;
  float rmin2 = point0.Perp2(), rmax2 = point1.Perp2();
  if (rmin2>=getRMax2() || rmax2<=getRMin2()) {
    return rval;
  }
  Ray ray(point0,point1);
  return rval;
}



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

