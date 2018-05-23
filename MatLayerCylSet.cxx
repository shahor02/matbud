// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file MatLayerCylSet.cxx
/// \brief Implementation of the wrapper for the set of cylindrical material layers 

#include "MatLayerCylSet.h"
#include "Ray.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "CommonConstants/MathConstants.h"
#include "DetectorsBase/GeometryManager.h"
#include <FairLogger.h>
#include <TFile.h>


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
  // merge similar (whose relative budget does not differ within maxRelDiff) phi slices
  for (int i=getNLayers();i--;) mLayers[i].optimizePhiSlices(maxRelDiff);
}

//________________________________________________________________________________
void MatLayerCylSet::populateFromTGeo(int ntrPerCel)
{
  ///< populate layers, using ntrPerCell test tracks per cell
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
  Ray ray(point0,point1);
  short lmin,lmax; // get innermost and outermost relevant layer
  if (!getLayersRange(ray,lmin,lmax)) {
    return rval;
  }
  short lrID = lmax;
  std::array<Ray::CrossPar,2> tpar;
  float tAcc = 0.f; // accumulated t over layers
  while (lrID>=lmin) { // go from outside to inside
    const auto& lr = getLayer(lrID);
    int nc = ray.crossLayer(lr);
    for (int ic=nc;ic--;) {
      auto &cross = ray.getCrossParams(ic); // tmax,tmin of crossing the layer
      auto phi0 = ray.getPhi(cross.first), phi1 = ray.getPhi(cross.second), dPhi = phi0 - phi1;
      auto phiID = lr.getPhiSliceID( phi0 ), phiIDLast = lr.getPhiSliceID( phi1 );
      // account for eventual wrapping around 0
      if (dPhi>0.f) { 
	if (dPhi>o2::constants::math::PI) { // wraps around phi=0
	  phiIDLast += lr.getNPhiSlices();
	}
      }
      else {
	if (dPhi<-o2::constants::math::PI) { // wraps around phi=0
	  phiID += lr.getNPhiSlices();
	}
      }
      int stepPhiID = phiID > phiIDLast ? -1 : 1; 
      bool checkMorePhi = true;
      auto tStartPhi = cross.first, tEndPhi = 0.f;
      do {
	// get the path in the current phi slice
	if (phiID==phiIDLast) {
	  tEndPhi = cross.second;
	  checkMorePhi = false;
	}
	else { // last phi slice still not reached
	  tEndPhi = ray.crossRadial( lr, (stepPhiID>0 ? phiID+1 : phiID) % lr.getNPhiSlices() );
	}
	auto zID = lr.getZBinID(ray.getZ(tStartPhi));
	auto zIDLast = lr.getZBinID(ray.getZ(tEndPhi));
	// check if Zbins are crossed
	printf("-- Zdiff (%3d : %3d) mode: t: %+e %+e\n",zID,zIDLast,tStartPhi,tEndPhi);
	if (zID!=zIDLast) {
	  auto stepZID = zID<zIDLast ? 1:-1;
	  bool checkMoreZ = true;
	  auto tStartZ = tStartPhi, tEndZ = 0.f;
	  do {
	    if (zID==zIDLast) {
	      tEndZ = tEndPhi;
	      checkMoreZ = false;
	    }
	    else {
	      tEndZ = ray.crossZ( lr.getZBinMin( stepZID>0 ? zID+1 : zID) );
	    }
	    auto pos0 = ray.getPos(tStartZ);
	    auto pos1 = ray.getPos(tEndZ);
	    printf("Lr#%3d / cross#%d : account %f<t<%f at phiSlice %d | Zbin: %3d (%3d) |[%+e %+e +%e]:[%+e %+e %+e]\n",
		   lrID,ic,tEndZ,tStartZ,phiID%lr.getNPhiSlices(),zID,zIDLast,
		   pos0.X(),pos0.Y(),pos0.Z(), pos1.X(),pos1.Y(),pos1.Z() );	    
	    tStartZ = tEndZ;
	    zID += stepZID;
	  } while(checkMoreZ);
	}
	else {
	  auto pos0 = ray.getPos(tStartPhi);
	  auto pos1 = ray.getPos(tEndPhi);

	  printf("Lr#%3d / cross#%d : account %f<t<%f at phiSlice %d | Zbin: %3d ----- |[%+e %+e +%e]:[%+e %+e %+e]\n",
		 lrID,ic,tEndPhi,tStartPhi,phiID%lr.getNPhiSlices(),zID,
		 pos0.X(),pos0.Y(),pos0.Z(), pos1.X(),pos1.Y(),pos1.Z() );
	}
	//
	tStartPhi = tEndPhi;
	phiID += stepPhiID;
	
      } while(checkMorePhi);
    }    
    lrID--;
  } // loop over layers
  
  return rval;
}

//_________________________________________________________________________________________________
bool MatLayerCylSet::getLayersRange(const Ray& ray, short& lmin,short& lmax) const
{
  // get range of layers corresponding to rmin/rmax
  //  
  // tmp,slow, TODO
  lmin = lmax = -1;
  float rmin2,rmax2;
  ray.getMinMaxR2(rmin2,rmax2);
  
  if (rmin2>=getRMax2() || rmax2<=getRMin2()) {
    return false;
  }

  for (lmin=0;lmin<getNLayers();lmin++) {
    const auto & lr =  getLayer(lmin);
    if (lr.getRMax2()>rmin2) {
      break;
    }
  }
  for (lmax=lmin;lmax<getNLayers();lmax++) {
    const auto & lr =  getLayer(lmin);
    if (lr.getRMax2()>rmax2-Ray::Tiny) {
      break;
    }
  }
  if (lmax==getNLayers()) lmax = getNLayers()-1;
  return true;
}
