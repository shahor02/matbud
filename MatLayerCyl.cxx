// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file MatLayerCyl.cxx
/// \brief Implementation of single cylindrical material layer 

#include "MatLayerCyl.h"
#include "Ray.h"
#include "MathUtils/Utils.h"
#include "CommonConstants/MathConstants.h"
#include "DetectorsBase/GeometryManager.h"

#include <array>
#include <FairLogger.h>

//________________________________________________________________________________
MatLayerCyl::MatLayerCyl(float rMin,float rMax,float zHalfSpan, float dzMin,float drphiMin)
{
  // main constructor
  if (dzMin<0.001f) dzMin = 0.001f;
  if (drphiMin<0.001f) drphiMin = 0.001f;
  float peri = (rMax+rMin)*o2::constants::math::PI;
  int nz = 2*zHalfSpan/dzMin, nphi = peri/drphiMin;
  initSegmentation(rMin,rMax,zHalfSpan, nz>0 ? nz:1, nphi>0 ? nphi:1);
}

//________________________________________________________________________________
void MatLayerCyl::initSegmentation(float rMin,float rMax,float zHalfSpan,short nz,short nphi)
{
  // recalculate aux parameters
  mRMin = rMin;
  mRMax = rMax;
  mZHalf = zHalfSpan;
  mNZBins = nz;
  mNPhiBins = nphi;
  
  assert(mRMin<mRMax);
  assert(mNZBins>0);
  assert(mNPhiBins>0);

  mRMin2 = mRMin*mRMin;
  mRMax2 = mRMax*mRMax;
  
  mDZ = (mZHalf + mZHalf) / mNZBins;
  mDZInv = 1.f/mDZ;

  mDPhi = o2::constants::math::TwoPI / mNPhiBins;
  mDPhiInv = 1.f/mDPhi;
  //
  mNPhiSlices = mNPhiBins;
  mPhiBin2Slice.resize(mNPhiBins);
  mSliceCosSin.resize(mNPhiBins);
  int nCells = int(mNZBins)*mNPhiSlices;
  mCells.resize(nCells);
  for (int i=mNPhiSlices;i--;) {
    mPhiBin2Slice[i] = i; // fill with trivial mapping
    mSliceCosSin[i].first = std::cos( getPhiBinMin(i) );
    mSliceCosSin[i].second = std::sin( getPhiBinMin(i) );
  }
}


//________________________________________________________________________________
short MatLayerCyl::getNPhiBinsInSlice(short iSlice, short &binMin, short &binMax) const
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
  /// populate layer with info extracted from TGeometry, using ntrPerCell test tracks per cell
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
  /// populate cell with info extracted from TGeometry, using ntrPerCell test tracks per cell
  
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
      mSliceCosSin[newSl] = mSliceCosSin[slMin];
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
  mSliceCosSin.resize(mNPhiSlices);
}
  
//________________________________________________________________________________
void MatLayerCyl::print(bool data) const
{
  ///< print layer data
  float szkb = float(getSize())/1024;
  printf("Cyl.Layer %.3f<r<%.3f %+.3f<Z<%+.3f | Nphi: %5d (%d slices) Nz: %5d Size: %.3f KB\n",
	 mRMin,mRMax,getZMin(),getZMax(),mNPhiBins, getNPhiSlices() ,mNZBins, szkb);
  if (!data) {
    return;
  }
  for (int ip=0;ip<getNPhiSlices();ip++) {
    short ib0,ib1;
    short nb = getNPhiBinsInSlice(ip,ib0,ib1);
    printf("phi slice: %d (%d bins %d-%d %.4f:%.4f) ... [iz/<rho>/<x/x0>] \n",ip,nb,ib0,ib1,mDPhi*ib0,mDPhi*(ib1+1));
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
