// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file MatLayerCyl.h
/// \brief Declarations for single cylindrical material layer class 

#ifndef ALICEO2_MATLAYERCYL_H
#define ALICEO2_MATLAYERCYL_H

#include "Rtypes.h"


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
  enum RangeStatus : short {Below=-1,Within=0,Above=1};
  
  MatLayerCyl() = default;
  MatLayerCyl(float rMin,float rMax,float zHalfSpan, float dzMin,float drphiMin);
  ~MatLayerCyl() = default;
  
  // ----------------------- segmentation 
  void  initSegmentation(float rMin,float rMax,float zHalfSpan,short nz,short nphi);
  float getRMin()        const {return mRMin;}
  float getRMax()        const {return mRMax;}
  float getZHalf()       const {return mZHalf;}
  float getZMin()        const {return -getZHalf();}
  float getZMax()        const {return getZHalf();}
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
  RangeStatus isZOutside(float z) const {return z<getZMin() ? Below : (z>getZMax() ? Above : Within); }
  int getZBinIDChecked(float z) const {return z<getZMin() || z>getZMax() ? -1 : int((z-getZMin())*mDZInv);}
  int getZBinID(float z)        const {return int((z-getZMin())*mDZInv);}

  // lower boundary of Z slice
  float getZBinMin(short id)  const {return getZMin() + id*mDZ;}

  // upper boundary of Z slice
  float getZBinMax(short id)  const {return getZMin() + (id+1)*mDZ;}

  // ---------------------- Phi slice manipulation (0:2pi convention, no check is done)
  // convert Phi (in 0:2pi convention) to PhiBinID
  short getPhiBinID(float phi)  const { return short(phi*mDPhiInv);}
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
  float mZHalf = 0.f;       ///< Z half span
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



#endif
