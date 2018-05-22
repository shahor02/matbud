// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file MatLayerCylSet.h
/// \brief Declarations for the wrapper for the set of cylindrical material layers 

#ifndef ALICEO2_MATLAYERCYLSET_H
#define ALICEO2_MATLAYERCYLSET_H

#include "Rtypes.h"
#include "MatLayerCyl.h"
#include "MathUtils/Cartesian3D.h"
#include <vector>

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
  const std::vector<MatLayerCyl>& getLayers() const {return mLayers;}
  bool  getLayersRange(float x0, float y0, float x1, float y1, short& lmin,short& lmax) const;
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
  std::vector<MatLayerCyl> mLayers; ///< set of cylinrical layers

  ClassDefNV(MatLayerCylSet,1);
};

#endif
