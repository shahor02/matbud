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
/// \brief Implementation of ray between start-end points for material budget estimate

#include "Ray.h"


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
int Ray::crossLayer(const MatLayerCyl& lr)
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



