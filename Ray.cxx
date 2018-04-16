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
int Ray::Cross(const MatLayerCyl& lr)
{
  ///< calculate parameters of crossing the layer
  std::array<CrossPar,2> cross;
  int nc = crossLayerR(lr.getRMin2(),lr.getRMax2(),cross);
  if (!nc) {
    return 0;
  }
  auto p0 = ray.getPos(cross[0].first); // high t
  auto p1 = ray.getPos(cross[0].second); // low t
  if (validateZRange(cross[0],p0,p1,lr)
  
}


