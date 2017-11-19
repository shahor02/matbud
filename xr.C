#include "AliExternalTrackParam.h"
#include "TMath.h"

inline double calcDyDx(double x2R, double sn, double csInv, double cs2Inv)
{
  // expansion of (s+(s+x2r))/( sqrt(1-s^2) + sqrt(1.-(s+x2r)^2) ) up to 3d power in x2r
  // = s/(1-s^2)^(1/2) + x/2/(1-s^2)^(3/2)+s*x^2/2/(1-s^2)^(5/2)+(4s^2+1)/8/(1-s^2)^(7/2)*x^3
  double x2Rcs2Inv = x2R*cs2Inv/2.;
  return csInv*(sn + x2Rcs2Inv*(1.+2.*x2Rcs2Inv*(sn+x2Rcs2Inv*(4.*sn*sn+1))));
}

inline double calcDyDxDD(double x2R, double sn, double csInv, double cs2Inv)
{
  // expansion of (s+(s+x2r))/( sqrt(1-s^2) + sqrt(1.-(s+x2r)^2) ) up to 3d power in x2r
  // = s/(1-s^2)^(1/2) + x/2/(1-s^2)^(3/2)+s*x^2/2/(1-s^2)^(5/2)+(4s^2+1)/8/(1-s^2)^(7/2)*x^3
  double x2Rcs2Inv = x2R*cs2Inv;
  return 0.5*csInv*cs2Inv*(1. + x2Rcs2Inv*(2.*sn+3./4.*x2Rcs2Inv*(4.*sn*sn+1)));
}


double xr(double r, const AliExternalTrackParam* trc)
{
  const double* par = trc->GetParameter();
  double sn = par[2];
  double cs2Inv = 1./((1.-sn)*(1+sn));
  double csInv = TMath::Sqrt(cs2Inv);
  double crv = trc->GetC(5.);
  double dx = r-trc->GetX();

  const double kMinEps = 1e-4; 
  for (int i=0;i<20;i++) {
    double dYdX = calcDyDx(dx*crv, sn, csInv, cs2Inv );
    double ddYdX = calcDyDxDD(dx*crv, sn, csInv, cs2Inv );
    double xUpd = trc->GetX()+dx;
    double yUpd = trc->GetY()+dx*dYdX;
    double F =  xUpd*xUpd + yUpd*yUpd - r*r;
    double dFdx = 2.*(xUpd + yUpd*(dYdX+dx*ddYdX*crv));
    double eps = -F/dFdx;
    dx += eps;
    if (TMath::Abs(eps)<kMinEps) break;
  }
  return trc->GetX()+dx;
}
