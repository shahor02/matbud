
MatLayerCylSet st;

void tst(float rmin=0,float rmax=50,int nr=50, float zmin=-50,float zmax=50,int nz=50,int nphi=30)
{
  o2::Base::GeometryManager::loadGeometry("~/tmp/match/pbpb_hijing/O2geometry.root");
  float dr = (rmax-rmin)/nr;
  float r=rmin;
  for (int ir=0;ir<rmax;ir++) {
    st.getLayers().emplace_back(r,r+dr, zmin,zmax, nz, nphi);
    r += dr;
  }
  st.populateFromTGeo(5);
}
