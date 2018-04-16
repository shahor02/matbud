
MatLayerCylSet st;

void tst(float rmin=0,float rmax=50,int nr=50, float zmax=50, float dz=1.f, float drphi=0.5f)
{
  gSystem->Load("libCommonUtils.so");
  gSystem->Load("libDetectorsBase.so");
  gSystem->Load("libMatBud.so");

  o2::Base::GeometryManager::loadGeometry("~/tmp/match/pbpb_hijing/O2geometry.root");
  float dr = (rmax-rmin)/nr;
  float r=rmin;
  for (int ir=0;ir<nr;ir++) {
    st.addLayer(r,r+dr, zmax, dz, drphi);
    r += dr;
  }
  st.populateFromTGeo(5);

  TFile stf("mb.root","recreate");
  stf.WriteObjectAny(&st, st.Class(), "matBud");
  stf.Close();
}
