#include "MatLayerCylSet.h"
#include "MatLayerCyl.h"
#include "DetectorsBase/GeometryManager.h"
#include <vector>
#include <TFile.h>
#include <TSystem.h>


MatLayerCylSet mbLUT;

struct LrData {
  float rMin = 0.f;
  float rMax = 0.f;
  float zHalf = 0.f;
  float dZMin = 999.f; // min Z bin
  float dRPhiMin = 999.f; // min r*phi bin
  
  LrData(float rMn=0.f,float rMx=0.f, float zHlf=0.f, float dzMn=9999.f, float drphMn=9999. ) :
    rMin(rMn), rMax(rMx), zHalf(zHlf), dZMin(dzMn), dRPhiMin(drphMn) {}
};

std::vector<LrData> lrData;
void configLayers();




void buildMat(std::string outName = "barrel", int nTst=10)
{
  gSystem->Exec("o2sim -n 0 -m PIPE ITS TPC TRD TOF EMCAL PHOS FIT > genGeom.log");
  
  //  o2::Base::GeometryManager::loadGeometry("~/tmp/match/pbpb_hijing/O2geometry.root");
  o2::Base::GeometryManager::loadGeometry("./O2geometry.root");

  configLayers();

  for (auto l : lrData) {
    if (l.rMin<75) continue;
    mbLUT.addLayer(l.rMin, l.rMax, l.zHalf, l.dZMin, l.dRPhiMin);
  }

  mbLUT.populateFromTGeo(nTst);

  TFile stf( (outName + "_raw.root").data(), "recreate");
  stf.WriteObjectAny(&mbLUT, mbLUT.Class(), "MatBud");
  stf.Close();

  mbLUT.dumpToTree( (outName + "_tree.root").data() );
  
  mbLUT.optimizePhiSlices();
  TFile stfo( (outName + "_opt.root").data(), "recreate");
  stfo.WriteObjectAny(&mbLUT, mbLUT.Class(), "MatBud");
  stfo.Close();
  
}

//_______________________________________________________________________
void configLayers()
{
  const float kToler = 1e-3;
  float drStep = 0.f, zSpanH = 0.f, zBin=0.f, rphiBin=0.f, phiBin=0.f;
  
  //                           rMin    rMax   zHalf
  lrData.emplace_back( LrData(  0.0f,  1.8f,  30.f ) );

  // beam pipe
  lrData.emplace_back( LrData(  lrData.back().rMax,  1.9f,  30.f ) );

  // ITS Inner Barrel
  drStep = 0.1;
  zSpanH = 17.;
  rphiBin = 0.2;// 0.1;
  zBin = 0.5;
  do {
    lrData.emplace_back( LrData(  lrData.back().rMax, lrData.back().rMax+drStep, zSpanH,  zBin, rphiBin) );
  } while (lrData.back().rMax<5.2+kToler);

  // air space between Inner and Middle Barrels
  lrData.emplace_back( LrData(  lrData.back().rMax, 18.0, zSpanH ) );

  // ITS Middle Barrel
  drStep = 0.2;
  zSpanH = 50.;
  rphiBin = 0.5;
  zBin = 0.5;
  do {
    lrData.emplace_back( LrData(  lrData.back().rMax, lrData.back().rMax+drStep, zSpanH,  zBin, rphiBin) );
  } while (lrData.back().rMax<29.+kToler);
  
  // air space between Middle and Outer Barrels
  zSpanH = 80.f;
  lrData.emplace_back( LrData(  lrData.back().rMax, 33.5, zSpanH ) );

  // ITS Outer barrel
  drStep = 0.5;
  zSpanH = 80.;
  rphiBin = 1.;
  zBin = 1.;
  do {
    lrData.emplace_back( LrData(  lrData.back().rMax, lrData.back().rMax+drStep, zSpanH,  zBin, rphiBin) );
  } while (lrData.back().rMax<45.5+kToler);

  // air space between Outer Barrel and shell
  zSpanH = 100.f;
  lrData.emplace_back( LrData(  lrData.back().rMax, 59.5, zSpanH ) );

  // Shell
  drStep = 0.5;
  zSpanH = 100.;
  rphiBin = 1.;
  zBin = 1.;
  do {
    lrData.emplace_back( LrData(  lrData.back().rMax, lrData.back().rMax+drStep, zSpanH,  zBin, rphiBin) );
  } while (lrData.back().rMax<63.+kToler);
  
  // air space between Shell and TPC
  zSpanH = 250.f;
  lrData.emplace_back( LrData(  lrData.back().rMax, 76, zSpanH ) );

  // TPC inner vessel
  // up to r = 78.5
  zSpanH = 250.f;
  rphiBin = 1.;
  zBin = 25.;
  lrData.emplace_back( LrData(  lrData.back().rMax, 78.5, zSpanH,  zBin, rphiBin) );

  //
  zSpanH = 250.f;
  rphiBin = 2.;
  zBin = 2;
  lrData.emplace_back( LrData(  lrData.back().rMax, 84.5, zSpanH,  zBin, rphiBin) );
 
  // TPC drum
  zSpanH = 250.f;
  lrData.emplace_back( LrData(  lrData.back().rMax, 250.0, zSpanH ) ); 

  // TPC outer vessel
  zSpanH = 247.f; // ignore large lumps of material at |z|>247
  rphiBin = 2.;
  zBin = 3.;
  lrData.emplace_back( LrData(  lrData.back().rMax, 258., zSpanH,  zBin, rphiBin) );

  zSpanH = 247.f; // ignore large lumps of material at |z|>247
  rphiBin = 2.;
  zBin = 999.; // no segmentation in Z
  lrData.emplace_back( LrData(  lrData.back().rMax, 280., zSpanH,  zBin, rphiBin) );

  drStep = 1;
  zBin = 5;
  rphiBin = 5.;
  do {
    zSpanH = lrData.back().rMax;
    lrData.emplace_back( LrData(  lrData.back().rMax, lrData.back().rMax+drStep, zSpanH,  zBin, rphiBin) );
  } while (lrData.back().rMax<400);  
}
