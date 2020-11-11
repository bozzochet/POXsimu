void drawtgeometry(){

  gSystem->Load("libGeom");

  TGeoManager::Import("plugins/libTestGeometry.vgm.root");

  // TObjArray* oav = gGeoManager->GetListOfVolumes();
  // for (int ii=0; ii<oav->GetEntries(); ii++) {
  //   oav->At(ii)->Dump();
  // }

  TCanvas* c1 = new TCanvas();
  TGeoVolume* world = gGeoManager->GetVolume("World");
  world->Draw();

  TCanvas* c2 = new TCanvas();
  TGeoVolume* magnet = gGeoManager->GetVolume("magnet");
  magnet->Draw();
  printf("%p\n", magnet->GetField());
  
}
