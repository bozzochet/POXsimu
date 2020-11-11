void drawtgeometry(){

  gSystem->Load("libGeom");

  TGeoManager::Import("plugins/libTestGeometry.vgm.root");

  // TObjArray* oav = gGeoManager->GetListOfVolumes();
  // for (int ii=0; ii<oav->GetEntries(); ii++) {
  //   oav->At(ii)->Dump();
  // }
  
  TGeoVolume* world = gGeoManager->GetVolume("World");

  world->Draw();
  
}
