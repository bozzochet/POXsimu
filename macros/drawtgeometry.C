void drawtgeometry(){

  gSystem->Load("libGeom");

  TGeoManager::Import("plugins/libTestGeometry.vgm.root");
  
  TBrowser *b=new TBrowser();
  std::cout<<"///////////////////////////////////////////\n";
  TGeoVolume *magnet = gGeoManager->GetVolume("magnet");
  TGeoUniformMagField *magField = new TGeoUniformMagField(0.,0.,0.05);
  magnet->SetField(magField);
  TObjArray* oav = gGeoManager->GetListOfVolumes();
  for (int ii=0; ii<oav->GetEntries(); ii++) {
    oav->At(ii)->Dump();
  }
  TGeoVolume* world = gGeoManager->GetVolume("World");

  world->Draw();
  
}
