{
  TString inclpath = gSystem->Getenv("GGS_SYS");
  inclpath.Append("/lib");
  inclpath.Prepend(".include ");
  gROOT->ProcessLine(inclpath.Data());
  gSystem->Load("libGGSReader");
  gSystem->Load("libGGSDataObjects");
}
