/*
 * SimpleAnalysis.C
 *
 *  Created on: 23 feb 2017
 *      Author: Nicola Mori
 */

/* Analysis ROOT script for the output file created by SimpleRun.mac. */

/* ********** INSTRUCTIONS ***********
 * Before loading this script in the Root shell you must load the dictionaries
 * for the GGS analysis classes, e.g.:
 *
 *   gROOT->ProcessLine(".include /path/to/GGSinstallation/lib);
 *   gSystem->Load("libGGSReader.so");
 *   gSystem->Load("libGGSDataObjects.so");
 *
 */

// GGS headers
#include "utils/GGSSmartLog.h"
/*
double
CosAngle(const TVector3& l, const TVector3& r)
{
  const double magnitudeA = l.GetMag();
  const double magnitudeB = r.GetMag();
    return (magnitudeA && magnitudeB) ?
      (l * r) / (magnitudeA * magnitudeB) : 1;
  }


double
Angle(const TVector3& left, const TVector3& right)
{
  //const double cosAngle = CosAngle(left, right);
  //  return (cosAngle >= 1) ? 0 : acos(cosAngle);
  // DV: this is more accurate for small and large angles
  const Vector a = 1/left.GetMag() * left;
    const Vector b = 1/right.GetMag() * right;
    const double d2 = (a - b).GetMag2();
    if (d2 <= 2)
      return 2 * std::asin(0.5 * std::sqrt(d2));
    else
      return 2 * std::acos(0.5 * (a + b).GetMag());
}

*/

// Event vars-----------------------
int nInt = 0;
int iEv = 0;
int nHits = 0;
int ppHit = -1;
double eLastZ = -1.;
double pLastZ = -1.;
  
const Int_t nMaxHits = 100;

//  Double_t zPath[nMaxHits];
Double_t gEne=-9999; 
Double_t zPath=-9999;
//  Double_t xCoord[nMaxHits];
//  Double_t yCoord[nMaxHits];
//  Double_t zCoord[nMaxHits];
//  Double_t xMom[nMaxHits];
//Double_t yMom[nMaxHits];
//Double_t zMom[nMaxHits];
Double_t eeDep[nMaxHits]={0.};
Double_t peDep[nMaxHits]={0.};
int chX[nMaxHits]={0};
int chY[nMaxHits]={0}; 

Int_t hPart[nMaxHits]={0};
Int_t hVol[nMaxHits]={0};
Double_t hVolZ[nMaxHits]={0.};
    
Double_t exCoord[nMaxHits]={0.};
Double_t eyCoord[nMaxHits]={0.};
Double_t ezCoord[nMaxHits]={0.};
Double_t eexxCoord=-9999.;
Double_t eexyCoord=-9999.;
Double_t eexzCoord=-9999.;
//  Double_t eexxCoord[nMaxHits];
//  Double_t eexyCoord[nMaxHits];
//  Double_t eexzCoord[nMaxHits];  
Double_t exMom[nMaxHits]={0.};
Double_t eyMom[nMaxHits]={0.};
Double_t ezMom[nMaxHits]={0.};
//Double_t eEne[nMaxHits];
Double_t eEne=-9999.;
int echX[nMaxHits]={0};
int echY[nMaxHits]={0}; 

Double_t pxCoord[nMaxHits]={0.};
Double_t pyCoord[nMaxHits]={0.};
Double_t pzCoord[nMaxHits]={0.};
//  Double_t pexxCoord[nMaxHits];
//  Double_t pexyCoord[nMaxHits];
//  Double_t pexzCoord[nMaxHits];
Double_t pexxCoord=-9999.;
Double_t pexyCoord=-9999.;
Double_t pexzCoord=-9999.;
Double_t pxMom[nMaxHits]={0.};
Double_t pyMom[nMaxHits]={0.};
Double_t pzMom[nMaxHits]={0.};
Double_t pEne=-9999.;
int pchX[nMaxHits]={0};
int pchY[nMaxHits]={0};
//----------------------------------

void CleanEvent(int nMaxHits){

  nInt = 0;
  iEv = 0;
  nHits = 0;
  ppHit = -1;
  eLastZ = -1.;
  pLastZ = -1.;

  //  Double_t zPath[nMaxHits];
  gEne=-9999; 
  zPath=-9999;
  //  Double_t xCoord[nMaxHits];
  //  Double_t yCoord[nMaxHits];
  //  Double_t zCoord[nMaxHits];
  //  Double_t xMom[nMaxHits];
  //Double_t yMom[nMaxHits];
  //Double_t zMom[nMaxHits];
  std::fill_n(eeDep, nMaxHits, 0.);
  std::fill_n(peDep, nMaxHits, 0.);
  std::fill_n(chX, nMaxHits, 0);
  std::fill_n(chY, nMaxHits, 0);
  
  std::fill_n(hPart, nMaxHits, 0);
  std::fill_n(hVol, nMaxHits, 0);
  std::fill_n(hVolZ, nMaxHits, 0.);
    
  std::fill_n(exCoord, nMaxHits, 0.);
  std::fill_n(eyCoord, nMaxHits, 0.);
  std::fill_n(ezCoord, nMaxHits, 0.);
  eexxCoord=-9999.;
  eexyCoord=-9999.;
  eexzCoord=-9999.;
  //  Double_t eexxCoord[nMaxHits];
  //  Double_t eexyCoord[nMaxHits];
  //  Double_t eexzCoord[nMaxHits];  
  std::fill_n(exMom, nMaxHits, 0.);
  std::fill_n(eyMom, nMaxHits, 0.);
  std::fill_n(ezMom, nMaxHits, 0.);
  //Double_t eEne[nMaxHits];
  eEne=-9999.;
  std::fill_n(echX, nMaxHits, 0);
  std::fill_n(echY, nMaxHits, 0);

  std::fill_n(pxCoord, nMaxHits, 0.);
  std::fill_n(pyCoord, nMaxHits, 0.);
  std::fill_n(pzCoord, nMaxHits, 0.);
  //  Double_t pexxCoord[nMaxHits];
  //  Double_t pexyCoord[nMaxHits];
  //  Double_t pexzCoord[nMaxHits];
  pexxCoord=-9999.;
  pexyCoord=-9999.;
  pexzCoord=-9999.;
  std::fill_n(pxMom, nMaxHits, 0.);
  std::fill_n(pyMom, nMaxHits, 0.);
  std::fill_n(pzMom, nMaxHits, 0.);
  pEne=-9999.;
  std::fill_n(pchX, nMaxHits, 0);
  std::fill_n(pchY, nMaxHits, 0);
  
}

void SimpleAnalysis(TString inputFileName, TString outputFileName) {
  static const std::string routineName("simpleanalysis");

  GGSSmartLog::verboseLevel = GGSSmartLog::INFO; // Print only INFO messages or more important

  COUT(INFO) << "Begin analysis" << ENDL;

  // Create the reader container and open the data file
  GGSTRootReader reader;
  if (!(reader.Open(inputFileName))) {
    COUT(ERROR) << "Cannot open input file " << inputFileName << ENDL;
    return 1;
  }

  // Print some information about the simulation and the geometry
  const GGSTSimInfo *simInfo = reader.GetSimInfo();
  if (simInfo) {
    std::cout << "*** Simulation informations:\n";
    std::cout << "Geant4 version     : " << simInfo->G4Version << "\n";
    std::cout << "GGS version        : " << simInfo->GGSVersion << "\n";
    //std::cout << "Geometry version   : " << simInfo->geometryVersion << "\n";
  }
  else {
    std::cout << "*** No simulation information are available\n";
  }  
  const GGSTGeoParams *geoParams = reader.GetGeoParams();
  if (geoParams) {
    std::cout << "*** Geometry parameters:\n";
    std::cout << "Target Layers Number:    " << geoParams->GetIntGeoParam("targetLayerNo") << "\n";
    std::cout << "\"New\" Sensor X Dimension:    " << geoParams->GetRealGeoParam("tileX") << " cm\n";
    std::cout << "\"New\" Sensor Y Dimension:    " << geoParams->GetRealGeoParam("tileY") << " cm\n";  
    std::cout << "\"New\" Sensor Thickness:    " << geoParams->GetRealGeoParam("tileThickness") << " mm\n";
    std::cout << "\"New\" Sensor S Pitch:    " << geoParams->GetRealGeoParam("tileSPitch") << " mm\n";
    std::cout << "\"New\" Sensor K Pitch:    " << geoParams->GetRealGeoParam("tileKPitch") << " mm\n";
    std::cout << "AMS Sensor X Dimension:    " << geoParams->GetRealGeoParam("AMStileX") << " cm\n";
    std::cout << "AMS Sensor Y Dimension:    " << geoParams->GetRealGeoParam("AMStileY") << " cm\n";  
    std::cout << "AMS Sensor Thickness:    " << geoParams->GetRealGeoParam("AMStileThickness") << " mm\n";
    std::cout << "AMS Sensor S Pitch:    " << geoParams->GetRealGeoParam("AMStileSPitch") << " mm\n";
    std::cout << "AMS Sensor K Pitch:    " << geoParams->GetRealGeoParam("AMStileKPitch") << " mm\n";
    std::cout << "DAMPE Sensor X Dimension:    " << geoParams->GetRealGeoParam("DAMPEtileX") << " cm\n";
    std::cout << "DAMPE Sensor Y Dimension:    " << geoParams->GetRealGeoParam("DAMPEtileY") << " cm\n";  
    std::cout << "DAMPE Sensor Thickness:    " << geoParams->GetRealGeoParam("DAMPEtileThickness") << " mm\n";
    std::cout << "DAMPE Sensor S Pitch:    " << geoParams->GetRealGeoParam("DAMPEtileSPitch") << " mm\n";
    std::cout << "DAMPE Sensor K Pitch:    " << geoParams->GetRealGeoParam("DAMPEtileKPitch") << " mm\n";
    std::cout << "Spectrometer Layers Number:    " << geoParams->GetIntGeoParam("smLayerNo") << "\n";
    std::cout << "Mag Field Value:    " << geoParams->GetRealGeoParam("magFieldVal") << " Tesla\n";
    std::cout << "Mag Field Z:    " << geoParams->GetRealGeoParam("magVolZ") << " cm\n";
  }
  std::cout << std::endl;

  bool DB=false;  
  // GEOMETRYPARS TO GET FROM READER

  // double mside = 7.04; // cm sensor measurement side length
  //double pitch = 0.0110; // cm sensor pitch

  //MD: FIX ME
  double msidex = geoParams->GetRealGeoParam("tileX"); // cm sensor measurement side length
  double msidey = geoParams->GetRealGeoParam("tileY"); // cm sensor measurement side length
  double pitch = geoParams->GetRealGeoParam("tileSPitch")/10.; // cm sensor pitch

  // double offsetx = geoParams->GetRealGeoParam("layersOffsetX"); // cm offset X
  // double offsety = geoParams->GetRealGeoParam("layersOffsetY"); // cm offset Y

  //----------------------------------------------------------------------------------------------------------------------------------
  // Create the output file
  TFile *outFile = TFile::Open(outputFileName, "RECREATE");
  if (!outFile || outFile->IsZombie()) {
    COUT(ERROR) << "Cannot create output file " << outputFileName << ENDL;
    return 1;
  }

  // Create and retrieve the hits sub-reader
  GGSTHitsReader *hReader = reader.GetReader<GGSTHitsReader>();
  // Set which hit detectors are to be read
  // hReader->SetDetector("siLayerMiniLog", kTRUE); // The name is the same of the sensitive logical volume name in the simulation
  // hReader->SetDetector("siLayerShortLog", kTRUE); // The name is the same of the sensitive logical volume name in the simulation
  // hReader->SetDetector("siLayerLongLog", kTRUE); // The name is the same of the sensitive logical volume name in the simulation
  hReader->SetDetector("siAMSTileLog", kTRUE); // The name is the same of the sensitive logical volume name in the simulation
  hReader->SetDetector("siDAMPETileLog", kTRUE); // The name is the same of the sensitive logical volume name in the simulation
  hReader->SetDetector("siTileLog", kTRUE); // The name is the same of the sensitive logical volume name in the simulation

  //   hReader->SetDetector("siStripLog", kTRUE); // The name is the same of the sensitive logical volume name in the simulation
   
  // Retrieve the MC truth sub-reader
  GGSTMCTruthReader *mcReader = reader.GetReader<GGSTMCTruthReader>();

  // Create and retrieve the hadronic interaction sub-reader
  GGSTHadrIntReader *hadrReader = reader.GetReader<GGSTHadrIntReader>();


  // Prepare the histograms
  gStyle->SetOptFit(1);
  
  //  TH1F *nHitHisto = new TH1F("nHitHisto", "NHits;NHit; Counts", 30, 0, 30); // nhits histo

  //  TH1F *eDepHisto = new TH1F("eDepHisto", "Total energy deposit;E (MeV); Counts", 1000, 0, 10); // Total energy release

  //  TH1F *ePDepHisto = new TH1F("ePDepHisto", "particle energy deposit;E (MeV); Counts", 1000, 0, 10); // Particle energy release
    
  //  TH1F *dyHisto = new TH1F("dyHisto", "Y bending of particle [mm]", 2000, -100., 100.); // y bending

  //  TH1F *dbHisto = new TH1F("dbHisto", "delta bending [mm]", 2000, -100., 100.); // theo - geant4 sim

  //  TH1F *dthHisto = new TH1F("dthHisto", "delta theta_{rec-true} [deg]", 2000, 0., 50.); // theo - geant4 sim

  //     TH1F *dthVHisto = new TH1F("dthVHisto", "track space angle_{rec-true} [deg]", 2000, 0., 50.); // theo - geant4 sim
	
     //  TH2F *xyHisto = new TH2F("xyHisto", "YX coordinates of particle [mm]", 1000, -10, 10,1000, -10,10); // YX primary
  
     // TH1F *zIntHisto = new TH1F("zIntHisto", "Z coordinate of interaction point [cm]", 1000, 0, 100); // Interaction point
     //TH2F *xyIntHisto = new TH2F("xyIntHisto", "YX coordinates of interaction point [cm]", 1000, -10, 10,1000, -10,10); // YX Interaction point


  int nLayers = geoParams->GetIntGeoParam("targetLayerNo")+geoParams->GetIntGeoParam("smLayerNo");


  Int_t xyAlign[nLayers];
  string sAlign=geoParams->GetStringGeoParam("layerAlignment");
  for (int isn=0; isn<nLayers; isn++) {
    xyAlign[isn]=sAlign[isn]-'0';
    //    printf("xyAlign[%d]=%d\n", isn, xyAlign[isn]);
  }
  /*
  // manual alignment for older files
  Int_t xyAlign[14]={0,1,0,1,0,1,0,1,0,1,1,1,1,1};
  //Int_t xyAlign[24]={0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1,1,1,1};
  */
  int nEvts = 0;
  
  //------------------------------------------------------

  TTree *runTree =  new TTree("runTree","tree of mcrun");
  runTree->Branch("nLayers",&nLayers,"nLayers/I");
  runTree->Branch("xyAlign",&xyAlign,"xyAlign[nLayers]/I");
  runTree->Branch("nEvts",&nEvts,"nEvts/I");
  runTree->Branch("gEne",&gEne,"gEne/D");
    
  TTree *hitTree =  new TTree("hitTree","tree of layer hits");
  hitTree->Branch("evID",&iEv,"evID/I");
  hitTree->Branch("nHits",&nHits,"nHits/I");
  hitTree->Branch("ppHit",&ppHit,"ppHit/I");
  hitTree->Branch("eLastZ",&eLastZ,"eLastZ/D");
  hitTree->Branch("pLastZ",&pLastZ,"pLastZ/D");

  // primary gamma and pprod energies
  //  hitTree->Branch("zPath",&zPath,"zPath[nHits]/D");
  hitTree->Branch("zPath",&zPath,"zPath/D");  
  hitTree->Branch("eEne",&eEne,"eEne/D");
  hitTree->Branch("pEne",&pEne,"pEne/D");
  
  //hitTree->Branch("xCoord",&xCoord,"xCoord[nHits]/D");
  //hitTree->Branch("yCoord",&yCoord,"yCoord[nHits]/D");
  //hitTree->Branch("zCoord",&zCoord,"zCoord[nHits]/D");  
  //hitTree->Branch("xMom",&xMom,"xMom[nHits]/D");
  //hitTree->Branch("yMom",&yMom,"yMom[nHits]/D");
  //hitTree->Branch("zMom",&zMom,"zMom[nHits]/D");

  // hit
  hitTree->Branch("hPart",&hPart,"hPart[nHits]/I");
  hitTree->Branch("hVol",&hVol,"hVol[nHits]/I");  
  hitTree->Branch("hVolZ",&hVolZ,"hVolZ[nHits]/D");  
  
  // electron
  hitTree->Branch("exCoord",&exCoord,"exCoord[nHits]/D");
  hitTree->Branch("eyCoord",&eyCoord,"eyCoord[nHits]/D");
  hitTree->Branch("ezCoord",&ezCoord,"ezCoord[nHits]/D");  
  hitTree->Branch("eexxCoord",&eexxCoord,"eexxCoord/D");
  hitTree->Branch("eexyCoord",&eexyCoord,"eexyCoord/D");
  hitTree->Branch("eexzCoord",&eexzCoord,"eexzCoord/D");  
  hitTree->Branch("exMom",&exMom,"exMom[nHits]/D");
  hitTree->Branch("eyMom",&eyMom,"eyMom[nHits]/D");
  hitTree->Branch("ezMom",&ezMom,"ezMom[nHits]/D");
  //  hitTree->Branch("eEne",&eEne,"eEne[nHits]/D");
  hitTree->Branch("eeDep",&eeDep,"eeDep[nHits]/D");
  
  // positron
  hitTree->Branch("pxCoord",&pxCoord,"pxCoord[nHits]/D");
  hitTree->Branch("pyCoord",&pyCoord,"pyCoord[nHits]/D");
  hitTree->Branch("pzCoord",&pzCoord,"pzCoord[nHits]/D");
  hitTree->Branch("pexxCoord",&pexxCoord,"pexxCoord/D");
  hitTree->Branch("pexyCoord",&pexyCoord,"pexyCoord/D");
  hitTree->Branch("pexzCoord",&pexzCoord,"pexzCoord/D");  
  hitTree->Branch("pxMom",&pxMom,"pxMom[nHits]/D");
  hitTree->Branch("pyMom",&pyMom,"pyMom[nHits]/D");
  hitTree->Branch("pzMom",&pzMom,"pzMom[nHits]/D");
  //  hitTree->Branch("pEne",&pEne,"pEne[nHits]/D");
  hitTree->Branch("peDep",&peDep,"peDep[nHits]/D");
    
  hitTree->Branch("chX",&chX,"chX[nHits]/I");
  hitTree->Branch("chY",&chY,"chY[nHits]/I");
  hitTree->Branch("echX",&echX,"echX[nHits]/I");
  hitTree->Branch("echY",&echY,"echY[nHits]/I");
  hitTree->Branch("pchX",&pchX,"pchX[nHits]/I");
  hitTree->Branch("pchY",&pchY,"pchY[nHits]/I");

  double mom=0.;
  int eno=0;
  int pno=0; 

  nEvts=reader.GetEntries();
  
  COUT(INFO) << "Begin loop over " << nEvts << " events" << ENDL;
  
  for (iEv = 0; iEv < nEvts; iEv++) {
  //for (int iEv = 0; iEv < 33; iEv++) {

    //    CleanEvent(nMaxHits);
    
    reader.GetEntry(iEv); // Reads all the data objects whose sub-readers have already been created
    
    // Retrieve inelastic interaction information
    GGSTHadrIntInfo *intInfo = hadrReader->GetInelastic(); // Get info about the inelastic interaction of the primary particle
    // Check if the inelastic interaction actually happened for this event
    if (intInfo) {
      //      zIntHisto->Fill(intInfo->GetInteractionPoint()[2]);
      //xyIntHisto->Fill(intInfo->GetInteractionPoint()[0],intInfo->GetInteractionPoint()[1]);
      std::cout<<"EVT "<<iEv<<" x,y,z [cm]: "<<intInfo->GetInteractionPoint()[0]<<" "<<intInfo->GetInteractionPoint()[1]<<" "<<intInfo->GetInteractionPoint()[2]<<std::endl;
      nInt++;
    }

    // Compute total energy release
    GGSTIntHit* thisHit;
    GGSTPartHit* thisPHit;

    std::vector<GGSTIntHit*> vTIntHit;

    nHits=0;
    // int nMHits = hReader->GetNHits("siLayerMiniLog"); //Number of hit siLayers for current event
    // int nSHits = hReader->GetNHits("siLayerShortLog"); //Number of hit siLayers for current event
    // int nHits = hReader->GetNHits("siStripLog"); //Number of hit siLayers for current event
    {
      int _nHits = hReader->GetNHits("siTileLog"); //Number of hit siLayers for current event
      for (int iHit = 0; iHit < _nHits; iHit++) {
	thisHit = hReader->GetHit("siTileLog", iHit);
	vTIntHit.push_back(thisHit);
      }
      nHits += _nHits;
    }
    {
      int _nHits = hReader->GetNHits("siAMSTileLog"); //Number of hit siLayers for current event
      for (int iHit = 0; iHit < _nHits; iHit++) {
	thisHit = hReader->GetHit("siAMSTileLog", iHit);
	vTIntHit.push_back(thisHit);
      }
      nHits += _nHits;
    }
    {
      int _nHits = hReader->GetNHits("siDAMPETileLog"); //Number of hit siLayers for current event
      for (int iHit = 0; iHit < _nHits; iHit++) {
	thisHit = hReader->GetHit("siDAMPETileLog", iHit);
	vTIntHit.push_back(thisHit);
      }
      nHits += _nHits;
    }
    std::cout<<"EVT "<<iEv<<" NHITS "<<nHits<<std::endl;
    //nHitHisto->Fill(nHits);
    
    float totEDep = 0.;   
   
    // efficiency calc
    bool pprod=false; 

    // true Y bending
    //double mom=1.;  // mom in GV (not geant4 ene
    if (hReader->GetHit("siTileLog",0)) {
      mom=hReader->GetHit("siTileLog",0)->GetPartHit(0)->entranceMomentum[2];
    }
    else {
      mom=0.0;
    }
    
    //double magf=1.; // mag field [tesla]
    double magf=geoParams->GetRealGeoParam("magFieldVal");
    double magz=geoParams->GetRealGeoParam("magVolZ")/100.; // mag z size [m]
    double ro=0.;
    double dYtrue=0.;
    double thtrue=0.;
    double thcenter=0.;
    double thtrueV=0.;
    
    if(magf){
      ro=mom/magf/0.3;  // radius of curvature [m]
      dYtrue = ro*(1-sqrt(1-(magz/ro)*(magz/ro)));
      thtrue = atan2(dYtrue,magz);
      thcenter = atan2(magz,(ro-dYtrue));
    }
 
    TVector3 vzeta(0.,0.,mom);
    //    TVector3 vtrue(0.,dYtrue,mom);
    TVector3 vtrue(0.,mom*sin(thcenter),mom*cos(thcenter));
    TVector3 vrec;

    thtrueV= vtrue.Angle(vzeta);

    if(DB)
      std::cout<<"DY true: "<<mom<<" "<<ro<<" "<<dYtrue*100<<" TH "<<thtrue*180/TMath::Pi()<<" THC "<<thcenter*180/TMath::Pi()<<" THV "<<thtrueV*180/TMath::Pi()<<std::endl;
    
    // Hits loop
    double yfront=0.;
    double yend1=0.;
    double yend2=0.;
    double threc=0.;

    for (int iHit = 0; iHit < nHits; iHit++) {
      hPart[iHit]=0;	
      thisHit = vTIntHit.at(iHit);
      int nPHits= thisHit->GetNPartHits();
      
      if(nPHits>=2||pprod)
	
	std::cout<<"******* HIT: "<< iHit <<" --->>>PARTHITS: "<<nPHits<<" detind: "<<hReader->GetDetectorIndex("siTileLog")<<" "<<thisHit->GetVolumeName()<<" "<<thisHit->GetVolumeID()<<" pos:"<<thisHit->GetVolumePosition()[2]<<std::endl;

      hVol[iHit]=thisHit->GetVolumeID();
      hVolZ[iHit]=thisHit->GetVolumePosition()[2];
      
      
      if (nPHits>=3)
	std::cout<<"*********************  EVT: "<<iEv<<" HIT: "<< iHit <<" --->>>3INT"<< std::endl;

      double zp=-9999;
      double ex=-9999;
      double ey=-9999;
      double ez=-9999;

      for (int iPHit =0; iPHit<nPHits; iPHit++){
	thisPHit  = thisHit->GetPartHit(iPHit);

	if (nPHits>=2||pprod)
	  thisPHit->DumpHit();
	
	if(thisPHit->parentID==0&&thisPHit->particlePdg==22){ /// primary particle hit values 
	  //zPath[iHit]=thisPHit->pathLength;
	  //xCoord[iHit]=thisPHit->entrancePoint[0];
	  //yCoord[iHit]=thisPHit->entrancePoint[1];
	  //zCoord[iHit]=thisPHit->entrancePoint[2];
	  //xMom[iHit]=thisPHit->entranceMomentum[0];
	  //yMom[iHit]=thisPHit->entranceMomentum[1];
	  //zMom[iHit]=thisPHit->entranceMomentum[2];
	  chX[iHit]=int((thisPHit->entrancePoint[0]+msidex/2.)/pitch);
	  chY[iHit]=int((thisPHit->entrancePoint[1]+msidey/2.)/pitch);
	  //	  ePDepHisto->Fill(1e3*thisPHit->eDep);
	  zp=thisPHit->pathLength;
	}  //// primary particle values


	if(thisPHit->parentID==1&&thisPHit->particlePdg==11){ /// electron
	  //eEne[iHit]=thisPHit->entranceEnergy;
	  exCoord[iHit]=thisPHit->entrancePoint[0];
	  eyCoord[iHit]=thisPHit->entrancePoint[1];
	  ezCoord[iHit]=thisPHit->entrancePoint[2];
	  //eexxCoord[iHit]=thisPHit->exitPoint[0];
	  //eexyCoord[iHit]=thisPHit->exitPoint[1];
	  //eexzCoord[iHit]=thisPHit->exitPoint[2];
	  ex=thisPHit->exitPoint[0];
	  ey=thisPHit->exitPoint[1];
	  ez=thisPHit->exitPoint[2];
	  eEne=thisPHit->entranceEnergy;
	  exMom[iHit]=thisPHit->entranceMomentum[0];
	  eyMom[iHit]=thisPHit->entranceMomentum[1];
	  ezMom[iHit]=thisPHit->entranceMomentum[2];
	  eeDep[iHit]=1e3*thisPHit->eDep;
	  echX[iHit]=int((exCoord[iHit]+msidex/2.)/pitch);
	  echY[iHit]=int((eyCoord[iHit]+msidey/2.)/pitch);
	  eno++;
	  hPart[iHit] +=1;
	}  //// electron values


	if(thisPHit->parentID==1&&thisPHit->particlePdg==-11){ /// positron
	  //pEne[iHit]=thisPHit->entranceEnergy;
	  pxCoord[iHit]=thisPHit->entrancePoint[0];
	  pyCoord[iHit]=thisPHit->entrancePoint[1];
	  pzCoord[iHit]=thisPHit->entrancePoint[2];
	  //pexxCoord[iHit]=thisPHit->exitPoint[0];
	  //pexyCoord[iHit]=thisPHit->exitPoint[1];
	  //pexzCoord[iHit]=thisPHit->exitPoint[2];
	  pxMom[iHit]=thisPHit->entranceMomentum[0];
	  pyMom[iHit]=thisPHit->entranceMomentum[1];
	  pzMom[iHit]=thisPHit->entranceMomentum[2];
	  peDep[iHit]=1e3*thisPHit->eDep;
	  pchX[iHit]=int((pxCoord[iHit]+msidex/2.)/pitch);
	  pchY[iHit]=int((pyCoord[iHit]+msidey/2.)/pitch);
	  pno++;
	  hPart[iHit] +=2;
	  if(DB)
	    std::cout<<"*********************  EVT: "<<iEv<<" HIT: "<< iHit <<" --->>> PPROD POSITRON PHIT: "<< iPHit<<" eno "<<eno<<" pno "<<pno<< std::endl;	  

	  if (!pprod){
	  pprod=true;
	  ppHit=iHit;
	  pexxCoord=thisPHit->exitPoint[0];
	  pexyCoord=thisPHit->exitPoint[1];
	  pexzCoord=thisPHit->exitPoint[2];
	  pEne=thisPHit->entranceEnergy;

	  gEne=mom;
	  zPath=zp;
	  eexxCoord=ex;
	  eexyCoord=ey;
	  eexzCoord=ez;
	  }	  
	} //// positron values

      } // loop on particle hits

      
      //eDep[iHit]=1e3*thisHit->eDep;
      totEDep += thisHit->eDep;
      if(DB)
	std::cout<<"*********************  EVT: "<<iEv<<" HIT: "<< iHit <<" nPART: "<< hPart[iHit]<< std::endl;
      
    } // loop on hits

    if(ppHit>=0){
      // layer ch part multiplicity
      Int_t lpmult=0;
      int fhz=nHits;
      eLastZ=hVolZ[ppHit];
      pLastZ=hVolZ[ppHit];
      double tmplz=hVolZ[ppHit];
      for (int ihh=ppHit+1;ihh<nHits;ihh++){

	lpmult=hPart[ihh];    	  
	if(lpmult==0)
	  continue;
	
	if (lpmult==3){
	  eLastZ=hVolZ[ihh];
	  pLastZ=hVolZ[ihh];
	  if(DB)
	    std::cout<<"EVT: "<<iEv<<" HIT: "<< ihh <<" VZ: "<<hVolZ[ihh]<<" HPART: "<<hPart[ihh]<<" LPMULT: "<<lpmult<<" EPLAST: "<< eLastZ<<" "<<pLastZ<< std::endl;
	  continue;
	}else{
	  
	  if(hPart[ihh]==1)
	    eLastZ=hVolZ[ihh];
	  if(hPart[ihh]==2)
	    pLastZ=hVolZ[ihh];
	  if(DB)
	    std::cout<<"***EVT: "<<iEv<<" HIT: "<< ihh <<" VZ: "<<hVolZ[ihh]<<" HPART: "<<hPart[ihh]<<" LPMULT: "<<lpmult<<" EPLAST: "<< eLastZ<<" "<<pLastZ<<std::endl;
	  
	  double hvz=hVolZ[ihh];
	  
	  for (int jhh=ihh+1;jhh<nHits;jhh++){
	    
	    if (hVolZ[jhh]==hvz&&hPart[jhh]!=0){
	      lpmult+=hPart[jhh];
	      if(hPart[jhh]==1)
		eLastZ=hVolZ[jhh];
	      if(hPart[jhh]==2)
		pLastZ=hVolZ[jhh];	    
	      fhz=jhh;
	      if(DB)
		std::cout<<"---> EVT: "<<iEv<<" HIT: "<< fhz <<" VZ: "<<hVolZ[jhh]<<" HPART: "<<hPart[jhh]<<" LPMULT: "<<lpmult<<" RECEPLASTZ: "<< eLastZ<<" "<<pLastZ<< std::endl;
	    }else{
	      continue;
	    }
	  }

	}
      }
    }//IF PPROD
    
    if(!intInfo&&pprod)
      {

	hitTree->Fill();
	
	//	  std::cout<<"LAST: "<<yend1<<" -> substitution: "<<threcV*180/TMath::Pi()<<" "<<lastZ<<" ="<<tan(threcV)*(60.-lastZ)<<std::endl;
	//        yend1=tan(threcV)*(60.-lastZ);

	// dbHisto->Fill(10*(yend1-yfront)-1e3*dYtrue);

	//	dthHisto->Fill((threc-thtrue)*180/TMath::Pi());
	//dthVHisto->Fill((vrec.Angle(vtrue))*180/TMath::Pi());
	//eDepHisto->Fill(1e3*totEDep);
      //      xyHisto->Fill(x,y);
      //      std::cout<<"EVT "<<iEv<<" EDep [MeV]: "<<1e3*totEDep<<std::endl;
      }
    
    //    GGSTMCTruthInfo *mcTruthInfo = mcReader-> ??; 
    // MC reader
    
  } // EV loop
  
  COUT(INFO) << "Event loop finished" << ENDL;
  
  runTree->Fill();
  
  // Save histograms
  outFile->cd();
  
  //  nHitHisto->Write();
  //dyHisto->Write();
  //dbHisto->Write();
  /*  
  TF1 * f1 = new TF1("f1","gaus",-8.5,8.5);
  dbHisto->Fit("f1","R");
  Double_t p0 = f1->GetParameter(0);
  Double_t p1 = f1->GetParameter(1);
  Double_t p2 = f1->GetParameter(2);
  //  delete f1;
  std::cout<<"DBHISTOPARS "<<mom<<" "<<dbHisto->GetMean()<<" "<<dbHisto->GetRMS()<<" "<<p1<<" "<<p2<<std::endl;



  TCanvas * cth=new TCanvas("cth","cth",600,400);
  Int_t nq = 100;
  Double_t xq[nq];  // position where to compute the quantiles in [0,1]
  Double_t yq[nq];  // array to contain the quantiles
  for (Int_t i=0;i<nq;i++) xq[i] = Float_t(i+1)/nq;
  dthHisto->GetQuantiles(nq,yq,xq);
  TLine * qline = new TLine(yq[68],0.,yq[68],dthHisto->GetBinContent(dthHisto->FindBin(yq[68])));
  qline->SetLineColor(2);
  qline->SetLineWidth(2);
  qline->SetLineStyle(2);
  dthHisto->Draw();
  qline->Draw("same");
  std::cout<<"DTH res "<<dthHisto->GetMean()<<" "<<dthHisto->GetRMS()<<" 68%: "<<yq[68]<<std::endl;
  dthHisto->Write();


  TCanvas * cthV=new TCanvas("cthV","cthV",600,400);
  Double_t xqV[nq];  // position where to compute the quantiles in [0,1]
  Double_t yqV[nq];  // array to contain the quantiles
  for (Int_t i=0;i<nq;i++) xqV[i] = Float_t(i+1)/nq;
  dthVHisto->GetQuantiles(nq,yqV,xqV);
  TLine * qlineV = new TLine(yqV[68],0.,yqV[68],dthVHisto->GetBinContent(dthVHisto->FindBin(yqV[68])));
  qlineV->SetLineColor(2);
  qlineV->SetLineWidth(2);
  qlineV->SetLineStyle(2);
  dthVHisto->Draw();
  qlineV->Draw("same");
  std::cout<<"DTHV res "<<dthVHisto->GetMean()<<" "<<dthVHisto->GetRMS()<<" 68%: "<<yqV[68]<<std::endl;
  dthVHisto->Write();

  eDepHisto->Write();
  ePDepHisto->Write();
  zIntHisto->Write();
  xyIntHisto->Write();
  */
  runTree->Write();  
  hitTree->Write();  
  outFile->Close();
  delete outFile;
  COUT(INFO) << "Analysis finished" << ENDL;

}
