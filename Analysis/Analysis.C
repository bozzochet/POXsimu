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
//#include "/home/vivi/POXsoft/NEW/GGSSoftware-install/include/utils/GGSSmartLog.h"
//#include "/Users/claudio/App/GGS/install/include/utils/GGSSmartLog.h"
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

#define nMaxHits 100
#define nMaxTotalHits 500
//#define nStrips 4096 // for now hardcoded but let's see

//int pri=0;
int pri=1;

// Event vars-----------------------
int nInt = 0;
int iEv = 0;
int nHits = 0;//number of TIntHits (essentially the number of firing volumes)
int nHitsData = 0;//
int nTotalHits = 0;//total number of TPathHits (~ nHits*<nPartHits>)
int ppHit = -1;
int nroSt = 0;
int nHS=0;
// double eLastZ = -1.;
// double pLastZ = -1.;

// Int_t hPart[nMaxHits]={0};


Int_t hVol[nMaxTotalHits]={0};
Double_t hVolZ[nMaxTotalHits]={0.};

// //  Double_t zPath[nMaxTotalHits];
Double_t gEne=-9999;
Double_t gX=-9999;
Double_t gY=-9999; 
Double_t zPath=-9999;

Double_t xCoord[nMaxTotalHits]={0.};
Double_t yCoord[nMaxTotalHits]={0.};
Double_t zCoord[nMaxTotalHits]={0.};
Double_t eDep[nMaxTotalHits]={0.};
Int_t PDG[nMaxTotalHits]={0};
Int_t TrID[nMaxTotalHits]={-1};
Int_t ParID[nMaxTotalHits]={-1};
Double_t xMom[nMaxTotalHits]={-1};
Double_t yMom[nMaxTotalHits]={-1};
Double_t zMom[nMaxTotalHits]={-1};
Double_t eEne[nMaxTotalHits]={-1};

Int_t chXY[nMaxTotalHits]={0};
//Int_t chY[nMaxTotalHits]={0}; 

Int_t hitStrips[nMaxTotalHits]={0};
Int_t simStrips[nMaxTotalHits]={0};

//

//std::pair<int,double> simData[nMaxHits*nStrips];
//Double_t simData[nMaxHits][nStrips];
/*
size_t N = 20;
size_t M = 20;
std::vector<int> normal;
normal.resize(N * M);

for (size_t i = 0; i < N; ++i)
    for (size_t j = 0; j < M; ++j)
        normal[i + j * N] = j;
*/

// nlayers*ntiles * nstrips

// hardcoded charge fractions approx capacitive coupling - width 3 strips
int cpw=3;
double chcp[7]={.04*.04*.04,.04*.04,.04,0.9167,.04,.04*.04,.04*.04*.04};

// vector of hit channels
std::vector<int> hitChan;  
// vector of hit deposits
//std::vector<Double_t> simData(nMaxHits*nStrips,0);
std::vector<double> hitDep;

// vector of hit channels with CC
std::vector<Int_t> simChan;  
// vector of hit deposits with CC
//std::vector<Double_t> simData(nMaxHits*nStrips,0);
std::vector<Double_t> simDep;


//----------------------------------

void CleanEvent(){

  nInt = 0;
  //  iEv = 0;//not initialize this!!!
  nHits = 0;
  nTotalHits = 0;
  ppHit = -1;
  nroSt = 0;
  nHS = 0;
  // eLastZ = -1.;
  // pLastZ = -1.;

  // std::fill_n(hPart, nMaxHits, 0);
  std::fill_n(hVol, nMaxHits, 0);
  std::fill_n(hVolZ, nMaxHits, 0.);
  
  // //  Double_t zPath[nMaxTotalHits];
  gEne=-9999;
  gX=-9999;
  gY=-9999; 
  // zPath=-9999;
  std::fill_n(xCoord, nMaxTotalHits, -9999.9);
  std::fill_n(yCoord, nMaxTotalHits, -9999.9);
  std::fill_n(zCoord, nMaxTotalHits, -9999.9);
  std::fill_n(eDep, nMaxTotalHits, -9999.9);
  std::fill_n(PDG, nMaxTotalHits, 0);
  std::fill_n(TrID, nMaxTotalHits, -1);
  std::fill_n(ParID, nMaxTotalHits, -1);
  std::fill_n(xMom, nMaxTotalHits, 0.);
  std::fill_n(yMom, nMaxTotalHits, 0.);
  std::fill_n(zMom, nMaxTotalHits, 0.);
  std::fill_n(eEne, nMaxTotalHits, 0.);

  std::fill_n(chXY, nMaxTotalHits, 0);
  //std::fill_n(chY, nMaxTotalHits, 0);

  std::fill_n(hitStrips, nMaxTotalHits, 0);
  std::fill_n(simStrips, nMaxTotalHits, 0);

  hitChan.clear();
  hitDep.clear();
  simChan.clear();
  simDep.clear();
  
  return;
}

//---------------------------------- MAKE LADDERCONF FILE -------------

// now only setting pitch
void ConfLadder(int ntdr, double spi, double kpi){
  /* #JINF TDR SPITCH (mm) KPITCH (mm) SRESO (mm) KRESO(mm) MultiplicyFlip SMirror KMirror BondType SHiThresh SLoThresh KHiThresh KLoThresh
0 0 0.220 0.220 0.010 0.010 0 0 0 0 3.5 1.0 3.5 1.0
0 1 0.220 0.220 0.010 0.010 0 0 0 0 3.5 1.0 3.5 1.0
	
  */
  ofstream fout;
  fout.open ("ladderconf_mc.dat");
  int dummy=0;
  double sRes=0.010;
  double kRes=0.010;
  double sHiT=3.5;
  double sLoT=1.0;
  double kHiT=3.5;
  double kLoT=1.0;
  
  fout<<"#JINF TDR SPITCH (mm) KPITCH (mm) SRESO (mm) KRESO(mm) MultiplicyFlip SMirror KMirror BondType SHiThresh SLoThresh KHiThresh KLoThresh"<<endl; 

  for (int il=0;il<ntdr;il++){
    fout<<dummy<<" "<<il<<" "<<showpoint<<setprecision(3)<<spi<<" "<<kpi<<" "<<setprecision(2)<<sRes<<" "<<kRes<<" "<<dummy<<" "<<dummy<<" "<<dummy<<" "<<dummy<<" "<<sHiT<<" "<<sLoT<<" "<<kHiT<<" "<<kLoT<<endl;

  }
  
  fout.close ();
  
}

//---------------------------------- MAKE ALIGNMENT.DAT FILE -------------

// now only setting Z of layers
void AlignmentMC(int ntdr, double * alz){
  // #JINF TDR S (mm) K (mm) Z (mm) MultiplicyFlip
//0 0 225.0 225.0 -299.925 0 
//0 1 225.0 225.0 -269.925 0 
//...
//0 23 225.0 225.0 750.075 0 	
//
  
  ofstream alfout;
  alfout.open ("alignment_mc.dat");
  int jinf=0;
  double als=225.0;
  double alk=225.0;
  int mflip=0;
  
  alfout<<"#JINF TDR S (mm) K (mm) Z (mm) MultiplicyFlip"<<endl; 

  for (int il=0;il<ntdr;il++){
    alfout<<jinf<<" "<<il<<" "<<showpoint<<setprecision(4)<<als<<" "<<alk<<" "<<setprecision(6)<<alz[il]<<" "<<mflip<<endl;
  }
  
  alfout.close ();
  
}

//cout << fixed << showpoint;
//    cout << setprecision(2);
// --------------------------------------------


void SimpleAnalysis(TString inputFileName, TString outputFileName,bool full50=1) {
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
    //    std::cout << "Geometry version   : " << simInfo->geometryVersion << "\n";
    std::cout << "Full version : " << full50 << "\n";
  }
  else {
    std::cout << "*** No simulation information are available\n";
  }  
  const GGSTGeoParams *geoParams = reader.GetGeoParams();
  /*
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
  std::cout << std::endl;*/
  
  bool DB=false;  
  // GEOMETRYPARS TO GET FROM READER
  
  // double mside = 7.04; // cm sensor measurement side length
  //double pitch = 0.0110; // cm sensor pitch
  
  //MD: FIX ME
  double msidex = geoParams->GetRealGeoParam("tileX"); // cm sensor measurement side length
  double msidey = geoParams->GetRealGeoParam("tileY"); // cm sensor measurement side length
  double spitch = geoParams->GetRealGeoParam("tileSPitch")/10.; // cm sensor pitch
  double kpitch = geoParams->GetRealGeoParam("tileKPitch")/10.; // cm sensor pitch
  //double kpitch = geoParams->GetRealGeoParam("tileSPitch")/10.; // cm sensor pitch
  double msidexy=0.;

  
  ///add insensitive border in geoParams
  //  double mborderx= msidex-pitch*nstr;
  
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
  //Set which hit detectors are to be read
  //hReader->SetDetector("siLayerMiniLog", kTRUE); // The name is the same of the sensitive logical volume name in the simulation
  //hReader->SetDetector("siLayerShortLog", kTRUE); // The name is the same of the sensitive logical volume name in the simulation
  //hReader->SetDetector("siLayerLongLog", kTRUE); // The name is the same of the sensitive logical volume name in the simulation
  if(!full50){
    hReader->SetDetector("siAMSTileLog", kTRUE); // The name is the same of the sensitive logical volume name in the simulation
    hReader->SetDetector("siDAMPETileLog", kTRUE); // The name is the same of the sensitive logical volume name in the simulation
  }
  hReader->SetDetector("siTileLog", kTRUE); // The name is the same of the sensitive logical volume name in the simulation

  //   hReader->SetDetector("siStripLog", kTRUE); // The name is the same of the sensitive logical volume name in the simulation
   
  // Retrieve the MC truth sub-reader
  GGSTMCTruthReader *mcReader = reader.GetReader<GGSTMCTruthReader>();

  // Create and retrieve the hadronic interaction sub-reader
  GGSTHadrIntReader *hadrReader = reader.GetReader<GGSTHadrIntReader>();


  // Prepare the histograms
  // gStyle->SetOptFit(1);
  
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

  ConfLadder(nLayers, geoParams->GetRealGeoParam("tileSPitch"), geoParams->GetRealGeoParam("tileKPitch"));
    
  Double_t alza[nLayers];

  // ALIGNMENT of AMS-like tiles
  // ladders in horizontal position always
  // 0: sensors K vertical strips -> measure horizontal coord (X)
  // 1: sensors S horizontal strips -> measure vertical coord (Y) 
  
  // manual alignment for older files
  //Int_t xyAlign[14]={0,1,0,1,0,1,0,1,0,1,1,1,1,1};
  //Int_t xyAlign[20]={0,1,0,1,0,1,0,1,0,1,1,1,1,1,0,0,0,0,0,0};
  //Int_t xyAlign[24]={0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1,1,1,1};
  Int_t xyAlign[nLayers];  

  //if(!full50){
  if(1){
    string sAlign=geoParams->GetStringGeoParam("layerAlignment");
    for (int isn=0; isn<nLayers; isn++) {
      xyAlign[isn]=sAlign[isn]-'0';
      //    printf("xyAlign[%d]=%d\n", isn, xyAlign[isn]);
    }
  }
  
  int nEvts = 0;
  
  //------------------------------------------------------

  TTree *runTree =  new TTree("runTree","tree of mcrun");
  runTree->Branch("nLayers",&nLayers,"nLayers/I");
  runTree->Branch("xyAlign",&xyAlign,"xyAlign[nLayers]/I");
  runTree->Branch("zPos",&alza,"zPos[nLayers]/D");
  runTree->Branch("nEvts",&nEvts,"nEvts/I");
  runTree->Branch("gEne",&gEne,"gEne/D");
  
  TTree *hitTree =  new TTree("hitTree","tree of layer hits");
  
  hitTree->Branch("evID",&iEv,"evID/I");
  hitTree->Branch("nHits",&nHits,"nHits/I");
  //nHitsData=nHits*nStrips;
  //hitTree->Branch("nHitsData",&nHitsData,"nHitsData/I");
  hitTree->Branch("nTotalHits",&nTotalHits,"nTotalHits/I");
  hitTree->Branch("ppHit",&ppHit,"ppHit/I");
  runTree->Branch("gX",&gX,"gX/D");
  runTree->Branch("gY",&gY,"gY/D");
  // hitTree->Branch("eLastZ",&eLastZ,"eLastZ/D");
  // hitTree->Branch("pLastZ",&pLastZ,"pLastZ/D");

  // hit
  // hitTree->Branch("hPart",&hPart,"hPart[nHits]/I");
  hitTree->Branch("hVol",&hVol,"hVol[nTotalHits]/I");  
  hitTree->Branch("hVolZ",&hVolZ,"hVolZ[nTotalHits]/D");  
  
  hitTree->Branch("xCoord",&xCoord,"xCoord[nTotalHits]/D");
  hitTree->Branch("yCoord",&yCoord,"yCoord[nTotalHits]/D");
  hitTree->Branch("zCoord",&zCoord,"zCoord[nTotalHits]/D");  
  hitTree->Branch("eDep",&eDep,"eDep[nTotalHits]/D");
  hitTree->Branch("PDG",&PDG,"PDG[nTotalHits]/I");
  hitTree->Branch("TrID",&TrID,"TrID[nTotalHits]/I");
  hitTree->Branch("ParID",&ParID,"ParID[nTotalHits]/I");
  hitTree->Branch("xMom",&xMom,"xMom[nTotalHits]/D");
  hitTree->Branch("yMom",&yMom,"yMom[nTotalHits]/D");
  hitTree->Branch("zMom",&zMom,"zMom[nTotalHits]/D");
  hitTree->Branch("eEne",&eEne,"eEne[nTotalHits]/D");
  hitTree->Branch("chXY",&chXY,"chXY[nTotalHits]/I");
  //  hitTree->Branch("chY",&chY,"chY[nTotalHits]/I");
  hitTree->Branch("hitStrips",&hitStrips,"hitStrips[nTotalHits]/I");
  hitTree->Branch("simStrips",&simStrips,"simStrips[nTotalHits]/I");
  //  hitTree->Branch("simData",&simData,"simData[nHits*nStrips]/D");
  hitTree->Branch("hitChan",&hitChan);
  hitTree->Branch("hitDep",&hitDep);
  hitTree->Branch("simChan",&simChan);
  hitTree->Branch("simDep",&simDep);
  
  double mom=0.;
  int eno=0;
  int pno=0; 

  nEvts=reader.GetEntries();
  
  COUT(INFO) << "Begin loop over " << nEvts << " events" << ENDL;
  //  int dones=2;  
  for (iEv = 0; iEv < nEvts; iEv++) {
    //    for (int iEv = 0; iEv < 33; iEv++) {
    
    CleanEvent();
    
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
	thisHit = static_cast<GGSTIntHit*>(hReader->GetHit("siTileLog", iHit));
	vTIntHit.push_back(thisHit);
      }
      nHits += _nHits;
    }
    if(!full50){
      {
        int _nHits = hReader->GetNHits("siAMSTileLog"); //Number of hit siLayers for current event
        for (int iHit = 0; iHit < _nHits; iHit++) {
          thisHit = static_cast<GGSTIntHit*>(hReader->GetHit("siAMSTileLog", iHit));
          vTIntHit.push_back(thisHit);
        }
        nHits += _nHits;
      }
      {
        int _nHits = hReader->GetNHits("siDAMPETileLog"); //Number of hit siLayers for current event
        for (int iHit = 0; iHit < _nHits; iHit++) {
          thisHit = static_cast<GGSTIntHit*>(hReader->GetHit("siDAMPETileLog", iHit));
          vTIntHit.push_back(thisHit);
        }
        nHits += _nHits;
      }
    }
    
    
    std::cout<<"EVT "<<iEv<<" NHITS "<<nHits<<std::endl;
    //nHitHisto->Fill(nHits);
    
    float totEDep = 0.;   
   
    // // efficiency calc
    bool pprod=false; 

    // true Y bending
    //double mom=1.;  // mom in GV (not geant4 ene)
    /*    if (hReader->GetHit("siTileLog",0)) {
	  mom=hReader->GetHit("siTileLog",0)->GetPartHit(0)->entranceMomentum[2];
	  }
	  else {
	  mom=0.0;
	  }
    */
    // //double magf=1.; // mag field [tesla]
    // double magf=geoParams->GetRealGeoParam("magFieldVal");
    // double magz=geoParams->GetRealGeoParam("magVolZ")/100.; // mag z size [m]
    // double ro=0.;
    // double dYtrue=0.;
    // double thtrue=0.;
    // double thcenter=0.;
    // double thtrueV=0.;
    
    //if(magf){
    //  ro=mom/magf/0.3;  // radius of curvature [m]
    //  dYtrue = ro*(1-sqrt(1-(magz/ro)*(magz/ro)));
    //  thtrue = atan2(dYtrue,magz);
    //  thcenter = atan2(magz,(ro-dYtrue));
    //}
 
    // TVector3 vzeta(0.,0.,mom);
    // //    TVector3 vtrue(0.,dYtrue,mom);
    // TVector3 vtrue(0.,mom*sin(thcenter),mom*cos(thcenter));
    // TVector3 vrec;

    // thtrueV= vtrue.Angle(vzeta);

    // if(DB)
    //   std::cout<<"DY true: "<<mom<<" "<<ro<<" "<<dYtrue*100<<" TH "<<thtrue*180/TMath::Pi()<<" THC "<<thcenter*180/TMath::Pi()<<" THV "<<thtrueV*180/TMath::Pi()<<std::endl;
    
    // // Hits loop
    // double yfront=0.;
    // double yend1=0.;
    // double yend2=0.;
    // double threc=0.;
    float VolPos=-9999;
    int layerID=-1;
    for (int iHit = 0; iHit < nHits; iHit++) {//una hit e' semplicemente un volume logico che si e' acceso      
      //      hPart[iHit]=0;	
      thisHit = vTIntHit.at(iHit);
      
      if (thisHit->GetVolumePosition()[2]!=VolPos){
        VolPos=thisHit->GetVolumePosition()[2];
        layerID++;
	if(iEv==0)
	  alza[layerID]=VolPos*10.;
      }
      
      int nPHits= thisHit->GetNPartHits();
      
      if(nPHits>=2||pprod) {
	std::cout<<"******* HIT: "<< iHit <<" --->>>PARTHITS: "<<nPHits<<" detind: "<<hReader->GetDetectorIndex("siTileLog")<<" "<<thisHit->GetVolumeName()<<" "<<thisHit->GetVolumeID()<<" pos:"<<thisHit->GetVolumePosition()[2]<<std::endl;
      }

      //      hVol[iHit]=thisHit->GetVolumeID();
      //      hVolZ[iHit]=thisHit->GetVolumePosition()[2];
      
      for (int iPHit =0; iPHit<nPHits; iPHit++){//queste sono le vere hit particella per particella
	nHS=0;        
        thisPHit = static_cast<GGSTPartHit*>(thisHit->GetPartHit(iPHit));
	
        if (nPHits>=2||pprod)
          thisPHit->Dump();
        
	//        if(gEne<0 && thisPHit->parentID==0 && thisPHit->particlePdg==22){ /// primary particle hit values
	        if(gEne<0 && thisPHit->parentID==0){ /// primary particle hit values 
          gEne=1e3*thisPHit->entranceEnergy;
	  gX=thisPHit->entrancePoint[0];
	  gY=thisPHit->entrancePoint[1];
        }
	
      	if(thisPHit->particlePdg!=22){  /// secondary particles (electron and positron)
          nTotalHits++;
	  
          hVolZ[nTotalHits-1]=thisHit->GetVolumePosition()[2];
          //Mappatura volumi
          hVol[nTotalHits-1]=layerID;
          
          xCoord[nTotalHits-1]=(thisPHit->entrancePoint[0]+thisPHit->exitPoint[0])/2;
          yCoord[nTotalHits-1]=(thisPHit->entrancePoint[1]+thisPHit->exitPoint[1])/2;
          zCoord[nTotalHits-1]=(thisPHit->entrancePoint[2]+thisPHit->exitPoint[2])/2;
	  
          eDep[nTotalHits-1]=1e3*thisPHit->eDep;
	  
          xMom[nTotalHits-1]=1e3*thisPHit->entranceMomentum[0];
          yMom[nTotalHits-1]=1e3*thisPHit->entranceMomentum[1];
          zMom[nTotalHits-1]=1e3*thisPHit->entranceMomentum[2];
          eEne[nTotalHits-1]=1e3*thisPHit->entranceEnergy;
          PDG[nTotalHits-1]=thisPHit->particlePdg;
          TrID[nTotalHits-1]=thisPHit->trackID;
          ParID[nTotalHits-1]=thisPHit->parentID;
	  
          //chY[nTotalHits-1]=int((yCoord[nTotalHits-1]+msidey/2.)/ropitch);
	  //int xych=xyAlign[hVol[nTotalHits-1]]?chY[nTotalHits-1]:chX[nTotalHits-1];
	  
	  //	  nroSt =xyAlign[hVol[nTotalHits-1]]?int(msidey/ropitch):int(msidex/ropitch);
	  nroSt =xyAlign[hVol[nTotalHits-1]]?int(msidey/spitch):int(msidex/kpitch); // Read-out strips are 4545 if side=1(Y<-->S) 2403 if side=0 (X<-->K) 
	  
	  
	  /// checking for readout strip channels 
	  int ichxy=-1;
	  int ochxy=-1;
	  double ixy=0.;
	  double oxy=0.;
	  double ropitch=0.;

	  if(xyAlign[hVol[nTotalHits-1]]){
	    ropitch=spitch;
	    ixy=thisPHit->entrancePoint[1];
	    oxy=thisPHit->exitPoint[1];
	    chXY[nTotalHits-1]=int(((yCoord[nTotalHits-1])+msidey/2.)/ropitch);
	    ichxy=int((ixy+msidey/2.)/ropitch);
	    ochxy=int((oxy+msidey/2.)/ropitch);
	    msidexy=msidey;
	  }else{
	    ropitch=kpitch;
	    ixy=thisPHit->entrancePoint[0];
	    oxy=thisPHit->exitPoint[0];
	    chXY[nTotalHits-1]=int(((xCoord[nTotalHits-1])+msidex/2.)/ropitch); 
	    ichxy=int((ixy+msidex/2.)/ropitch);
	    ochxy=int((oxy+msidex/2.)/ropitch);
	    msidexy=msidex;
	  }
	    
	  double impitch=0.5*ropitch; /// implant pitch
	  
	  /// digital reso	  
	  /// check for inclined tracks crossing multiple strips
	  /// fill all strips with appropriate eDep and include capacitive coupling (3 left and 3 right)
	  
	  hitStrips[nTotalHits-1]=abs(ochxy-ichxy)+1;
	  simStrips[nTotalHits-1]=2*cpw;
	  	  
	 	  
	  bool iro=false;
	  bool oro=false;
	  
	  double itailfr=0.;
	  double iheadfr=0.;
	  double otailfr=0.;
	  double oheadfr=0.;
	  
	  double trl=abs(ixy-oxy);
	  
	  double dich = (ixy+msidexy/2.)/ropitch-ichxy;
	  if(dich<=0.5){
	    // inch isreadout
	      iro=true;
	      itailfr=(ixy+msidexy/2.-ichxy*ropitch)/trl; 
	      iheadfr=((ichxy+0.5)*ropitch-ixy-msidexy/2.)/trl;
	  }else{
	    itailfr=(ixy+msidexy/2.-(ichxy+0.5)*ropitch)/trl;
	    iheadfr=((ichxy+1)*ropitch-ixy-msidexy/2.)/trl;
	    /// if(ochxy<ichxy)
	    if(oxy<ixy) 
	      hitStrips[nTotalHits-1]+=1; /// in ch is tail and not ro-> add right strip (ichxy+1)
	  }
	  
	  double doch = (oxy+msidexy/2.)/ropitch-ochxy;
	  if(doch<=0.5){
	    // outch is readout
	    oro=true;
	    otailfr=(oxy+msidexy/2.-ochxy*ropitch)/trl; 
	    oheadfr=((ochxy+0.5)*ropitch-oxy-msidexy/2.)/trl;
	  }else{
	    otailfr=(oxy+msidexy/2.-(ochxy+0.5)*ropitch)/trl; 
	    oheadfr=((ochxy+1)*ropitch-oxy-msidexy/2.)/trl;
	    /// if(ochxy>ichxy)
	    if(oxy>ixy) 
	      hitStrips[nTotalHits-1]+=1; /// out ch is tail and not ro -> add right strip (ochxy+1)
	  }
	  
	  
	  int strn=0;
	  int swn=0;
	  /// fully contained strip Dep frac
	  double  sfrac = 0.5*ropitch/trl; // strip implant = 0.5*ropitch
	  double iDep[hitStrips[nTotalHits-1]];
	  for(int i=0;i<hitStrips[nTotalHits-1];i++){
	    iDep[i]=0.;
	  }
	    

	  if (ochxy==ichxy){
	    if(iro){
	      iDep[0]=eDep[nTotalHits-1]*(oro?1.:(iheadfr+otailfr/2.));
	      if(!oro)
		iDep[1]=eDep[nTotalHits-1]*(otailfr/2.);
	    
	    if(pri)
	      cout<<"------------>SINGLE IRO "<<(oro?"ORO: ":": ")<<iEv<<" "<<hVol[nTotalHits-1]<<" AL (1=Y, 0=X) "<<xyAlign[hVol[nTotalHits-1]]<<" PDG "<<PDG[nTotalHits-1]<<" PAR "<<ParID[nTotalHits-1]<<" TrID "<<TrID[nTotalHits]<<" CH "<<chXY[nTotalHits-1]<<" DCH I/O: "<<dich<<" "<<doch<<" EDEP "<<eDep[nTotalHits-1]<<" iDep:"<<iDep[0]<<" "<<(oro?iDep[0]:iDep[1])<<endl;

	    }
	    else{
	      iDep[0]=eDep[nTotalHits-1]*(oro?(itailfr/2. + oheadfr):0.5);
	      iDep[1]=eDep[nTotalHits-1]*(oro?itailfr/2.:0.5);

	      cout<<"------------>SINGLE INRO "<<(oro?"ORO: ":": ")<<iEv<<" "<<hVol[nTotalHits-1]<<" AL (1=Y, 0=X) "<<xyAlign[hVol[nTotalHits-1]]<<" PDG "<<PDG[nTotalHits-1]<<" PAR "<<ParID[nTotalHits-1]<<" TrID "<<TrID[nTotalHits]<<" CHS "<<ichxy<<" "<<ichxy+1<<" HSTR "<< hitStrips[nTotalHits-1]<<" DCH I/O: "<<dich<<" "<<doch<<" EDEP "<<eDep[nTotalHits-1]<<" iDep:"<<iDep[0]<<" "<<iDep[1]<<endl;
	    }		
	    
	  }
	  else{ // ochxy!=ichxy
	    
	    
	    if(ochxy>ichxy){ /// in ch is head out ch is tail	    
	    //hitChan.push_back(ichxy);
	    
	      if(oro){
		strn=hitStrips[nTotalHits-1]-1;
		iDep[strn]=eDep[nTotalHits-1]*(otailfr + sfrac/2.);
	      }else{
		strn=hitStrips[nTotalHits-1]-2; // as hitStrips has been already incremented
		iDep[strn]=eDep[nTotalHits-1]*(otailfr/2. + sfrac*3./2.);
		iDep[hitStrips[nTotalHits-1]-1]=eDep[nTotalHits-1]*otailfr/2.;
	      }
	      
	      if(strn>1){ // check for vertical
		for (int hs=1;hs<strn;hs++){
		  //hitChan.push_back(ichxy+hs);
		  iDep[hs]=eDep[nTotalHits-1]*2*sfrac;  	  
		}	      
	      }
	     
	      //hitChan.push_back(ochxy);
	      //if(!oro)
	      //hitChan.push_back(ochxy+1);
	      
	      if(iro){
		iDep[0]=eDep[nTotalHits-1]*(iheadfr+sfrac/2.);
	      }else{
		iDep[0]=eDep[nTotalHits-1]*iheadfr/2.;
		iDep[1]+=(eDep[nTotalHits-1]*(iheadfr/2. - sfrac/2.));
	      }
	      
	      // if (eDep[nTotalHits-1]==0)
	      //	hitDep.insert(hitDep.end(),hitStrips[nTotalHits-1],0.);
	      //else
	      //	hitDep.insert(hitDep.end(),iDep,iDep+hitStrips[nTotalHits-1]);
	      
	    }else{ /// ochxy<ichxy
	      //hitChan.push_back(ochxy);
	      
	      if(iro){
		strn=hitStrips[nTotalHits-1]-1;
		iDep[strn]=eDep[nTotalHits-1]*(itailfr+sfrac/2.);
	      }
	      else{
		strn=hitStrips[nTotalHits-1]-2; // as hitStrips has been already incremented
		iDep[strn]=eDep[nTotalHits-1]*(itailfr/2. + sfrac*3/2.);
		iDep[hitStrips[nTotalHits-1]-1]=eDep[nTotalHits-1]*itailfr/2.;
	      }
	      
	      if(strn>1){ // check for vertical
		for (int hs=1;hs<strn;hs++){
		  //hitChan.push_back(ichxy+hs);
		  iDep[hs]=eDep[nTotalHits-1]*2*sfrac;  	  
		}	      
	      }
	      //hitChan.push_back(ichxy);
	      //	      if(!iro)
	      //		hitChan.push_back(ichxy+1);
	      
	      if(oro){
		iDep[0]=eDep[nTotalHits-1]*(oheadfr+sfrac/2.);
	      }else{
		iDep[0]=eDep[nTotalHits-1]*iheadfr/2.;
		iDep[1]+=(eDep[nTotalHits-1]*(iheadfr/2.-sfrac/2.));
	      }
	      
	    }/// ochxy<ichxy
	    
	    if(pri)
	      cout<<"------------>INCLINED: "<<iEv<<" "<<hVol[nTotalHits-1]<<" AL (1=Y, 0=X) "<<xyAlign[hVol[nTotalHits-1]]<<" PDG "<<PDG[nTotalHits-1]<<" PAR "<<ParID[nTotalHits-1]<<" TrID "<<TrID[nTotalHits]<<" CoG "<<chXY[nTotalHits-1]<<" HSTR "<<hitStrips[nTotalHits-1]<<" TRL "<<trl<<" EDEP "<<eDep[nTotalHits-1]<<"  CHX IN "<<ichxy<<" OUT "<<ochxy<<"  HFR "<<(ichxy<ochxy?iheadfr:oheadfr)<<" TFR "<<(ichxy<ochxy?otailfr:itailfr)<<endl;

	  }/// ichxy!=ochxy

	  //cout<<"INCL HIT VALS: "<<hitStrips[nTotalHits-1]<<" ";
	  for (int id=0; id<hitStrips[nTotalHits-1];id++){
	    //cout<<(ochxy>ichxy?ichxy:ochxy)+id<<":"<<iDep[id]<<" ";
	    //hitChan.push_back((ochxy>ichxy?ichxy:ochxy)+id);
	    hitChan.push_back((oxy>ixy?ichxy:ochxy)+id);
	    hitDep.push_back(iDep[id]);
	  }
	      //cout<<endl;
	    
	    
	  /// now apply capacitive coupling 
	  
	  double cciDep[hitStrips[nTotalHits-1]+2*cpw];
	  for(int i=0;i<hitStrips[nTotalHits-1]+2*cpw;i++)
	    cciDep[i]=0.;
	    
	  for (int hs=0;hs<hitStrips[nTotalHits-1];hs++){
	    simStrips[nTotalHits-1]+=1;
	    for (int icc=0;icc<2*cpw+1;icc++){
	      cciDep[icc+hs]+=chcp[icc]*iDep[hs];
	    }
	  }
	  //	    cout<<endl;
	  if(pri)
	    cout<<"**** SIM VALS: "<<simStrips[nTotalHits-1]<<" ";
	  for (int id=0; id<hitStrips[nTotalHits-1]+2*cpw;id++){
	    if(pri)
	      cout<<(oxy>ixy?ichxy:ochxy)-cpw+id<<":"<<cciDep[id]<<" ";
		
	    if((oxy>ixy?ichxy:ochxy)-cpw+id>=0&&(oxy>ixy?ichxy:ochxy)-cpw+id<=nroSt){
	      simChan.push_back((oxy>ixy?ichxy:ochxy)-cpw+id);
	      simDep.push_back(cciDep[id]);
	    }else{
	      simStrips[nTotalHits-1]-=1;
	      if(pri)
		cout<<"SENSOR BORDER REACHED: "<<chXY[nTotalHits-1]-cpw+id<<" "<<cciDep[id]<<" LOST SIGNAL"<<endl;
	    }
	  }
	  cout<<endl;
	  
	  if(thisPHit->parentID==1&&thisPHit->particlePdg==-11){ /// positron
	    if (!pprod){
	      pprod=true;
	      ppHit=iHit;
	    }
	  }
	}
	
      } // loop on particle hits
      
      
      // //eDep[nTotalHits-1]=1e3*thisHit->eDep;
      // totEDep += thisHit->eDep;
      // if(DB)
      // 	std::cout<<"*********************  EVT: "<<iEv<<" HIT: "<< iHit <<" nPART: "<< hPart[iHit]<< std::endl;
      
    }// loop on hits
    

    // if(ppHit>=0){
    //   // layer ch part multiplicity
    //   Int_t lpmult=0;
    //   int fhz=nHits;
    //   eLastZ=hVolZ[ppHit];
    //   pLastZ=hVolZ[ppHit];
    //   double tmplz=hVolZ[ppHit];
    //   for (int ihh=ppHit+1;ihh<nHits;ihh++){

    // 	lpmult=hPart[ihh];    	  
    // 	if(lpmult==0)
    // 	  continue;
	
    // 	if (lpmult==3){
    // 	  eLastZ=hVolZ[ihh];
    // 	  pLastZ=hVolZ[ihh];
    // 	  if(DB)
    // 	    std::cout<<"EVT: "<<iEv<<" HIT: "<< ihh <<" VZ: "<<hVolZ[ihh]<<" HPART: "<<hPart[ihh]<<" LPMULT: "<<lpmult<<" EPLAST: "<< eLastZ<<" "<<pLastZ<< std::endl;
    // 	  continue;
    // 	}else{
	  
    // 	  if(hPart[ihh]==1)
    // 	    eLastZ=hVolZ[ihh];
    // 	  if(hPart[ihh]==2)
    // 	    pLastZ=hVolZ[ihh];
    // 	  if(DB)
    // 	    std::cout<<"***EVT: "<<iEv<<" HIT: "<< ihh <<" VZ: "<<hVolZ[ihh]<<" HPART: "<<hPart[ihh]<<" LPMULT: "<<lpmult<<" EPLAST: "<< eLastZ<<" "<<pLastZ<<std::endl;
	  
    // 	  double hvz=hVolZ[ihh];
	  
    // 	  for (int jhh=ihh+1;jhh<nHits;jhh++){
	    
    // 	    if (hVolZ[jhh]==hvz&&hPart[jhh]!=0){
    // 	      lpmult+=hPart[jhh];
    // 	      if(hPart[jhh]==1)
    // 		eLastZ=hVolZ[jhh];
    // 	      if(hPart[jhh]==2)
    // 		pLastZ=hVolZ[jhh];	    
    // 	      fhz=jhh;
    // 	      if(DB)
    // 		std::cout<<"---> EVT: "<<iEv<<" HIT: "<< fhz <<" VZ: "<<hVolZ[jhh]<<" HPART: "<<hPart[jhh]<<" LPMULT: "<<lpmult<<" RECEPLASTZ: "<< eLastZ<<" "<<pLastZ<< std::endl;
    // 	    }else{
    // 	      continue;
    // 	    }
    // 	  }

    // 	}
    //   }
    // }//IF PPROD

    if (!intInfo&&pprod) {
      
      cout<<"********* SUMMARY EVENT DATA: "<<iEv<<" nTotalHits "<<nTotalHits<<" SIMD: "<<simDep.size()<<endl;  
      hitTree->Fill();
      int inx=0;
      for(int in=0;in<nTotalHits;in++){
	
	cout<<"LAYER "<<hVol[in]<<" PDG "<<PDG[in]<<" ParID "<<ParID[in]<<": "<<endl;
	for (int id=0;id<simStrips[in];id++)
	  {
	    cout <<simChan[inx]<<" "<<simDep[inx]<<endl;
	    inx++;
	  }
	cout<<endl;
      }
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
  AlignmentMC(nLayers,alza);  
  COUT(INFO) << "Analysis finished" << ENDL;
  
  return;
}
