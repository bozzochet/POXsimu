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

// Event vars-----------------------
int nInt = 0;
int iEv = 0;
int nHits = 0;//number of TIntHits (essentially the number of firing volumes)
int nTotalHits = 0;//total number of TPathHits (~ nHits*<nPartHits>)
int ppHit = -1;
// double eLastZ = -1.;
// double pLastZ = -1.;

// Int_t hPart[nMaxHits]={0};
Int_t hVol[nMaxTotalHits]={0};
Double_t hVolZ[nMaxTotalHits]={0.};

// //  Double_t zPath[nMaxTotalHits];
Double_t gEne=-9999; 
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

Int_t chX[nMaxTotalHits]={0};
Int_t chY[nMaxTotalHits]={0}; 

//----------------------------------

void CleanEvent(){

  nInt = 0;
  //  iEv = 0;//not initialize this!!!
  nHits = 0;
  nTotalHits = 0;
  ppHit = -1;
  // eLastZ = -1.;
  // pLastZ = -1.;

  // std::fill_n(hPart, nMaxHits, 0);
  std::fill_n(hVol, nMaxHits, 0);
  std::fill_n(hVolZ, nMaxHits, 0.);
  
  // //  Double_t zPath[nMaxTotalHits];
  gEne=-9999; 
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

  std::fill_n(chX, nMaxTotalHits, 0);
  std::fill_n(chY, nMaxTotalHits, 0);

  return;
}

void SimpleAnalysis(TString inputFileName, TString outputFileName,bool full50=0) {
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
  // manual alignment for older files
  Int_t xyAlign[20]={0,1,0,1,0,1,0,1,0,1,1,1,1,1,0,0,0,0,0,0};
  //Int_t xyAlign[24]={0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1,1,1,1};

  if(!full50){
  //  Int_t xyAlign[nLayers];
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
  runTree->Branch("nEvts",&nEvts,"nEvts/I");
  runTree->Branch("gEne",&gEne,"gEne/D");
  
  TTree *hitTree =  new TTree("hitTree","tree of layer hits");

  
  hitTree->Branch("evID",&iEv,"evID/I");
  hitTree->Branch("nHits",&nHits,"nHits/I");
  hitTree->Branch("nTotalHits",&nTotalHits,"nTotalHits/I");
  hitTree->Branch("ppHit",&ppHit,"ppHit/I");
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
  hitTree->Branch("chX",&chX,"chX[nTotalHits]/I");
  hitTree->Branch("chY",&chY,"chY[nTotalHits]/I");

  double mom=0.;
  int eno=0;
  int pno=0; 

  nEvts=reader.GetEntries();
  
  COUT(INFO) << "Begin loop over " << nEvts << " events" << ENDL;
  
  for (iEv = 0; iEv < nEvts; iEv++) {
    //for (int iEv = 0; iEv < 33; iEv++) {

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
	      thisHit = hReader->GetHit("siTileLog", iHit);
	      vTIntHit.push_back(thisHit);
      }
      nHits += _nHits;
    }
    if(!full50){
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
    }
    std::cout<<"EVT "<<iEv<<" NHITS "<<nHits<<std::endl;
    //nHitHisto->Fill(nHits);
    
    float totEDep = 0.;   
   
    // // efficiency calc
    bool pprod=false; 

    // true Y bending
    //double mom=1.;  // mom in GV (not geant4 ene)
    if (hReader->GetHit("siTileLog",0)) {
      mom=hReader->GetHit("siTileLog",0)->GetPartHit(0)->entranceMomentum[2];
    }
    else {
      mom=0.0;
    }
    
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
      }

      int nPHits= thisHit->GetNPartHits();
      
      if(nPHits>=2||pprod) {
	      std::cout<<"******* HIT: "<< iHit <<" --->>>PARTHITS: "<<nPHits<<" detind: "<<hReader->GetDetectorIndex("siTileLog")<<" "<<thisHit->GetVolumeName()<<" "<<thisHit->GetVolumeID()<<" pos:"<<thisHit->GetVolumePosition()[2]<<std::endl;
      }

      //      hVol[iHit]=thisHit->GetVolumeID();
      //      hVolZ[iHit]=thisHit->GetVolumePosition()[2];
      
      if (nPHits>=3) {
	      std::cout<<"*********************  EVT: "<<iEv<<" HIT: "<< iHit <<" --->>>3INT"<< std::endl;
      }

      for (int iPHit =0; iPHit<nPHits; iPHit++){//queste sono le vere hit particella per particella
        
        thisPHit = thisHit->GetPartHit(iPHit);

        if (nPHits>=2||pprod)
          thisPHit->DumpHit();
        
        if(gEne<0 && thisPHit->parentID==0 && thisPHit->particlePdg==22){ /// primary particle hit values 
          gEne=1e3*thisPHit->entranceEnergy;
        }

      	if(thisPHit->particlePdg!=22){  /// secondary particles (electron and positron)
          nTotalHits++;

          hVolZ[nTotalHits-1]=thisHit->GetVolumePosition()[2];
          //Mappatura volumi
          hVol[nTotalHits-1]=layerID;
          
          xCoord[nTotalHits-1]=(thisPHit->entrancePoint[0]+thisPHit->exitPoint[0])/2;
          yCoord[nTotalHits-1]=(thisPHit->entrancePoint[1]+thisPHit->exitPoint[1])/2;
          //	  zCoord[nTotalHits-1]=(thisPHit->entrancePoint[2]+thisPHit->exitPoint[2])/2;
          zCoord[nTotalHits-1]=thisHit->GetVolumePosition()[2];
          eDep[nTotalHits-1]=1e3*thisPHit->eDep;
          xMom[nTotalHits-1]=1e3*thisPHit->entranceMomentum[0];
          yMom[nTotalHits-1]=1e3*thisPHit->entranceMomentum[1];
          zMom[nTotalHits-1]=1e3*thisPHit->entranceMomentum[2];
          eEne[nTotalHits-1]=1e3*thisPHit->entranceEnergy;
          PDG[nTotalHits-1]=thisPHit->particlePdg;
          TrID[nTotalHits-1]=thisPHit->trackID;
          ParID[nTotalHits-1]=thisPHit->parentID;

          chX[nTotalHits-1]=int((xCoord[nTotalHits-1]+msidex/2.)/pitch);
          chY[nTotalHits-1]=int((yCoord[nTotalHits-1]+msidey/2.)/pitch);

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
      
    } // loop on hits

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
  
  return;
}
