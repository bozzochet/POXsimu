#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>
#include "ConstField.h"

#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <TGeoMaterialInterface.h>

#include <EventDisplay.h>

#include <HelixTrackModel.h>
#include <MeasurementCreator.h>

#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TRandom.h>
#include <TVector3.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TVirtualMagField.h>

#include <vector>

/*
struct event{
  double Z;
  double XorY;
  int PDG;
  int TrID;
  int ParID;
  TLorentzVector Mom;
  void operator=(event h){
    Z=h.Z;
    XorY=h.XorY;
    Mom=h.Mom;
    PDG=h.PDG;
    TrID=h.TrID;
    ParID=h.ParID;
    return;
  }
};
*/

int main() {

  // reading the input root file ----------------------
  TString inputFileName="anaOut.root";
  TFile *inFile=new TFile(inputFileName,"READ");
  
  TTree *hitTree=(TTree*)inFile->Get("hitTree");
  
  //Declaration of leaves types
  Int_t           evID;
  Int_t           nHits;
  Int_t           nTotalHits;
  Int_t           ppHit;
  Int_t           hVol[12];
  Double_t        hVolZ[12];
  Double_t        xCoord[12];
  Double_t        yCoord[12];
  Double_t        zCoord[12];
  Double_t        eDep[12];
  Int_t           PDG[12];
  Int_t           TrID[12];
  Int_t           ParID[12];
  Double_t        xMom[12];
  Double_t        yMom[12];
  Double_t        zMom[12];
  Double_t        eEne[12];
  Int_t           chXY[12];
  Int_t           hitStrips[12];
  Int_t           simStrips[12];
  // std::vector<int>     hitChan;
  // std::vector<double>  hitDep;
  // std::vector<int>     simChan;
  // std::vector<double>  simDep;
  
  // Set branch addresses.
  hitTree->SetBranchAddress("evID",&evID);
  hitTree->SetBranchAddress("nHits",&nHits);
  hitTree->SetBranchAddress("nTotalHits",&nTotalHits);
  hitTree->SetBranchAddress("ppHit",&ppHit);
  hitTree->SetBranchAddress("hVol",hVol);
  hitTree->SetBranchAddress("hVolZ",hVolZ);
  hitTree->SetBranchAddress("xCoord",xCoord);
  hitTree->SetBranchAddress("yCoord",yCoord);
  hitTree->SetBranchAddress("zCoord",zCoord);
  hitTree->SetBranchAddress("eDep",eDep);
  hitTree->SetBranchAddress("PDG",PDG);
  hitTree->SetBranchAddress("TrID",TrID);
  hitTree->SetBranchAddress("ParID",ParID);
  hitTree->SetBranchAddress("xMom",xMom);
  hitTree->SetBranchAddress("yMom",yMom);
  hitTree->SetBranchAddress("zMom",zMom);
  hitTree->SetBranchAddress("eEne",eEne);
  hitTree->SetBranchAddress("chXY",chXY);
  hitTree->SetBranchAddress("hitStrips",hitStrips);
  hitTree->SetBranchAddress("simStrips",simStrips);
  // hitTree->SetBranchAddress("hitChan",&hitChan);
  // hitTree->SetBranchAddress("hitDep",&hitDep);
  // hitTree->SetBranchAddress("simChan",&simChan);
  // hitTree->SetBranchAddress("simDep",&simDep);
    
  Long64_t nentries = hitTree->GetEntries();
  std::cout<<"got "<<nentries<<" events"<<std::endl;
  Long64_t nbytes = 0;    

  //--------------------------------------------------------
  
  // init MeasurementCreator
  genfit::MeasurementCreator measurementCreator;
  // init geometry and mag. field
  //  gSystem->Load("libGeom");
  //  TGeoManager::Import("~/Documents/c++/thesis/POXsimu_build/plugins/libTestGeometry.vgm.root");
  TGeoManager::Import("plugins/libTestGeometry.vgm.root");
  
  // TGeoVolume *magnet = gGeoManager->GetVolume("magnet");
  // magnet->SetField(new genfit::ConstField(0.,0.,0.5));//0.5 kGauss = 0.05T
  // just to see the Geometry
  //  TGeoVolume *world = gGeoManager->GetVolume("World");
  //  world->Draw();
  //  magnetic field on all the volumes, not ok for the simulation, only for testing
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,15.)); // 15 kGauss

  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  //  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

  // main loop (loops on the events)
  for (Long64_t iEvent=0; iEvent<nentries; iEvent++) {
    nbytes += hitTree->GetEntry(iEvent);
      
    // true start values                                      
    TVector3 pos(0, 0, 0);
    pos.SetX(xCoord[0]);
    pos.SetY(yCoord[0]);
    pos.SetZ(zCoord[0]);
    TVector3 mom(1.,0,0);
    mom.SetPhi(TMath::ATan2(xMom[0], yMom[0]));
    mom.SetMag(TMath::Sqrt( TMath::Power(xMom[0],2) + TMath::Power(yMom[0],2) + TMath::Power(zMom[0],2)));
    mom.SetTheta(TMath::ACos(zMom[0]/mom.Mag()));
    
    // helix track model    
    // get the charge of the particle
    const double charge = TDatabasePDG::Instance()->GetParticle(PDG[0])->Charge()/(3.);
    // TDatabasePDG::Instance()->GetParticle(PDG[0])->Dump();
    // printf("%d %f\n", PDG[0], TDatabasePDG::Instance()->GetParticle(PDG[0])->Charge());
    genfit::HelixTrackModel* helix = new genfit::HelixTrackModel(pos, mom, charge);
    measurementCreator.setTrackModel(helix);

    unsigned int nMeasurements = gRandom->Uniform(5, 15);

    // smeared start values                                                                                                                                                            
    const bool smearPosMom = true;     // init the Reps with smeared pos and mom         
    const double posSmear = 0.1;     // cm                                               
    const double momSmear = 3. /180.*TMath::Pi();     // rad                             
    const double momMagSmear = 0.1;   // relative                                          
    TVector3 posM(pos);
    TVector3 momM(mom);
    if (smearPosMom) {
      posM.SetX(gRandom->Gaus(posM.X(),posSmear));
      posM.SetY(gRandom->Gaus(posM.Y(),posSmear));
      posM.SetZ(gRandom->Gaus(posM.Z(),posSmear));

      momM.SetPhi(gRandom->Gaus(mom.Phi(),momSmear));
      momM.SetTheta(gRandom->Gaus(mom.Theta(),momSmear));
      momM.SetMag(gRandom->Gaus(mom.Mag(), momMagSmear*mom.Mag()));
    }

    // approximate covariance
    TMatrixDSym covM(6);
    double resolution = 0.01;
    for (int i = 0; i < 3; ++i)
      covM(i,i) = resolution*resolution;
    for (int i = 3; i < 6; ++i)
      covM(i,i) = TMath::Power(resolution / nMeasurements / sqrt(3), 2);

    // trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(PDG[0]);

    // smeared start state
    genfit::MeasuredStateOnPlane stateSmeared(rep);
    stateSmeared.setPosMomCov(posM, momM, covM);

    // create track
    TVectorD seedState(6);
    TMatrixDSym seedCov(6);
    stateSmeared.get6DStateCov(seedState, seedCov);
    genfit::Track fitTrack(rep, seedState, seedCov);

    // create random measurement types
    std::vector<genfit::eMeasurementType> measurementTypes;
    for (unsigned int i = 0; i < nMeasurements; ++i)
      measurementTypes.push_back(genfit::eMeasurementType(gRandom->Uniform(8)));


    // create smeared measurements and add to track
    try{
      for (unsigned int i=0; i<measurementTypes.size(); ++i){
        std::vector<genfit::AbsMeasurement*> measurements = measurementCreator.create(measurementTypes[i], i*5.);
        fitTrack.insertPoint(new genfit::TrackPoint(measurements, &fitTrack));
      }
    }
    catch(genfit::Exception& e){
      std::cerr<<"Exception, next track"<<std::endl;
      std::cerr << e.what();
      continue;
    }

    //check
    fitTrack.checkConsistency();

    // do the fit
    fitter->processTrack(&fitTrack);

    //check
    fitTrack.checkConsistency();


    // if (iEvent < 1000) {
    //   // add track to event display
    //   display->addEvent(&fitTrack);
    // }

  }// end loop over events

  delete fitter;

  // open event display
  //  display->open();

}

