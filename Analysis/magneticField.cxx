#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>

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
#include <vector>

#include "TDatabasePDG.h"
#include <TMath.h>

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
const int nHitsMax = 1000;
int evID=-9999;
int nHits=0;
int nTotalHits=0;
Double_t xCoord[nHitsMax]={-9999.};
Double_t yCoord[nHitsMax]={-9999.};
Double_t zCoord[nHitsMax]={-9999.};
Double_t xMom[nHitsMax]={-9999.};
Double_t yMom[nHitsMax]={-9999.};
Double_t zMom[nHitsMax]={-9999.};
Double_t eDep[nHitsMax]={-9999.};
Double_t eEne[nHitsMax]={-9999.};
int PDG[nHitsMax]={0};
int TrID[nHitsMax]={-1};
int ParID[nHitsMax]={-1};
int nEvs=0;
void readFile(){
  
  TString inputFileName="anaOut.root";
  TFile *inFile=new TFile(inputFileName,"READ");

  TTree *T=(TTree*)inFile->Get("hitTree");

  T->SetBranchAddress( "evID", &evID );
  T->SetBranchAddress( "nHits", &nHits );
  T->SetBranchAddress( "nTotalHits", &nTotalHits );
  T->SetBranchAddress( "xCoord", &xCoord );
  T->SetBranchAddress( "yCoord", &yCoord );
  T->SetBranchAddress( "zCoord", &zCoord );
  T->SetBranchAddress( "xMom", &xMom );
  T->SetBranchAddress( "yMom", &yMom );
  T->SetBranchAddress( "zMom", &zMom );
  T->SetBranchAddress( "eEne", &eEne );
  T->SetBranchAddress( "eDep", &eDep );
  T->SetBranchAddress( "PDG", &PDG );
  T->SetBranchAddress( "TrID", &TrID );
  T->SetBranchAddress( "ParID", &ParID );

  nEvs=T->GetEntries();
  std::cout<<"got "<<nEvs<<" muons"<<std::endl;
}

int main() {
  // reading the root file
  // init MeasurementCreator
  genfit::MeasurementCreator measurementCreator;
  // init geometry and mag. field
  TGeoManager("Geometry", "Geane geometry"); 
  TGeoManager::Import("~/Documents/c++/thesis/POXsimu_build/plugins()/libTestGeometry.vgm.root");
  
  TGeoVolume *magnet = gGeoManager->GetVolume("magnet");
  TGeoUniformMagField *magField = new TGeoUniformMagField(0.,0.,0.05);
  magnet->SetField(magField);
  //to fix like prof. Duranti said
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  //look if you can use it sourcing them
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

  // main loop (loops on the events)
  for (unsigned int iEvent=0; iEvent < nEvs; iEvent++){
    // true start values                                      
    TVector3 pos(0, 0, 0);
    pos.SetX(xCoord[iEvent]);
    pos.SetY(yCoord[iEvent]);
    pos.SetZ(zCoord[iEvent]);
    TVector3 mom(1.,0,0);
    mom.SetPhi(TMath::ATan(xMom[iEvent] / yMom[iEvent]));
    mom.SetMag(TMath::Sqrt( TMath::Power(xMom[iEvent],2) + TMath::Power(yMom[iEvent],2) + TMath::Power(zMom[iEvent],2)));
    mom.SetTheta(TMath::ACos(zMom[iEvent] / mom.Mag() ));
    
    // helix track model    
    //get the charge of the particle
    const double charge = TDatabasePDG::Instance()->GetParticle(PDG[iEvent])->Charge()/(3.);
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
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(PDG[iEvent]);

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


    if (iEvent < 1000) {
      // add track to event display
      display->addEvent(&fitTrack);
    }



  }// end loop over events

  delete fitter;

  // open event display
  display->open();

}

