#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackCand.h>
#include <TrackPoint.h>

#include <MeasurementProducer.h>
#include <MeasurementFactory.h>

#include "mySpacepointDetectorHit.h"
#include "mySpacepointMeasurement.h"

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

#include <TDatabasePDG.h>
#include <TGeoMaterialInterface.h>
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
#include <TClonesArray.h>


class field: public genfit::AbsBField{
private:
  double field_[3];
  TGeoVolume* magnet;
public:
  field(double b1, double b2, double b3, TGeoVolume* &mag) {
    magnet = mag;
    field_[0]=b1;
    field_[1]=b2;
    field_[2]=b3;
  };
  ~field(){};
  inline TVector3 get(const TVector3& v) const {
    double x = v.x();
    double y = v.y();
    double z = v.z();
    double Bx;
    double By;
    double Bz;
    get(x, y, z, Bx, By, Bz);
    return TVector3(Bx, By, Bz);
  }
  inline void get(const double& x, const double& y, const double& z, double& Bx, double& By, double& Bz)  const { 
    double Point[3] = {x, y, z};
    double Bfield[3];
    if(magnet->Contains(Point)){
      Bfield[0]=field_[0];
      Bfield[1]=field_[1];
      Bfield[2]=field_[2];
    }else{
      Bfield[0]=0.;
      Bfield[1]=0.;
      Bfield[2]=0.;
    }
    
    Bx = Bfield[0];
    By = Bfield[1];
    Bz = Bfield[2]; 
  }
};


int main() {
  // reading the input root file ----------------------
  TString inputFileName="anaOut.root";
  TFile *inFile=new TFile(inputFileName,"READ");
  
  TTree *hitTree=(TTree*)inFile->Get("hitTree");

  const int MaxTotalHits=500;
  
  //Declaration of leaves types
  Int_t           evID;
  Int_t           nHits;
  Int_t           nTotalHits;
  Int_t           ppHit;
  Int_t           hVol[MaxTotalHits];
  Double_t        hVolZ[MaxTotalHits];
  Double_t        xCoord[MaxTotalHits];
  Double_t        yCoord[MaxTotalHits];
  Double_t        zCoord[MaxTotalHits];
  Double_t        eDep[MaxTotalHits];
  Int_t           PDG[MaxTotalHits];
  Int_t           TrID[MaxTotalHits];
  Int_t           ParID[MaxTotalHits];
  Double_t        xMom[MaxTotalHits];
  Double_t        yMom[MaxTotalHits];
  Double_t        zMom[MaxTotalHits];
  Double_t        eEne[MaxTotalHits];
  Int_t           chXY[MaxTotalHits];
  Int_t           hitStrips[MaxTotalHits];
  Int_t           simStrips[MaxTotalHits];
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
  TGeoManager::Import("plugins/libTestGeometry.vgm.root");
  
  TGeoVolume *magnet = gGeoManager->GetVolume("magnet");
  
  genfit::FieldManager::getInstance()->init(new field(0., 0., 0.5, magnet));//0.5 kGauss = 0.05T   
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());


  // init event display
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();


  // init fitter
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();


  TClonesArray myDetectorHitArray("genfit::mySpacepointDetectorHit");

  // init the factory
  int myDetId(1);
  genfit::MeasurementFactory<genfit::AbsMeasurement> factory;
  genfit::MeasurementProducer<genfit::mySpacepointDetectorHit, genfit::mySpacepointMeasurement> myProducer(&myDetectorHitArray);
  factory.addProducer(myDetId, &myProducer);


  // main loop
  for (unsigned int iEvent=0; iEvent<100; ++iEvent){

    myDetectorHitArray.Clear();

    //TrackCand
    genfit::TrackCand myCand;

    // true start values
    TVector3 pos(0, 0, 0);
    TVector3 mom(1.,0,0);
    mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
    mom.SetTheta(gRandom->Uniform(0.4*TMath::Pi(),0.6*TMath::Pi()));
    mom.SetMag(gRandom->Uniform(0.2, 1.));


    // helix track model
    const int pdg = 13;               // particle pdg code
    const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/(3.);
    genfit::HelixTrackModel* helix = new genfit::HelixTrackModel(pos, mom, charge);
    measurementCreator.setTrackModel(helix);


    unsigned int nMeasurements = gRandom->Uniform(5, 15);

    // covariance
    double resolution = 0.01;
    TMatrixDSym cov(3);
    for (int i = 0; i < 3; ++i)
      cov(i,i) = resolution*resolution;

    for (unsigned int i=0; i<nMeasurements; ++i) {
      // "simulate" the detector
      TVector3 currentPos = helix->getPos(i*2.);
      currentPos.SetX(gRandom->Gaus(currentPos.X(), resolution));
      currentPos.SetY(gRandom->Gaus(currentPos.Y(), resolution));
      currentPos.SetZ(gRandom->Gaus(currentPos.Z(), resolution));

      // Fill the TClonesArray and the TrackCand
      // In a real experiment, you detector code would deliver mySpacepointDetectorHits and fill the TClonesArray.
      // The patternRecognition would create the TrackCand.
      new(myDetectorHitArray[i]) genfit::mySpacepointDetectorHit(currentPos, cov);
      myCand.addHit(myDetId, i);
    }


    // smeared start values (would come from the pattern recognition)
    const bool smearPosMom = true;   // init the Reps with smeared pos and mom
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

    // initial guess for cov
    TMatrixDSym covSeed(6);
    for (int i = 0; i < 3; ++i)
      covSeed(i,i) = resolution*resolution;
    for (int i = 3; i < 6; ++i)
      covSeed(i,i) = pow(resolution / nMeasurements / sqrt(3), 2);


    // set start values and pdg to cand
    myCand.setPosMomSeedAndPdgCode(posM, momM, pdg);
    myCand.setCovSeed(covSeed);


    // create track
    genfit::Track fitTrack(myCand, factory, new genfit::RKTrackRep(pdg));


    // do the fit
    try{
      fitter->processTrack(&fitTrack);
    }
    catch(genfit::Exception& e){
      std::cerr << e.what();
      std::cerr << "Exception, next track" << std::endl;
      continue;
    }

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


