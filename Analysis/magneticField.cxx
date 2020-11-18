#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackCand.h>
#include <TrackPoint.h>
#include "ConstField.h"

#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <TGeoMaterialInterface.h>
#include <PlanarMeasurement.h>
#include <MeasurementFactory.h>
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
#include <TClonesArray.h>

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
  // gSystem->Load("libGeom");
  TGeoManager::Import("plugins/libTestGeometry.vgm.root");
  
  TGeoVolume *magnet = gGeoManager->GetVolume("magnet");
  // magnetic field on all the volumes, not ok for the simulation, only for testing
  
  genfit::FieldManager::getInstance()->init(new field(0., 0., 0.5, magnet ));//0.5 kGauss = 0.05T   
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

  genfit::MeasurementFactory<genfit::PlanarMeasurement> *measFact = new genfit::MeasurementFactory<genfit::PlanarMeasurement>();
  
  // main loop (loops on the events)  
  for(Long64_t iEvent=0; iEvent<nentries; iEvent++) {
    nbytes += hitTree->GetEntry(iEvent);

    
    // start values
    
    // declaration of the variables
    bool filled = false;
    int PDGid;
    TVector3 pos, mom;
    int nMeasurements = 0,index = 0;
    genfit::TrackCandHit* hit;
    //std::vector<genfit::eMeasurementType> measurementTypes;

    // what is it???
    TClonesArray *data = new TClonesArray();

    // loop on the measurements for each event
    for(int s=0;s<12;s++){
      if((xCoord[s] != 0)&&(yCoord[s] != 0)){ // hit is on both x and y
	if(!filled){
	  pos.SetX(xCoord[s]);
	  pos.SetY(yCoord[s]);
	  pos.SetZ(zCoord[s]);
	  mom.SetPhi(TMath::ATan2(xMom[s], yMom[s]));
	  mom.SetMag(TMath::Sqrt( TMath::Power(xMom[s],2) + TMath::Power(yMom[s],2) + TMath::Power(zMom[s],2)));
	  mom.SetTheta(TMath::ACos(zMom[s]/mom.Mag()));
	  index=s;
	  PDGid=PDG[s];
	  // uncomment for further information on the particle
	       // TDatabasePDG::Instance()->GetParticle(PDG[s])->Dump();
	       // printf("%d %f\n", PDG[s], TDatabasePDG::Instance()->GetParticle(PDG[s])->Charge());
    
	  filled=true;
	}
	// filling the data cluster to create the measurement
	/*hit = new genfit::TrackCandHit(evID,s,hVol[s],0.);
	genfit::MeasurementProducer<genfit::TrackCandHit,genfit::PlanarMeasurement> *prod = new genfit::MeasurementProducer<genfit::TrackCandHit,genfit::PlanarMeasurement>(data);
	measFact->addProducer(prod);
	*/
	nMeasurements++;
      }
    }

    // getting the charge of the particle
    const double charge = TDatabasePDG::Instance()->GetParticle(PDGid)->Charge()/(3.);
	  
    // helix track model    
    // get the charge of the particle
    genfit::HelixTrackModel* helix = new genfit::HelixTrackModel(pos, mom, charge);
    measurementCreator.setTrackModel(helix);
    
    // approximate covariance
    TMatrixDSym covM(6);
    double resolution = 0.01;
    for (int i = 0; i < 3; ++i)
      covM(i,i) = resolution*resolution;
    for (int i = 3; i < 6; ++i)
      covM(i,i) = TMath::Power(resolution / nMeasurements / sqrt(3), 2);

    // trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(PDGid);
      
    // smeared start state
    genfit::MeasuredStateOnPlane stateSmeared(rep);
    stateSmeared.setPosMomCov(pos, mom, covM);

    // create track
    TVectorD seedState(6);
    TMatrixDSym seedCov(6);
    stateSmeared.get6DStateCov(seedState, seedCov);
    genfit::Track fitTrack(rep, seedState, seedCov);

    genfit::TrackCand trackHits;// = new genfit::TrackCand();
    trackHits.addHit(hit);
    
    try{

      for (unsigned int i=0; i<nMeasurements; ++i){
      std::vector<genfit::PlanarMeasurement*> measurements = measFact->createMany(trackHits);
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

