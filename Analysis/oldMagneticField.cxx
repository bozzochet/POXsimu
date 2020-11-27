#include <vector>

#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackCand.h>
#include <TrackPoint.h>
#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <PlanarMeasurement.h>
#include <MeasurementFactory.h>
#include <EventDisplay.h>
#include <HelixTrackModel.h>
#include <MeasurementCreator.h>
#include <mySpacepointDetectorHit.h>
#include <mySpacepointMeasurement.h>

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

  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

  TClonesArray data("genfit::mySpacepointDetectorHit", 500);//size massima
  //genfit::MeasurementFactory<genfit::mySpacepointMeasurement>* measFact = new genfit::MeasurementFactory<genfit::mySpacepointMeasurement>();
     genfit::MeasurementFactory<genfit::AbsMeasurement>* measFact = new genfit::MeasurementFactory<genfit::AbsMeasurement>();
  genfit::AbsMeasurementProducer<genfit::AbsMeasurement>* prod = new genfit::MeasurementProducer<genfit::mySpacepointDetectorHit, genfit::mySpacepointMeasurement>(&data);
  for (int ii=0; ii<19; ii++) {//ogni piano è un detector e quind un producer
    //measFact->addProducer(ii, (genfit::AbsMeasurementProducer<genfit::mySpacepointMeasurement>*)prod);
        measFact->addProducer(ii, prod);
  }
    
  // main loop (loops on the events)  
  for(Long64_t iEvent=0; iEvent<nentries; iEvent++) {
    nbytes += hitTree->GetEntry(iEvent);

    genfit::TrackCand trackHits;// = new genfit::TrackCand();
    
    // declaration of the variables
    unsigned int nMeasurements = 0;
    
    TMatrixDSym covM(3);
    TVector3 posTemp;
    
    // approximate covariance
    double resolution = 0.01;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
	if (i==j) covM(i, j) = resolution*resolution;
	else covM(i, j) = 0.0;
      }
    }
    //    std::cout<<"n: "<<nTotalHits<<std::endl;
    // loop on the measurements for each event
    for(int s=0;s<nTotalHits;s++){
      if(eDep[s]>0.01){// if the energy is deposited on the layer then the hit is detectable
	posTemp.SetX(xCoord[s]);
	posTemp.SetY(yCoord[s]);
	posTemp.SetZ(zCoord[s]);
	// uncomment for further information on the particle
	  // TDatabasePDG::Instance()->GetParticle(PDG[s])->Dump();
	  // printf("%d %f\n", PDG[s], TDatabasePDG::Instance()->GetParticle(PDG[s])->Charge());
	//	std::cout<<"iEvent: "<<iEvent<<"s: "<<s<<"layerID (da come c'è scritto su ):"<<hVol[s]<<std::endl;
	trackHits.addHit(hVol[s], s);
	new(data[s]) genfit::mySpacepointDetectorHit(posTemp, covM);
	nMeasurements++;
      }
    }
    //genfit::Track fitTrack(trackHits, (*((genfit::MeasurementFactory<genfit::AbsMeasurement> *)measFact)));
    genfit::Track fitTrack(trackHits, *measFact);
    trackHits.Print();
    fitTrack.Print();
    // checkConsistency is already called in the constructor of track
    // do the fit
    fitter->processTrackWithRep(&fitTrack,new genfit::RKTrackRep());
    
    //check
    fitTrack.checkConsistency();
    
    if (iEvent < 1000) {
      // add track to event display
      //display->addEvent(&fitTrack);
    }

    data.Delete();
  }// end loop over events
  
  // open event display
  display->open();

  delete fitter;
  delete measFact;
  delete prod;
  
}

