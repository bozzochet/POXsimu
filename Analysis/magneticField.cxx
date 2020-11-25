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
#include <mySpacepointMeasurement.h>
#include <MeasurementFactory.h>
#include <EventDisplay.h>
#include <HelixTrackModel.h>
#include <MeasurementCreator.h>
#include <mySpacepointDetectorHit.h>

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
  
  //Declaration of leaves types
  Int_t           evID;
  Int_t           nHits;
  Int_t           nTotalHits;
  Int_t           ppHit;
  Int_t           hVol[19];
  Double_t        hVolZ[19];
  Double_t        xCoord[19];
  Double_t        yCoord[19];
  Double_t        zCoord[19];
  Double_t        eDep[19];
  Int_t           PDG[19];
  Int_t           TrID[19];
  Int_t           ParID[19];
  Double_t        xMom[19];
  Double_t        yMom[19];
  Double_t        zMom[19];
  Double_t        eEne[19];
  Int_t           chXY[19];
  Int_t           hitStrips[19];
  Int_t           simStrips[19];
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
  
  genfit::FieldManager::getInstance()->init(new field(0., 0., 0.5, magnet ));//0.5 kGauss = 0.05T   
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

  genfit::MeasurementFactory<genfit::mySpacepointMeasurement> *measFact = new genfit::MeasurementFactory<genfit::mySpacepointMeasurement>();

  // main loop (loops on the events)  
  for(Long64_t iEvent=0; iEvent<nentries; iEvent++) {
    nbytes += hitTree->GetEntry(iEvent);
    
    // declaration of the variables
    bool filled = false;
    std::vector<TVector3> pos, mom;
    unsigned int nMeasurements = 0;
    int index = 0;
    genfit::TrackCandHit* hit;
    
    TClonesArray *data = new TClonesArray("genfit::mySpacepointDetectorHit", 19);//size massima.

    
    // con la classe AbsMeasurementProducer non da errore, ma così non è stata inizializzata con un TClonesArray che contiene le posizioni, come le prende? se lo creo non con la classe astratta da errore il metodo per addarlo al MeasurementFactory, se non uso il factory, ma solo il metodo produce(int index,const TrackCandHit *hit) di MeasurementProducer viene creata la Measurement con la classe astratta, quindi da errore dopo nel passarla alla traccia che poi va fittata
    genfit::AbsMeasurementProducer<genfit::mySpacepointMeasurement> *prod;// =new genfit::MeasurementProducer<genfit::TrackCandHit,genfit::mySpacepointMeasurement>(data);

    int detID=5;
    TMatrixDSym covM(6);
    genfit::TrackCand trackHits;// = new genfit::TrackCand();
    TVector3 posTemp,momTemp;

    // approximate covariance
    double resolution = 0.01;
    for (int i = 0; i < 3; ++i)
      covM(i,i) = resolution*resolution;
    for (int i = 3; i < 6; ++i)
      covM(i,i) = TMath::Power(resolution / nMeasurements / sqrt(3), 2);

    // loop on the measurements for each event
    for(int s=0;s<19;s++){
      if(eDep[s]>0.01){// if the energy is deposited on the layer then the hit is detectable
	if(!filled){
	  index=s;
	  filled=true;
	}
	posTemp.SetX(xCoord[s]);
	posTemp.SetY(yCoord[s]);
	posTemp.SetZ(zCoord[s]);
	// non passando il momento in coordinate polari non da l'errore "TVector3 can't be stretched" ma dice comunque che il momento è zero
	
	momTemp.SetPhi(TMath::ATan2(yMom[s],xMom[s]));
	std::cout<<"phi: "<<momTemp.Phi()<<std::endl;
	momTemp.SetMag(TMath::Sqrt( TMath::Power(xMom[s],2) + TMath::Power(yMom[s],2) + TMath::Power(zMom[s],2)));
	momTemp.SetTheta(TMath::ACos(zMom[s]/momTemp.Mag()));
	pos.push_back(posTemp);
	mom.push_back(momTemp);
	// uncomment for further information on the particle
	      // TDatabasePDG::Instance()->GetParticle(PDG[s])->Dump();
	      // printf("%d %f\n", PDG[s], TDatabasePDG::Instance()->GetParticle(PDG[s])->Charge());
    
	// filling the data cluster to create the measurement	
	hit = new genfit::TrackCandHit(s,nMeasurements,evID,0.); // detID is s because is the index that loops on the layers
	                                                         // hitID is nMeasurements because is the counter on the detectable hits
	trackHits.addHit(hit);
	//	data->AddAt(new genfit::mySpacepointDetectorHit(posTemp,covM),s);
	measFact->addProducer(s,prod);	
	nMeasurements++;
      }
    }
    
    // trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep();//c'era il pdg ID può essere anche inizializzato senza  //stesso commento di sopra sulla carica...

    //    prod = new genfit::MeasurementProducer<genfit::TrackCandHit,genfit::mySpacepointMeasurement>(data);

    // smeared start state
    genfit::MeasuredStateOnPlane stateSmeared(rep);
    stateSmeared.setPosMomCov(pos[index], mom[index], covM);

    // create track
    TVectorD seedState(6);
    TMatrixDSym seedCov(6);
    stateSmeared.get6DStateCov(seedState, seedCov);
    genfit::Track fitTrack(rep,seedState,seedCov);//trackHits,measFact,rep);

    
    std::vector<genfit::mySpacepointMeasurement*> measurements = measFact->createMany(trackHits);
    
    try{
      for (unsigned int i=0; i<nMeasurements; ++i){
	fitTrack.insertPoint(new genfit::TrackPoint(measurements[i], &fitTrack));
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
    delete prod;
    delete hit;
    delete rep;
    delete data;
  }// end loop over events
  
  delete fitter;
  delete measFact;
  // open event display
  display->open();

}

