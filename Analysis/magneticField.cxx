#include <vector>
#include <fstream>

#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitter.h>
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
#include <RectangularFinitePlane.h>

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
#include <TApplication.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TText.h>

#define _TC_(L) try { L; } catch(genfit::Exception& e){std::cerr<<"Exception, next track"<<std::endl; std::cerr<< e.what(); continue;}

#define NODISPLAY

FILE *muteOut(fopen("/dev/null", "w"));
std::ofstream fout("/dev/null");

FILE *stdoutSave = nullptr;
FILE *stderrSave = nullptr;
std::streambuf *cout_sbuf = nullptr;
std::streambuf *cerr_sbuf = nullptr;

bool isMuted = false;

void MuteOutput() {
  if (!std::exchange(isMuted, true)) {
    stdoutSave = stdout;
    stdout = muteOut;
    stderrSave = stderr;
    stderr = muteOut;

    // cout_sbuf = std::cout.rdbuf();
    // cerr_sbuf = std::cerr.rdbuf();
    // std::cout.rdbuf(fout.rdbuf());
    // std::cerr.rdbuf(fout.rdbuf());
    cout_sbuf = genfit::printOut.rdbuf();
    cerr_sbuf = genfit::errorOut.rdbuf();
    genfit::printOut.rdbuf(fout.rdbuf());
    genfit::errorOut.rdbuf(fout.rdbuf());
  }
}

void UnmuteOutput() {
  if (std::exchange(isMuted, false)) {
    if (stdoutSave) {
      stdout = stdoutSave;
      stdoutSave = nullptr;
    }
    if (stderrSave) {
      stderr = stderrSave;
      stderrSave = nullptr;
    }

    if (cout_sbuf) {
      //      std::cout.rdbuf(cout_sbuf);
      genfit::printOut.rdbuf(cout_sbuf);
      cout_sbuf = nullptr;
    }
    if (cerr_sbuf) {
      //      std::cerr.rdbuf(cerr_sbuf);
      genfit::errorOut.rdbuf(cerr_sbuf);
      cerr_sbuf = nullptr;
    }
  }
}

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


int main(int argc, char* argv[]){
  TApplication* app = new TApplication("app", &argc, argv);

  TCanvas *cchi = new TCanvas();
  TCanvas* cmom = new TCanvas();
  TCanvas* cinvmom = new TCanvas();
  TCanvas* creso = new TCanvas();
  TH1D* hchi = new TH1D("chi","chi",100,-1,18);
  TH1D* hinvmom = new TH1D("1.0/mom","invmom",1000,-1.0, 1.0);
  TH1D* hmom = new TH1D("mom","mom",10000,-100, 200);
  TH1D* hreso = new TH1D("reso","reso",1000,-10.0,10.0);//,-0.005,0.0005);
  
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

  // init geometry and mag. field
  TGeoManager::Import("plugins/libTestGeometry.vgm.root");
  TGeoVolume *magnet = gGeoManager->GetVolume("magnet");

  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  
  genfit::FieldManager::getInstance()->init(new field(0., 0.5, 0.0, magnet));//0.5 kGauss = 0.05T
  //  genfit::FieldManager::getInstance()->init(new field(0., 10., 0., magnet));//10.0 kGauss = 1T   

#ifndef NODISPLAY
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();
#endif
  
  //  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitter();
  
  // particle pdg code; muon hypothesis
  const int pdg = 13;  
      
  TMatrixDSym covM(2);
  // approximate covariance
  double resolution = TMath::Power(10,-5);//10 um;
  covM.UnitMatrix();
  covM *= resolution*resolution;

  TVectorD hitCoords(2);
  std::vector<genfit::PlanarMeasurement*> v_m;
  
  // main loop (loops on the events)  
  for(Long64_t iEvent=0; iEvent<nentries; iEvent++) {
    nbytes += hitTree->GetEntry(iEvent);
    //    printf("********+ %lld ********\n", iEvent);
    if(
       // iEvent==1941 ||
       // iEvent==3221 ||
       // iEvent==7440 ||
       // iEvent==2573 ||
       // iEvent==5964 ||
       // iEvent==4836 ||
       // iEvent==9491 ||
       // iEvent==6143 ||
       // iEvent==6763 ||
       // iEvent== 282 ||
       // iEvent==2438 ||
       // iEvent==3771 ||
       // iEvent==4157 ||
       // iEvent== 613 ||
       // iEvent==6326 ||
       // iEvent==6786 ||
       // iEvent==1593 ||
       // iEvent==2570 ||
       // iEvent== 867 ||
       // iEvent==1574 ||
       // iEvent==1912 ||
       // iEvent==3725 || //MD
       // iEvent== 612 || //MD
       // iEvent== 131 || //MD
       // iEvent==6857 || //MD
       // iEvent==8711 || //MD
       // iEvent==9168 || //MD
       nTotalHits<10 || //MD 
       0) continue;

    std::cout<<"event: "<<iEvent<<std::endl;
    
    // trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

    // start values for the fit, e.g. from pattern recognition
    //    TVector3 pos(0, 0, 0);
    TVector3 pos(xCoord[0], yCoord[0], zCoord[0]);

    //    TVector3 mom(xMom[0], yMom[0], zMom[0]);//ma deve essere giusto quello sotto
    TVector3 mom(xMom[0]/1000.0, yMom[0]/1000.0, zMom[0]/1000.0);
    //    mom.Print();

    // create track
    genfit::Track fitTrack;
    //    genfit::Track fitTrack(rep, pos, mom);
    
    // declaration of the variables
    TVector3 posTruth;
    // std::cout<<"n: "<<nTotalHits<<std::endl;
    // loop on the measurements for each event
    for(int s=0;s<nTotalHits;s++){
      if(eDep[s]>0.01){// if the energy is deposited on the layer then the hit is detectable
	posTruth.SetX(xCoord[s]);
	posTruth.SetY(yCoord[s]);
	posTruth.SetZ(zCoord[s]);
	//	posTruth.Print();
	//	posTruth.SetZ(hVolZ[s]);//MD: no? sembrano sempre uguali
	//	printf("%f %f\n", zCoord[s], hVolZ[s]);
	// uncomment for further information on the particle
	  // TDatabasePDG::Instance()->GetParticle(PDG[s])->Dump();
	  // printf("%d %f\n", PDG[s], TDatabasePDG::Instance()->GetParticle(PDG[s])->Charge());
	//	std::cout<<"iEvent: "<<iEvent<<"s: "<<s<<"layerID (da come c'è scritto su ):"<<hVol[s]<<std::endl;
	hitCoords[0] = posTruth.X();
	hitCoords[1] = posTruth.Y();

	//maybe here we have to create the planes and not only the planarmeasurement? Something like:
	// 	      genfit::AbsFinitePlane* recta = new RectangularFinitePlane( -4.047, 4.012, -3.3918, -1.47 );
	// genfit::SharedPlanePtr detectorplane (new genfit::DetPlane( origin_, TVector3(0,0,1), recta));

	//in the FOOT sw they do:
	// PlanarMeasurement* hit = new PlanarMeasurement(planarCoords, planarCov, m_detectorID_map["VT"], iHit, nullptr );
	// hit->setPlane(m_detectorPlanes[clus->GetPlaneNumber()], clus->GetPlaneNumber());
	// m_hitCollectionToFit_dataLike[ track_ID ].push_back( hit );
	// where
	// map <int, vector<AbsMeasurement*> > m_hitCollectionToFit_dataLike;

	//	for (unsigned int g = 0; g < m_hitCollectionToFit_dataLike[iTrack].size(); ++g){
	//	  fitTrack_->insertMeasurement( m_hitCollectionToFit_dataLike[iTrack].at(g) );
	//	}

	//SharedPlanePtr (that is a #typedef for a shared pointer to a DetPlane)
	// has
	// isInActive() to ask if a position is in the active area
	// the FOOT SW has also
	// isInActiveX() and isInActiveY() but I do not find in the GenFit Doxygen...
	
	genfit::PlanarMeasurement* measurement = NULL;
	_TC_(measurement = new genfit::PlanarMeasurement(hitCoords, covM, hVol[s], s, nullptr));
	_TC_(measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0, posTruth.Z()), TVector3(1,0,0), TVector3(0,1,0),nullptr)), s));
	genfit::TrackPoint* tp = new genfit::TrackPoint(measurement, &fitTrack);
	_TC_(tp->setSortingParameter(posTruth.Z()));
	_TC_(fitTrack.insertPoint(tp));
      }
    }
    bool sort_changed = false;
    sort_changed = fitTrack.sort();
    //    printf("%d\n", sort_changed);
    //    fitTrack.Print();
    _TC_(fitTrack.checkConsistency());
      
    // do the fit
    MuteOutput();
    _TC_(fitter->processTrack(&fitTrack));
    UnmuteOutput();
    
    // print fit result
    //    fitTrack.getFittedState().Print();
    //    fitTrack.checkConsistency();

    double chisq = fitTrack.getFitStatus()->getChi2();
    double ndf = fitTrack.getFitStatus()->getNdf();
    double chisq_red = chisq/ndf;
    //    printf("%f %f %d\n", chisq, ndf, nTotalHits);
    if (/*chisq>0.1 &&*/ ndf==33) {
      TVector3 posFitted, momFitted;
      TMatrixDSym covFitted;
      double rhoMomTruth = TMath::Sqrt(TMath::Power(xMom[0],2)+TMath::Power(yMom[0],2)+TMath::Power(zMom[0],2))/1000.0;
      //    momFitted =  TVector3( (fitTrack.getFittedState().get6DState())[3], (fitTrack.getFittedState().get6DState())[4], (fitTrack.getFittedState().get6DState())[5] );
      momFitted = fitTrack.getFittedState().getMom();
      //      momFitted.Print();
      //    fitTrack.getFittedState().getPosMomCov(posFitted, momFitted, covFitted);
      hchi->Fill(chisq_red);
      //      double rhoMomFitted = -1.0*fitTrack.getFittedState().getCharge()*TMath::Sqrt(TMath::Power(momFitted.X(),2)+TMath::Power(momFitted.Y(),2)+ TMath::Power(momFitted.Z(),2) );
      double rhoMomFitted = TMath::Sqrt(TMath::Power(momFitted.X(),2)+TMath::Power(momFitted.Y(),2)+ TMath::Power(momFitted.Z(),2) );
      double val= (1.0 - (rhoMomTruth/rhoMomFitted));
      hinvmom->Fill(1.0/rhoMomFitted);
      hmom->Fill(rhoMomFitted);
      hreso->Fill(val);
      //    std::cout<<"val "<<rhoMomTruth<<"      "<<rhoMomFitted<<std::endl;
      //    printf("val = %f, Ptruth = %f, Pfitted = %f, Charge = %f\n", val, rhoMomTruth, rhoMomFitted, fitTrack.getFitStatus()->getCharge());
    }
    
#ifndef NODISPLAY
    if (iEvent<1000){
      // add track to event display
      display->addEvent(&fitTrack);
    }
    else {
      break;
    }
#endif
    //    std::cout<<"Pval"<<fitter->getPVal(fitTrack, rep,-1)<<std::endl;
  }// end loop over events

#ifndef NODISPLAY
  // open event display
  display->open();
#endif

  cmom->cd();
  //  hMom3->Fit("gaus","","",-1000,0.005);
  hmom->Draw();
  
  cinvmom->cd();
  //  hMom4->Fit("gaus","","",-1000,0.005);
  hinvmom->Draw();
  
  creso->cd();
  hreso->Fit("gaus","","",-1000,0.005);
  hreso->Draw();
  TText T;
  T.SetTextFont(42);
  T.SetTextAlign(21);
  T.DrawTextNDC(.5,.95,"(1/momTruth - 1/momFitted)/(1/MomTruth)");
  
  cchi->cd();
  hchi->Draw();

  app->Run();
  
  delete fitter;
  delete cchi;
  delete cmom;
  delete cinvmom;
  delete creso;

  return 0;
}

