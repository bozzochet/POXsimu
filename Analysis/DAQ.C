using namespace std;

int evID=-9999;

const int nHitsMax=50;
int nHits=0;
int nTotalHits=0;
Double_t xCoord[nHitsMax]={0.};
Double_t yCoord[nHitsMax]={0.};
Double_t zCoord[nHitsMax]={0.};
Double_t eDep[nHitsMax]={0.};

///////////////////////////////////////////////////////////////////

void TrackFinding( vector<vector<double>> hcoord, vector<vector<double>> &dir, string name, bool kDraw=true){

  int divisions= gStyle->GetNumberContours();
  gStyle->SetNumberContours(8);
  
  static int iter=0;
  iter++;

  printf("%s, iter=%d\n", name.c_str(), iter);
  
  int n_search=2;               //number of track to search
  if (iter>n_search) return;
  int nHits=hcoord.size();

  TCanvas* c2 = new TCanvas(Form("%s_%d", name.c_str(), iter), Form("%s %d", name.c_str(), iter));
  
  //this 2D histogram will contain the Hough transform
  //  int binth=150; double minth=-TMath::Pi(), maxth=TMath::Pi();  //binning of theta
  int binth=100; double minth=-2, maxth=2;  //binning of theta
  int binr=100; double minr=-3, maxr=3;     //binning of r
  //double scartoth=0;//double(maxth-minth)/double(2*binth);
  //double scartor=0;//double(maxr-minr)/double(2*binr);
  TH2D* h = new TH2D(Form("h_%s_%d", name.c_str(), iter), "hough", binth, minth, maxth, binr, minr, maxr);

  //Generate the Hough transform
  double
    m,  //la pendenza 
    r,  //la distanza tra retta e origine
    th; //l'angolo tra r e l'asse z
  for (int i=0;i<nHits;i++){
    for (int j=i+1;j<nHits;j++){
      if(hcoord[i][0]!=hcoord[j][0]){
	m=(hcoord[i][1]-hcoord[j][1])/(hcoord[i][0]-hcoord[j][0]);
	th=atan(m);
	//	th=atan2((hcoord[i][1]-hcoord[j][1]), (hcoord[i][0]-hcoord[j][0]));
	r=cos(th)*hcoord[i][1]-sin(th)*hcoord[i][0];
	//	printf("%f %f\n", th, r);
	h->Fill(th,r);
      }
    }
  }

  //Search for peaks in the Hough transform
  TSpectrum2* s=new TSpectrum2();
  int n_found=s->Search(h, 1, "nobackground", .5);
  //  s->Print();
  double* thpeaks=s->GetPositionX();
  double* rpeaks=s->GetPositionY();

  if (kDraw) {
    //Plot the transform
    h->SetXTitle("angolo");
    h->SetYTitle("distanza");
    //    h->Draw("lego2");
    h->Draw("colz");
    TString filename=Form("%s_%d%s", name.c_str(), iter, ".png");
    c2->SaveAs(filename);
  }
  if (c2 && !kDraw) delete c2;
  if (h && !kDraw) delete h;

  gStyle->SetNumberContours(divisions);
  
  cout<<s->GetNPeaks()<<endl;
  if(s->GetNPeaks()==0) {
    iter--;
    return;
  }
  // if (kDraw) {
  //   TPolyMarker* pm = new TPolyMarker(s->GetNPeaks(), thpeaks, rpeaks);
  //   pm->Draw("same");
  // }

  //Add the best track (x0,theta) to the OUTPUT vector
  vector<double> best_dir(2);
  best_dir[0]=thpeaks[0]; //angle
  best_dir[1]= rpeaks[0] / cos(best_dir[0]);  //q
  dir.push_back(best_dir);
  
  //remove fitted points
  vector<vector<double>> newcoord;
  int z0=10000;int i0=0;double sum=0;double conta=0;
  for (int i=0;i<nHits;i++){
    double xexp=best_dir[1]+tan(best_dir[0])*hcoord[i][0];
    double chi=pow(hcoord[i][1]-xexp,2);
    if(chi>1*abs(xexp)){
      newcoord.push_back(hcoord[i]);
    }
    else {
      sum=sum+chi;conta++;
    }
    if(hcoord[i][0]<z0){
      z0=hcoord[i][1];
      i0=i;
    }
  }
  newcoord.push_back(hcoord[i0]);

  //Add chi squared to the output vector
  //chisq.push_back(sum/conta);
  
  TrackFinding(newcoord,dir,name, kDraw);

  iter--;
  return;
}

//////////////////////////////////////////////////////////////

void PlotHits(vector<vector<Double_t>> hcoord,vector<vector<Double_t>> dir,string name,int ev, int evID, bool kDraw=true){
  
  int n=hcoord.size();

  //create array of x and array of z for the plot
  Double_t x[n],z[n];
  for(int i=0;i<n;i++){
    z[i]=hcoord[i][0];
    x[i]=hcoord[i][1];
  }

  int n_found=dir.size();
  //  printf("%s) n dir found %d\n", name.c_str(), n_found);
  TF1* f[n_found];
  for(int i=0;i<n_found;i++){
    f[i]=new TF1("f1","[0]*x+[1]",-30,30);
    f[i]->SetParameter(0,tan(dir[i][0]));
    f[i]->SetParameter(1,dir[i][1]);
  }

  TGraph *g=new TGraph(n,z,x);
  TCanvas *c=new TCanvas(Form("2D_%s_%s", name.c_str(), "Hits2D"),"Hits2D");

  //plot the hits
  g->SetMarkerStyle(3);
  g->Draw("AP");
  for(int i=0;i<n_found;i++)
    f[i]->Draw("SAME");
  TString filename = Form("%s_%d_%d%s", name.c_str(), ev, evID, ".png");
  c->SaveAs(filename);

  if (!kDraw) delete c;
}

///////////////////////////////////////////////////////////////////////

void PlotHits3D(Double_t x[],Double_t y[],Double_t z[],int n,string name,int ev, int evID, bool kDraw=true){

  TGraph2D *g=new TGraph2D(nTotalHits,x,y,z);
  TCanvas *c=new TCanvas(Form("3D_%s_%s", name.c_str(), "Hits3D"),"Hits3D");

  //plot the hits
  g->SetMarkerStyle(3);
  g->Draw("AP");
  TString filename = Form("%s_%d_%d%s", name.c_str(), ev, evID, ".png");
  c->SaveAs(filename);
  
  if (!kDraw) delete c;
}

/////////////////////////////////////////////////////////////

void clear(){
  nHits=0;
  nTotalHits=0;
  fill_n(xCoord, nHitsMax, 0.);
  fill_n(yCoord, nHitsMax, 0.);
  fill_n(zCoord, nHitsMax, 0.);
  fill_n(eDep, nHitsMax, 0.);
}

//////////////////////////////////////////////////////////////

void DAQ(int reqEv=-9999){

  system("mkdir -p ./Immagini/");
  
  TString inputFileName="anaOut.root";

  TFile *inFile=new TFile(inputFileName,"READ");

  TTree *T=(TTree*)inFile->Get("hitTree");

  T->SetBranchAddress( "evID",&evID );
  T->SetBranchAddress( "nHits",&nHits );
  T->SetBranchAddress( "nTotalHits",&nTotalHits );
  T->SetBranchAddress( "xCoord",&xCoord );
  T->SetBranchAddress( "yCoord",&yCoord );
  T->SetBranchAddress( "zCoord",&zCoord );
  T->SetBranchAddress( "eDep",&eDep );

  int nEvs=T->GetEntries();
  cout<<"Si sono verificate "<<nEvs<<" conversioni"<<endl;

  //This vectors of vectors will contain the hits' x and y coordinates
  vector<vector<Double_t>> zx;
  vector<vector<Double_t>> zy;

  vector<vector<Double_t>> traccex;
  vector<vector<Double_t>> traccey;

  bool last=false;

  //The for loop run through every event in the tree
  int startEv=0;
  int stopEv=nEvs-1;
  if (reqEv>=0) {
    startEv=reqEv;
    stopEv=reqEv;
  }
  for (int Ev=startEv; Ev<=stopEv; Ev++){
    
    if (Ev==stopEv) last=true;
    
    clear();
    zx.clear();
    zy.clear();
    
    traccex.clear();
    traccey.clear();
    
    T->GetEntry(Ev);

    if (reqEv>=0) PlotHits3D(zCoord,xCoord,yCoord,nTotalHits,"Immagini/simu",Ev, evID, last);
    
    //this loop fill the vectors of vectors from the tree branches
    for (int i=0;i<nTotalHits;i++){
      if(eDep[i]>.01 && zCoord[i]<30){    //we're taking into consideration only the detectable hits in the tracker
	vector<Double_t> hitx(2);
	hitx[0]=zCoord[i];
	hitx[1]=xCoord[i];
	zx.push_back(hitx);
	
	vector<Double_t> hity(2);
	hity[0]=zCoord[i];
	hity[1]=yCoord[i];
	zy.push_back(hity);
      }
    }
    if (zx.size()<3 && zy.size()<3) continue;
    TrackFinding(zx,traccex, "Immagini/houghZX", last);
    TrackFinding(zy,traccey, "Immagini/houghZY", last);
    PlotHits(zx,traccex,"Immagini/simuZX",Ev, evID, last);
    PlotHits(zy,traccey,"Immagini/simuZY",Ev, evID, last);
  }

  return;
}
