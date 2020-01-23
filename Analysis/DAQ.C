#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "TRandom.h"
#include "TLorentzVector.h"

using namespace std;

struct Hit2D {
  double Z;
  double XorY;
  int PDG;
  int TrID;
  int ParID;
  TLorentzVector Mom;
  void operator=(Hit2D h){
    Z=h.Z;
    XorY=h.XorY;
    Mom=h.Mom;
    PDG=h.PDG;
    TrID=h.TrID;
    ParID=h.ParID;
    return;
  }
};

typedef vector<Hit2D> Hit2DColl;

struct Dir {
  double Th;
  double R;
  Hit2DColl voters;
};
typedef vector<Dir> DirColl;

int evID=-9999;

const int nHitsMax=100;
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

//Hough transform parameter
int binth=150; double minth=-TMath::Pi()/2, maxth=TMath::Pi()/2;  //binning of theta
int binr=100; double minr=-3, maxr=3;     //binning of r

///////////////////////////////////////////////////////////////////

void Refine(Hit2DColl &hits, Dir &bdir, Hit2DColl &renegade, int iter=0){
  
  //the function recorsively calls itself studing the hits at different z values each time
  if(iter>45) return;   //exit condition for the loop 
  double pos=15-iter;   //current z position for the study
  //this vector collects the hits we keep for the next iteration
  Hit2DColl newhits;
  //variables for the search of the minimum distance
  double min=9999;
  int jmin=-1;
  
  //search for hits in the currently studied layer
  int s=hits.size();
  for(int j=0;j<s;j++){
   // cout<<j<<"/"<<s<<endl;
    if (hits[j].Z<pos && hits[j].Z>pos-1){    //layer belonging condition
      double xexp = bdir.R/cos(bdir.Th) + tan(bdir.Th)*hits[j].Z;
      double d=TMath::Abs(hits[j].XorY-xexp);

      if(d<min){      //between the coplanar hits we look for the nearest to the expected point
        if(jmin>-1) { //we collect as renegade all the hits at non-minimum distance
	        renegade.push_back(hits[jmin]);
	      }   
        min=d;
        jmin=j;
      }
      else {   //we collect as renegade all the hits at non-minimum distance
        renegade.push_back(hits[j]);
      }
    }
    else {
      newhits.push_back(hits[j]);    //we keep all the other hits
    }
    //cout<<newhits.size()<<endl;
  }
  if(jmin<0){   //if there weren't any hit in the current layer the iteration continue without any news
    Refine(hits,bdir,renegade,iter+1);
  }
  else {
    newhits.push_back(hits[jmin]);     //the nearest is the only coplanar we keep
    bdir.voters.push_back(hits[jmin]);
    
    //compute a new best dir as the average theta and r of the newhits
    int nnew=newhits.size(); int conta=0;
    double sumr=0,sumth=0;
    double m,th,r;
    for(int i=0;i<nnew;i++){
      for (int j=i+1;j<nnew;j++){
        if(newhits[i].Z!=newhits[j].Z){
	        m=(newhits[i].XorY-newhits[j].XorY)/(newhits[i].Z-newhits[j].Z);
	        th=atan(m);
	        r=cos(th)*newhits[i].XorY-sin(th)*newhits[i].Z;
	  
          sumr+=r;
          sumth+=th;
          conta++;
        }
      }
    }
    bdir.Th=sumth/conta;
    bdir.R=(sumr/conta);
    
    Refine(newhits,bdir,renegade,iter+1);
  }

  return;
}

///////////////////////////////////////////////////////////////////

void TrackFinding(Hit2DColl hcoord, DirColl &dir, string directory, string name, bool vert, bool kDraw=true,int iter=1){
  
  int divisions= gStyle->GetNumberContours();
  gStyle->SetNumberContours(8);
  
  //static int iter=0;
  //iter++;
  printf("%s, iter=%d\n", name.c_str(), iter);
  
  int n_search=2;              //number of track to search
  if (iter>n_search) {
    printf("Two tracks already found...\n");
    //iter--;
    return;
  }
  int _nHits=hcoord.size();

  TCanvas* c2 = new TCanvas(Form("%s_%d", name.c_str(), iter), Form("%s_%d", name.c_str(), iter));
  //this 2D histogram will contain the Hough transform
  TH2D* h = new TH2D(Form("h_%s_%d", name.c_str(), iter), "Hough Trasform", binth, minth, maxth, binr, minr, maxr);
  h->GetXaxis()->SetTitle("#theta (rad)");
  h->GetYaxis()->SetTitle("r (cm)");

  //Generate the Hough transform
  double
    m,  //pendenza 
    r,  //distanza tra retta e origine
    th; //angolo tra retta e asse z
  for (int i=0;i<_nHits;i++){
    for (int j=i+1;j<_nHits;j++){
      //      if(hcoord[i].Z!=hcoord[j].Z){
      if(fabs(hcoord[i].Z-hcoord[j].Z)>1.0e-1){
      	m=(hcoord[i].XorY-hcoord[j].XorY)/(hcoord[i].Z-hcoord[j].Z);
	      double atanm=atan(m);
	      //double atan2m=atan2((hcoord[i].XorY-hcoord[j].XorY), (hcoord[i].Z-hcoord[j].Z));
	      th=atanm;
	      r=cos(th)*hcoord[i].XorY-sin(th)*hcoord[i].Z;
	      //	printf("atan=%f atan2=%f, m=%f r=%f\n", atanm, atan2m, m, r);
	      h->Fill(th,r);
      }
    }
  }

  //Get the position of the maximum
  int max_bin=h->GetMaximumBin();
  int xmax,ymax,zmax;
  h->GetBinXYZ(max_bin,xmax,ymax,zmax);

  //Add the best track (x0,theta) to the OUTPUT vector
  Dir best_dir;
  best_dir.Th = h->GetXaxis()->GetBinCenter(xmax); //angle
  best_dir.R =  h->GetYaxis()->GetBinCenter(ymax); //distance

  if (kDraw) {
    //Plot the transform
    h->Draw("colz");
    TString filename=Form("%s/%s_%d%s", directory.c_str(), name.c_str(), iter, ".png");
    c2->SaveAs(filename);
  }
  if (c2 && !kDraw) delete c2;
  if (h && !kDraw) delete h;
  
  gStyle->SetNumberContours(divisions);
  
  Hit2DColl newcoord;  //this vector will contain the hits that do not fit the best dir
  Refine(hcoord, best_dir, newcoord);
  dir.push_back(best_dir);
  if(vert) newcoord.push_back(hcoord[0]);

  TrackFinding(newcoord, dir, directory, name, vert, kDraw,iter+1);

  ///iter--;
  return;
}

//////////////////////////////////////////////////////////////

void Join(DirColl d, DirColl dv, Hit2DColl &h1, Hit2DColl &h2){
//  cout<<dv[0].voters.size()<<" vs "<<dv[1].voters.size()<<endl;

  Hit2D jv;
  jv = dv[0].voters.back();
  int j[2]={0,0};
  for(int i=0;i<2;i++){
    while(d[i].voters[j[i]].Z<=jv.Z){
      j[i]++;
    }
    j[i]--;
  }
  
  double min=9999;
  int imin=-1;int imax=-1;
  for (int i=0;i<2;i++){
    double dist=TMath::Abs(jv.XorY - d[i].voters[j[i]].XorY);
  //  cout<<"dist "<<dist<<endl;
    if(dist<min){
      imax=imin;
      min=dist;
      imin=i;
    }
  }
  imax=1-imin;
  //cout<<"imin"<<imin<<endl;
  h1 = dv[0].voters;
  int n0=dv[0].voters.size();
  int nmin=d[imin].voters.size();
  if(n0<nmin){
    for(int i=n0;i<nmin;i++){
        h1.push_back(d[imin].voters[i]);
    }
  }
  
  h2 = dv[1].voters; 
  int n1=dv[1].voters.size();
  int nmax=d[imax].voters.size();
  if(n1<nmax){
    for(int i=n1; i<nmax;i++){
        h2.push_back(d[imax].voters[i]);
    }
  }
}

/////////////////////////////////////////////////////////////

void PlotHits(Hit2DColl hcoord, DirColl dir, DirColl dirf, string directory, string name, int ev, int _evID, bool kDraw=true){
  int n=hcoord.size();
/*
  TGraph* eg = new TGraph();
  eg->SetMarkerColor(4);
  eg->SetMarkerStyle(24);

  TGraph* pg = new TGraph();
  pg->SetMarkerColor(2);
  pg->SetMarkerStyle(24);
*/
  TGraph* g = new TGraph();
  g->SetMarkerStyle(24);

  for(int i=0;i<n;i++){
    /*if(hcoord[i].ParID==1 && hcoord[i].PDG==11){
    eg->SetPoint(eg->GetN(),hcoord[i].Z, hcoord[i].XorY);
    }
    else if(hcoord[i].ParID==1 && hcoord[i].PDG==-11){
    pg->SetPoint(pg->GetN(),hcoord[i].Z, hcoord[i].XorY);
    }
    else {*/
    g->SetPoint(g->GetN(),hcoord[i].Z, hcoord[i].XorY);
    //}
  }

  TCanvas* c = new TCanvas(Form("Hits2D_%s", name.c_str()), Form("Hits2D_%s", name.c_str()));
/*  TMultiGraph* mg=new TMultiGraph;
  mg->Add(eg);
  mg->Add(pg);
  mg->Add(g);
*/
  int n_found=dir.size();
  //  printf("%s) n dir found %d\n", name.c_str(), n_found);
  TF1* f[n_found];
  for(int i=0;i<n_found;i++){
    f[i] = new TF1("f1", "[0]*x+[1]", -30, 30);
    f[i]->SetParameter(0,tan(dir[i].Th));
    f[i]->SetParameter(1,dir[i].R/cos(dir[i].Th));
    f[i]->SetLineColor(kBlue+2);
  }
  
  int n_foundf=dirf.size();
  //  printf("%s) n dir found %d\n", name.c_str(), n_found);
  TF1* ff[n_foundf];
  for(int i=0;i<n_foundf;i++){
    ff[i] = new TF1("f1", "[0]*x+[1]", -30, 30);
    ff[i]->SetParameter(0,tan(dirf[i].Th));
    ff[i]->SetParameter(1,dirf[i].R/cos(dirf[i].Th));
    ff[i]->SetLineColor(kGreen+2);
  }

  //plot the hits
  g->GetXaxis()->SetTitle("z");
  g->Draw("AP");
  for(int i=0;i<n_found;i++) f[i]->Draw("SAME");
  for(int i=0;i<n_foundf;i++) ff[i]->Draw("SAME");
  
  TString filename = Form("%s/%s_%d_%d.png", directory.c_str(), name.c_str(), ev, _evID);
  c->SaveAs(filename);
  
  if (!kDraw) {
    for(int i=0; i<n_found; i++){
      delete f[i];
      delete ff[i];
    }
    delete c;
  }

  return;
}

///////////////////////////////////////////////////////////////////////

void PlotHitsL(Hit2DColl hcoord,Hit2DColl hcoord1,Hit2DColl hcoord2,DirColl dir, DirColl dirf, string directory, string name, int ev, int _evID, bool kDraw=true){

  //hit da elettroni
  TGraph* eg = new TGraph();
  eg->SetMarkerColor(4);
  eg->SetMarkerStyle(20);
  //hit da positroni
  TGraph* pg = new TGraph();
  pg->SetMarkerColor(2);
  pg->SetMarkerStyle(20);
  //hit da particelle terziarie
  TGraph* g = new TGraph();
  g->SetMarkerStyle(24);
  g->SetLineWidth(0);

  int n=hcoord.size();
  for(int i=0;i<n;i++){
    if(hcoord[i].ParID==1 && hcoord[i].PDG==11){
      eg->SetPoint(eg->GetN(),hcoord[i].Z, hcoord[i].XorY);
    }
    else if(hcoord[i].ParID==1 && hcoord[i].PDG==-11){
    pg->SetPoint(pg->GetN(),hcoord[i].Z, hcoord[i].XorY);
    }
    else {
    g->SetPoint(g->GetN(),hcoord[i].Z, hcoord[i].XorY);
    }
  }
  //set di hits 1
  TGraph* g1 = new TGraph();
  g1->SetLineColor(4);
  g1->SetLineWidth(2);
  g1->SetLineStyle(10);
  g1->SetMarkerStyle(1);
  //set di hits 2
  TGraph* g2 = new TGraph();
  g2->SetLineColor(2);
  g2->SetLineWidth(2);
  g2->SetLineStyle(10);
  g2->SetMarkerStyle(1);
 
  int n1=hcoord1.size();
  for(int i=0;i<n1;i++){
    g1->SetPoint(g1->GetN(),hcoord1[i].Z, hcoord1[i].XorY);
  }
  int n2=hcoord2.size();
  for(int i=0;i<n2;i++){
    g2->SetPoint(g2->GetN(),hcoord2[i].Z, hcoord2[i].XorY);
  }  

  TCanvas* c = new TCanvas(Form("Hits2D_%s", name.c_str()), Form("Hits2D_%s", name.c_str()));
  TMultiGraph* mg=new TMultiGraph();
  mg->Add(eg);
  mg->Add(pg);
  mg->Add(g1);
  mg->Add(g2);
  if(g->GetN()>0)
    mg->Add(g);

  int n_foundf=dirf.size();
  //  printf("%s) n dir found %d\n", name.c_str(), n_found);
  TF1* ff[n_foundf];
  for(int i=0;i<n_foundf;i++){
    ff[i] = new TF1("f1", "[0]*x+[1]", -30, 30);
    ff[i]->SetParameter(0,tan(dirf[i].Th));
    ff[i]->SetParameter(1,dirf[i].R/cos(dirf[i].Th));
    ff[i]->SetLineStyle(3);
    ff[i]->SetLineColor(kGreen+2);
  }

  //plot the hits
  mg->GetXaxis()->SetTitle("z");
  mg->Draw("APL");
  for(int i=0;i<n_foundf;i++) ff[i]->Draw("SAME");
  
  TString filename = Form("%s/%s_%d_%d.png", directory.c_str(), name.c_str(), ev, _evID);
  c->SaveAs(filename);
  
  if (!kDraw) {
    for(int i=0; i<n_foundf; i++){
      delete ff[i];
    }
    delete c;
  }

  return;
}

///////////////////////////////////////////////////////////////////////

double CountCorrect(Hit2DColl h){
  double contae=0,contap=0;
  double n=h.size();
  for (int i=0;i<n;i++){
    if(h[i].ParID==1 && h[i].PDG==11) contae++;
    else if(h[i].ParID==1 && h[i].PDG==-11) contap++;
  }
  if (contae > contap) return contae/n;
  else return contap/n;
}

bool Tracked(Hit2DColl h){
  if(h[0].PDG==h.back().PDG && h[0].ParID==h.back().ParID) return 1;
  else return 0;
}

/////////////////////////////////////////////////////////////

void ordina(Hit2DColl &zx){
  int n=zx.size();
  
  Hit2D box;
  bool ordinato=0;
  
  while(!ordinato){
    ordinato=1;
    for(int i=0;i<n-1;i++){
      if(zx[i].Z>zx[i+1].Z){
        ordinato=0;
        box=zx[i];
        zx[i]=zx[i+1];
        zx[i+1]=box;
      }
    }
  }
}

/////////////////////////////////////////////////////////////

void raddrizza(DirColl &d){
  for(int i=0;i<2;i++){
    double m=(d[i].voters[1].XorY-d[i].voters[0].XorY)/(d[i].voters[1].Z-d[i].voters[0].Z);
    d[i].Th = atan(m);
    d[i].R = d[i].voters[0].XorY*cos(d[i].Th)-d[i].voters[0].Z*sin(d[i].Th);

    int n=d[i].voters.size();
    if(n>2){
      for(int j=2;j<n;j++){
        d[i].voters.pop_back();
      }
    }
  }
}

/////////////////////////////////////////////////////////////

void clear(){
  nHits=0;
  nTotalHits=0;
  fill_n(xCoord, nHitsMax, -9999.);
  fill_n(yCoord, nHitsMax, -9999.);
  fill_n(zCoord, nHitsMax, -9999.);
  fill_n(xMom, nHitsMax, -9999.);
  fill_n(yMom, nHitsMax, -9999.);
  fill_n(zMom, nHitsMax, -9999.);
  fill_n(eEne, nHitsMax, -9999.);
  fill_n(eDep, nHitsMax, -9999.);
  fill_n(PDG, nHitsMax, 0);
  fill_n(TrID, nHitsMax, -1);
  fill_n(ParID, nHitsMax, -1);
}

//////////////////////////////////////////////////////////////

void DAQ(int reqEv=-9999) {

  //system("mkdir -p ./Immagini/");
  //system("mkdir -p ./Immagini2/");

  TRandom *rnd=new TRandom();
    
  TString inputFileName="anaOut2.root";
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

  int nEvs=T->GetEntries();
  cout<<"Si sono verificate "<<nEvs<<" conversioni"<<endl;

  //This vectors of vectors will contain the hits' x and y coordinates
  Hit2DColl zx, zy;
  DirColl TrackX, TrackY;

  //This vector will contain just the first three hits from the conversion
  Hit2DColl zxV, zyV;
  DirColl TrackVX, TrackVY;

  auto* PSF=new TH2D("","",80,-.03,.03, 80,-.03,.03);
  auto* Correct=new TH1D("","",7,0.3,1);
  double momx=0,momy=0;

  int contacontati=0;
  int contatotali=0;

  int contatracciati=0;

  bool last=false;

  //The for loop run through every event in the tree
  int startEv=0;
  int stopEv=nEvs-1;
  if (reqEv>=0) {
    startEv=reqEv;
    stopEv=reqEv;
  }
  for (int Ev=startEv; Ev<=stopEv; Ev++){
    cout<<"siamo al "<<Ev<<endl;
    if (Ev==stopEv) last=true;
    contatotali++;

    clear();
    zx.clear();
    zy.clear();
    
    TrackX.clear();
    TrackY.clear();

    zxV.clear();
    zyV.clear();
    
    TrackVX.clear();
    TrackVY.clear();
    
    T->GetEntry(Ev);

    //this loop fill the vectors of vectors from the tree branches
    for (int i=0;i<nTotalHits;i++){
      if(eDep[i]>0.01 && zCoord[i]<30){    //we're taking into consideration only the detectable hits in the tracker
      	Hit2D hitx;
      	hitx.Z=zCoord[i];
      	hitx.XorY=xCoord[i];
        hitx.Mom.SetXYZT(xMom[i],yMom[i],zMom[i],eEne[i]);
        hitx.PDG=PDG[i];
        hitx.TrID=TrID[i];
        hitx.ParID=ParID[i];
	      zx.push_back(hitx);

	      Hit2D hity;
      	hity.Z=zCoord[i];
	      hity.XorY=yCoord[i];
        hity.Mom.SetXYZT(xMom[i],yMom[i],zMom[i],eEne[i]);
        hity.PDG=PDG[i];
        hity.TrID=TrID[i];
        hity.ParID=ParID[i];
	      zy.push_back(hity);
      }
    }

    ordina(zx);
    ordina(zy);

    int opt_n=3;
    bool vertice=zx[0].Z!=zx[1].Z;
    if (!vertice) opt_n=6;
    //    printf("%d\n", opt_n);

    int n=zx.size();
    if (n<opt_n) continue;
    contacontati++;

    for(int i=0;i<opt_n;i++){
      zxV.push_back(zx[i]);
      zyV.push_back(zy[i]);
    }

    TrackFinding(zx, TrackX, "Immagini", "houghZX", vertice, last);
    TrackFinding(zy, TrackY, "Immagini", "houghZY", vertice, last);
    
    TrackFinding(zxV, TrackVX, "Immagini", "houghZXf", vertice, 0);
    TrackFinding(zyV, TrackVY, "Immagini", "houghZYf", vertice, 0);



    //These vectors will contain the 2 subset of hits corresponding to each particle according to the join function
    Hit2DColl SSetX1,SSetX2,SSetY1,SSetY2;
    for(int i=0;i<2;i++){
      ordina(TrackX[i].voters);
      ordina(TrackVX[i].voters);
      ordina(TrackY[i].voters);
      ordina(TrackVY[i].voters);
    }

    raddrizza(TrackVX);
    raddrizza(TrackVY);

    PlotHits(zx,TrackX, TrackVX, "Immagini", "simuZX", Ev, evID, last);
    PlotHits(zy,TrackY, TrackVY, "Immagini", "simuZY", Ev, evID, last);

    Join(TrackX, TrackVX, SSetX1, SSetX2);
    Join(TrackY, TrackVY, SSetY1, SSetY2);

    if(reqEv>=0){
      //PlotHitsL(zx,SSetX1,SSetX2,TrackX, TrackVX, "Immagini2", "simuZX", Ev, evID, last);
      //PlotHitsL(zy,SSetY1,SSetY2,TrackY, TrackVY, "Immagini2", "simuZY", Ev, evID, last);
    } 
    else {
      double PX1=SSetX1.back().Mom[3];
      double PX2=SSetX2.back().Mom[3]; 
      double PY1=SSetY1.back().Mom[3];
      double PY2=SSetY2.back().Mom[3];

      double PX1r= 1/(rnd -> Gaus(1/PX1, 1/(.3*PX1)));
      double PX2r= 1/(rnd -> Gaus(1/PX2, 1/(.3*PX2)));
      double PY1r= 1/(rnd -> Gaus(1/PY1, 1/(.3*PY1)));
      double PY2r= 1/(rnd -> Gaus(1/PY2, 1/(.3*PY2)));


      momx= (( 1/(rnd -> Gaus(1/(SSetX1.back().Mom[0]), 1/(0.3*SSetX1.back().Mom[0]))) )*sin(TrackVX[0].Th) +
            ( 1/(rnd -> Gaus(1/(SSetX2.back().Mom[0]), 1/(0.3*SSetX2.back().Mom[0]))) )*sin(TrackVX[1].Th))/(PX1r+PX2r);

      momy= (( 1/(rnd -> Gaus(1/(SSetY1.back().Mom[0]), 1/(0.3*SSetY1.back().Mom[0]))) )*sin(TrackVY[0].Th) +
            ( 1/(rnd -> Gaus(1/(SSetY2.back().Mom[0]), 1/(0.3*SSetY2.back().Mom[0]))) )*sin(TrackVY[1].Th))/(PY1r+PY2r);
      cout<<momx<<", "<<momy<<endl;
      PSF->Fill(momx,momy);

      Correct->Fill(CountCorrect(SSetX1));
      Correct->Fill(CountCorrect(SSetX2));
      Correct->Fill(CountCorrect(SSetY1));
      Correct->Fill(CountCorrect(SSetY2));

      if(Tracked(SSetX1)) contatracciati++;
      if(Tracked(SSetX2)) contatracciati++;
      if(Tracked(SSetY1)) contatracciati++;
      if(Tracked(SSetY2)) contatracciati++;
    }
  }
  if(reqEv<0){
    auto canv1=new TCanvas();
    PSF->Draw("colz");
    canv1->SaveAs("psf.png");
    canv1->SaveAs("psf.c");

    auto canv2=new TCanvas();
    Correct->Draw();
    canv2->SaveAs("corrette.png");
    canv2->SaveAs("corrette.c");
    cout<<contacontati<<" / "<<contatotali<<" eventi con almeno opt_n tracce"<<endl;
    cout<<contatracciati<<" / "<<4*contacontati<<" eventi tracciati correttamente allo spettrometro"<<endl;
  }
  delete rnd;
  return;
}

int main(int argc, char* argv[]) {
  TApplication *myapp=new TApplication("myapp",0,0);

 // int par=atoi(argv[1]);
  //cout<<par<<endl;
  DAQ(136);

  myapp->Run();

  return 0;
}