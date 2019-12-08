using namespace std;

int evID=-9999;

const int nHitsMax=50;
int nHits=0;
int nTotalHits=0;
Double_t xCoord[nHitsMax]={0.};
Double_t yCoord[nHitsMax]={0.};
Double_t zCoord[nHitsMax]={0.};
Double_t eDep[nHitsMax]={0.};

//Hough transform parameter
int binth=150; double minth=-TMath::Pi()/2, maxth=TMath::Pi()/2;  //binning of theta
int binr=100; double minr=-3, maxr=3;     //binning of r
double scartoth=2*(maxth-minth)/double(binth);
double scartor=2*(maxr-minr)/double(binr);

///////////////////////////////////////////////////////////////////

void Refine(vector<vector<double>> &hits, vector<double> &bdir, vector<vector<double>> &renegade, int iter=0){
  //the function recorsively calls itself studing the hits at different z values each time
  if(iter>45) return;   //exit condition for the loop 
  double pos=15-iter;   //current z position for the study
  //this vector collects the hits we keep for the next iteration
  vector<vector<double>> newhits;
  //variables for the search of the minimum distance
  double min=9999;
  int jmin=-1;

  //search for hits in the currently studied layer
  int s=hits.size();
  for(int j=0;j<s;j++){
   // cout<<j<<"/"<<s<<endl;
    if (hits[j][0]<pos && hits[j][0]>pos-1){    //layer belonging condition
      double xexp = bdir[1]/cos(bdir[0]) + tan(bdir[0])*hits[j][0];
      double d=TMath::Abs(hits[j][1]-xexp);

      if(d<min){      //between the coplanar hits we look for the nearest to the expected point
        if(jmin>-1) {renegade.push_back(hits[jmin]);}   //we collect as renegade all the hits at non-minimum distance
        min=d;
        jmin=j;
      } else {
        renegade.push_back(hits[j]);   //we collect as renegade all the hits at non-minimum distance
      }
    }
    else {
      newhits.push_back(hits[j]);    //we keep all the other hits
    }
    //cout<<newhits.size()<<endl;
  }
  if(jmin<0){   //if there weren't any hit in the current layer the iteration continue without any news
    Refine(hits,bdir,renegade,iter+1);
  }  else {
    newhits.push_back(hits[jmin]);     //the nearest is the only coplanar we keep

    //compute a new best dir as the average theta and r of the newhits
    int nnew=newhits.size();int conta=0;
    double sumr=0,sumth=0;
    double m,th,r;
    for(int i=0;i<nnew;i++){
      for (int j=i+1;j<nnew;j++){
        if(newhits[i][0]!=newhits[j][0]){
        	m=(newhits[i][1]-newhits[j][1])/(newhits[i][0]-newhits[j][0]);
	        th=atan(m);
	        r=cos(th)*newhits[i][1]-sin(th)*newhits[i][0];

          sumr+=r;
          sumth+=th;
          conta++;
        }
      }
    }
    bdir[0]=sumth/conta;
    bdir[1]=(sumr/conta);

    Refine(newhits,bdir,renegade,iter+1);
  }
}

///////////////////////////////////////////////////////////////////

void TrackFinding( vector<vector<double>> hcoord, vector<vector<double>> &dir, string name, bool vert, bool kDraw=true,int iter=1){

  int divisions= gStyle->GetNumberContours();
  gStyle->SetNumberContours(8);
  
  //static int iter=0;
  //iter++;

  //printf("%s, iter=%d\n", name.c_str(), iter);

  int n_search=2;               //number of track to search
  if (iter>n_search) {
   // iter--;
    return;
  }
  int nHits=hcoord.size();

  TCanvas* c2 = new TCanvas(Form("%s_%d", name.c_str(), iter), "Hough Transform");
  //this 2D histogram will contain the Hough transform
  TH2D* h = new TH2D(Form("h_%s_%d", name.c_str(), iter), "hough", binth, minth, maxth, binr, minr, maxr);

  //Generate the Hough transform
  double m, r, th; //pendenza, distanza dall'origine e angolo tra retta e asse z
  for (int i=0;i<nHits;i++){
    for (int j=i+1;j<nHits;j++){
      if(hcoord[i][0]!=hcoord[j][0]){
      	m=(hcoord[i][1]-hcoord[j][1])/(hcoord[i][0]-hcoord[j][0]);
	      th=atan(m);
	      r=cos(th)*hcoord[i][1]-sin(th)*hcoord[i][0];
	  //	printf("%f %f\n", th, r);
	      h->Fill(th,r);
      }
    }
  }

  //Get the position of the maximum
  int max_bin=h->GetMaximumBin();
  int xmax,ymax,zmax;
  h->GetBinXYZ(max_bin,xmax,ymax,zmax);

  //Add the best track (x0,theta) to the OUTPUT vector
  vector<double> best_dir(2);
  best_dir[0]=h->GetXaxis()->GetBinCenter(xmax); //angle
  best_dir[1]=h->GetYaxis()->GetBinCenter(ymax);  //distance

  if (kDraw) {
    //Plot the transform
    h->SetXTitle("angolo");
    h->SetYTitle("distanza");
    h->Draw("colz");
    //    h->DrawCopy("colz");
    TString filename=Form("%s_%d%s", name.c_str(), iter, ".png");
    c2->SaveAs(filename);
  }
  if (c2 && !kDraw) delete c2;
  if (h && !kDraw) delete h;

  gStyle->SetNumberContours(divisions);

  vector<vector<double>> newcoord;  //this vector will contain the hits that do not fit the best dir
  Refine(hcoord,best_dir,newcoord);
  dir.push_back(best_dir);
  if(vert) newcoord.push_back(hcoord[0]);
  TrackFinding(newcoord,dir,name,vert,kDraw,iter+1);

 // iter--;
  return;
}

//////////////////////////////////////////////////////////////

void PlotHits(vector<vector<Double_t>> hcoord,vector<vector<Double_t>> dir,vector<vector<Double_t>> dirf,string name,int ev, int evID, bool kDraw=true){
  
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
    f[i]->SetParameter(1,dir[i][1]/cos(dir[i][0]));
  }
  int n_foundf=dirf.size();
  //  printf("%s) n dir found %d\n", name.c_str(), n_found);
  TF1* ff[n_foundf];
  for(int i=0;i<n_foundf;i++){
    ff[i]=new TF1("f1","[0]*x+[1]",-30,30);
    ff[i]->SetParameter(0,tan(dirf[i][0]));
    ff[i]->SetParameter(1,dirf[i][1]/cos(dirf[i][0]));
    ff[i]->SetLineColor(3);
  }

  TGraph *g=new TGraph(n,z,x);
  TCanvas *c=new TCanvas(Form("2D_%s_%s", name.c_str(), "Hits2D"),"Hits2D");

  //plot the hits
  g->SetMarkerStyle(24);
  g->GetXaxis()->SetTitle("z");
  g->Draw("AP");
  for(int i=0;i<n_found;i++) f[i]->Draw("SAME");
  for(int i=0;i<n_foundf;i++) ff[i]->Draw("SAME");
  TString filename = Form("%s_%d_%d%s", name.c_str(), ev, evID, ".png");
  c->SaveAs(filename);

  if (!kDraw) delete c;
}

///////////////////////////////////////////////////////////////////////

void ordina(vector<vector<double>> &zx){
  int n=zx.size();
  double box;
  bool ordinato=0;
  while(!ordinato){
    ordinato=1;
    for(int i=0;i<n-1;i++){
      if(zx[i][0]>zx[i+1][0]){
        ordinato=0;
        box=zx[i][0];
        zx[i][0]=zx[i+1][0];
        zx[i+1][0]=box;

        box=zx[i][1];
        zx[i][1]=zx[i+1][1];
        zx[i+1][1]=box;
      }
    }
  }
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

void DAQ(int reqEv=-9999) {

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
  vector<vector<Double_t>> zx, zy;
  vector<vector<Double_t>> traccex, traccey;

  //This vector will contain just the first three hits from the conversion
  vector<vector<Double_t>> zxf, zyf;
  vector<vector<Double_t>> traccexf, tracceyf;

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

    zxf.clear();
    zyf.clear();
    
    traccexf.clear();
    tracceyf.clear();
    
    T->GetEntry(Ev);

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
    int n=zx.size();
    if (n<3) continue;

    ordina(zx);
    ordina(zy);

    int opt_n=3;
    bool vertice=zx[0][0]!=zx[1][0];
    if(!vertice) opt_n=6;

    for(int i=0;i<opt_n;i++){
	    zxf.push_back(zx[i]);
	    zyf.push_back(zy[i]);
    }

    TrackFinding(zx,traccex, "Immagini/houghZX",vertice, last);
    TrackFinding(zy,traccey, "Immagini/houghZY",vertice, last);

    TrackFinding(zxf,traccexf, "Immagini/houghZXf",vertice, last);
    TrackFinding(zyf,tracceyf, "Immagini/houghZYf",vertice, last);

    PlotHits(zx,traccex,traccexf,"Immagini/simuZX",Ev, evID, last);
    PlotHits(zy,traccey,tracceyf,"Immagini/simuZY",Ev, evID, last);
  }

  return;
}