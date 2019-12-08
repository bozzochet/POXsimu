using namespace std;

struct Hit2D {
  double Z;
  double XorY;
};
typedef vector<Hit2D> Hit2DColl;

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

///////////////////////////////////////////////////////////////////

void Refine(Hit2DColl &hits, Hit2D &bdir, Hit2DColl &renegade, int iter=0){
  
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
      double xexp = bdir.Z/cos(bdir.Z) + tan(bdir.Z)*hits[j].Z;
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
    bdir.Z=sumth/conta;
    bdir.XorY=(sumr/conta);
    
    Refine(newhits,bdir,renegade,iter+1);
  }

  return;
}

///////////////////////////////////////////////////////////////////

void TrackFinding(Hit2DColl hcoord, Hit2DColl &dir, string directory, string name, bool vert, bool kDraw=true){
  
  int divisions= gStyle->GetNumberContours();
  gStyle->SetNumberContours(8);
  
  static int iter=0;
  iter++;
  printf("%s, iter=%d\n", name.c_str(), iter);
  
  int n_search=2;               //number of track to search
  if (iter>n_search) {
    printf("Two tracks already found...\n");
    iter--;
    return;
  }
  int _nHits=hcoord.size();

  TCanvas* c2 = new TCanvas(Form("%s_%d", name.c_str(), iter), Form("%s_%d", name.c_str(), iter));
  //this 2D histogram will contain the Hough transform
  TH2D* h = new TH2D(Form("h_%s_%d", name.c_str(), iter), "HoughTrasform:#theta (rad):r (cm)", binth, minth, maxth, binr, minr, maxr);

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
	double atan2m=atan2((hcoord[i].XorY-hcoord[j].XorY), (hcoord[i].Z-hcoord[j].Z));
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
  Hit2D best_dir;
  best_dir.Z=h->GetXaxis()->GetBinCenter(xmax); //angle
  best_dir.XorY=h->GetYaxis()->GetBinCenter(ymax);  //distance

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
  Refine(hcoord,best_dir,newcoord);
  dir.push_back(best_dir);
  if(vert) newcoord.push_back(hcoord[0]);

  TrackFinding(newcoord, dir, directory, name, vert, kDraw);

  iter--;
  return;
}

//////////////////////////////////////////////////////////////

void PlotHits(Hit2DColl hcoord, Hit2DColl dir, Hit2DColl dirf, string directory, string name, int ev, int _evID, bool kDraw=true){
  
  int n=hcoord.size();

  //create array of x and array of z for the plot
  double* x;
  double* z;
  x = new double[n];
  z = new double[n];
  for(int i=0;i<n;i++){
    z[i]=hcoord[i].Z;
    x[i]=hcoord[i].XorY;
  }

  int n_found=dir.size();
  //  printf("%s) n dir found %d\n", name.c_str(), n_found);
  TF1* f[n_found];
  for(int i=0;i<n_found;i++){
    f[i] = new TF1("f1", "[0]*x+[1]", -30, 30);
    f[i]->SetParameter(0,tan(dir[i].Z));
    f[i]->SetParameter(1,dir[i].XorY/cos(dir[i].Z));
    f[i]->SetLineColor(kBlue+2);
  }
  
  int n_foundf=dirf.size();
  //  printf("%s) n dir found %d\n", name.c_str(), n_found);
  TF1* ff[n_foundf];
  for(int i=0;i<n_foundf;i++){
    ff[i] = new TF1("f1", "[0]*x+[1]", -30, 30);
    ff[i]->SetParameter(0,tan(dirf[i].Z));
    ff[i]->SetParameter(1,dirf[i].XorY/cos(dirf[i].Z));
    ff[i]->SetLineColor(kGreen+2);
  }

  TGraph* g = new TGraph(n,z,x);
  TCanvas* c = new TCanvas(Form("Hits2D_%s", name.c_str()), Form("Hits2D_%s", name.c_str()));

  //plot the hits
  g->SetMarkerStyle(24);
  g->GetXaxis()->SetTitle("z");
  g->Draw("AP");
  for(int i=0;i<n_found;i++) f[i]->Draw("SAME");
  for(int i=0;i<n_foundf;i++) ff[i]->Draw("SAME");
  
  TString filename = Form("%s/%s_%d_%d.png", directory.c_str(), name.c_str(), ev, _evID);
  c->SaveAs(filename);
  
  if (!kDraw) {
    for(int i=0; i<n_found; i++){
      delete f[i];
    }
    delete c;
  }

  return;
}

///////////////////////////////////////////////////////////////////////

void ordina(Hit2DColl &zx){
  
  int n=zx.size();
  
  double box;
  bool ordinato=0;
  
  while(!ordinato){
    ordinato=1;//come cazzo fara' a funzionare...
    for(int i=0;i<n-1;i++){
      if(zx[i].Z>zx[i+1].Z){
        ordinato=0;
        box=zx[i].Z;
        zx[i].Z=zx[i+1].Z;
        zx[i+1].Z=box;
	
        box=zx[i].XorY;
        zx[i].XorY=zx[i+1].XorY;
        zx[i+1].XorY=box;
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

  T->SetBranchAddress( "evID", &evID );
  T->SetBranchAddress( "nHits", &nHits );
  T->SetBranchAddress( "nTotalHits", &nTotalHits );
  T->SetBranchAddress( "xCoord", &xCoord );
  T->SetBranchAddress( "yCoord", &yCoord );
  T->SetBranchAddress( "zCoord", &zCoord );
  T->SetBranchAddress( "eDep", &eDep );

  int nEvs=T->GetEntries();
  cout<<"Si sono verificate "<<nEvs<<" conversioni"<<endl;

  //This vectors of vectors will contain the hits' x and y coordinates
  Hit2DColl zx, zy;
  Hit2DColl traccex, traccey;

  //This vector will contain just the first three hits from the conversion
  Hit2DColl zxf, zyf;
  Hit2DColl traccexf, tracceyf;

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
      	Hit2D hitx;
      	hitx.Z=zCoord[i];
	hitx.XorY=xCoord[i];
	zx.push_back(hitx);
	
	Hit2D hity;
	hity.Z=zCoord[i];
	hity.XorY=yCoord[i];
	zy.push_back(hity);
      }
    }

    int n=zx.size();
    if (n<3) continue;

    ordina(zx);
    ordina(zy);

    int opt_n=3;
    bool vertice=zx[0].Z!=zx[1].Z;
    if (!vertice) opt_n=6;
    //    printf("%d\n", opt_n);
    
    for(int i=0;i<opt_n;i++){
      zxf.push_back(zx[i]);
      zyf.push_back(zy[i]);
    }

    TrackFinding(zx, traccex, "Immagini", "houghZX", vertice, last);
    TrackFinding(zy, traccey, "Immagini", "houghZY", vertice, last);
    
    TrackFinding(zxf, traccexf, "Immagini", "houghZXf", vertice, last);
    TrackFinding(zyf, tracceyf, "Immagini", "houghZYf", vertice, last);
    
    PlotHits(zx,traccex, traccexf, "Immagini", "simuZX", Ev, evID, last);
    PlotHits(zy,traccey, tracceyf, "Immagini", "simuZY", Ev, evID, last);
  }

  return;
}
