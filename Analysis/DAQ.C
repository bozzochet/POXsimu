using namespace std;


    //devo definire una serie di variabili su cui rovesciare le info che devo cavare dal file di output
    const int nHitsMax=50;

    int nHits=0;
    Double_t xCoord[nHitsMax]={0.};
    Double_t yCoord[nHitsMax]={0.};
    Double_t zCoord[nHitsMax]={0.};
    Double_t eDep[nHitsMax]={0.};

//////////////////////////////////////////////////////////////

void PlotHits(vector<vector<Double_t>> hcoord,string name,int ev){

    int nHits=hcoord.size();

    //create array of x and array of z for the plot
    Double_t x[nHits],z[nHits];
    for(int i=0;i<nHits;i++){
        x[i]=hcoord[i][0];
        z[i]=hcoord[i][1];
    }

    TGraph *g=new TGraph(nHits,z,x);
    TCanvas *c=new TCanvas("c","Hits");

    //plot the hits
    g->SetMarkerStyle(3);
    g->Draw("AP");
    stringstream ss;
    ss<<ev;
    name=name+ss.str()+".png";
    c->SaveAs(name.c_str());
    //delete c;

}

void PlotHits3D(Double_t x[],Double_t y[],Double_t z[],int n,string name,int ev){

    TGraph2D *g=new TGraph2D(nHits,x,y,z);
    TCanvas *c=new TCanvas("c","Hits");

    //plot the hits
    g->SetMarkerStyle(3);
    g->Draw("AP");
    stringstream ss;
    ss<<ev;
    name=name+ss.str()+".png";
    c->SaveAs(name.c_str());
    //delete c;

}

/////////////////////////////////////////////////////////////

void clear(){
    nHits=0;
    fill_n(xCoord, nHitsMax, 0.);
    fill_n(yCoord, nHitsMax, 0.);
    fill_n(zCoord, nHitsMax, 0.);
    fill_n(eDep, nHitsMax, 0.);
}

//////////////////////////////////////////////////////////////

void DAQ(int robo=9999){
    TString inputFileName="anaOut.root";

    TFile *inFile=new TFile(inputFileName,"READ");

    TTree *T=(TTree*)inFile->Get("hitTree");

    T->SetBranchAddress( "nHits",&nHits );
    T->SetBranchAddress( "xCoord",&xCoord );
    T->SetBranchAddress( "yCoord",&yCoord );
    T->SetBranchAddress( "zCoord",&zCoord );
    T->SetBranchAddress( "eDep",&eDep );

    int nEvs=T->GetEntries();
    cout<<"Si sono verificate "<<nEvs<<" conversioni"<<endl;

    //This vectors of vectors will contain the hits' x and y coordinates
    vector<vector<Double_t>> xz;
    vector<vector<Double_t>> yz;

    if (robo==9999){
        //The for loop run through every event in the tree
        for (int Ev=0;Ev<nEvs;Ev++){
            clear();
            xz.clear();
            yz.clear();

            T->GetEntry(Ev);

            //this loop fill the vectors of vectors from the tree branches
            for (int i=0;i<nHits;i++){
                if(eDep[i]>.01 && zCoord[i]<60){
                    vector<Double_t> hitx(2);
                    hitx[0]=xCoord[i];
                    hitx[1]=zCoord[i];
                    xz.push_back(hitx);

                    vector<Double_t> hity(2);
                    hity[0]=yCoord[i];
                    hity[1]=zCoord[i];
                    yz.push_back(hity);
                }
            }
        PlotHits(xz,"simuXZ",Ev);
        PlotHits(yz,"simuYZ",Ev);
        }
    }else{
        clear();
        xz.clear();
        yz.clear();

        T->GetEntry(robo);

        PlotHits3D(zCoord,xCoord,yCoord,nHits,"simu",robo);
    }

}