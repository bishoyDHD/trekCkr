/*****************************************************************
 *
 *             WARNING!! ACHTUNG!! WARNING!! ACHTUNG!!
 *  (Need calibration parameter file in order to run this plugin)
 *
 * ***************************************************************/

#include <Det_ClusterCsI.h>
#include <mn2CsIfn.h>
#include<cmath>
mn2CsIfn minu2;
Det_ClusterCsI::Det_ClusterCsI(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p){
  // Set defaults for various options
  treeCali=0;
  h1Mnft[12][2][2][16]=NULL;
  h1Diff[12][2][2][16]=NULL;
  h1Fits[12][2][2][16]=NULL;
  pi0Etot=NULL;
  h1clust=NULL; h1sclus=NULL;
  E2g=NULL; h2Ene=NULL; h2ang=NULL; h2deg=NULL;
  h2theta=NULL; h2phi=NULL; h2Ang=NULL;
  E_cut=NULL, cosTheta=NULL;
  h1Mpi0=NULL;
  h1Mpi02=NULL;
  h1pi0px=NULL; h1pi0py=NULL; h1pi0pz=NULL;
  h1vertpx=NULL; h1vertpy=NULL; h1vertpz=NULL;
  vertOp=NULL;
  h2corrAng=NULL;
  //paramFile.open("kpi2evenlist.txt");
  //parfile.open("calibPar.txt");
  std::cout<<"....checking this shit \n";
  resetH=false;
  notfire=false;
};

Det_ClusterCsI::~Det_ClusterCsI(){
  delete s;
  delete treeClus;
  //outFile.close();
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<16;iModule++){
          delete h1Fits[iClock][iFB][iUD][iModule];
          delete h1Mnft[iClock][iFB][iUD][iModule];
          delete h1Diff[iClock][iFB][iUD][iModule];
        }
      }
    }
  }
  cout<<"  Exiting fitting program \n";
};

// this function template is used to obtain max std::vector value
// then iterate in decending order from max->min
typedef vector<double> ve;
extern ve indexph;
ve indexph;
std::vector<double> *phval, *clusth, *clusphi;
std::size_t get_nthIndex(ve, std::size_t k){
  std::vector<std::size_t> indexes(indexph.size());
  std::iota(indexes.begin(), indexes.end(), 0);

  std::nth_element(indexes.begin(), indexes.begin() + k, indexes.end(),
    [&](int lhs, int rhs){
      return indexph[lhs] > indexph[rhs];
    }
  );
  return indexes[k];
}

// Method to read in external files such as calibration par's
void Det_ClusterCsI::readFiles(){
  string fileName="calibPar.txt";
  double calib;
  int iclock, ifb, iud, imodule;
  parfile.open(fileName);
  if(parfile){
    for(int iClock=0;iClock<12;iClock++){
      for(int iFB=0;iFB<2;iFB++){
        for(int iUD=0;iUD<2;iUD++){
          for(int iModule=0;iModule<16;iModule++){
            parfile>>iclock>>ifb>>iud>>imodule>>calib;
            calibpar[iclock][ifb][iud][imodule]=calib;
          }
        }
      }
    }
    parfile.close();
  }else{
    std::cerr<<" ****** Error opening file ''"<<fileName<<".'' Please make sure you have the file \n";
    std::abort();
  }
}

double Det_ClusterCsI::round(double var){
  // type cast to int in order to obtain integer var
  // then divide by 100 or the value is set to correct decimal points
  // Useful upto 1e4 higher than this is basically a no-no
  double value=(int)(var*100+.5);
  return (double)value/100;
}
Long_t Det_ClusterCsI::angleCsI(int id, int module, int channel, int yy, int zz){
  //std::cout<<"\n --- "<<id<<"-"<<module<<"-"<<channel<<"-"<<yy<<"-"<<zz<<" ---\n";
  phiCsI[module][channel]=yy;
  thetaCsI[module][channel]=zz;
  std::cout<<"\n --- thetaCsI["<<module<<"]["<<channel<<"] = "<<yy<<" <---> "<<thetaCsI[module][channel]<<" ---\n";

  return 0;
}

// utilize this as an initialization method along with histogram definitions
Long_t Det_ClusterCsI::histos(){
  pi0Etot=dH1("pi0Etot", " Total Energy of #pi^0",62.5,0.,0.5);
  E2g=dH1("E2g", " Total Energy of 2#gamma",82.5,0.,0.9);
  h1Mpi0= dH1("Mpi0", " Invariant mass of #pi^{0}",62.5,0.,0.5);
  h1Mpi02=dH1("Mpi02", " Invariant mass M^{2} of #pi^{0}",62.5,-.2,0.2);
  h1pi0px=dH1("h1pi0px", " #gamma momentum direction p_{x}",62.5, -0.4, 0.4);
  h1pi0py=dH1("h1pi0py", " #gamma momentum direction p_{y}",62.5, -0.4, 0.4);
  h1pi0pz=dH1("h1pi0pz", " #gamma momentum direction p_{z}",62.5, -0.4, 0.4);
  h1vertpx=dH1("h1vertpx", " #gamma momentum direction p_{x}",62.5, -0.4, 0.4);
  h1vertpy=dH1("h1vertpy", " #gamma momentum direction p_{y}",62.5, -0.4, 0.4);
  h1vertpz=dH1("h1vertpz", " #gamma momentum direction p_{z}",62.5, -0.4, 0.4);
  h2ang=dH2("h2ang", "Angular distribution #phi vs #theta", 24.0,0,M_PI, 48.0,0.,2.*M_PI);
  h2deg=dH2("h2deg", "Angular distribution #phi vs #theta deg", 24.0,0.,180., 48.0,0.,360.);
  h2Ang=dH2("h2Ang", "Mapped Angular dist #phi vs #theta deg", 24.0,0.,180., 48.0,0.,360.);
  h2Ene=dH2("h2Ene","E_{tot}(#pi^{+} + #pi^{0}) vs. E_{tot}(2#gamma + #pi^{0})", 62.5,0.,1.,62.5,0.,.7);
  h2theta=dH2("h2th","#theta(izz) comparison", 48.0,0.,180.,48.,0.,180.);
  h2phi=dH2("h2ph","#phi(iyy) comparison", 48.0,0.,360.,48.,0.,360.);
  //h2corrAng=dH2("h2cAng", "cos(#theta_{#pi^{+}}) vs. cos(#theta_{#gamma#gamma})",25.,-1.,1.,25,-1.,1.);
  E_cut=dH1("E_cut", "#pi^{0} total energy", 62.5, 0., 0.25);
  cosTheta=dH1("cosTheta", "opening angle for 2 #gamma's", 25., 0., 100.);
  vertOp=dH1("vertOp", "Opening angle between #pi^{+} and #pi^{0}", 50.,-1.1,1.1);
  h1clust=dH1("h1clust", "Cluster multiplicity", 12, 0., 12.);
  h1sclus=dH1("h1sclust", "Single cluster multiplicity", 12, 0., 12.);
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<16;iModule++){
          std::ostringstream name1, name2, name3;
          name1<<"stat_"; name2<<"Mnfit_"; name3<<"Diff_";
          name1<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name2<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name3<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          h1Fits[iClock][iFB][iUD][iModule]=new TH1D(name1.str().c_str(),"stat",250,0,250);
          h1Mnft[iClock][iFB][iUD][iModule]=new TH1D(name2.str().c_str(),"stat",250,0,250);
          h1Diff[iClock][iFB][iUD][iModule]=new TH1D(name3.str().c_str(),"stat",250,0,250);
          calibpar[12][2][2][16]=0;
        }
      }
    }
  }
  // parameter specific histos
  s = new TSpectrum(4);
  for(int n=0; n<4; n++){
    px[n]=0; py[n]=0; pz[n]=0; // initialize array
  }
  std::ostringstream title;
  title<<"CsI(Tl) clusters";
  //h2clus=dH2("h2clust",title.str().c_str(), 30,-15,15,50,0,50);
  h2clus=dH2("h2clust",title.str().c_str(), 24,0.0,180,48,0,360.); //<-- Mapped angles
  h2clus->GetXaxis()->SetTitle("#theta [deg]");
  h2clus->GetYaxis()->SetTitle("#phi [deg]");
  c1= new TCanvas("c1","",900,800);
  c1->cd();
  c1->Update();
  readFiles();
  return 0;
}

Long_t Det_ClusterCsI::startup(){
  getBranchObject("vf48",(TObject **) &treeRaw);
  getBranchObject("testBranch",(TObject **) &tracktree);
  getBranchObject("RawBeamInfo",(TObject **) &treeBeam);
  treeClus=new CRTClusterCsI();
  makeBranch("treeClus",(TObject **) &treeClus);
  gStyle->SetOptStat(0);
  //outFile.open("kpi2evtList.dat");
  string fileN="mapping.txt";
  //mapfile.open(fileN);
  //mapfile<<"module"<<"\t"<<"channel"<<"\t"<<"iyy"<<"\t"<<"izz"<<"\t"<<"TKO"<<"\t"<<"label"<<"\t"<<"U/D"<<"\t"<<"F/B"<<"\t"<<"theta(izz)"<<"\t"<<"theta(Bishoy)"<<std::endl;

  return 0;
}

void Det_ClusterCsI::initVar(){
  int dummy=-1000;
  // init 2 gamma variables
  treeClus->g1Px=dummy;       treeClus->g2Px=dummy;
  treeClus->g1Py=dummy;       treeClus->g2Py=dummy;
  treeClus->g1Pz=dummy;       treeClus->g2Pz=dummy;
  // init pion vars
  treeClus->piPpx=dummy;       treeClus->pi0px=dummy;
  treeClus->piPpy=dummy;       treeClus->pi0py=dummy;
  treeClus->piPpz=dummy;       treeClus->pi0pz=dummy;
  treeClus->piPpi0=dummy;      treeClus->piP2g=dummy;
  treeClus->ggCosTheta=dummy;  treeClus->piCosTheta=dummy;
  // WaveID and cluster var
  treeClus->waveID=dummy;
  treeClus->dubP_1=dummy;
  treeClus->channel=dummy;
  // initialize some vector variables 
  treeClus->clusterM=dummy;
  treeClus->ClustCrys=dummy;
  treeClus->Ncrys=dummy;
  treeClus->thetaE=dummy;
  treeClus->phiE=dummy;

  treeClus->E_pi0=dummy;
  treeClus->M_pi0=dummy;
  treeClus->pi0M2=dummy;
  treeClus->M_k=dummy;
  treeClus->kM2=dummy;
  treeClus->g1E=dummy;
  treeClus->g2E=dummy;
}
Long_t Det_ClusterCsI::process(){
  phval=new vector<double>();
  clusth=new vector<double>();
  clusphi=new vector<double>();
  int iHist=0;
  clus_csi=false;
  idCrys.clear(),  indexph.clear();
  typeAB.clear(),  gud.clear(),      gno.clear();
  csThet.clear(),  csPhi.clear();
  csiph.clear();   phval->clear();   csiClus.clear();
  clusth->clear(); clusphi->clear(); clusEne.clear();
  singleEne.clear();
  clusThetaE.clear(); clusPhiE.clear();
  int evtNum=treeRaw->eventNo;
  treeClus->evtNo=evtNum;
  initVar(); // make sure variables are initialized V events
  //std::cout<<"\n\n event number is :"<<treeRaw->eventNo<<"\n\n";
  //std::cout<<"\n\n event number2 is:"<<tracktree->evtNum<<"\n\n";
  /*
  if(ppip>=0.195 && ppip<=0.215){
    outFile<<tracktree->evtNum+1<<std::endl;
    std::cout<<"  >>>>>>>>>>>>>>>>>> "<<ppip<< "<<>> "<<tracktree->evtNum<<std::endl;
  }*/
  ppip=tracktree->pVertpi0;
  if(ppip<0.195 || ppip>0.215) goto exitLoop;
  //if(resetH)
    //h2clus->Reset(); //need to reset stats in cluster event viewer
    //std::cout<<"  >>>>>>>>>>>>>>>>>> "<<ppip<<std::endl;
  for(UInt_t i=0;i<treeRaw->nChannel;i++){ // loop over fired crystals
    char* p=(char*)&(treeRaw->nameModule[i]);
    int moduleName=(p[3]-'0')*10+(p[2]-'0')-1;
    //cout<< " Index clock: "<<indexClock<<endl;
    string nameModule;
    nameModule+=(*p);
    p++;
    nameModule+=*p;
    p++;
    nameModule+=*p;
    p++;
    nameModule+=*p;
    string nameCsI;
    p=(char*)&(treeRaw->nameCsI[i]);
    int indexClock=(p[3]-'0')*10+(p[2]-'0')-1;
    p+=3;
    nameCsI+=(*p);
    p--;
    nameCsI+=*p;
    p--;
    nameCsI+=*p;
    p--;
    nameCsI+=*p;
  
    //cout<< "  ***** THIS IS A TEST! TESTING! TESTING! "<<nameCsI<<"\t"<<nameModule<<endl;
    int indexModule=treeRaw->indexCsI[i]-1;
    int indexFB=0;
    if(p[1]=='b' || p[1]=='B') indexFB=1;
    int indexUD=0;
    ud=1;
    if(p[0]=='d' || p[0]=='D'){
      ud=0;
      indexUD=1;
    }
    /*
    if(treeRaw->indexCsI[i]==16 && indexClock==0 && indexFB==0 && indexUD==0){
      //if(indexClock==0 && indexFB==0 && indexUD==0){
        for(UInt_t iData=0;iData<treeRaw->nSample[i];iData++){
          //h1time[indexClock][indexFB][indexUD][indexModule]->SetBinContent(iData+1,treeRaw->data[i][iData]);
        }
      //}
    }*/
    // Looking at signal modules: ignore timing and reference modules
    if(!(treeRaw->indexCsI[i]==16 && indexFB==0 && indexUD==0) ||
                    !(indexClock==0 || indexClock==2 || indexClock==4 ||
                            indexClock==6 || indexClock==8 || indexClock==10)){
      for(UInt_t iData=0;iData<treeRaw->nSample[i];iData++){
        h1Fits[indexClock][indexFB][indexUD][indexModule]->SetBinContent(iData+1,treeRaw->data[i][iData]);
      }
      #pragma omp parallel num_threads(8)
      xpos.clear();
      val.clear();
      // Utilize fitting method
      x1=h1Fits[indexClock][indexFB][indexUD][indexModule]->
      	GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
      calib=calibpar[indexClock][indexFB][indexUD][indexModule];
      calib=calib/1000.; // convert MeV to GeV
      //FIXME this is makeshift solution to be corrected later
      if(treeRaw->indexCsI[i]==16) calib=.22/1000.;
      if(x1>=58 && x1<=68){
        if(treeRaw->nChannel>=7){ // Start by checking how many CsI crystals have fired
	  //if(treeRaw->indexCsI[i]==16){
            std::cout<< "\n\n ****************************************** "<<endl;
            std::cout<< " Index clock: "<<indexClock<<endl;
            std::cout<< " Gap config FB is  : " <<p[1]<<endl;
            std::cout<< " Gap config UD is  : " <<p[0]<<endl;
	    std::cout<< " TKO Module Name   : "<<nameModule<<std::endl;
            std::cout<< " Channel number is : " <<treeRaw->indexChannel[i]<<endl;
            std::cout<< " size of nSample is  : " <<treeRaw->nSample[i]<<endl;
            std::cout<< "\n   ** Event Number ------> "<<treeRaw->eventNo<<"\n\n";
	  //}
	  // convert TKO module to module number here
	  if(nameModule=="TK34") moduleNo=1;
	  if(nameModule=="TK32") moduleNo=2;
	  if(nameModule=="TK36") moduleNo=3;
	  if(nameModule=="TK37") moduleNo=4;
	  if(nameModule=="TK38") moduleNo=5;
	  if(nameModule=="TK08") moduleNo=6;
	  if(nameModule=="TK50") moduleNo=7;
	  if(nameModule=="TK09") moduleNo=8;
	  if(nameModule=="TK54") moduleNo=9;
	  if(nameModule=="TK31") moduleNo=10;
	  if(nameModule=="TK04") moduleNo=11;
	  if(nameModule=="TK45") moduleNo=12;
	  if(nameModule=="TK33") moduleNo=13;
	  if(nameModule=="TK39") moduleNo=14;
	  if(nameModule=="TK40") moduleNo=15;
	  if(nameModule=="TK41") moduleNo=16;

	  std::cout<< " CsI module Nmber is : "<<moduleNo<<std::endl;
	  std::cout<<" ---- Checking XML Parameter --> "<<thetaCsI[moduleNo-1][treeRaw->indexChannel[i]]<<", "<<phiCsI[moduleNo-1][treeRaw->indexChannel[i]]<<std::endl;
          double mn = h1Fits[indexClock][indexFB][indexUD][indexModule]->
          	GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMinimumBin());
          y1 = h1Fits[indexClock][indexFB][indexUD][indexModule]->
          	GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(x1));
          double bl = h1Fits[indexClock][indexFB][indexUD][indexModule]->
          	GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(mn));
          for(int ivar=1; ivar<=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetNbinsX(); ivar+=1){
            xpos.push_back(ivar);
            val.push_back(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinContent(ivar));
          }
          if(y1 == 1023){
            double x;
            for(int ovr=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin(); ovr<250; ovr +=1){
              while(y1 == 1023){
                x = ovr++;
                y1 = h1Fits[indexClock][indexFB][indexUD][indexModule]->
          	      GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(x));
              }
            }
	    std::cout<< " Okay this thing is smart x =" << x << endl;
	    TF1* f1=new TF1("wave1",(minu2.overrangemodel()).c_str(), 0.5, x1);
            // fill parameters for the fit function(s)
            for(int ivar=0; ivar<10; ivar+=1){
              f1->SetParameter(ivar, minu2.par(ivar));
              f1->SetParLimits(ivar,minu2.parmin(ivar),minu2.parlim(ivar));
            }
            nfound=s->Search(h1Fits[indexClock][indexFB][indexUD][indexModule], 2,"",0.10);
            double *xpeaks=NULL;
            xpeaks=s->GetPositionX();
            double posX[2];
            for(int ivar=0; ivar<nfound; ivar++){
              double a=xpeaks[ivar];
              int bin=1+Int_t(a+.5);
              posX[ivar]=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(bin);
            }
            sort(xpeaks,xpeaks+nfound);
            //std::cout<<"\n ---->  Gotta check the hell outta this shit \n";
            double rtime=(xx1+xx2)/2;
            f1->SetParameter(0,y1+1023*2.8);
            f1->SetParLimits(0,y1-61.7,y1+1723.7);
            f1->SetParameter(1,rtime+20.1);
            f1->SetParLimits(1,rtime-261.7,rtime+271.7);
            f1->SetParameter(8,bl);
            f1->SetParLimits(8,bl-61.7,bl+171.7);
            h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f1); //, "0");
            f1chi2=f1->GetChisquare();
            xx1=(int)x1; xx2=(int)x; ymax=y1;
            ovrfn ffcn(xpos, xx1, xx2, val, ymax);
            // Create wrapper for minimizer
            MnUserParameters upar2;
            std::vector<double> par(10), err(10);
            par.clear(); err.clear();
            //std::cout<< "  ----> Testing left and right x-limits: "<<xx1<< ", "<<xx2<<endl;
            for(int n=0; n<10;n+=1){
              upar2.Add(minu2.nameL(n).c_str(), f1->GetParameter(n), 0.1);
            }/*
            ovrfn ffcn(xpos, xx1, xx2, val, ymax);
            // Create wrapper for minimizer 
	    std::cout<<"  ....Oops just making sure this thing works! \n";
            MnUserParameters upar2;
            std::vector<double> par(10), err(10);
            par.clear(); err.clear();
	    std::cout<< "  ----> Testing left and right x-limits: "<<xx1<< ", "<<xx2<<endl;
            for(int n=0; n<10;n+=1){
              upar2.Add((minu2.nameL(n)).c_str(), ovrpar[n], 0.1);
            }*/
            // create Migrad minimizer
            MnMigrad migrad(ffcn, upar2);
            FunctionMinimum min = migrad(180,1e-2);
            std::cout<<"minimum: "<<min<<std::endl;
            std::vector<double> param;
            param.clear();
            for(int imn2=0; imn2<10; imn2+=1){
              param.push_back(migrad.Value((minu2.nameL(imn2)).c_str()));
              //std::cout<< "  par["<<i<<"] value --> ["<<param[i]<<"] \n";
            }
            for(int imn2=1; imn2<+h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetNbinsX()+1; imn2+=1){
              double x=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(imn2);
              double yv=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinContent(imn2);
              double mnfit=minu2.overR(x, param);
              double res=100*(yv-mnfit)/yv;
              h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(imn2, mnfit);
              if(imn2>=xx1 && imn2<=xx2){
                h1Diff[indexClock][indexFB][indexUD][indexModule]->SetBinContent(imn2, 0);
              }else{
                h1Diff[indexClock][indexFB][indexUD][indexModule]->SetBinContent(imn2, res);
              }
            }
            double mnx,mny,max,may;
            max=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                  GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
            if(max>=60 && max<=65){
	      if(!clus_csi)
                clus_csi=true;
              if(!resetH) resetH=true;
              clock=indexClock;
	      //std::cout<< "\n\n  ======> "<<clock<<" <========\n\n";
              fb=indexFB;
              module=indexModule;
              mnx=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                    GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMinimumBin());
              may=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                    GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(max));
              mny=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                    GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(mnx));
              //Integrating fitting function for more accurate gain calibration
              TAxis* axis=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetXaxis();
              int xmin=0, xmax=250;
              int bmin=axis->FindBin(xmin);
              int bmax=axis->FindBin(xmax);
              int integral=h1Mnft[indexClock][indexFB][indexUD][indexModule]->Integral(bmin,bmax);
              integral-=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
        		GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
              integral-=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
        		GetBinContent(bmax)*(axis->GetBinLowEdge(bmin)-bmax)/axis->GetBinWidth(bmax);
              double area=integral-(mny*250);
              calInt=area;
	      std::cout<<" =====> Value for integral1 is:" <<area<<endl;
              double diff=(may-mny)*calib; //apply E-calib
              int imod=0;
	      int igap=4*indexClock+2, igapU=4*indexClock+3;
              int csimod=(-1*treeRaw->indexCsI[i])+1;  //default
              if(p[0]=='d' || p[0]=='D'){
                igap=4*indexClock+2;
                csimod=(-1*treeRaw->indexCsI[i])+1;
                if(treeRaw->indexCsI[i]>=10){
                  igap=4*indexClock+1;
                  csimod=(-1*treeRaw->indexCsI[i])+7;
                  if(treeRaw->indexCsI[i]==16) igap=4*indexClock+1.5;
                }
                if(p[1]=='b' || p[1]=='B'){
                  csimod=treeRaw->indexCsI[i];
                  if(treeRaw->indexCsI[i]>=10){
                    igap=4*indexClock+1;
                    csimod=(treeRaw->indexCsI[i])-6;
                    if(treeRaw->indexCsI[i]==16) igap=4*indexClock+1.5;
                  }
	        }	
	      }
              if(p[0]=='u' || p[0]=='U'){
                igap=4*indexClock+3;
                if(treeRaw->indexCsI[i]>=10){
                  igap=4*(indexClock+1);
                  csimod=(-1*treeRaw->indexCsI[i])+7;
                  if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
                }
		//FIXME
                if(p[1]=='b' || p[1]=='B'){
                  csimod=treeRaw->indexCsI[i];
                  if(treeRaw->indexCsI[i]>=10){
                    igap=4*(indexClock+1);
                    csimod=(treeRaw->indexCsI[i])-6;
                    if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
	          }
                }
	      }
              indexph.push_back(diff);
              phval->push_back(diff);
              int tAB=0;
              if((p[0]=='u' || p[0]=='U') && (indexModule>=9)) tAB=1;
	      std::cout<< " reading out value for Up and TypeB Crystal: "<<tAB<<endl;
              if((p[0]=='d' || p[0]=='D') && (indexModule<=8)) tAB=1;
              idCrys.push_back(indexModule), gfb.push_back(fb);
              typeAB.push_back(tAB), gno.push_back(indexClock), gud.push_back(ud);
              double csitheta=theta[fb][indexModule], csiphi=phi[indexClock][ud][tAB];
	      // handeling phi for channel 16: Crystal type 10
	      if(treeRaw->indexCsI[i]==16){
	        csiphi=indexClock*30.;
		if(p[0]=='u' || p[0]=='U') csiphi=indexClock*30.+15.;
	      }
	      int thetaIndex=thetaCsI[moduleNo-1][treeRaw->indexChannel[i]];
	      int phiIndex=phiCsI[moduleNo-1][treeRaw->indexChannel[i]];
	      //otheta=mapTheta[thetaIndex];
	      //otheta=mapTheta[20-thetaIndex];
	      //reading from Osaka map as is: otheta and ophi angles
	      double cosTheta=crysZ[20-thetaIndex]/std::sqrt(std::pow(crysZ[20-thetaIndex],2)+
			      std::pow(crysr[20-thetaIndex],2));
	      double acosTheta=TMath::ACos(cosTheta);
	      otheta=TMath::RadToDeg()*acosTheta; // convert to deg.
	      otheta=round(otheta); // make sure to obtain 2 dp
	      mapPhi=2*180*(phiIndex-0.5)/48.; // from mapping init file
	      ophi=mapPhi; // convert to global scheme: gap3=+xaxis
	      //ophi=90.-mapPhi; // convert to global scheme: gap3=+xaxis
	      //if(ophi < 0.) ophi=360.+ophi; // loop around the clock
              cout<< " *** Regular Angle "<<csitheta<<"  "<<csiphi<<endl;
              cout<< " ***  Osaka Angles "<<otheta<<"  "<<ophi<<endl;
	      h2Ang->Fill(otheta,ophi);
	      h2theta->Fill(otheta,csitheta);
	      h2phi->Fill(ophi,csiphi);
              auto angles=std::make_pair(otheta,ophi);
              csThet.push_back(otheta), csPhi.push_back(ophi);
              //thetaPhi.push_back(angles);
              csiph[angles]=diff;
              csiClus[angles]=true;
              //h2clus->Fill(csimod,igap,diff);
              h2clus->Fill(otheta,ophi,diff);
              // writing to check mapping file
              //mapfile<<moduleNo<<"\t"<<treeRaw->indexChannel[i]+1<<"\t"<<phiCsI[moduleNo-1][treeRaw->indexChannel[i]]<<"\t"<<thetaCsI[moduleNo-1][treeRaw->indexChannel[i]]<<"\t"<<nameModule<<"\t"<<indexClock+1<<"\t"<<p[0]<<"\t"<<p[1]<<std::endl;
	      treeClus->waveID=5;
            }// <--- Use this to get rid of double and single fitting functions * /
          } //<-- end of overrange if loop
	  if(y1<1023){
            //cout<< "  Size of x is:  "<<xpos.size()<<endl;
            xx1=xpos.size(); xx2=0; ymax=y1;
            nfound=s->Search(h1Fits[indexClock][indexFB][indexUD][indexModule], 2,"",0.10);
            if(nfound>=3){
	      goto exitLoop;
              double *xpeaks=s->GetPositionX();
              sort(xpeaks,xpeaks+nfound);
	      if(nfound==4){
	        if(xpeaks[3]-xpeaks[2]<=25) nfound==3;
	      }
	      int parV=13;
              std::cout<<"\n ------- Within Signal Loop Event number is:  "<<treeRaw->eventNo<<" -------\n\n";
	      std::string pileUp=minu2.triplemodel();
	      if(nfound==4)
                pileUp=minu2.quadruplemodel();
              TF1* f1=new TF1("f1",pileUp.c_str(),0.0,250);
              for(int nnp=0; nnp<13; nnp+=1){
                f1->SetParameter(nnp,minu2.par(nnp));
                f1->SetParLimits(nnp,minu2.parmin(nnp),minu2.parlim(nnp));
              }
              double posX[4];
              double valY[4];
	      for(int nval=0; nval<4; nval++){
	        posX[nval]=0;
		valY[nval]=0;
	      }
              for(int ivar=0; ivar<nfound; ivar++){
                double a=xpeaks[ivar];
                int bin=1+Int_t(a+.5);
                posX[ivar]=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(bin);
              }
	      for(int valy=0; valy<nfound; valy++)
                valY[valy]=h1Fits[indexClock][indexFB][indexUD][indexModule]->
          	    GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(xpeaks[valy]));
              f1->SetParameter(0,valY[0]);
              f1->SetParLimits(0,valY[0]-11.7,valY[0]+71.7);
              //f1->SetParLimits(0,y1-61.7,y1+971.7);
              f1->SetParameter(1,xpeaks[0]-5.);
              f1->SetParLimits(1,xpeaks[0]-25.7,xpeaks[0]+10.7);
              f1->SetParameter(8,bl);
              f1->SetParLimits(8,bl-61.7,bl+171.7);
              f1->SetParameter(9,xpeaks[1]-1.9);
              f1->SetParLimits(9,xpeaks[1]-31.7,xpeaks[1]+21.7);
              f1->SetParameter(10,valY[1]);
              f1->SetParLimits(10,valY[1]-61.7,valY[1]+17.7);
              f1->FixParameter(12,xpeaks[2]-15.1);
              //f1->SetParLimits(12,xpeaks[2]-61.7,xpeaks[2]+71.7);
              f1->SetParameter(11,valY[2]);
              f1->SetParLimits(11,valY[2]-41.7,valY[2]+25.77);
	      if(nfound==4){
	        f1->SetParameter(13,valY[3]);
                f1->SetParLimits(13,valY[3]-61.7,valY[3]+10.7);
                f1->SetParameter(14,xpeaks[3]-9.1);
                f1->SetParLimits(14,xpeaks[3]-21.7,xpeaks[3]+11.7);
	        parV=15;
	      }
              h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f1); //,"0");
      	      f1chi2=f1->GetChisquare();
              std::cout<<" ****** Checking the position of the peaks: "<<xpeaks[0]<<", "<<xpeaks[1];
	      std::cout<<", "<<xpeaks[2]<<", "<<xpeaks[3]<<std::endl;
      	      //std::cout<<" **************\n ******** Chi2 for F1 fit "<<f1->GetChisquare()<<endl;
              Minuit2Minimizer* mnu2=new Minuit2Minimizer("Minuit2");
              // Create wrapper for minimizer
              fitfn3 ffcn1(xpos, xx1, xx2, val, ymax);
              fitfn4 ffcn2(xpos, xx1, xx2, val, ymax);
              std::vector<double> param;
              std::vector<double> parm(15), err(15);
              param.clear(); parm.clear(); err.clear();
              MnUserParameters upar;
              for(int n=0; n<parV;n+=1){
                upar.Add(minu2.nameL(n).c_str(), f1->GetParameter(n),0.1); //,parmin(n), parlim(n), 0.1);
                //upar.Add(nameL(n).c_str(), f1->GetParameter(n), parmin(n), parlim(n), 1e-3);
              }
              // create Migrad minimizer
	      if(nfound==3){
                MnMigrad migrad(ffcn1, upar);
                //FunctionMinimum min = migrad();  //6000,1e-9);
                FunctionMinimum min = migrad(90,1e-12);
                std::cout<<"minimum: "<<min<<std::endl;
                //MnHesse hesse;
                for(int ivar=0; ivar<parV; ivar+=1){
                  param.push_back(migrad.Value(minu2.nameL(ivar).c_str()));
                  //std::cout<< "  par["<<i<<"] value --> ["<<param[ivar]<<"] \n";
                }
	      }
	      if(nfound==4){
	        MnMigrad migrad(ffcn2, upar);
                //FunctionMinimum min = migrad();  //6000,1e-9);
                FunctionMinimum min = migrad(90,1e-12);
                std::cout<<"minimum: "<<min<<std::endl;
                //MnHesse hesse;
                for(int ivar=0; ivar<parV; ivar+=1){
                  param.push_back(migrad.Value(minu2.nameL(ivar).c_str()));
                  //std::cout<< "  par["<<i<<"] value --> ["<<param[ivar]<<"] \n";
                }
	      }
              double mnfit=0;
              for(int ivar=0; ivar<h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetNbinsX()+1; ivar+=1){
                double x=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(ivar);
                double yv=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinContent(ivar);
                mnfit=minu2.model3(x, param);
		if(nfound==4)
                  mnfit=minu2.model4(x, param);
                h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(ivar, mnfit);
                double res=100*(yv-mnfit)/yv;
                h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(ivar, mnfit);
                h1Diff[indexClock][indexFB][indexUD][indexModule]->SetBinContent(ivar, res);
              }
              double mnx,mny,max,may;
              max=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
          	    GetBinLowEdge(h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
              if(max>=60 && max<=65){
	        treeClus->waveID=nfound;
		//goto exitLoop;
	        if(!clus_csi)
                  clus_csi=true;
		if(!resetH) resetH=true;
                mnx=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinLowEdge(h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetMinimumBin());
                may=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(h1Mnft[indexClock][indexFB][indexUD][indexModule]->FindBin(max));
                mny=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(h1Mnft[indexClock][indexFB][indexUD][indexModule]->FindBin(mnx));
                TAxis* axis=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetXaxis();
                int xmin=0, xmax=250;
                int bmin=axis->FindBin(xmin);
                int bmax=axis->FindBin(xmax);
                int integral=h1Mnft[indexClock][indexFB][indexUD][indexModule]->Integral(bmin,bmax);
                integral-=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
                integral-=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(bmax)*(axis->GetBinLowEdge(bmin)-bmax)/axis->GetBinWidth(bmax);
                double diff=(may-mny)*calib;
		double area=integral-mny*250;
        	clock=indexClock;
        	fb=indexFB;
        	ud=indexUD;
        	module=indexModule;
                int imod=0, igap=4*indexClock+2, igapU=4*indexClock+3;
                int csimod=(-1*treeRaw->indexCsI[i])+1;  //default
	        //std::cout<< "\n\n  ======> "<<clock<<" peaks>=3 <========\n\n";
                if(p[0]=='d' || p[0]=='D'){
                  igap=4*indexClock+2;
                  csimod=(-1*treeRaw->indexCsI[i])+1;
                  if(treeRaw->indexCsI[i]>=10){
                    igap=4*indexClock+1;
                    csimod=(-1*treeRaw->indexCsI[i])+7;
                    if(treeRaw->indexCsI[i]==16) igap=4*indexClock+1.5;
                  }
                  if(p[1]=='b' || p[1]=='B'){
                    csimod=treeRaw->indexCsI[i];
                    if(treeRaw->indexCsI[i]>=10){
                      igap=4*indexClock+1;
                      csimod=(treeRaw->indexCsI[i])-6;
                      if(treeRaw->indexCsI[i]==16) igap=4*indexClock+1.5;
                    }
	          }	
		}
                if(p[0]=='u' || p[0]=='U'){
                  igap=4*indexClock+3;
                  csimod=(-1*treeRaw->indexCsI[i])+1;
                  if(treeRaw->indexCsI[i]>=10){
                    igap=4*(indexClock+1);
                    csimod=(-1*treeRaw->indexCsI[i])+7;
                    if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
                  }
	          //FIXME
                  if(p[1]=='b' || p[1]=='B'){
                    csimod=treeRaw->indexCsI[i];
                    if(treeRaw->indexCsI[i]>=10){
                      igap=4*(indexClock+1);
                      csimod=(treeRaw->indexCsI[i])-6;
                      if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
	            }
                  }
	        }
                indexph.push_back(diff);
                phval->push_back(diff);
                int tAB=0;
                if((p[0]=='u' || p[0]=='U') && (indexModule>=9)) tAB=1;
	        std::cout<< " reading out value for Up and TypeB Crystal: "<<tAB<<endl;
                if((p[0]=='d' || p[0]=='D') && (indexModule<=8)) tAB=1;
                idCrys.push_back(indexModule), gfb.push_back(fb);
                typeAB.push_back(tAB), gno.push_back(indexClock), gud.push_back(ud);
                double csitheta=theta[fb][indexModule], csiphi=phi[indexClock][ud][tAB];
	        // handeling phi for channel 16: Crystal type 10
	        if(treeRaw->indexCsI[i]==16){
	          csiphi=indexClock*30.;
	          if(p[0]=='u' || p[0]=='U') csiphi=indexClock*30.+15.;
	        }
	        int thetaIndex=thetaCsI[moduleNo-1][treeRaw->indexChannel[i]];
	        int phiIndex=phiCsI[moduleNo-1][treeRaw->indexChannel[i]];
	        //otheta=mapTheta[thetaIndex];
	        //otheta=mapTheta[20-thetaIndex];
	        //reading from Osaka map as is: otheta and ophi angles
	        double cosTheta=crysZ[20-thetaIndex]/std::sqrt(std::pow(crysZ[20-thetaIndex],2)+
	          	      std::pow(crysr[20-thetaIndex],2));
	        double acosTheta=TMath::ACos(cosTheta);
	        otheta=TMath::RadToDeg()*acosTheta; // convert to deg.
	        otheta=round(otheta); // make sure to obtain 2 dp
	        mapPhi=2*180*(phiIndex-0.5)/48.;
	        ophi=mapPhi; // convert to global scheme: gap3=+xaxis
	        //ophi=90.-mapPhi; // convert to global scheme: gap3=+xaxis
	        //if(ophi < 0.) ophi=360.+ophi; // loop around the clock
                cout<< " *** Regular Angle "<<csitheta<<"  "<<csiphi<<endl;
                cout<< " ***  Osaka Angles "<<otheta<<"  "<<ophi<<endl;
	        h2Ang->Fill(otheta,ophi);
	        h2theta->Fill(otheta,csitheta);
	        h2phi->Fill(ophi,csiphi);
                auto angles=std::make_pair(otheta,ophi);
                csThet.push_back(otheta), csPhi.push_back(ophi);
		// writing to check mapping file
                //mapfile<<moduleNo<<"\t"<<treeRaw->indexChannel[i]+1<<"\t"<<phiCsI[moduleNo-1][treeRaw->indexChannel[i]]<<"\t"<<thetaCsI[moduleNo-1][treeRaw->indexChannel[i]]<<"\t"<<nameModule<<"\t"<<indexClock+1<<"\t"<<p[0]<<"\t"<<p[1]<<std::endl;
                //thetaPhi.push_back(angles);
                csiph[angles]=diff;
                csiClus[angles]=true;
                //cout<< " *********** theta "<<csitheta<<"  "<<csiphi<<endl;
                //h2clus->Fill(csimod,igap,diff);
                h2clus->Fill(otheta,ophi,diff);
		delete f1;
              } // <-- End of K+ decay time if loop
            } //<-- Use to get rid of 3 peaks functions here * /
            if(nfound==2){
	      //goto exitLoop;
	      std::string dubFit = minu2.doublemodel();
              TF1* f1=new TF1("f1",dubFit.c_str(),0.0,250);
              for(int n=0; n<15; n+=1){
                f1->SetParameter(n,minu2.par(n));
                f1->SetParLimits(n,minu2.parmin(n),minu2.parlim(n));
              }
              double *xpeaks=s->GetPositionX();
              double posX[2];
              for(int ipks=0; ipks<nfound; ipks++){
                double a=xpeaks[ipks];
                int bin=1+Int_t(a+.5);
                posX[ipks]=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(bin);
              }
              sort(xpeaks,xpeaks+nfound);
	      std::cout<<" ****** Checking the position of the peaks: "<<xpeaks[0]<<", "<<xpeaks[1]<<endl;
              double yp2=h1Fits[indexClock][indexFB][indexUD][indexModule]->
          	    GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(xpeaks[1]));
              f1->SetParameter(0,y1);
              f1->SetParLimits(0,y1-6.7,y1+7.7);
              f1->SetParameter(1,xpeaks[0]-10.1);
              f1->SetParLimits(1,xpeaks[0]-26.7,xpeaks[0]+27.7);
              f1->SetParameter(8,bl);
              f1->SetParLimits(8,bl-161.7,bl+11.7);
              f1->SetParameter(9,xpeaks[1]-15.1);
              f1->SetParLimits(9,xpeaks[1]-61.7,xpeaks[1]+17.7);
              f1->SetParameter(10,yp2);
              f1->SetParLimits(10,yp2-7.7,yp2+7.7);
              h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f1); //,"0");
              // Create wrapper for minimizer 
              fitfn2 ffcn1(xpos, xx1, xx2, val, ymax);
              std::vector<double> param;
              std::vector<double> parm(15), err(15);
              param.clear(); parm.clear(); err.clear();
              MnUserParameters upar;
              for(int n=0; n<11;n+=1){
                upar.Add((minu2.nameL(n)).c_str(), f1->GetParameter(n),1e-3); //,parmin(n), parlim(n), 0.1);
                //upar.Add(nameL(n).c_str(), f1->GetParameter(n), parmin(n), parlim(n), 1e-3);
              }
              // create Migrad minimizer
              MnMigrad migrad(ffcn1, upar);
              //FunctionMinimum min = migrad();  //6000,1e-9);
              FunctionMinimum min = migrad(10,1e-2);
              std::cout<<"minimum: "<<min<<std::endl;
              //MnHesse hesse;
              for(int imn2=0; imn2<11; imn2+=1){
                param.push_back(migrad.Value((minu2.nameL(imn2)).c_str()));
                //std::cout<< "  par["<<i<<"] value --> ["<<param[i]<<"] \n";
              }
              for(int imn2=0; imn2<h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetNbinsX()+1; imn2+=1){
                double x=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(imn2);
                double yv=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinContent(imn2);
                double mnfit=minu2.model2(x, param);
                h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(imn2, mnfit);
                double res=100*(yv-mnfit)/yv;
                h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(imn2, mnfit);
                h1Diff[indexClock][indexFB][indexUD][indexModule]->SetBinContent(imn2, res);
              }
              double mnx,mny,max,may;
              max=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
          	    GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
              if(max>=60 && max<=65){
	        treeClus->waveID=nfound;
	        if(xpeaks[0] < 60){
	          treeClus->waveID=nfound;
	          treeClus->dubP_1=1;
	        }
	        if(!clus_csi)
                  clus_csi=true;
		if(!resetH) resetH=true;
                mnx=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMinimumBin());
                may=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(max));
                mny=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(mnx));
                double diff=(may-mny)*calib;
                indexph.push_back(diff);
                phval->push_back(diff);
                clock=indexClock;
                fb=indexFB;
                module=indexModule;
      	        idCrys.push_back(indexModule), gfb.push_back(fb);
      	        int tAB=0;
      	        if((p[0]=='u' || p[0]=='U') && (indexModule>=9)) tAB=1;
		std::cout<< " reading out value for Up and TypeB Crystal: "<<tAB<<endl;
      	        if((p[0]=='d' || p[0]=='D') && (indexModule<=8)) tAB=1;
      	        typeAB.push_back(tAB), gno.push_back(indexClock), gud.push_back(ud);
      	        int imod=0, igap=4*indexClock+2, igapU=4*indexClock+3;
                int csimod=(-1*treeRaw->indexCsI[i])+1;  //default
      	        double csitheta=theta[fb][indexModule], csiphi=phi[indexClock][ud][tAB];
	        // handeling phi for channel 16: Crystal type 10
	        if(treeRaw->indexCsI[i]==16){
	          csiphi=indexClock*30.;
	          if(p[0]=='u' || p[0]=='U') csiphi=indexClock*30.+15.;
	        }
	        int thetaIndex=thetaCsI[moduleNo-1][treeRaw->indexChannel[i]];
	        int phiIndex=phiCsI[moduleNo-1][treeRaw->indexChannel[i]];
	        //otheta=mapTheta[thetaIndex];
	        //otheta=mapTheta[20-thetaIndex];
	        //ophi=2*180*(phiIndex-0.5)/48.;
	        //reading from Osaka map as is: otheta and ophi angles
	        double cosTheta=crysZ[20-thetaIndex]/std::sqrt(std::pow(crysZ[20-thetaIndex],2)+
	          	      std::pow(crysr[20-thetaIndex],2));
	        double acosTheta=TMath::ACos(cosTheta);
	        otheta=TMath::RadToDeg()*acosTheta; // convert to deg.
	        otheta=round(otheta); // make sure to obtain 2 dp
	        mapPhi=2*180*(phiIndex-0.5)/48.;
	        ophi=mapPhi; // convert to global scheme: gap3=+xaxis
	        //ophi=90.-mapPhi; // convert to global scheme: gap3=+xaxis
	        //if(ophi < 0.) ophi=360.+ophi; // loop around the clock
                cout<< " *** Regular Angle "<<csitheta<<"  "<<csiphi<<endl;
                cout<< " ***  Osaka Angles "<<otheta<<"  "<<ophi<<endl;
	        h2Ang->Fill(otheta,ophi);
	        h2theta->Fill(otheta,csitheta);
	        h2phi->Fill(ophi,csiphi);
      	        auto angles=std::make_pair(otheta,ophi);
      	        //thetaPhi.push_back(angles);
      	        csThet.push_back(otheta), csPhi.push_back(ophi);
      	        csiph[angles]=diff;
      	        csiClus[angles]=true;
		// writing to check mapping file
                //mapfile<<moduleNo<<"\t"<<treeRaw->indexChannel[i]+1<<"\t"<<phiCsI[moduleNo-1][treeRaw->indexChannel[i]]<<"\t"<<thetaCsI[moduleNo-1][treeRaw->indexChannel[i]]<<"\t"<<nameModule<<"\t"<<indexClock+1<<"\t"<<p[0]<<"\t"<<p[1]<<std::endl;
      	        //cout<< " *********** theta "<<csitheta<<"  "<<csiphi<<endl;
	        //std::cout<< "\n\n  ======> "<<clock<<" peaks>=2 <========\n\n";
                if(p[0]=='d' || p[0]=='D'){
                  igap=4*indexClock+2;
                  csimod=(-1*treeRaw->indexCsI[i])+1;
      	          if(treeRaw->indexCsI[i]>=10){
      	            igap=4*indexClock+1;
                    csimod=(-1*treeRaw->indexCsI[i])+7;
      	            if(treeRaw->indexCsI[i]==16) igap=4*indexClock+1.5;
      	          }
                  if(p[1]=='b' || p[1]=='B'){
      	            csimod=treeRaw->indexCsI[i];
      	            if(treeRaw->indexCsI[i]>=10){
      	              igap=4*indexClock+1;
      	              csimod=(treeRaw->indexCsI[i])-6;
      	              if(treeRaw->indexCsI[i]==16) igap=4*indexClock+1.5;
      	            }
      	          } 
		}
                if(p[0]=='u' || p[0]=='U'){
      	          igap=4*indexClock+3;
                  csimod=(-1*treeRaw->indexCsI[i])+1;
      	          if(treeRaw->indexCsI[i]>=10){
      	            igap=4*(indexClock+1);
      	            csimod=(-1*treeRaw->indexCsI[i])+7;
      	            if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
      	          }
		  //FIXME
      	          if(p[1]=='b' || p[1]=='B'){
      	            csimod=treeRaw->indexCsI[i];
      	            if(treeRaw->indexCsI[i]>=10){
                      igap=4*(indexClock+1);
                      csimod=(treeRaw->indexCsI[i])-6;
      	              if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
                    }
      	          }
      	        } 
                //h2clus->Fill(csimod,igap,diff);
                h2clus->Fill(otheta,ophi,diff);
              }
	      delete f1;
            } //<-- Use to get rid of 2 peaks functions here * /
            if(nfound==1){
              TF1* f1=new TF1("f1",(minu2.singlemodel()).c_str(),1.0,250);
              for(int n=0; n<9; n+=1){
                f1->SetParameter(n,minu2.par(n));
                f1->SetParLimits(n,minu2.parmin(n),minu2.parlim(n));
              }
              f1->SetParameter(0,y1);
              f1->SetParLimits(0,y1-61.7,y1+971.7);
              f1->SetParameter(1,x1);
              f1->SetParLimits(1,x1-261.7,x1+571.7);
              f1->SetParameter(8,bl);
              f1->SetParLimits(8,bl-161.7,bl+171.7);
              h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f1,"0");
              Minuit2Minimizer* mnu2=new Minuit2Minimizer("Minuit2");
              // Create wrapper for minimizer 
              //TMinuit* gmin = new TMinuit(8);
              fitfn ffcn1(xpos, xx1, xx2, val, ymax);
              std::vector<double> param;
              std::vector<double> parm(15), err(15);
              param.clear(); parm.clear(); err.clear();
              MnUserParameters upar;
              for(int n=0; n<9;n+=1){
                upar.Add((minu2.nameL(n)).c_str(), f1->GetParameter(n),1e-3); //,parmin(n), parlim(n), 0.1);
              }
	      std::cout<<"  Okay we have cleared the loop!"<<endl;
              // create Migrad minimizer
              //MnStrategy mnstra = MnStrategy(2);
              MnMigrad migrad(ffcn1, upar);
              //MnMigrad migrad(ffcn1, parm, err);
              //FunctionMinimum min = migrad();  //6000,1e-9);
              //try{
              FunctionMinimum min = migrad(40,1e-3);
                //throw logic_error("Assertion `s0.IsValid()' failed.");
              std::cout<<"minimum: "<<min<<std::endl;
              /*}
              catch(const std::logic_error & e){
                std::cerr<<" Found assertion error \n";
              }*/
              //MnHesse hesse;
              //FunctionMinimum min1 = migrad(3000,1e-9);
              for(int imn2=0; imn2<9; imn2+=1){
                param.push_back(migrad.Value((minu2.nameL(imn2)).c_str()));
                //std::cout<< "  par["<<i<<"] value --> ["<<param[i]<<"] \n";
              }
              for(int imn2=0; imn2<h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetNbinsX()+1; imn2+=1){
                double x=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(imn2);
                double yv=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinContent(imn2);
                double mnfit=minu2.model(x, param);
                double res=100*(yv-mnfit)/yv;
                h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(imn2, mnfit);
                h1Diff[indexClock][indexFB][indexUD][indexModule]->SetBinContent(imn2, res);
              }
              double mnx,mny,max,may;
              max=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
          	    GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
              if(max>=60 && max<=65){
	        if(!clus_csi)
                  clus_csi=true;
		if(!resetH) resetH=true;
                clock=indexClock;
                fb=indexFB;
                module=indexModule;
                mnx=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                    GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMinimumBin());
                may=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                    GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(max));
                mny=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                    GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(mnx));
                //Integrating fitting function for more accurate gain calibration
                TAxis* axis=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetXaxis();
                int xmin=0, xmax=250;
                int bmin=axis->FindBin(xmin);
                int bmax=axis->FindBin(xmax);
                int integral=h1Mnft[indexClock][indexFB][indexUD][indexModule]->Integral(bmin,bmax);
                integral-=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
                integral-=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(bmax)*(axis->GetBinLowEdge(bmin)-bmax)/axis->GetBinWidth(bmax);
                double area=integral-mny*250;
		std::cout<<" =====> Value for integral1 is:" <<area<<endl;
                calInt=area;
                double diff=(may-mny)*calib;
      	        indexph.push_back(diff);
      	        phval->push_back(diff);
      	        idCrys.push_back(indexModule), gfb.push_back(fb);
      	        int tAB=0;
      	        if((p[0]=='u' || p[0]=='U') && (indexModule>=9)) tAB=1;
      	        if((p[0]=='d' || p[0]=='D') && (indexModule<=8)) tAB=1;
      	        typeAB.push_back(tAB), gno.push_back(indexClock), gud.push_back(ud);
      	        double csitheta=theta[fb][indexModule], csiphi=phi[indexClock][ud][tAB];
	        // handeling phi for channel 16: Crystal type 10
	        if(treeRaw->indexCsI[i]==16){
	          csiphi=phiCh16[indexClock];
	          if(p[0]=='u' || p[0]=='U') csiphi=phiCh16[indexClock]+15.;
	        }
                cout<< " *********** theta, phi "<<csitheta<<"  "<<csiphi<<endl;
	        int thetaIndex=thetaCsI[moduleNo-1][treeRaw->indexChannel[i]];
	        int phiIndex=phiCsI[moduleNo-1][treeRaw->indexChannel[i]];
	        //otheta=mapTheta[thetaIndex];
	        //otheta=mapTheta[20-thetaIndex];
	        //reading from Osaka map as is: otheta and ophi angles
	        double cosTheta=crysZ[20-thetaIndex]/std::sqrt(std::pow(crysZ[20-thetaIndex],2)+
	          	      std::pow(crysr[20-thetaIndex],2));
	        double acosTheta=TMath::ACos(cosTheta);
	        otheta=TMath::RadToDeg()*acosTheta; // convert to deg.
	        otheta=round(otheta); // make sure to obtain 2 dp
	        mapPhi=2*180*(phiIndex-0.5)/48.;
	        ophi=mapPhi; // convert to global scheme: gap3=+xaxis
	        //ophi=2*180*(phiIndex-0.5)/48.;
	        //mapPhi=2*180*(phiIndex-0.5)/48.;
	        //ophi=90.-mapPhi; // convert to global scheme: gap3=+xaxis
	        //if(ophi < 0.) ophi=360.+ophi; // loop around the clock
                cout<< " *** Regular Angle "<<csitheta<<"  "<<csiphi<<endl;
                cout<< " ***  Osaka Angles "<<otheta<<"  "<<ophi<<endl;
      	        auto angles=std::make_pair(otheta,ophi);
      	        //thetaPhi.push_back(angles);
      	        csThet.push_back(otheta), csPhi.push_back(ophi);
		//h2ang->Fill(csitheta*M_PI/180,csiphi*M_PI/180);
		//h2deg->Fill(csitheta,csiphi);
      	        csiph[angles]=diff;
      	        //csiph[thetaPhi]=diff;
      	        csiClus[angles]=true;
	        h2Ang->Fill(otheta,ophi);
	        h2theta->Fill(otheta,csitheta);
	        h2phi->Fill(ophi,csiphi);
		// writing to check mapping file
                //mapfile<<moduleNo<<"\t"<<treeRaw->indexChannel[i]+1<<"\t"<<phiCsI[moduleNo-1][treeRaw->indexChannel[i]]<<"\t"<<thetaCsI[moduleNo-1][treeRaw->indexChannel[i]]<<"\t"<<nameModule<<"\t"<<indexClock+1<<"\t"<<p[0]<<"\t"<<p[1]<<std::endl;
      	        int imod=0, igap=4*indexClock+2;  // igapU=4*indexClock+3;
                int csimod=(-1*treeRaw->indexCsI[i])+1;  //default
		//indexClock=0;
                if(p[0]=='d' || p[0]=='D'){
      	          igap=4*(indexClock)+2;
                  csimod=(-1*treeRaw->indexCsI[i])+1;  //default
      	          if(treeRaw->indexCsI[i]>=10){
      	            igap=4*indexClock+1;
                    csimod=(-1*treeRaw->indexCsI[i])+7;
      	            if(treeRaw->indexCsI[i]==16) igap=4*indexClock+1.5;
      	          }
                  if(p[1]=='b' || p[1]=='B'){
      	            csimod=treeRaw->indexCsI[i];
      	            if(treeRaw->indexCsI[i]>=10){
      	              igap=4*indexClock+1;
      	              csimod=(treeRaw->indexCsI[i])-6;
      	              if(treeRaw->indexCsI[i]==16) igap=4*indexClock+1.5;
      	            }
      	          } 
		}
                if(p[0]=='u' || p[0]=='U'){
      	          igap=4*indexClock+3;
                  csimod=(-1*treeRaw->indexCsI[i])+1;  //default
      	          if(treeRaw->indexCsI[i]>=10){
      	            igap=4*(indexClock+1);
      	            csimod=(-1*treeRaw->indexCsI[i])+7;
      	            if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
      	          }
		  //FIXME
      	          if(p[1]=='b' || p[1]=='B'){
      	            csimod=treeRaw->indexCsI[i];
      	          igap=4*indexClock+3;
      	            if(treeRaw->indexCsI[i]>=10){
                      igap=4*(indexClock+1);
                      csimod=(treeRaw->indexCsI[i])-6;
      	              if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
	            }
      	          }
      	        } 
                h2clus->Fill(otheta,ophi,diff);
                //h2clus->Fill(csimod,igap,diff);
	        treeClus->waveID=nfound;
              } // <--- Use this to get rid of double and single fitting functions * /
            }
          }
        } // <--- End of No. CsI crystals that fired* /
      }  // <--- End of prelim. timing cut if loop
      //pCali->Fill();
      c1->cd();
      h2clus->GetZaxis()->SetRangeUser(0., .5);
      h2clus->Draw("colz");
      // Create a drawing of CsI outline using anlges
      for(int i=0;i<11;i++){
        hbox1[i]=new TLine(67.5,30*i+7.5,112.5,30*i+7.5); //lower gap edge in phi
        vbox1[i]=new TLine(67.5,30*i+22.5,67.5,30*i+37.5); // create left gap edge in theta
        hbox1[i+11]=new TLine(67.5,30*i+22.5,112.5,30*i+22.5); // upper gap edge in phi
        vbox1[i+11]=new TLine(112.5,30*i+22.5,112.5,30*i+37.5);
        hline1[i]=new TLine(15.0,30*i+15.0,165.5,30*i+15.0);
        hline1[i+11]=new TLine(15.0,30*i+15.0,165.5,30*i+15.0);
        hline2[i]=new TLine(15.0,30*i,67.5,30*i);
        hline2[i+11]=new TLine(15.0,30*i,67.5,30*i);
        hline3[i]=new TLine(112.5,30*i,165.5,30*i);
        hline3[i+11]=new TLine(112.5,30*i,165.5,30*i);
      }
      hline2[22]=new TLine(15.0,330.0,67.5,330.0);
      hline2[23]=new TLine(112.5,330.0,165.5,330.0);
      hline2[24]=new TLine(15.0,352.5,67.5,352.5);
      hline2[25]=new TLine(112.5,352.5,165.5,352.5);
      hbox2[0]=new TLine(67.5,337.5,112.5,337.5);
      hbox2[1]=new TLine(67.5,352.5,112.5,352.5);
      vbox2[0]=new TLine(67.5,0.0,67.5,7.5);
      vbox2[1]=new TLine(112.5,0.0,112.5,7.5);
      vbox2[2]=new TLine(67.5,352.5,67.5,360.0);
      vbox2[3]=new TLine(112.5,352.5,112.5,360.0);
      hline1[22]=new TLine(15.0,345.0,165.5,345.0);
      for(int n=0;n<22;n++){
        hbox1[n]->Draw("l");
        vbox1[n]->Draw("l");
        hline1[n]->Draw("l");
        hline2[n]->Draw("l");
        hline3[n]->Draw("l");
      }
      for(int n=0;n<2;n++){
        hbox2[n]->Draw("l");
      }
      for(int n=0;n<4;n++){
        vbox2[n]->Draw("l");
      }
      hline1[22]->Draw("l");
      hline2[22]->Draw("l");
      hline2[23]->Draw("l");
      vbox2[0]->Draw("l");
      c1->Modified();
      c1->Update();
      empty();
    } // <--- End of if loop
  } // <--- End of number of crystals that fired loop
  //treeRaw->Clear();

  if(clus_csi){
    std::cout<<"\n ***************************************************************************\n";
	  std::cout<<" -------------->  || "<<tracktree->pVertpi0<<" || <------------------\n";
    for(UInt_t cr=0;cr<phval->size();cr++){
      std::cout<<" angles of corr. pulse heigh["<<cr<<"] = "<<(*phval)[cr]<<endl;
    }/*
    for(auto itrr=csiClus.begin();itrr!=csiClus.end();itrr ++){
      std::cout<<"  cluster bool has the following: "<<itrr->second<<std::endl;
    }*/
    double ntheta, nphi;
    int numOfClus=0, numOfsingleClus=0;
    for(std::size_t mm=0; mm !=indexph.size(); mm++){
      const auto index=get_nthIndex(indexph, mm);
      std::cout<<"  the greater index --> "<<index <<std::endl;
              //<< " with value "<<indexph[index]<<std::endl;
      //std::nth_element(begin(hcrys), begin(hcrys)+ii, end(hcrys));
      ntheta=csThet[index], nphi=csPhi[index];
      //erpair.clear();
      auto tppair=std::make_pair(ntheta,nphi);
      std::cout<<"   ==>  theta and phi "<<ntheta<<"  "<<nphi<<endl;
      auto angP1=std::make_pair(ntheta+7.5,nphi);      auto angP2=std::make_pair(ntheta-7.5,nphi);
      auto angP3=std::make_pair(ntheta,nphi+7.5);      auto angP4=std::make_pair(ntheta,nphi-7.5);
      auto angP5=std::make_pair(ntheta+7.5,nphi+7.5);  auto angP6=std::make_pair(ntheta-7.5,nphi-7.5);
      auto angP7=std::make_pair(ntheta-7.5,nphi+7.5);  auto angP8=std::make_pair(ntheta+7.5,nphi-7.5);
      // accounting for the complicated case of Crystal No. 10 (index=16)
      // This is strange because they are 2X larger in phi
      std::pair<double,double> angE1, angE2, angE3, angE4, angE5;
      std::pair<double,double> angE6, angE7, angE8, angE9, angE10;
      // from chan16 -> chan15: Back
      if(ntheta==161.25){
        angE1=std::make_pair(ntheta-7.5,nphi-7.50);  angE2=std::make_pair(ntheta-7.5,nphi+15.0);
        angE3=std::make_pair(ntheta-7.5,nphi);       angE4=std::make_pair(ntheta-7.5,nphi+7.50);
        // only phi angle of chan16 changes
        angE9=std::make_pair(ntheta, nphi+15);       angE10=std::make_pair(ntheta,nphi-15);
      }
      if(ntheta==18.75){// from chan16 -> chan15: Front
        angE5=std::make_pair(ntheta+7.5,nphi-7.50);  angE6=std::make_pair(ntheta+7.5,nphi+15.0);
        angE7=std::make_pair(ntheta+7.5,nphi);       angE8=std::make_pair(ntheta+7.5,nphi+7.5);
        // only phi angle of chan16 changes
        angE9=std::make_pair(ntheta, nphi+15);       angE10=std::make_pair(ntheta,nphi-15);
      }
      // handeling of going from chan15->chan16: Back
      if((ntheta+7.5)==161.25){
        angE1=std::make_pair(ntheta+7.5,nphi-7.50);  angE2=std::make_pair(ntheta+7.5,nphi+7.50);
        angE3=std::make_pair(ntheta+7.5,nphi);       //angE4=std::make_pair(ntheta+7.5,nphi+3.75);
      }
      // from chan15 -> chan16: Front
      if((ntheta-7.5)==18.75){
        angE5=std::make_pair(ntheta-7.5,nphi-7.50);  angE6=std::make_pair(ntheta-7.5,nphi+7.50);
        angE7=std::make_pair(ntheta-7.5,nphi);       //angE8=std::make_pair(ntheta-7.5,nphi+3.75);
      }
      // dealing with the revolution case
      std::pair<double,double> angR1, angR2, angR3, angR4, angR5, angR6;
      if(nphi==348.75 || nphi==356.25){
        // only phi angle of chan16 changes
        angR1=std::make_pair(ntheta, 3.75);       angR2=std::make_pair(ntheta+7.5,3.75);
        angR3=std::make_pair(ntheta-7.5, 3.75);   //angR4=std::make_pair(ntheta+7.5,3.75);
      }
      if(nphi==3.75){
        // only phi angle of chan16 changes
        angR1=std::make_pair(ntheta, 348.75);       angR2=std::make_pair(ntheta,356.25);
        angR3=std::make_pair(ntheta-7.5, 348.75);   angR4=std::make_pair(ntheta+7.5,356.25);
        angR5=std::make_pair(ntheta-7.5, 356.25);   angR6=std::make_pair(ntheta+7.5,348.75);
      }
      clusCrys=0;
      Eclus=0.;
      thetaE=0; phiE=0;
      if(csiph[angP1] > 0 &&  csiph[angP2] > 0 &&  csiph[angP3] > 0 && csiph[angP4] > 0 &&
        csiph[angP5] > 0 && csiph[angP6] > 0 && csiph[angP7] > 0 &&
        csiph[angP8] > 0){
        std::cout<<"  Total number of cluster crystals is 8 \n";
      }else if(csiph[angP1]>0 || csiph[angP2]>0 || csiph[angP3]>0 || csiph[angP4]>0 || 
        csiph[angP5]>0 || csiph[angP6]>0 || csiph[angP7]>0 || csiph[angP8]>0 ||
	csiph[angE1]>0 || csiph[angE2]>0 || csiph[angE3]>0 || csiph[angE4]>0 ||
	csiph[angE5]>0 || csiph[angE6]>0 || csiph[angE7]>0 || csiph[angE8]>0 ||
	csiph[angE9]>0 || csiph[angE10]>0 || csiph[angR1]>0 || csiph[angR2]>0){
        //clusCrys=clusCrys+1;
        if(csiph[angP1]>0){
          std::cout<<" This crystal Cluster pulse-height P1: "<<csiph[angP1];
	  std::cout<<" ["<<std::get<0>(angP1)<<", "<<std::get<1>(angP1)<<"] \n";
          if(csiClus[angP1]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angP1];
	    thetaE=thetaE+csiph[angP1]*(std::get<0>(angP1));
	    phiE  =phiE  +csiph[angP1]*(std::get<1>(angP1));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angP1]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angP2]>0){
          std::cout<<" This crystal Cluster finder in pair loop P2: "<<csiph[angP2];
	  std::cout<<" ["<<std::get<0>(angP2)<<", "<<std::get<1>(angP2)<<"] \n";
          if(csiClus[angP2]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angP2];
	    thetaE=thetaE+csiph[angP2]*(std::get<0>(angP2));
	    phiE  =phiE  +csiph[angP2]*(std::get<1>(angP2));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angP2]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angP3]>0){
          std::cout<<" This crystal Cluster finder in pair loop P3: "<<csiph[angP3];
	  std::cout<<" ["<<std::get<0>(angP3)<<", "<<std::get<1>(angP3)<<"] \n";
          if(csiClus[angP3]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angP3];
	    thetaE=thetaE+csiph[angP3]*(std::get<0>(angP3));
	    phiE  =phiE  +csiph[angP3]*(std::get<1>(angP3));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angP3]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angP4]>0){
          std::cout<<" This crystal Cluster finder in pair loop P4: "<<csiph[angP4];
	  std::cout<<" ["<<std::get<0>(angP4)<<", "<<std::get<1>(angP4)<<"] \n";
          if(csiClus[angP4]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angP4];
	    thetaE=thetaE+csiph[angP4]*(std::get<0>(angP4));
	    phiE  =phiE  +csiph[angP4]*(std::get<1>(angP4));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angP4]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angP5]>0){
          std::cout<<" This crystal Cluster finder in pair loop P5: "<<csiph[angP5];
	  std::cout<<" ["<<std::get<0>(angP5)<<", "<<std::get<1>(angP5)<<"] \n";
          if(csiClus[angP5]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angP5];
	    thetaE=thetaE+csiph[angP5]*(std::get<0>(angP5));
	    phiE  =phiE  +csiph[angP5]*(std::get<1>(angP5));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angP5]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angP6]>0){
          std::cout<<" This crystal Cluster finder in pair loop P6: "<<csiph[angP6];
	  std::cout<<" ["<<std::get<0>(angP6)<<", "<<std::get<1>(angP6)<<"] \n";
          if(csiClus[angP6]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angP6];
	    thetaE=thetaE+csiph[angP6]*(std::get<0>(angP6));
	    phiE  =phiE  +csiph[angP6]*(std::get<1>(angP6));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angP6]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angP7]>0){
          std::cout<<" This crystal Cluster finder in pair loop P7: "<<csiph[angP7];
	  std::cout<<" ["<<std::get<0>(angP7)<<", "<<std::get<1>(angP7)<<"] \n";
          if(csiClus[angP7]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angP7];
	    thetaE=thetaE+csiph[angP7]*(std::get<0>(angP7));
	    phiE  =phiE  +csiph[angP7]*(std::get<1>(angP7));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angP7]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angP8]>0){
          std::cout<<" This crystal Cluster finder in pair loop P8: "<<csiph[angP8];
	  std::cout<<" ["<<std::get<0>(angP8)<<", "<<std::get<1>(angP8)<<"] \n";
          if(csiClus[angP8]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angP8];
	    thetaE=thetaE+csiph[angP8]*(std::get<0>(angP8));
	    phiE  =phiE  +csiph[angP8]*(std::get<1>(angP8));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angP8]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angE1]>0){
          std::cout<<" Edge effect Cluster finder in pair loop E1: "<<csiph[angE1];
	  std::cout<<" ["<<std::get<0>(angE1)<<", "<<std::get<1>(angE1)<<"] \n";
          if(csiClus[angE1]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angE1];
	    thetaE=thetaE+csiph[angE1]*(std::get<0>(angE1));
	    phiE  =phiE  +csiph[angE1]*(std::get<1>(angE1));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angE1]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angE2]>0){
          std::cout<<" Edge effect Cluster finder in pair loop E2: "<<csiph[angE2];
	  std::cout<<" ["<<std::get<0>(angE2)<<", "<<std::get<1>(angE2)<<"] \n";
          if(csiClus[angE2]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angE2];
	    thetaE=thetaE+csiph[angE2]*(std::get<0>(angE2));
	    phiE  =phiE  +csiph[angE2]*(std::get<1>(angE2));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angE2]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angE3]>0){
          std::cout<<" Edge effect Cluster finder in pair loop E3: "<<csiph[angE3];
	  std::cout<<" ["<<std::get<0>(angE3)<<", "<<std::get<1>(angE3)<<"] \n";
          if(csiClus[angE3]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angE3];
	    thetaE=thetaE+csiph[angE3]*(std::get<0>(angE3));
	    phiE  =phiE  +csiph[angE3]*(std::get<1>(angE3));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angE3]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
	// Got rid of angE4 and angE8 pairs because they are not needed at this time
	// Dongwi: 07.14.2019
        if(csiph[angE5]>0){
          std::cout<<" Edge effect Cluster finder in pair loop E5: "<<csiph[angE5];
	  std::cout<<" ["<<std::get<0>(angE5)<<", "<<std::get<1>(angE5)<<"] \n";
          if(csiClus[angE5]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angE5];
	    thetaE=thetaE+csiph[angE5]*(std::get<0>(angE5));
	    phiE  =phiE  +csiph[angE5]*(std::get<1>(angE5));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angE5]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angE6]>0){
          std::cout<<" Edge effect Cluster finder in pair loop E6: "<<csiph[angE6];
	  std::cout<<" ["<<std::get<0>(angE6)<<", "<<std::get<1>(angE6)<<"] \n";
          if(csiClus[angE6]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angE6];
	    thetaE=thetaE+csiph[angE6]*(std::get<0>(angE6));
	    phiE  =phiE  +csiph[angE6]*(std::get<1>(angE6));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angE6]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angE7]>0){
          std::cout<<" Edge effect Cluster finder in pair loop E7: "<<csiph[angE7];
	  std::cout<<" ["<<std::get<0>(angE7)<<", "<<std::get<1>(angE7)<<"] \n";
          if(csiClus[angE7]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angE7];
	    thetaE=thetaE+csiph[angE7]*(std::get<0>(angE7));
	    phiE  =phiE  +csiph[angE7]*(std::get<1>(angE7));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angE7]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angE9]>0){
          std::cout<<" Edge effect Cluster finder in pair loop E9: "<<csiph[angE9];
	  std::cout<<" ["<<std::get<0>(angE9)<<", "<<std::get<1>(angE9)<<"] \n";
          if(csiClus[angE9]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angE9];
	    thetaE=thetaE+csiph[angE9]*(std::get<0>(angE9));
	    phiE  =phiE  +csiph[angE9]*(std::get<1>(angE9));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angE9]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angE10]>0){
          std::cout<<" Edge effect Cluster finder in pair loop E10: "<<csiph[angE10];
	  std::cout<<" ["<<std::get<0>(angE10)<<", "<<std::get<1>(angE10)<<"] \n";
          if(csiClus[angE10]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angE10];
	    thetaE=thetaE+csiph[angE10]*(std::get<0>(angE10));
	    phiE  =phiE  +csiph[angE10]*(std::get<1>(angE10));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angE10]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
	// revolution around CsI(Tl) barrel
        if(csiph[angR1]>0){
          std::cout<<" Revolution case Cluster pulse-height R1: "<<csiph[angR1];
	  std::cout<<" ["<<std::get<0>(angR1)<<", "<<std::get<1>(angR1)<<"] \n";
          if(csiClus[angR1]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angR1];
	    thetaE=thetaE+csiph[angR1]*(std::get<0>(angR1));
	    phiE  =phiE  +csiph[angR1]*(std::get<1>(angR1));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angR1]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angR2]>0){
          std::cout<<" Revolution case Cluster finder in pair loop R2: "<<csiph[angR2];
	  std::cout<<" ["<<std::get<0>(angR2)<<", "<<std::get<1>(angR2)<<"] \n";
          if(csiClus[angR2]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angR2];
	    thetaE=thetaE+csiph[angR2]*(std::get<0>(angR2));
	    phiE  =phiE  +csiph[angR2]*(std::get<1>(angR2));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angR2]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angR3]>0){
          std::cout<<" Edge effect Cluster finder in pair loop R3: "<<csiph[angR3];
	  std::cout<<" ["<<std::get<0>(angR3)<<", "<<std::get<1>(angR3)<<"] \n";
          if(csiClus[angR3]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angR3];
	    thetaE=thetaE+csiph[angR3]*(std::get<0>(angR3));
	    phiE  =phiE  +csiph[angR3]*(std::get<1>(angR3));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angR3]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angR4]>0){
          std::cout<<" Edge effect Cluster finder in pair loop R4: "<<csiph[angR4];
	  std::cout<<" ["<<std::get<0>(angR4)<<", "<<std::get<1>(angR4)<<"] \n";
          if(csiClus[angR4]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angR4];
	    thetaE=thetaE+csiph[angR4]*(std::get<0>(angR4));
	    phiE  =phiE  +csiph[angR4]*(std::get<1>(angR4));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angR4]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angR5]>0){
          std::cout<<" Edge effect Cluster finder in pair loop R5: "<<csiph[angR5];
	  std::cout<<" ["<<std::get<0>(angR5)<<", "<<std::get<1>(angR5)<<"] \n";
          if(csiClus[angR5]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angR5];
	    thetaE=thetaE+csiph[angR5]*(std::get<0>(angR5));
	    phiE  =phiE  +csiph[angR5]*(std::get<1>(angR5));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angR5]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        if(csiph[angR6]>0){
          std::cout<<" Edge effect Cluster finder in pair loop R6: "<<csiph[angR6];
	  std::cout<<" ["<<std::get<0>(angR6)<<", "<<std::get<1>(angR6)<<"] \n";
          if(csiClus[angR6]){
            clusCrys=clusCrys+1;
	    Eclus=Eclus+csiph[angR6];
	    thetaE=thetaE+csiph[angR6]*(std::get<0>(angR6));
	    phiE  =phiE  +csiph[angR6]*(std::get<1>(angR6));
            std::cout<<" This crystal is now removed from the list: "<<std::endl;
            csiClus[angR6]=false;
          }else{
            std::cout<<" Already checked this crystal \n";
          }
        }
        std::cout<<" some crystals actually have hits "<<clusCrys<<std::endl;
      }else{
        clusCrys=0;//clusCrys+1;
        std::cout<<" Single cluster crystals here \n"; //<<clusCrys<<std::endl;
      }
      if(csiClus[tppair]){
        clusCrys=clusCrys+1;
	Eclus=Eclus+csiph[tppair];
	thetaE=thetaE+csiph[tppair]*(std::get<0>(tppair));
	phiE  =phiE  +csiph[tppair]*(std::get<1>(tppair));
	std::cout<<" -->  pulse-heignt for central crystal: "<<csiph[tppair];
	std::cout<<" ["<<std::get<0>(tppair)<<", "<<std::get<1>(tppair)<<"] \n";
	std::cout<<" >>>  Cluster energy is ------------->: "<<Eclus<<" [GeV]";
      }
      if(clusCrys>=2){
        numOfClus++;
	// perform energy-weighting and convert from deg-->rad
	rtheta=TMath::DegToRad()*(thetaE/Eclus);
	rphi=TMath::DegToRad()*(phiE/Eclus);
	// Fill the theta, phi distributions in rad
	treeClus->thetaE=rtheta;
	treeClus->phiE=rphi;
	clusEne.push_back(Eclus);
	treeClus->Ncrys=clusCrys;
	treeClus->ClustCrys=clusCrys;
	clusThetaE.push_back(rtheta);
	clusPhiE.push_back(rphi);
        h2ang->Fill(rtheta,rphi);
        h2deg->Fill(rtheta*180./M_PI,rphi*180/M_PI);
      }
      if(clusCrys==1){
        numOfsingleClus++;
        h1sclus->Fill(1);
	treeClus->Ncrys=1;
	singleEne.push_back(Eclus);
      }
      csiClus[tppair]=false; // mute central crystal
      cout<<"  Number of crystals is:  "<<clusCrys<<endl;
      std::cout<<"\n --------------------------------------------------------------------->\n\n";
    }
    for(UInt_t idc=0;idc<csThet.size();idc++){
      std::cout<<"    theta[m][n] "<<csThet[idc]<<endl;
      //cout<<"    pairs pHeig "<<csiph.find(thetaPhi)->second<<endl;
    }
    /***********************************************
     *   Vertex state vector info. for comparison
     * *********************************************/
    // momentum/Energy in GeV/c (c=1)
    ppip=tracktree->pVertpi0;
    // vertex momentum direction for pi+
    piPpx=ppip*(tracktree->nxVert);
    piPpy=ppip*(tracktree->nyVert);
    piPpz=ppip*(tracktree->nzVert);
    // vertex pi+ position
    pi0x=tracktree->xVert;
    pi0y=tracktree->yVert;
    pi0z=tracktree->zVert;
    T_pi0=std::sqrt(std::pow(ppip,2)+std::pow(M_pi0,2));//-M_pi0;
    pipEtot=std::sqrt(std::pow(ppip,2)+std::pow(M_piP,2));//-M_pi0;
    std::cout<<"\n ----------  pi0 E_tot = "<<T_pi0<<" ---------------\n";
      std::cout<<"\n  Checking cos(theta)s:      "<<std::cos(2*3.142)<<endl;
    if(numOfClus==2){
      // calculate 3-momentum direction for pi0: from (theta,phi) of 2*gamma
      g1px=clusEne[0]*std::sin(clusThetaE[0])*std::cos(clusPhiE[0]);
      g1py=clusEne[0]*std::sin(clusThetaE[0])*std::sin(clusPhiE[0]);
      g1pz=clusEne[0]*std::cos(clusThetaE[0]);
      g2px=clusEne[1]*std::sin(clusThetaE[1])*std::cos(clusPhiE[1]);
      g2py=clusEne[1]*std::sin(clusThetaE[1])*std::sin(clusPhiE[1]);
      g2pz=clusEne[1]*std::cos(clusThetaE[1]);
      // calculate pi0 invariant mass from above info.
      TLorentzVector gamma1;
      TLorentzVector gamma2;
      gamma1.SetPxPyPzE(g1px, g1py, g1pz, clusEne[0]);
      gamma2.SetPxPyPzE(g2px, g2py, g2pz, clusEne[1]);
      TLorentzVector pi0=gamma1+gamma2;
      // K+ Lorentz vector info. from pi+ and pi0
      TLorentzVector piPl;
      piPl.SetPxPyPzE(piPpx, piPpy, piPpz,pipEtot);
      TLorentzVector kaon=piPl+pi0;
      //ThreeVector for angular analysis
      TVector3 piPv(piPpx, piPpy,piPpz);
      TVector3 pi0v=-1*piPv; //(pi0.Px(), pi0.Py(),pi0.Pz());
      TVector3 gv1(g1px, g1py, g1pz);
      TVector3 gv2(g2px, g2py, g2pz);
      // Fill histos
      E2g->Fill(clusEne[0]+clusEne[1]);
      pi0Etot->Fill(T_pi0);
      h1Mpi0->Fill(pi0.M());
      h1Mpi02->Fill(pi0.M2());
      h1pi0px->Fill(piPpx);
      h1pi0py->Fill(piPpy);
      h1pi0pz->Fill(piPpz);
      //vertex info.
      h1vertpx->Fill(pi0.Px());
      h1vertpy->Fill(pi0.Py());
      h1vertpz->Fill(pi0.Pz());
      vertOp->Fill(std::cos(piPv.Angle(pi0v)));
      if(clusEne[0]+clusEne[1] >=.225 && clusEne[0]+clusEne[1]<=.252) E_cut->Fill(clusEne[0]+clusEne[1]);
      if((clusEne[0]-clusEne[1])/(clusEne[0]+clusEne[1])<=0.1) cosTheta->Fill(gv1.Angle(gv2)*(180./M_PI));
      h2Ene->Fill(clusEne[0]+clusEne[1]+pipEtot,pipEtot+T_pi0);
      // Fill tree var
      treeClus->E_pi0=clusEne[0]+clusEne[1];
      treeClus->g1Px=g1px;       treeClus->g2Px=g2px;      treeClus->pi0vx=pi0v.X();
      treeClus->g1Py=g1py;       treeClus->g2Py=g2py;      treeClus->pi0vy=pi0v.Y();
      treeClus->g1Pz=g1pz;       treeClus->g2Pz=g2pz;      treeClus->pi0vz=pi0v.Z();
      treeClus->piPpx=piPpx;     treeClus->pi0px=pi0.Px();
      treeClus->piPpy=piPpy;     treeClus->pi0py=pi0.Py();
      treeClus->piPpz=piPpz;     treeClus->pi0pz=pi0.Pz();
      treeClus->ggCosTheta=pi0.CosTheta();
      treeClus->piCosTheta=pi0.CosTheta();
      treeClus->piP2g=clusEne[0]+clusEne[1]+pipEtot;
      treeClus->piPpi0=T_pi0+pipEtot;
      treeClus->M_pi0=pi0.M();
      treeClus->pi0M2=pi0.M2();
      treeClus->M_k=kaon.M();
      treeClus->kM2=kaon.M2();
      treeClus->g1E=clusEne[0];
      treeClus->g2E=clusEne[1];
      std::cout<<"\n  piPecking total Cluster Energy:  "<<clusEne[0]+clusEne[1]<<endl;
      std::cout<<"\n  Angular1 checking (centriod)   ("<<clusThetaE[0]<<", "<<clusPhiE[0]<<")\n";
      std::cout<<"\n  Angular2 checking (centriod)   ("<<clusThetaE[1]<<", "<<clusPhiE[1]<<")\n";
      std::cout<<"\n  Checking pi0 InvMass:      "<<pi0.M()<<endl;
      std::cout<<"\n  Checking cos(theta):       "<<pi0.CosTheta()<<endl;
      std::cout<<"\n  Checking vertex opening    "<<piPv.Angle(pi0v)<<endl;
    }
    if(numOfClus>0){
      h1clust->Fill(numOfClus);
      treeClus->clusterM=numOfClus;
    }
    treeClus->channel=(numOfClus+numOfsingleClus);
    std::cout<<"\n\n  Number of clusters is   :  "<<numOfClus<<endl;
    std::cout<<"  Number of single clusters is:  "<<numOfsingleClus<<endl;
    std::cout<<" ***************************************************************************\n";
  }//<--- end cluster if loop* /
exitLoop:
  delete phval;
  delete clusth;
  delete clusphi;

  return 0;
}

void Det_ClusterCsI::empty(){
  for(int i=0;i<23;i++) delete hline1[i];
  for(int i=0;i<26;i++) delete hline2[i];
  for(int i=0;i<4;i++) delete vbox2[i];
  for(int i=0;i<2;i++){
    delete hbox2[i];
  }
  for(int i=0;i<22;i++){
    delete hbox1[i];
    delete vbox1[i];
    delete hline3[i];
  }
}

Long_t Det_ClusterCsI::finalize(){
  //mapfile.close();

  return 0;
}

Long_t Det_ClusterCsI::cmdline(char *cmd){
  //add cmdline hanling here

  return 0; // 0 = all ok
};

extern "C"{
  Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p){
    return (Plugin *) new Det_ClusterCsI(in,out,inf_,outf_,p);
  }
}
ClassImp(Det_ClusterCsI);
