#include <Det_CsI.h>
#include <mn2CsIfn.h>
#include "TMath.h"
#include "TF1.h"
#include<iostream>
#include<cmath>
mn2CsIfn mn2;
Det_CsI::Det_CsI(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p){
  // Set defaults for various options
  treeCali=0;
  h1Mnft[12][2][2][16]=NULL;
  h1Diff[12][2][2][16]=NULL;
  h1rate[12][2][2][16]=NULL; // normalized diff b/n fit and histogram
  h2Fits[12][2][2][16]=NULL;
  h1Fits[12][2][2][16]=NULL;
  h1Amps[12][2][2][16]=NULL;
  h1time[12][2][2][16]=NULL;
  std::cout<<" checking this shit \n";
  loopX=false, notfire=false;
  firedCsI=false; csiT=false;
  dummy=-1000;
  xpeaks=NULL;
  name<<"stat_"; name2<<"Mnfit_"; name3<<"F1fit_"; name4<<"Dhfit_"; name5<<"pHeight", name6<<"Diff";
  tname<<"time";
};

Det_CsI::~Det_CsI(){
  delete s;
  delete treeSing;
  cout<<"  Exiting fitting program \n";
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<16;iModule++){
          //delete h2Fits[iClock][iFB][iUD][iModule];
          delete h1time[iClock][iFB][iUD][iModule];
          delete h1Fits[iClock][iFB][iUD][iModule];
          delete h1Amps[iClock][iFB][iUD][iModule];
          delete h1Mnft[iClock][iFB][iUD][iModule];
          delete h1Diff[iClock][iFB][iUD][iModule];
        }
      }
    }
  }
};


Long_t Det_CsI::cmdline(char *cmd){
  //add cmdline hanling here

  return 0; // 0 = all ok
};
std::string singleFit = mn2.singlemodel();   std::string quadFit= mn2.quadruplemodel();
std::string doubleFit = mn2.doublemodel();   std::string overr  = mn2.overrangemodel();
std::string tripleFit = mn2.triplemodel();   std::string pileUp = mn2.triplemodel();
Long_t Det_CsI::histos(){
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<16;iModule++){
          //std::ostringstream name, name2, name3, name4, name5, name6, tname;
          name2<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name3<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name4<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name5<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name6<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          tname<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          h1time[iClock][iFB][iUD][iModule]=new TH1D(tname.str().c_str(),"stat",250,0,250);
          h1Amps[iClock][iFB][iUD][iModule]=new TH1D(name5.str().c_str(),"stat",250,0,250);
          h1Mnft[iClock][iFB][iUD][iModule]=new TH1D(name2.str().c_str(),"stat",250,0,250);
          h1Diff[iClock][iFB][iUD][iModule]=new TH1D(name6.str().c_str(),"stat",250,0,250);
        }
      }
    }
  }
  // parameter specific histos
  s = new TSpectrum(4);
  p0=dH1("p0","Parameter p0", 250, 0, 1100);
  p1=dH1("p1","Parameter p1", 250, 0, 250);
  p9=dH1("p9","Parameter p9", 250, 0, 300);
  p10=dH1("p10","Parameter p10", 200, 0, 800);
  h2clus=dH2("hclust","Clusters in the CsI(Tl)  ", 50,-25,25,27,0,54);
  h1Pamp=dH1("hpulse","Pulse height distribution", 250, 0, 1000);
  h1kmu2=dH1("kmu2DP","Pulse height distribution", 250, 0, 1000);
  h1Intg=dH1("Integr","Integrated pulse height distribution", 250, 0, 100000);
  //h1cali=dH1("Calibr","Integrated pulse height distribution", 250, 0, 1000);
  h1ped=dH1("Ped","Pedestals for the waveform ", 250, 0, 1000);

  return 0;
}

Long_t Det_CsI::startup(){
  getBranchObject("vf48",(TObject **) &treeRaw);
  getBranchObject("RawBeamInfo",(TObject **) &treeBeam);
  treeSing=new CRTSingleCsI();
  makeBranch("treeSing",(TObject **) &treeSing);
  gStyle->SetOptStat(0);

  return 0;
}

//Initialize storage variables here
void Det_CsI::initVar(){
  treeSing->thSing=dummy;
  treeSing->indexCsI=dummy;      treeSing->phiSing=dummy;
  treeSing->tpeak=dummy;         
  treeSing->trise=dummy;         treeSing->typeAB=dummy;
  treeSing->calInt=dummy;        treeSing->crysID=dummy;
  treeSing->csiArrange[0]=dummy; treeSing->fb=dummy;
  treeSing->csiArrange[1]=dummy; treeSing->ped=dummy;
  treeSing->clock=dummy;
  treeSing->ovrped=dummy;
  treeSing->ovrpH=dummy;
  treeSing->ud=dummy;           
  treeSing->phei=dummy;         
  treeSing->phdstr=dummy;
  //Single pulse                 //double pulse
  treeSing->sphei=dummy;         treeSing->kmu2=dummy;
  treeSing->sptime=dummy;        treeSing->dubPed=dummy;
  treeSing->sped=dummy;          treeSing->dubphei=dummy;
  treeSing->waveID=dummy;        treeSing->intKmu2=dummy;
  for(int i=0;i<3;i++){
    treeSing->tref[i]=dummy;
    treeSing->refpk[i]=dummy;
    treeSing->tcorr[i]=dummy;
    treeSing->refmn[i]=dummy;
    treeSing->rgaus[i]=dummy;
  }
}

//function to get timing from ref. module
std::string Det_CsI::refT(){
  char reft[1000], timing[1000];
  sprintf(reft,"1-exp(-(x-[1])/[2])");
  sprintf(timing,"[0]*(x-[1])/(%s)+[3]",reft);
  return timing;
}
Long_t Det_CsI::process(){
  initVar(); //initialize storage variables
  for(UInt_t i=0;i<treeRaw->nChannel;i++){
    char* p=(char*)&(treeRaw->nameModule[i]);
    int moduleName=(p[3]-'0')*10+(p[2]-'0')-1;
    std::string nameModule;
    nameModule+=(*p);
    p++;
    nameModule+=*p;
    p++;
    nameModule+=*p;
    p++;
    nameModule+=*p;
    std::string nameCsI;
    p=(char*)&(treeRaw->nameCsI[i]);
    int indexClock=(p[3]-'0')*10+(p[2]-'0')-1;
    //std::cout<< "\n Index clock 11: "<<indexClock<<endl;
    p+=3;
    nameCsI+=(*p);
    p--;
    nameCsI+=*p;
    p--;
    nameCsI+=*p;
    p--;
    nameCsI+=*p;

    //std::cout<< "  ***** THIS IS A TEST! TESTING! TESTING! "<<nameCsI<<"\t"<<nameModule<<endl;
    int indexModule=treeRaw->indexCsI[i]-1;
    int indexFB=0;
    if(p[1]=='b' || p[1]=='B') indexFB=1;
    int indexUD=0;
    if(p[0]=='d' || p[0]=='D') indexUD=1;
    // since this is a calibration plugin,
    // we should only consider single crystal hits
    if(treeRaw->nChannel>7) goto exitLoop;
    name<<indexClock<<"_"<<indexFB<<"_"<<indexUD<<"_"<<indexModule;
    h1Fits[indexClock][indexFB][indexUD][indexModule]=new TH1D(name.str().c_str(),"stat",250,0,250);
/*
    // reference timing from 3 modules
    // timing from all 3 modules will considered
    if((treeRaw->indexCsI[i]==16 && indexFB==0 && indexUD==0) && 
		    (indexClock==0 || indexClock==4 || indexClock==8)){
      name<<indexClock<<"_"<<indexFB<<"_"<<indexUD<<"_"<<indexModule;
      h1Fits[indexClock][indexFB][indexUD][indexModule]=new TH1D(name.str().c_str(),"stat",250,0,250);
      for(UInt_t iData=0;iData<treeRaw->nSample[i];iData++){
        h1Fits[indexClock][indexFB][indexUD][indexModule]->SetBinContent(iData+1,treeRaw->data[i][iData]);
      }
      x1=h1Fits[indexClock][indexFB][indexUD][indexModule]->
          GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
      // use switch statement for various timing modules
      // all 3 will be use and the timing analysis will be conducted accordingly
      lowRange=x1-6; upRange=x1+6;
      TF1* f1=new TF1("f1","gaus",lowRange,upRange);
      TF1* f2=new TF1("f2",refT().c_str(),0,50);
      switch(indexClock){
        case 0:
          f2->SetParameters(21.3,29.6,1.85,120);
          h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f2,"QR+");
          h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f1,"QR+");
          maxfn[0]=f1->GetMaximum();
          minfn[0]=f2->GetMinimum();
          cf50[0]=.5*maxfn[0]+.5*minfn[0];
          // Fill tree variables here:
          T_ref[0]=f2->GetX(cf50[0]);
          refgaus[0]=f1->GetX(maxfn[0]);
          treeSing->indexCsI=treeRaw->indexCsI[i];
          treeSing->phei=dummy;
          treeSing->rgaus[0]=refgaus[0];
          treeSing->csiArrange[0]=p[0];
          treeSing->csiArrange[1]=p[1];
          treeSing->clock=indexClock+1;
          treeSing->ud=indexUD;         treeSing->fb=indexFB;
          treeSing->tref[0]=T_ref[0];   treeSing->refpk[0]=maxfn[0];
          treeSing->refmn[0]=minfn[0];
          //std::cout<<"\n Event number is:  "<<treeRaw->eventNo<<std::endl; 
          //std::cout<< " Index clock: "<<indexClock<<std::endl;
          //std::cout<< " Gap config FB is  : " <<p[1]<<std::endl;
          //std::cout<< " Gap config UD is  : " <<p[0]<<std::endl;
          //std::cout<< " size of nChannel is : " <<treeRaw->indexCsI[i]-1<<std::endl;
          //std::cout<< " size of nSample is  : " <<treeRaw->nSample[i]<<std::endl;
          //std::cout<< " Chan No.: "<<treeRaw->nChannel<<std::endl;
          //std::cout<<"\n ---------------------------------------------------------\n";
          //std::cout<<" \n\n  ------>ref time and peak time(1): "<<T_ref[0]<<" "<<refgaus[0]<<"\n";
          //std::cout<<" \n\n  ------> CDF timing:  "<<(valx2-valx1)<<" \n\n";
          delete f1; delete f2;
          break;
        case 4:
          f2->SetParameters(21.3,29.6,1.85,120);
          h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f2,"QR+");
          h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f1,"QR+");
          maxfn[1]=f1->GetMaximum();
          minfn[1]=f2->GetMinimum();
          cf50[1]=.5*maxfn[1]+.5*minfn[1];
          // Fill tree variables here:
          T_ref[1]=f2->GetX(cf50[1]);
          refgaus[1]=f1->GetX(maxfn[1]);
          treeSing->indexCsI=treeRaw->indexCsI[i];
          treeSing->phei=dummy;
          treeSing->rgaus[1]=refgaus[1];
          treeSing->csiArrange[0]=p[0];
          treeSing->csiArrange[1]=p[1];
          treeSing->clock=indexClock+1;
          treeSing->ud=indexUD;         treeSing->fb=indexFB;
          treeSing->tref[1]=T_ref[1];   treeSing->refpk[1]=maxfn[1];
          treeSing->refmn[1]=minfn[1];
          //std::cout<< " \n Index clock: "<<indexClock<<std::endl;
          //std::cout<< " Gap config FB is  : " <<p[1]<<std::endl;
          //std::cout<< " Gap config UD is  : " <<p[0]<<std::endl;
          //std::cout<< " size of nChannel is : " <<treeRaw->indexCsI[i]-1<<std::endl;
          //std::cout<< " size of nSample is  : " <<treeRaw->nSample[i]<<std::endl;
          //std::cout<< " Chan No.: "<<treeRaw->nChannel<<std::endl;
          //std::cout<<"\n ---------------------------------------------------------\n";
          //std::cout<<" \n\n  ------>ref time and peak time (2): "<<T_ref[1]<<" "<<refgaus[1]<<"\n";
          //std::cout<<" \n\n  ------> CDF timing:  "<<(valx2-valx1)<<" \n\n";
          delete f1; delete f2;
          break;
        case 8:
          csiT=true;
          f2->SetParameters(21.3,29.6,1.85,120);
          h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f2,"QR+");
          h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f1,"QR+");
          maxfn[2]=f1->GetMaximum();
          minfn[2]=f2->GetMinimum();
          cf50[2]=.5*maxfn[2]+.5*minfn[2];
          T_ref[2]=f2->GetX(cf50[2]);
          refgaus[2]=f1->GetX(maxfn[2]);
          // Fill tree variables here:
          treeSing->indexCsI=treeRaw->indexCsI[i];
          treeSing->phei=dummy;
          treeSing->rgaus[2]=refgaus[2];
          treeSing->csiArrange[0]=p[0];
          treeSing->csiArrange[1]=p[1];
          treeSing->clock=indexClock+1;
          treeSing->ud=indexUD;         treeSing->fb=indexFB;
          treeSing->tref[2]=T_ref[2];   treeSing->refpk[2]=maxfn[2];
          treeSing->refmn[2]=minfn[2];
          //std::cout<< " \n Index clock: "<<indexClock<<std::endl;
          //std::cout<< " Gap config FB is  : " <<p[1]<<std::endl;
          //std::cout<< " Gap config UD is  : " <<p[0]<<std::endl;
          //std::cout<< " size of nChannel is : " <<treeRaw->indexCsI[i]-1<<std::endl;
          //std::cout<< " size of nSample is  : " <<treeRaw->nSample[i]<<std::endl;
          //std::cout<<"\n ---------------------------------------------------------\n";
          //std::cout<<" \n\n  ------>ref time and peak time (3): "<<T_ref[2]<<" "<<refgaus[2]<<"\n";
          ////std::cout<<" \n\n  ------> CDF timing:  "<<(valx2-valx1)<<" \n\n";
          delete f1; delete f2;
          break;
      }// end of switch statement
      delete h1Fits[indexClock][indexFB][indexUD][indexModule];
      //delete f1; delete f2;
    } // end of ref-time fired "if" loop
    
    if(firedCsI)
      goto jailbreak; */
    // Painlessly remove the both event-tag and timing modules from 
    // larger analysis --> hope for loss in performance
    if(!(treeRaw->indexCsI[i]==16 && indexFB==0 && indexUD==0) ||
		    !(indexClock==0 || indexClock==2 || indexClock==4 ||
			    indexClock==6 || indexClock==8 || indexClock==10)){
      for(UInt_t iData=0;iData<treeRaw->nSample[i];iData++){
        h1Fits[indexClock][indexFB][indexUD][indexModule]->SetBinContent(iData+1,treeRaw->data[i][iData]);
      }
      #pragma omp parallel num_threads(8)
      // Get x-bin corresponding to the max
      x1=h1Fits[indexClock][indexFB][indexUD][indexModule]->
      	GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
      if(x1>=50 && x1<=70){ // <-- prelim. timing cut if loop
        //std::cout<<"\n ------- Within Signal Loop Event number is:  "<<treeRaw->eventNo<<" -------\n\n";
        if(treeRaw->nChannel>=7){ // Start by checking how many CsI crystals have fired
          if(treeRaw->indexCsI[i]-1==15){
            std::cout<<"\n ------- Within Signal Loop Event number is:  "<<treeRaw->eventNo<<" -------\n\n";
            std::cout<< " \n Index clock: "<<indexClock<<std::endl;
            std::cout<< " Gap config FB is  : " <<p[1]<<std::endl;
            std::cout<< " Gap config UD is  : " <<p[0]<<std::endl;
            std::cout<< " size of nChannel is : " <<treeRaw->indexCsI[i]-1<<std::endl;
            std::cout<<"\n ---------------------------------------------------------\n";
            //std::cout<< " size of nSample is  : " <<treeRaw->nSample[i]<<std::endl;
            //std::cout<< " Chan No.: "<<treeRaw->nChannel<<std::endl;
          }
          xpos.clear();
          val.clear();
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
          //std::cout<<" The baseline is: "<<bl<<std::endl;
          if(y1 == 1023){
            double x;
            for(int ivar=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin(); ivar< 250; ivar+= 1){
              while(y1 == 1023){
                x = ivar++;
                y1 = h1Fits[indexClock][indexFB][indexUD][indexModule]->
          	      GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(x));
              }
            }
            std::cout<< "\n Okay this thing is smart x =" << x << endl;
            TF1* f1=new TF1("wave1",overr.c_str(), 0.5, x1);
            // fill parameters for the fit function(s)
            for(int ivar=0; ivar<10; ivar+=1){
              f1->SetParameter(ivar, mn2.par(ivar));
              f1->SetParLimits(ivar,mn2.parmin(ivar),mn2.parlim(ivar));
            }
	    nfound=s->Search(h1Fits[indexClock][indexFB][indexUD][indexModule], 2,"",0.10);
	    xpeaks=new double();
	    xpeaks=s->GetPositionX();
            double posX[2];
            for(int ivar=0; ivar<nfound; ivar++){
              double a=xpeaks[ivar];
              int bin=1+Int_t(a+.5);
              posX[ivar]=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(bin);
            }
            sort(xpeaks,xpeaks+nfound);
	    //std::cout<<"\n ---->  Gotta check the hell outta this shit \n";
            rtime=(xx1+xx2)/2;
            f1->SetParameter(0,y1+1023*2);
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
            par.clear(); err.clear();
            //std::cout<< "  ----> Testing left and right x-limits: "<<xx1<< ", "<<xx2<<endl;
            for(int n=0; n<10;n+=1){
              upar2.Add(mn2.nameL(n).c_str(), f1->GetParameter(n), 0.1);
            }
            // create Migrad minimizer
            MnMigrad migrad(ffcn, upar2);
            //FunctionMinimum min = migrad();
            FunctionMinimum min = migrad(180,1e-5);
            try{
              throw logic_error("Assertion `s0.IsValid()' failed.");
            }
            catch(const std::logic_error & e){
              std::cerr<<" Found assertion error \n";
            }
            std::cout<<"minimum: "<<min<<std::endl;
            MnHesse hesse;
            param.clear();
            for(int ivar=0; ivar<10; ivar+=1){
              param.push_back(migrad.Value(mn2.nameL(ivar).c_str()));
              //std::cout<< "  par["<<ivar<<"] value --> ["<<param[ivar]<<"] \n";
            }
            for(int ivar=1; ivar<+h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetNbinsX()+1; ivar+=1){
              double x=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(ivar);
              double yv=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinContent(ivar);
              double mnfit=mn2.overR(x, param);
              double res=100*(yv-mnfit)/yv;
              h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(ivar, mnfit);
              if(ivar>=xx1 && ivar<=xx2){
                h1Diff[indexClock][indexFB][indexUD][indexModule]->SetBinContent(ivar, 0);
              }else{
                h1Diff[indexClock][indexFB][indexUD][indexModule]->SetBinContent(ivar, res);
              }
            }
            max=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                  GetBinLowEdge(h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
            if(max>=60 && max<=65){
	      loopX=true; firedCsI=true;
              clock=indexClock;
              fb=indexFB;
              ud=indexUD;
              module=indexModule;
              mnx=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                    GetBinLowEdge(h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetMinimumBin());
              may=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                    GetBinContent(h1Mnft[indexClock][indexFB][indexUD][indexModule]->FindBin(max));
              mny=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                    GetBinContent(h1Mnft[indexClock][indexFB][indexUD][indexModule]->FindBin(mnx));
              //Integrating fitting function for more accurate gain calibration
              TAxis* axis=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetXaxis();
              bmin=axis->FindBin(xmin);
              bmax=axis->FindBin(xmax);
              integral=h1Mnft[indexClock][indexFB][indexUD][indexModule]->Integral(bmin,bmax);
              integral-=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
        		GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
              integral-=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
        		GetBinContent(bmax)*(axis->GetBinLowEdge(bmin)-bmax)/axis->GetBinWidth(bmax);
              area=integral-(mny*250);
              h1Intg->Fill(area);
              std::cout<<" =====> Value for integral1 is:" <<area<<"  "<<xpeaks[1]<<std::endl;
	      int tAB=0;
              if((p[0]=='u' || p[0]=='U') && (indexModule>=9)) tAB=1;
              std::cout<< " reading out value for Up and TypeB Crystal: "<<tAB<<endl;
              if((p[0]=='d' || p[0]=='D') && (indexModule<=8)) tAB=1;
              double csitheta=theta[indexFB][indexModule];
              double csiphi=phi[indexClock][ud][tAB];
              csThet.push_back(csitheta), csPhi.push_back(csiphi);
              treeSing->thSing=csitheta;
              treeSing->phiSing=csiphi;
              diff=may-mny; int imod=0, igap=4*indexClock+2;
	      //double tsigH=h1Mnft[indexClock][indexFB][indexUD][indexModule]->FindFirstBinAbove(.5*diff+mny);
              int csimod=(-1*treeRaw->indexCsI[i])+1;
              if(p[0]=='u' || p[0]=='U') igap=4*(indexClock+1);
              if(p[1]=='b' || p[1]=='B') csimod=treeRaw->indexCsI[i];
              h1Pamp->Fill(diff); h1ped->Fill(mny);      tpeak=max;
              h2clus->Fill(csimod,igap,diff);
              h1kmu2->Fill(diff);
              treeSing->indexCsI=treeRaw->indexCsI[i];
              treeSing->tpeak=max;
              treeSing->trise=param[1];
	      tsigL=h1Mnft[indexClock][indexFB][indexUD][indexModule]->FindFirstBinAbove(.5*diff+mny);
              treeSing->calInt=area;
              treeSing->csiArrange[0]=p[0];
              treeSing->csiArrange[1]=p[1];
              treeSing->clock=indexClock+1;
              treeSing->ovrped=mny;
              treeSing->ovrpH=diff;
              treeSing->waveID=4;
              treeSing->ud=indexUD;                  treeSing->fb=indexFB;
              treeSing->phei=diff;                   treeSing->ped=mny;
              treeSing->tcorr[0]=(tsigL-T_ref[0]);   
              treeSing->tcorr[1]=(tsigL-T_ref[1]);   
              treeSing->tcorr[2]=(tsigL-T_ref[2]);   
	      if(nfound==2){
                treeSing->ovrpLoc=xpeaks[1];
	      }
	      delete f1,xpeaks;
	      //if(loopX && csiT)
	        goto exitLoop;
            } // <---  End of K+ decay time if loop
          }// <--- Use this to get rid of double and single fitting functions * /
	  if(y1<1023){ // forgo if else for more efficient less nested code
            //std::cout<< "  Size of x is:  "<<xpos.size()<<endl;
            xx1=xpos.size(); xx2=0; ymax=y1;
            nfound=s->Search(h1Fits[indexClock][indexFB][indexUD][indexModule], 2,"",0.10);
            if(nfound==2){
              TF1* f1=new TF1("f1",doubleFit.c_str(),1.0,250);
              for(int n=0; n<15; n+=1){
                f1->SetParameter(n,mn2.par(n));
                f1->SetParLimits(n,mn2.parmin(n),mn2.parlim(n));
              }
	      xpeaks=new double();
              xpeaks=s->GetPositionX();
              double posX[2];
              for(int ivar=0; ivar<nfound; ivar++){
                double a=xpeaks[ivar];
                int bin=1+Int_t(a+.5);
                posX[ivar]=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(bin);
              }
              sort(xpeaks,xpeaks+nfound);
              //std::cout<<" ****** Checking the position of the peaks: "<<xpeaks[0]<<", "<<xpeaks[1]<<endl;
              double yp2=h1Fits[indexClock][indexFB][indexUD][indexModule]->
          	    GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(xpeaks[1]));
              f1->SetParameter(0,y1);
              f1->SetParLimits(0,y1-61.7,y1+971.7);
              f1->SetParameter(1,xpeaks[0]+.1);
              f1->SetParameter(8,bl);
              f1->SetParLimits(8,bl-61.7,bl+171.7);
              f1->SetParameter(9,xpeaks[1]+.1);
              f1->SetParLimits(1,xpeaks[0]-61.7,xpeaks[0]+71.7);
              f1->SetParLimits(9,xpeaks[1]-61.7,xpeaks[1]+71.7);
              f1->SetParameter(10,yp2);
              f1->SetParLimits(10,yp2-61.7,yp2+971.7);
      	      p0->Fill(f1->GetParameter(0));
      	      p1->Fill(f1->GetParameter(1));
      	      p9->Fill(f1->GetParameter(9));
      	      p10->Fill(f1->GetParameter(10));
              h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f1,"0");
      	      f1chi2=f1->GetChisquare();
      	      //std::cout<<" **************\n ******** Chi2 for F1 fit "<<f1->GetChisquare()<<endl;
              // Create wrapper for minimizer
              fitfn2 ffcn1(xpos, xx1, xx2, val, ymax);
              param.clear(); parm.clear(); err.clear();
              MnUserParameters upar;
              for(int n=0; n<15;n+=1){
                upar.Add(mn2.nameL(n).c_str(), f1->GetParameter(n),1e-3); //,parmin(n), parlim(n), 0.1);
                //upar.Add(nameL(n).c_str(), f1->GetParameter(n), parmin(n), parlim(n), 1e-3);
              }
              for(int n=0; n<15;n+=1){
                upar.Add(mn2.nameL(n).c_str(), mn2.par(n), 1e-3);
              }
              // create Migrad minimizer
              MnMigrad migrad(ffcn1, upar);
              //FunctionMinimum min = migrad();  //6000,1e-9);
              FunctionMinimum min = migrad(180,1e-6);
              std::cout<<"minimum: "<<min<<std::endl;
              //MnHesse hesse;
              for(int ivar=0; ivar<15; ivar+=1){
                param.push_back(migrad.Value(mn2.nameL(ivar).c_str()));
                //std::cout<< "  par["<<i<<"] value --> ["<<param[ivar]<<"] \n";
              }
              for(int ivar=0; ivar<h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetNbinsX()+1; ivar+=1){
                double x=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(ivar);
                double yv=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinContent(ivar);
                double mnfit=mn2.model2(x, param);
                h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(ivar, mnfit);
                double res=100*(yv-mnfit)/yv;
                h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(ivar, mnfit);
                h1Diff[indexClock][indexFB][indexUD][indexModule]->SetBinContent(ivar, res);
              }
              max=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
          	    GetBinLowEdge(h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
              if(max>=60 && max<=65){
		loopX=true; firedCsI=true;
                mnx=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinLowEdge(h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetMinimumBin());
                may=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(h1Mnft[indexClock][indexFB][indexUD][indexModule]->FindBin(max));
                mny=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(h1Mnft[indexClock][indexFB][indexUD][indexModule]->FindBin(mnx));
                TAxis* axis=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetXaxis();
                bmin=axis->FindBin(xmin);
                bmax=axis->FindBin(xmax);
                integral=h1Mnft[indexClock][indexFB][indexUD][indexModule]->Integral(bmin,bmax);
                integral-=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
                integral-=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(bmax)*(axis->GetBinLowEdge(bmin)-bmax)/axis->GetBinWidth(bmax);
                diff=may-mny;
		area=integral-mny*250;
        	clock=indexClock;
        	fb=indexFB;
        	ud=indexUD;
        	module=indexModule;
                int imod=0, igap=4*indexClock+2;
                int csimod=(-1*treeRaw->indexCsI[i])+1;
                if(p[0]=='u' || p[0]=='U') igap=4*(indexClock+1);
                if(p[1]=='b' || p[1]=='B') csimod=treeRaw->indexCsI[i];
                h2clus->Fill(csimod,igap,diff);
                h1kmu2->Fill(diff); h1ped->Fill(mny);
		int tAB=0;
                if((p[0]=='u' || p[0]=='U') && (indexModule>=9)) tAB=1;
                std::cout<< " reading out value for Up and TypeB Crystal: "<<tAB<<endl;
                if((p[0]=='d' || p[0]=='D') && (indexModule<=8)) tAB=1;
                double csitheta=theta[indexFB][indexModule];
                double csiphi=phi[indexClock][ud][tAB];
                csThet.push_back(csitheta), csPhi.push_back(csiphi);
		//Filling the TTree here:
                treeSing->indexCsI=treeRaw->indexCsI[i];
                treeSing->tpeak=max;
                treeSing->trise=param[1];
	        tsigL=h1Mnft[indexClock][indexFB][indexUD][indexModule]->FindFirstBinAbove(.5*diff+mny);
                treeSing->kmu2=diff;          
                treeSing->phei=diff;          treeSing->ped=mny;
                treeSing->dubPed=mny;         
                treeSing->calInt=area;
                treeSing->intKmu2=area;
                treeSing->csiArrange[0]=p[0];
                treeSing->csiArrange[1]=p[1];
                treeSing->clock=indexClock+1;
                treeSing->thSing=csitheta;
                treeSing->phiSing=csiphi;
                treeSing->dubphei=xpeaks[1];
                treeSing->waveID=2;
                treeSing->ud=indexUD;         treeSing->fb=indexFB;
                treeSing->tcorr[0]=(tsigL-T_ref[0]);
                treeSing->tcorr[1]=(tsigL-T_ref[1]);
                treeSing->tcorr[2]=(tsigL-T_ref[2]);
                treeSing->tcorr[2]=(tsigL-T_ref[2]);
		std::cout<<" --->rise time: "<<param[1]<<"\n";
		delete f1,xpeaks;
	        //if(loopX && csiT)
	          goto exitLoop;
              } // <-- End of K+ decay time if loop
            } //<-- Use to get rid of 2 peaks functions here * /
	    if(nfound==1){
              TF1* f1=new TF1("f1",singleFit.c_str(),1,250);
              for(int n=0; n<9; n+=1){
                f1->SetParameter(n,mn2.par(n));
                f1->SetParLimits(n,mn2.parmin(n),mn2.parlim(n));
              }
              f1->SetParameter(0,y1);
              f1->SetParLimits(0,y1-61.7,y1+971.7);
              f1->SetParameter(1,x1);
              f1->SetParLimits(1,x1-261.7,x1+571.7);
              f1->SetParameter(8,bl);
              f1->SetParLimits(8,bl-161.7,bl+171.7);
              h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f1);//,"0");
      	      f1chi2=f1->GetChisquare();
              // Create wrapper for minimizer
              fitfn ffcn1(xpos, xx1, xx2, val, ymax);
              param.clear(); parm.clear(); err.clear();
              MnUserParameters upar;
              for(int n=0; n<9;n+=1){
                upar.Add(mn2.nameL(n).c_str(), f1->GetParameter(n),1e-3); //,parmin(n), parlim(n), 0.1);
              }
              std::cout<<"  Okay we have cleared the loop!"<<endl;
              // create Migrad minimizer
              MnMigrad migrad(ffcn1, upar);
              FunctionMinimum min = migrad(40,1e-4);
              std::cout<<"minimum: "<<min<<std::endl;
              //FunctionMinimum min1 = migrad(3000,1e-9);
              for(int ivar=0; ivar<9; ivar+=1){
                param.push_back(migrad.Value(mn2.nameL(ivar).c_str()));
                //std::cout<< "  par["<<i<<"] value --> ["<<param[ivar]<<"] \n";
              }
      	      bool dpval=false;
              for(int ivar=0; ivar<h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetNbinsX()+1; ivar+=1){
                double x=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(ivar);
                double yv=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinContent(ivar);
                double mnfit=mn2.model(x, param);
                double res=100*(yv-mnfit)/yv;
                h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(ivar, mnfit);
                h1Diff[indexClock][indexFB][indexUD][indexModule]->SetBinContent(ivar, res);
              }
              max=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
          	    GetBinLowEdge(h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
              if(max>=60 && max<=65){
		loopX=true; firedCsI=true;
                clock=indexClock;
                fb=indexFB;
                ud=indexUD;
                module=indexModule;
                mnx=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                    GetBinLowEdge(h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetMinimumBin());
                may=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                    GetBinContent(h1Mnft[indexClock][indexFB][indexUD][indexModule]->FindBin(max));
                mny=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                    GetBinContent(h1Mnft[indexClock][indexFB][indexUD][indexModule]->FindBin(mnx));
                //Integrating fitting function for more accurate gain calibration
                TAxis* axis=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetXaxis();
                bmin=axis->FindBin(xmin);
                bmax=axis->FindBin(xmax);
                integral=h1Mnft[indexClock][indexFB][indexUD][indexModule]->Integral(bmin,bmax);
                integral-=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
                integral-=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(bmax)*(axis->GetBinLowEdge(bmin)-bmax)/axis->GetBinWidth(bmax);
                area=integral-mny*250;
                std::cout<<" =====> Value for integral1 is:" <<area<<endl;
                h1Intg->Fill(area);  
		int tAB=0;
                if((p[0]=='u' || p[0]=='U') && (indexModule>=9)) tAB=1;
                std::cout<< " reading out value for Up and TypeB Crystal: "<<tAB<<endl;
                if((p[0]=='d' || p[0]=='D') && (indexModule<=8)) tAB=1;
                double csitheta=theta[indexFB][indexModule];
                double csiphi=phi[indexClock][ud][tAB];
                csThet.push_back(csitheta), csPhi.push_back(csiphi);
                diff=may-mny;
                int imod=0, igap=4*indexClock+2;
                int csimod=(-1*treeRaw->indexCsI[i])+1;
                if(p[0]=='u' || p[0]=='U') igap=4*(indexClock+1);
                if(p[1]=='b' || p[1]=='B') csimod=treeRaw->indexCsI[i];
                h2clus->Fill(csimod,igap,diff);
                h1Pamp->Fill(diff); h1ped->Fill(mny);
		// fill TTree here:
                treeSing->thSing=csitheta;
                treeSing->phiSing=csiphi;
                treeSing->indexCsI=treeRaw->indexCsI[i];
                treeSing->tpeak=max;
                treeSing->trise=param[1];
	        tsigL=h1Mnft[indexClock][indexFB][indexUD][indexModule]->FindFirstBinAbove(.5*diff+mny);
                treeSing->calInt=area;
                treeSing->csiArrange[0]=p[0];
                treeSing->csiArrange[1]=p[1];
                treeSing->clock=indexClock+1;
                treeSing->phei=diff;          treeSing->ped=mny;
		treeSing->sphei=diff;         treeSing->sptime=max;
                treeSing->ud=indexUD;         treeSing->fb=indexFB;
                treeSing->waveID=1;
                treeSing->tcorr[0]=(tsigL-T_ref[0]);
                treeSing->tcorr[1]=(tsigL-T_ref[1]);
                treeSing->tcorr[2]=(tsigL-T_ref[2]);
		std::cout<<" --->rise time: "<<param[1]<<"\n";
		delete f1;
	        //if(loopX && csiT)
	          goto exitLoop;
              } // <--- Use this to get rid of double and single fitting functions * /
            } // <-- End of nfound==1, single peak
          }
        } // <--- End of No. CsI crystals that fired
      } // <--- End of timing cut if loop
    jailbreak:
      firedCsI=true;
    } // <--- End of if loop
    delete h1Fits[indexClock][indexFB][indexUD][indexModule];
  } // <--- End of nChannel for loop
  exitLoop:
    loopX=false;
    firedCsI=false; csiT=false;

  return 0;
}

Long_t Det_CsI::finalize(){
  return 0;
}

extern "C"{
  Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p){
    return (Plugin *) new Det_CsI(in,out,inf_,outf_,p);
  }
}
ClassImp(Det_CsI);
