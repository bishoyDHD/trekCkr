#include <Det_CsI.h>
//#include <function.h>
#include <mn2CsIfn.h>
#include<iostream>
#include<cmath>

double gaanKak=800815;
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
};

Det_CsI::~Det_CsI(){
  delete s;
};


Long_t Det_CsI::cmdline(char *cmd){
  //add cmdline hanling here

  return 0; // 0 = all ok
};
string singleFit = singlemodel();   string quadFit   = quadruplemodel();
string doubleFit = doublemodel();   string overr     = overrangemodel();
string tripleFit = triplemodel();
Long_t Det_CsI::histos(){
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<16;iModule++){
          std::ostringstream name, name2, name3, name4, name5, name6, tname;
          name<<"stat_"; name2<<"Mnfit_"; name3<<"F1fit_"; name4<<"Dhfit_"; name5<<"pHeight", name6<<"Diff";
          tname<<"time";
          name<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name2<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name3<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name4<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name5<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name6<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          tname<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          h1time[iClock][iFB][iUD][iModule]=dH1(tname.str().c_str(),"stat",250,0,250);
          h1Fits[iClock][iFB][iUD][iModule]=dH1(name.str().c_str(),"stat",250,0,250);
          h1Amps[iClock][iFB][iUD][iModule]=dH1(name5.str().c_str(),"stat",250,0,250);
          h1Mnft[iClock][iFB][iUD][iModule]=dH1(name2.str().c_str(),"stat",250,0,250);
          h1Diff[iClock][iFB][iUD][iModule]=dH1(name6.str().c_str(),"stat",250,0,250);
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
  h1cali=dH1("Calibr","Integrated pulse height distribution", 250, 0, 1000);
  h1ped=dH1("Ped","Pedestals for the waveform ", 250, 0, 1000);
  return 0;
}

Long_t Det_CsI::startup(){
  getBranchObject("vf48",(TObject **) &treeRaw);
  getBranchObject("RawBeamInfo",(TObject **) &treeBeam);
  gStyle->SetOptStat(0);

  return 0;
}

Long_t Det_CsI::process(){
  //std::cout<<" ---> Baisically the number of cluster: "<<treeRaw->nChannel<<std::endl;
  if(treeRaw->nChannel==7){ // Start by checking how many CsI crystals have fired
    for(UInt_t i=0;i<treeRaw->nChannel;i++){
      char* p=(char*)&(treeRaw->nameModule[i]);
      int moduleName=(p[3]-'0')*10+(p[2]-'0')-1;
      //std::cout<< " Index clock: "<<indexClock<<endl;
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

      if(treeRaw->indexCsI[i]==16){
        //std::cout<< " CHAN 16:  Index clock: "<<indexClock<<std::endl;
        //std::cout<< " CHAN 16:  Gap config FB is  : " <<p[1]<<std::endl;
        //std::cout<< " CHAN 16:  Gap config UD is  : " <<p[0]<<std::endl;
        //std::cout<< " CHAN 16:  size of nChannel is : " <<treeRaw->indexCsI[i]-1<<std::endl;
        //std::cout<< " CHAN 16:  size of nSample is  : " <<treeRaw->nSample[i]<<std::endl;
        //std::cout<< " CHAN 16:  Chan No.: "<<treeRaw->nChannel<<std::endl;
	/*
        if((indexClock==0 || indexClock==4) && indexFB==0 && indexUD==0){
          for(UInt_t iData=0;iData<treeRaw->nSample[i];iData++){
            h1Fits[indexClock][indexFB][indexUD][indexModule]->SetBinContent(iData+1,treeRaw->data[i][iData]);
          }
          x1=h1Fits[indexClock][indexFB][indexUD][indexModule]->
              GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
          tref=x1;
          clock=indexClock;
          //std::cout<< "   ********************************************************* \n";
          //std::cout<< "            this is gap 16, ref time= "<<x1<<endl;
          //std::cout<< "   ********************************************************* \n";
        }*/
      }
      if((treeRaw->indexCsI[i]!=16) /*&& (indexClock==0 && indexFB==1 && indexUD==0)*/){
        //std::cout<< " Index clock: "<<indexClock<<std::endl;
        //std::cout<< " Gap config FB is  : " <<p[1]<<std::endl;
        //std::cout<< " Gap config UD is  : " <<p[0]<<std::endl;
        //std::cout<< " size of nChannel is : " <<treeRaw->indexCsI[i]-1<<std::endl;
        //std::cout<< " size of nSample is  : " <<treeRaw->nSample[i]<<std::endl;
        //std::cout<< " Chan No.: "<<treeRaw->nChannel<<std::endl;
        for(UInt_t iData=0;iData<treeRaw->nSample[i];iData++){
          h1Fits[indexClock][indexFB][indexUD][indexModule]->SetBinContent(iData+1,treeRaw->data[i][iData]);
        }
        #pragma omp parallel num_threads(8)
        xpos.clear();
        val.clear();
        // Utilize fitting method
        x1=h1Fits[indexClock][indexFB][indexUD][indexModule]->
        	GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
        double mn = h1Fits[indexClock][indexFB][indexUD][indexModule]->
        	GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMinimumBin());
        y1 = h1Fits[indexClock][indexFB][indexUD][indexModule]->
        	GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(x1));
        double bl = h1Fits[indexClock][indexFB][indexUD][indexModule]->
        	GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(mn));
	std::cout<<" The baseline is: "<<bl<<std::endl;
        for(int i=1; i<=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetNbinsX(); i+=1){
          xpos.push_back(i);
          val.push_back(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinContent(i));
        }
        //std::cout<< " iHist: "<<iHist<<endl;
        if(x1>=50 && x1<=70){
          if(y1 == 1023){
            double x;
            for(int i = h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin(); i < 250; i += 1){
              while(y1 == 1023){
                x = i++;
                y1 = h1Fits[indexClock][indexFB][indexUD][indexModule]->
          	      GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(x));
              }
            }
            std::cout<< " Okay this thing is smart x =" << x << endl;
            TF1* f1=new TF1("wave1",overr.c_str(), 0.5, x1);
            // fill parameters for the fit function(s)
            for(int i = 1; i < 10; i+=1){
              f1->SetParameter(i, par(i));
              f1->SetParLimits(i,parmin(i),parlim(i));
            }
            double rtime=(xx1+xx2)/2;
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
            std::vector<double> par(10), err(10);
            par.clear(); err.clear();
            std::cout<< "  ----> Testing left and right x-limits: "<<xx1<< ", "<<xx2<<endl;
            for(int n=0; n<10;n+=1){
              upar2.Add(nameL(n).c_str(), f1->GetParameter(n), 0.1);
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
            // Minos factory
            //MnMinos minos(ffcn, min);
            //FunctionMinimum min3 = migrad();
            //MnMigrad migrad(ffcn, upar2);
            //FunctionMinimum min2 = migrad(8000,1e-5);
            MnHesse hesse;
            std::vector<double> param;
            param.clear();
            for(int i=0; i<10; i+=1){
              param.push_back(migrad.Value(nameL(i).c_str()));
              //std::cout<< "  par["<<i<<"] value --> ["<<param[i]<<"] \n";
            }
            for(int i=1; i<+h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetNbinsX()+1; i+=1){
              double x=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(i);
              double yv=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinContent(i);
              double mnfit=overR(x, param);
              double res=100*(yv-mnfit)/yv;
              h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(i, mnfit);
              if(i>=xx1 && i<=xx2){
                h1Diff[indexClock][indexFB][indexUD][indexModule]->SetBinContent(i, 0);
              }else{
                h1Diff[indexClock][indexFB][indexUD][indexModule]->SetBinContent(i, res);
              }
            }
            double mnx,mny,max,may;
            max=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                  GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
            if(max>=60 && max<=65){
              clock=indexClock;
              fb=indexFB;
              ud=indexUD;
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
              h1Intg->Fill(area);
              calInt=area;
              std::cout<<" =====> Value for integral1 is:" <<area<<endl;
              double diff=may-mny; int imod=0, igap=4*indexClock+2;
              int csimod=(-1*treeRaw->indexCsI[i])+1;
              if(p[0]=='u' || p[0]=='U') igap=4*(indexClock+1);
              if(p[1]=='b' || p[1]=='B') csimod=treeRaw->indexCsI[i];
              h1Pamp->Fill(diff); h1ped->Fill(mny);      tpeak=max;
              h2clus->Fill(csimod,igap,diff);
              phei=diff;          ped=mny;
              kmu2=diff;          h1kmu2->Fill(diff);
            }else{
              phei=-100, fb=-100, ud=-100, module=-100; tpeak=-100;
            } // <--- Use this to get rid of double and single fitting functions * /
            //fileOut->cd();
            //h1Intg->Write();
            //h1Fits[indexClock][indexFB][indexUD][indexModule]->Write();
            //h1Mnft[indexClock][indexFB][indexUD][indexModule]->Write();
            //h1Diff[indexClock][indexFB][indexUD][indexModule]->Write();
          }
	  if(y1<1023){ // forgo if else for more efficient less nested code
            //std::cout<< "  Size of x is:  "<<xpos.size()<<endl;
            xx1=xpos.size(); xx2=0; ymax=y1;
            nfound=s->Search(h1Fits[indexClock][indexFB][indexUD][indexModule], 2,"",0.10);
            if(nfound==2){
              TF1* f1=new TF1("f1",doublemodel().c_str(),1.0,250);
              for(int n=0; n<15; n+=1){
                f1->SetParameter(n,par(n));
                f1->SetParLimits(n,parmin(n),parlim(n));
              }
              double *xpeaks=s->GetPositionX();
              double posX[2];
              for(int i=0; i<nfound; i++){
                double a=xpeaks[i];
                int bin=1+Int_t(a+.5);
                posX[i]=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(bin);
              }
              sort(xpeaks,xpeaks+nfound);
              std::cout<<" ****** Checking the position of the peaks: "<<xpeaks[0]<<", "<<xpeaks[1]<<endl;
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
      	      //std::cout<<" *******************************\n ******** Chi2 for F1 fit "<<f1->GetChisquare()<<endl;
              Minuit2Minimizer* mnu2=new Minuit2Minimizer("Minuit2");
              // Create wrapper for minimizer
              fitfn2 ffcn1(xpos, xx1, xx2, val, ymax);
              std::vector<double> param;
              std::vector<double> parm(15), err(15);
              param.clear(); parm.clear(); err.clear();
              MnUserParameters upar;
              for(int n=0; n<15;n+=1){
                upar.Add(nameL(n).c_str(), f1->GetParameter(n),1e-3); //,parmin(n), parlim(n), 0.1);
                //upar.Add(nameL(n).c_str(), f1->GetParameter(n), parmin(n), parlim(n), 1e-3);
              }
              for(int n=0; n<15;n+=1){
                upar.Add(nameL(n).c_str(), par(n), 1e-3);
              }
              // create Migrad minimizer
              MnMigrad migrad(ffcn1, upar);
              //FunctionMinimum min = migrad();  //6000,1e-9);
              FunctionMinimum min = migrad(180,1e-6);
              std::cout<<"minimum: "<<min<<std::endl;
              /*
              try{
                throw logic_error("Assertion `s0.IsValid()' failed.");
              }
              catch(const std::logic_error & e){
                std::cerr<<" Found assertion error \n";
              }*/
              //MnHesse hesse;
              for(int i=0; i<15; i+=1){
                param.push_back(migrad.Value(nameL(i).c_str()));
                //std::cout<< "  par["<<i<<"] value --> ["<<param[i]<<"] \n";
              }
              //hesse(ffcn1, min1);
              //std::cout<<"minimum after hesse: "<<min<<std::endl;
              //FunctionMinimum min1 = migrad();
              //std::cout<<"minimum1: "<<min1<<std::endl;
              for(int i=0; i<h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetNbinsX()+1; i+=1){
                double x=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(i);
                double yv=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinContent(i);
                double mnfit=model2(x, param);
                h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(i, mnfit);
                double res=100*(yv-mnfit)/yv;
                h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(i, mnfit);
                h1Diff[indexClock][indexFB][indexUD][indexModule]->SetBinContent(i, res);
              }
              double mnx,mny,max,may;
              max=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
          	    GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
              if(max>=60 && max<=65){
                mnx=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMinimumBin());
                may=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(max));
                mny=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(mnx));
                double diff=may-mny;
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
        	kmu2=diff;          ped=mny;
        	dubPed=mny;         tpeak=max;
        	if(indexClock==0 && indexFB==1 && indexUD==0){
        	  h1cali->Fill(diff);
        	  TSpectrum *sp = new TSpectrum(4);
                  Int_t nf = sp->Search(h1cali,1,"",0.10);
                  double lowRange,upRange;
        	  Double_t *xps = sp->GetPositionX();
                  std::sort(xps,xps+nf);
                  lowRange=h1cali->GetBinLowEdge(h1cali->GetMaximumBin()-xps[0]/2.0);
                  upRange=h1cali->GetBinLowEdge(h1cali->GetMaximumBin()+xps[0]*3.0/2.0);
                  h1cali->Fit("gaus","Q","",lowRange,upRange);
                }
              }else{
                kmu2=-100; calInt=-100;
        	module=-100;   tpeak=-100;
        	dubPed=-100;
              } //<-- Use to get rid of 2 peaks functions here * /
            }else if(nfound==1){
              TF1* f1=new TF1("f1",singlemodel().c_str(),1.0,250);
              for(int n=0; n<9; n+=1){
                f1->SetParameter(n,par(n));
                f1->SetParLimits(n,parmin(n),parlim(n));
              }
              f1->SetParameter(0,y1);
              f1->SetParLimits(0,y1-61.7,y1+971.7);
              f1->SetParameter(1,x1);
              f1->SetParLimits(1,x1-261.7,x1+571.7);
              f1->SetParameter(8,bl);
              f1->SetParLimits(8,bl-161.7,bl+171.7);
              h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f1,"0");
      	      f1chi2=f1->GetChisquare();
              Minuit2Minimizer* mnu2=new Minuit2Minimizer("Minuit2");
              // Create wrapper for minimizer
              //TMinuit* gmin = new TMinuit(8);
              fitfn ffcn1(xpos, xx1, xx2, val, ymax);
              std::vector<double> param;
              std::vector<double> parm(15), err(15);
              param.clear(); parm.clear(); err.clear();
              MnUserParameters upar;
              for(int n=0; n<9;n+=1){
                upar.Add(nameL(n).c_str(), f1->GetParameter(n),1e-3); //,parmin(n), parlim(n), 0.1);
                //upar.Add(nameL(n).c_str(), f1->GetParameter(n),parmin(n), parlim(n), 1e-3);
                //upar.SetLimits(n,parmin(n), parlim(n));
                //parm.push_back(f1->GetParameter(n)); err.push_back(0);
                //upar.SetLimits(nameL(n).c_str(),parmin(n), parlim(n));
              }
              std::cout<<"  Okay we have cleared the loop!"<<endl;
              // create Migrad minimizer
              //MnStrategy mnstra = MnStrategy(2);
              MnMigrad migrad(ffcn1, upar);
              //MnMigrad migrad(ffcn1, parm, err);
              //FunctionMinimum min = migrad();  //6000,1e-9);
              //try{
              FunctionMinimum min = migrad(40,1e-4);
                //throw logic_error("Assertion `s0.IsValid()' failed.");
              std::cout<<"minimum: "<<min<<std::endl;
              /*}
              catch(const std::logic_error & e){
                std::cerr<<" Found assertion error \n";
              }*/
              //MnHesse hesse;
              //FunctionMinimum min1 = migrad(3000,1e-9);
              for(int i=0; i<9; i+=1){
                param.push_back(migrad.Value(nameL(i).c_str()));
                //std::cout<< "  par["<<i<<"] value --> ["<<param[i]<<"] \n";
              }
              //hesse(ffcn1, min1);
              //std::cout<<"minimum after hesse: "<<min<<std::endl;
              //FunctionMinimum min1 = migrad();
              //std::cout<<"minimum1: "<<min1<<std::endl;
      	      bool dpval=false;
              for(int i=0; i<h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetNbinsX()+1; i+=1){
                double x=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(i);
                double yv=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinContent(i);
                double mnfit=model(x, param);
                double res=100*(yv-mnfit)/yv;
      	        /*if(x>=70){
      	        if(abs(res)>5){
      	          dpval=true;
      	          goto dpulse;
      	        }
      	      }*/
                h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(i, mnfit);
                h1Diff[indexClock][indexFB][indexUD][indexModule]->SetBinContent(i, res);
              }
              double mnx,mny,max,may;
              max=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
          	    GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
              if(max>=60 && max<=65){
                clock=indexClock;
                fb=indexFB;
                ud=indexUD;
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
                h1Intg->Fill(area);  tpeak=max;
                calInt=area;
                double diff=may-mny;
                int imod=0, igap=4*indexClock+2;
                int csimod=(-1*treeRaw->indexCsI[i])+1;
                if(p[0]=='u' || p[0]=='U') igap=4*(indexClock+1);
                if(p[1]=='b' || p[1]=='B') csimod=treeRaw->indexCsI[i];
                h2clus->Fill(csimod,igap,diff);
                h1Pamp->Fill(diff); h1ped->Fill(mny);
                phei=diff;          ped=mny;
              }else{
                phei=-100, fb=-100, ud=-100, module=-100; tpeak=-100;
              } // <--- Use this to get rid of double and single fitting functions * /
      	      dpulse:  // <-- On off chance that double pulse misdiagnosed as single pulse
      	        if(dpval){
      	          std::cout<<" *** Checking to make sure this is called at the right time. \n";
                  TF1* f2=new TF1("f2",doublemodel().c_str(),1.0,250);
                  for(int n=0; n<13; n+=1){
                    f2->SetParameter(n,par(n));
                    f2->SetParLimits(n,parmin(n),parlim(n));
                  }
                  f2->SetParameter(0,y1);
                  f2->SetParLimits(0,y1-61.7,y1+971.7);
                  f2->SetParameter(1,x1);
                  f2->SetParameter(8,bl);
                  f2->SetParLimits(8,bl-61.7,bl+171.7);
                  f2->SetParameter(9,100.1);
                  f2->SetParLimits(1, -21.7, 171.7);
                  f2->SetParLimits(9,21.7,371.7);
                  f2->SetParameter(10,y1*.7);
                  f2->SetParLimits(10,.7*y1-61.7,.7*y1+971.7);
                  h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f2,"0");
                  Minuit2Minimizer* mnu2=new Minuit2Minimizer("Minuit2");
                  // Create wrapper for minimizer
                  fitfn2 ffcn1(xpos, xx1, xx2, val, ymax);
                  std::vector<double> param;
                  std::vector<double> parm(15), err(15);
                  param.clear(); parm.clear(); err.clear();
                  MnUserParameters upar;
                  for(int n=0; n<13;n+=1){
                    upar.Add(nameL(n).c_str(), f2->GetParameter(n),1e-3); //,parmin(n), parlim(n), 0.1);
                  }
                  for(int n=0; n<13;n+=1){
                    upar.Add(nameL(n).c_str(), par(n), 1e-3);
                  }
                  // create Migrad minimizer
                  MnMigrad migrad(ffcn1, upar);
                  std::cout<<"minimum: "<<min<<std::endl;
                  //MnHesse hesse;
                  for(int i=0; i<13; i+=1){
                    param.push_back(migrad.Value(nameL(i).c_str()));
                    //std::cout<< "  par["<<i<<"] value --> ["<<param[i]<<"] \n";
                  }
                  for(int i=0; i<h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetNbinsX()+1; i+=1){
                    double x=h1Mnft[indexClock][indexFB][indexUD][indexModule]->GetBinCenter(i);
                    double yv=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinContent(i);
                    double mnfit=model2(x, param);
                    h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(i, mnfit);
                    double res=100*(yv-mnfit)/yv;
                    h1Mnft[indexClock][indexFB][indexUD][indexModule]->SetBinContent(i, mnfit);
                    h1Diff[indexClock][indexFB][indexUD][indexModule]->SetBinContent(i, res);
                  }
                  double mnx,mny,max,may;
                  max=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
          	        GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMaximumBin());
                  if(max>=60 && max<=65){
                    mnx=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                          GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMinimumBin());
                    may=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                          GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(max));
                    mny=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                          GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(mnx));
                    double diff=may-mny;
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
                    kmu2=diff;          ped=mny;
                    dubPed=mny;         tpeak=max;
                    if(indexClock==0 && indexFB==1 && indexUD==0){
                      h1cali->Fill(diff);
                      TSpectrum *sp = new TSpectrum(4);
                      Int_t nf = sp->Search(h1cali,1,"",0.10);
                      double lowRange,upRange;
                      Double_t *xps = sp->GetPositionX();
                      std::sort(xps,xps+nf);
                      lowRange=h1cali->GetBinLowEdge(h1cali->GetMaximumBin()-xps[0]/2.0);
                      upRange=h1cali->GetBinLowEdge(h1cali->GetMaximumBin()+xps[0]*3.0/2.0);
                      h1cali->Fit("gaus","Q","",lowRange,upRange);
        	    }
                  }else{
                    kmu2=-100; calInt=-100;
        	    module=-100;   tpeak=-100;
        	    dubPed=-100;
                  } //<-- Use to get rid of 2 peaks functions here * /
      	        } // <--- End dpulse block
            }else{
              kmu2=-100; phei=-100; calInt=-100; module=-100;
            }
          }
        }else{
          kmu2=-100; phei=-100; calInt=-100; module=-100; tpeak=-100;
        }
        //pCali->Fill();
        //if(iHist>=17) return;
      } // <--- End of if loop
    } // <--- End of nChannel loop
  } // <--- End of No. CsI crystals that fired

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
