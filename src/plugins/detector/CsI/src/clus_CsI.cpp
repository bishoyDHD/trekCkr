#include <Det_CsI.h>
#include <clus_var.h>
#include<iostream>
#include<cmath>

Long_t Det_CsI::histos_clus(){
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
  for(int i=0;i<11;i++){
    hbox1[i]=new TLine(-2.0,4*(i+1),4.00,4*(i+1));
    vbox1[i]=new TLine(-2.0,4*(i+1),-2.0,4*(i+1)+2);
    hbox1[i+11]=new TLine(-2.0,4*(i+1)+2,4.00,4*(i+1)+2);
    vbox1[i+11]=new TLine(4.0,4*(i+1),4.00,4*(i+1)+2);
    hline1[i]=new TLine(-9.0,4*(i+1)-1,11.0,4*(i+1)-1);
    hline1[i+11]=new TLine(-9.0,4*(i+1)-1,11.0,4*(i+1)-1);
    hline2[i]=new TLine(-9.0,4*(i+1)+1,-2.0,4*(i+1)+1);
    hline2[i+11]=new TLine(-9.0,4*(i+1)+1,-2.0,4*(i+1)+1);
    hline3[i]=new TLine(4.00,4*(i+1)+1,11.0,4*(i+1)+1);
    hline3[i+11]=new TLine(4.00,4*(i+1)+1,11.0,4*(i+1)+1);
  }
  hline2[22]=new TLine(-9.0,1,-2.0,1);
  hline2[23]=new TLine(4.00,1,11.0,1);
  hline2[24]=new TLine(-9.0,49,-2.0,49);
  hline2[25]=new TLine(4.00,49,11.0,49);
  hbox2[0]=new TLine(-2.0,2.0,4,2.0);
  hbox2[1]=new TLine(-2.0,48.,4,48.);
  vbox2[0]=new TLine(-2.0,1.0,-2,2.0);
  vbox2[1]=new TLine(4.00,1.0,4.,2.0);
  vbox2[2]=new TLine(-2.0,48.,-2,49.);
  vbox2[3]=new TLine(4.00,48.,4.,49.);
  hline1[22]=new TLine(-9.0,47,11.0,47);
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
    hline2[n+22]->Draw("l");
  }
  hline1[22]->Draw("l");
  return 0;
}

Long_t Det_CsI::startup_clus(){
  getBranchObject("vf48",(TObject **) &treeRaw);
  getBranchObject("RawBeamInfo",(TObject **) &treeBeam);
  gStyle->SetOptStat(0);

  return 0;
}

Long_t Det_CsI::process_clus(){
  phval=new vector<double>();
  clusth=new vector<double>();
  clusphi=new vector<double>();
  int iHist=0;
  clus_csi=false;
  idCrys.clear(),  indexph.clear();
  typeAB.clear(),  gud.clear(),     gno.clear();
  csThet.clear(),  csPhi.clear();
  csiph.clear();   phval->clear();   csiClus.clear();
  clusth->clear(); clusphi->clear();
  if(treeRaw->nChannel>=7){ // Start by checking how many CsI crystals have fired
    for(UInt_t i=0;i<treeRaw->nChannel;i++){
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
      
      if(treeRaw->indexCsI[i]==16){
        //if(indexClock==0 && indexFB==0 && indexUD==0){
          for(UInt_t iData=0;iData<treeRaw->nSample[i];iData++){
            h1time[indexClock][indexFB][indexUD][indexModule]->SetBinContent(iData+1,treeRaw->data[i][iData]);
          }
        //}
      }
      if((treeRaw->indexCsI[i]!=16) /*&& (indexClock==0 && indexFB==1 && indexUD==0)*/){
        //std::cout<< " Index clock: "<<indexClock<<endl;
        //std::cout<< " Gap config FB is  : " <<p[1]<<endl;
        //std::cout<< " Gap config UD is  : " <<p[0]<<endl;
        //std::cout<< " size of nChannel is : " <<treeRaw->indexCsI[i]-1<<endl;
        //std::cout<< " size of nSample is  : " <<treeRaw->nSample[i]<<endl;
        //std::cout<< " Chan No.: "<<treeRaw->nChannel<<endl;
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
        //cout<<" The baseline is: "<<bl<<endl;
        for(int i=1; i<=h1Fits[indexClock][indexFB][indexUD][indexModule]->GetNbinsX(); i+=1){
          xpos.push_back(i);
          val.push_back(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetBinContent(i));
        }
        if(x1>=57 && x1<=68){
          //cout<< " evtNo: "<<ev<<endl;
          //cluster finding algorithm: Finding neighbours!
          //gud=0, typeAB=0,gno=indexClock, fb=indexFB; crysID=indexModule;
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
            FunctionMinimum min = migrad(180,1e-2);
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
              clus_csi=true;
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
              double area=integral-(mny*250);
              h1Intg->Fill(area);
              calInt=area;
	      std::cout<<" =====> Value for integral1 is:" <<area<<endl;
              double diff=may-mny; 
              int imod=0, igap=4*indexClock+2, igapU=4*indexClock+3;
              int csimod=(-1*treeRaw->indexCsI[i])+1;  //default
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
              if(p[0]=='u' || p[0]=='U'){
                igap=4*indexClock+3;
                if(treeRaw->indexCsI[i]>=10){
                  igap=4*(indexClock+1);
                  csimod=(-1*treeRaw->indexCsI[i])+7;
                  if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
                }
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
              auto angles=std::make_pair(csitheta,csiphi);
              csThet.push_back(csitheta), csPhi.push_back(csiphi);
              //thetaPhi.push_back(angles);
              csiph[angles]=diff;
              csiClus[angles]=true;
              //cout<< " *********** theta "<<csitheta<<"  "<<csiphi<<endl;
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
          } //<-- end of overrange if loop
          if(y1<1023){
            //cout<< "  Size of x is:  "<<xpos.size()<<endl;
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
              h1Fits[indexClock][indexFB][indexUD][indexModule]->Fit(f1,"0");
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
              FunctionMinimum min = migrad(180,1e-2);
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
                clus_csi=true;
                mnx=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinLowEdge(h1Fits[indexClock][indexFB][indexUD][indexModule]->GetMinimumBin());
                may=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(max));
                mny=h1Mnft[indexClock][indexFB][indexUD][indexModule]->
                      GetBinContent(h1Fits[indexClock][indexFB][indexUD][indexModule]->FindBin(mnx));
                double diff=may-mny;
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
      	        auto angles=std::make_pair(csitheta,csiphi);
      	        //thetaPhi.push_back(angles);
      	        csThet.push_back(csitheta), csPhi.push_back(csiphi);
      	        csiph[angles]=diff;
      	        //csiph[thetaPhi]=diff;
      	        csiClus[angles]=true;
      	        //cout<< " *********** theta "<<csitheta<<"  "<<csiphi<<endl;
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
                if(p[0]=='u' || p[0]=='U'){
      	          igap=4*indexClock+3;
      	          if(treeRaw->indexCsI[i]>=10){
      	            igap=4*(indexClock+1);
      	            csimod=(-1*treeRaw->indexCsI[i])+7;
      	            if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
      	          }
      	          if(p[1]=='b' || p[1]=='B'){
      	            csimod=treeRaw->indexCsI[i];
      	            if(treeRaw->indexCsI[i]>=10){
                      igap=4*(indexClock+1);
                      csimod=(treeRaw->indexCsI[i])-6;
      	              if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
                    }
      	          }
      	        } 
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
            }
            if(nfound==1){
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
              FunctionMinimum min = migrad(40,1e-3);
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
                clus_csi=true;
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
                h1Intg->Fill(area);  tpeak=max;
                calInt=area;
                double diff=may-mny;
      	        indexph.push_back(diff);
      	        phval->push_back(diff);
      	        idCrys.push_back(indexModule), gfb.push_back(fb);
      	        int tAB=0;
      	        if((p[0]=='u' || p[0]=='U') && (indexModule>=9)) tAB=1;
      	        if((p[0]=='d' || p[0]=='D') && (indexModule<=8)) tAB=1;
      	        typeAB.push_back(tAB), gno.push_back(indexClock), gud.push_back(ud);
      	        double csitheta=theta[fb][indexModule], csiphi=phi[indexClock][ud][tAB];
      	        auto angles=std::make_pair(csitheta,csiphi);
      	        //thetaPhi.push_back(angles);
      	        csThet.push_back(csitheta), csPhi.push_back(csiphi);
      	        csiph[angles]=diff;
      	        //csiph[thetaPhi]=diff;
      	        csiClus[angles]=true;
      	        //cout<< " *********** theta "<<csitheta<<"  "<<csiphi<<endl;
      	        int imod=0, igap=4*indexClock+2, igapU=4*indexClock+3;
                int csimod=(-1*treeRaw->indexCsI[i])+1;  //default
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
                if(p[0]=='u' || p[0]=='U'){
      	          igap=4*indexClock+3;
      	          if(treeRaw->indexCsI[i]>=10){
      	            igap=4*(indexClock+1);
      	            csimod=(-1*treeRaw->indexCsI[i])+7;
      	            if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
      	          }
      	          if(p[1]=='b' || p[1]=='B'){
      	            csimod=treeRaw->indexCsI[i];
      	            if(treeRaw->indexCsI[i]>=10){
                      igap=4*(indexClock+1);
                      csimod=(treeRaw->indexCsI[i])-6;
      	              if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
	            }
      	          }
      	        } 
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
                    int imod=0, igap=4*indexClock+2, igapU=4*indexClock+3;
                    int csimod=(-1*treeRaw->indexCsI[i])+1;  //default
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
                    if(p[0]=='u' || p[0]=='U'){
                      igap=4*indexClock+3;
                      if(treeRaw->indexCsI[i]>=10){
                        igap=4*(indexClock+1);
                        csimod=(-1*treeRaw->indexCsI[i])+7;
                        if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
                      }
                      if(p[1]=='b' || p[1]=='B'){
                        csimod=treeRaw->indexCsI[i];
                        if(treeRaw->indexCsI[i]>=10){
                          igap=4*(indexClock+1);
                          csimod=(treeRaw->indexCsI[i])-6;
                          if(treeRaw->indexCsI[i]==16) igap=4*indexClock+3.5;
                        }
                      }
                    }
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
            } /*else{
              kmu2=-100; phei=-100; calInt=-100; module=-100;
            }*/
          }
        }else{
          kmu2=-100; phei=-100; calInt=-100; module=-100; tpeak=-100;
        }
        //pCali->Fill();
      } // <--- End of if loop
    } // <--- End of nChannel loop
    //treeRaw->Clear();
    if(clus_csi){
      std::cout<<"\n ***************************************************************************\n";
      for(UInt_t cr=0;cr<phval->size();cr++){
        //int idFB=gfb[cr]; int idCryst=idCrys[cr];
        //int gNo=gno[cr], gUD=gud[cr], tyAB=typeAB[cr];
        //cout<<" pulse height order is hcrys["<<cr<<"] = "<<theta[idFB][idCryst]<<" "<<phi[gNo][gUD][tyAB]<<" "<<hcrys[cr]<<": Gap No."<<gno[cr]<<endl;
        std::cout<<" angles of corr. pulse heigh["<<cr<<"] = "<<(*phval)[cr]<<endl;
      }
      //std::map<std::pair<double,double>,bool>::iterator itrr;
      for(auto itrr=csiClus.begin();itrr!=csiClus.end();itrr ++){
        std::cout<<"  cluster bool has the following: "<<itrr->second<<std::endl;
      }
      double ntheta, nphi;
      int ii=0;
      for(std::size_t mm=0; mm !=indexph.size(); mm++){
        const auto index=get_nthIndex(indexph, mm);
        std::cout<<"  the greater index --> "<<index
                << " with value "<<indexph[index]<<std::endl;
        //std::nth_element(begin(hcrys), begin(hcrys)+ii, end(hcrys));
        //cout<<"   checking value of nth element finder :"<<hcrys[ii]<<endl;
        //auto phval=std::find(begin(refpulse), end(refpulse),hcrys[ii]);
        //double loci=std::distance(begin(refpulse),phval);
        //cout<<" ****  checking the heck out of this| "<<itr->first[ii].first<<"\t"<<itr->first[ii].second<<"\t"<<itr->second<<endl;
        ntheta=csThet[index], nphi=csPhi[index];
        //erpair.clear();
        auto tppair=std::make_pair(ntheta,nphi);
	std::cout<<"   ==>  theta and phi "<<ntheta<<"  "<<nphi<<endl;
        //cout<<"  Pulse-height for given theta phi pair: "<< csiph[tpair]<<endl;
        auto angP1=std::make_pair(ntheta+7.5,nphi);      auto angP2=std::make_pair(ntheta-7.5,nphi);
        auto angP3=std::make_pair(ntheta,nphi+7.5);      auto angP4=std::make_pair(ntheta,nphi-7.5);
        auto angP5=std::make_pair(ntheta+7.5,nphi+7.5);  auto angP6=std::make_pair(ntheta-7.5,nphi-7.5);
        auto angP7=std::make_pair(ntheta-7.5,nphi+7.5);  auto angP8=std::make_pair(ntheta+7.5,nphi-7.5);
        clusCrys=0;
        if(csiph[angP1] > 0 &&  csiph[angP2] > 0 &&  csiph[angP3] > 0 && csiph[angP4] > 0 &&
          csiph[angP5] > 0 && csiph[angP6] > 0 && csiph[angP7] > 0 &&
          csiph[angP8] > 0){
          std::cout<<"  Total number of cluster crystals is 8 \n";
        }else if(csiph[angP1] > 0 || csiph[angP2] > 0 || csiph[angP3] > 0 || csiph[angP4] > 0 || 
          csiph[angP5] > 0 || csiph[angP6] > 0 || csiph[angP7] > 0 || 
          csiph[angP8] > 0){
          //clusCrys=clusCrys+1;
          if(csiph[angP1]>0){
            std::cout<<" This crystal Cluster pulse-height P1: "<<csiClus[angP1]<<std::endl;
            if(csiClus[angP1]){
              clusCrys=clusCrys+1;
              std::cout<<" This crystal is now removed from the list: "<<std::endl;
              csiClus[angP1]=false;
            }else{
              std::cout<<" Already checked this crystal \n";
            }
            //csiph.erase(erpair); //csiph.erase(itr);
            //erpair.clear();
          }
          if(csiph[angP2]>0){
            std::cout<<" This crystal Cluster finder in pair loop P2: "<<csiClus[angP2]<<std::endl;
            if(csiClus[angP2]){
              clusCrys=clusCrys+1;
              std::cout<<" This crystal is now removed from the list: "<<std::endl;
              csiClus[angP2]=false;
            }else{
              std::cout<<" Already checked this crystal \n";
            }
          }
          if(csiph[angP3]>0){
            std::cout<<" This crystal Cluster finder in pair loop P3: "<<csiClus[angP3]<<std::endl;
            if(csiClus[angP3]){
              clusCrys=clusCrys+1;
              std::cout<<" This crystal is now removed from the list: "<<std::endl;
              csiClus[angP3]=false;
            }else{
              std::cout<<" Already checked this crystal \n";
            }
            //csiph.erase(erpair); //csiph.erase(itr);
          }
          if(csiph[angP4]>0){
            std::cout<<" This crystal Cluster finder in pair loop P4: "<<csiClus[angP4]<<std::endl;
            if(csiClus[angP4]){
              clusCrys=clusCrys+1;
              std::cout<<" This crystal is now removed from the list: "<<std::endl;
              csiClus[angP4]=false;
            }else{
              std::cout<<" Already checked this crystal \n";
            }
          }
          if(csiph[angP5]>0){
            std::cout<<" This crystal Cluster finder in pair loop P5: "<<csiClus[angP5]<<std::endl;
            if(csiClus[angP5]){
              clusCrys=clusCrys+1;
              std::cout<<" This crystal is now removed from the list: "<<std::endl;
              csiClus[angP5]=false;
            }else{
              std::cout<<" Already checked this crystal \n";
            }
          }
          if(csiph[angP6]>0){
            std::cout<<" This crystal Cluster finder in pair loop P6: "<<csiClus[angP6]<<std::endl;
            if(csiClus[angP6]){
              clusCrys=clusCrys+1;
              std::cout<<" This crystal is now removed from the list: "<<std::endl;
              csiClus[angP6]=false;
            }else{
              std::cout<<" Already checked this crystal \n";
            }
          }
          if(csiph[angP7]>0){
            std::cout<<" This crystal Cluster finder in pair loop P7: "<<csiClus[angP7]<<std::endl;
            if(csiClus[angP7]){
              clusCrys=clusCrys+1;
              std::cout<<" This crystal is now removed from the list: "<<std::endl;
              csiClus[angP7]=false;
            }else{
              std::cout<<" Already checked this crystal \n";
            }
          }
          if(csiph[angP8]>0){
            std::cout<<" This crystal Cluster finder in pair loop P8: "<<csiClus[angP8]<<std::endl;
            if(csiClus[angP8]){
              clusCrys=clusCrys+1;
              std::cout<<" This crystal is now removed from the list: "<<std::endl;
              csiClus[angP8]=false;
            }else{
              std::cout<<" Already checked this crystal \n";
            }
          }
          std::cout<<" some crystals actually have hits "<<clusCrys<<std::endl;
        }else{
          clusCrys=clusCrys+1;
          std::cout<<" number of single cluster crystals is "<<clusCrys<<std::endl;
          /*if(std::find(thetaPhi.begin(), thetaPhi.end(), angP2) != thetaPhi.end()){
            std::cout<<"  Checking DeMorgan's laws in C++ \n";
          }*/
        }
        if(csiClus[tppair]){
        }
        csiClus[tppair]=false; // mute central crystal
        /*
        cout<<"  Number of crystals is:  "<<clusCrys<<endl;
        ii++;
        //itr++, ii++;*/
      }
      for(UInt_t idc=0;idc<csThet.size();idc++){
        std::cout<<"    theta[m][n] "<<csThet[idc]<<endl;
        //cout<<"    pairs pHeig "<<csiph.find(thetaPhi)->second<<endl;
      }
      /*
      for(UInt_t id=0;id<idCrys.size();id++){
        int gNo=gno[id], gUD=gud[id], tyAB=typeAB[id];
        cout<<"    phi[gNo][UD][type] "<<phi[gNo][gUD][tyAB]<<endl;
      }*/
      //cout<<"  ~~~~~~~~~~  making sure this works:  "<<csiph[thetaPhi]/*.find(thetaPhi)->second*/<<endl;
      std::cout<<" ***************************************************************************\n";
    }//<--- end cluster if loop
  } // <--- End of No. CsI crystals that fired
  else{
    //continue;
  } // <--- else to No. CsI crystals that fired if loop
  //clus_csi=false;
  delete phval;
  delete clusth;
  delete clusphi;
}

Long_t Det_CsI::finalize_clus(){
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
  return 0;
}
