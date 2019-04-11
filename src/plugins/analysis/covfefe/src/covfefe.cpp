#include <covfefe.h>
#include <iostream>
#include <cmath>
#include "TVector2.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"
#include <TCanvas.h>
#include <stdio.h>
#include <math.h>

covfefe::covfefe(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p){
  calibcsi=0;
  h1time[12][2][2][16]=NULL;
  h1cali[12][2][2][16]=NULL;
  calibHist=NULL;
  Ecorr=NULL;
  integHist=NULL;
  intEn=NULL;
  std::cout<<" checking this shit \n";
};

covfefe::~covfefe(){
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<15;iModule++){
          delete h1time[iClock][iFB][iUD][iModule];
          delete h1cali[iClock][iFB][iUD][iModule];
        }
      }
    }
  }
  delete calibHist;
  delete Ecorr;
  delete integHist;
  delete intEn;
};

Long_t covfefe::histos(){
  std::ostringstream nameCal, ename, nameInt, inEne;
  nameCal<<"CalibCsI";
  ename<<"E_corr";
  nameInt<<"integCsI";
  inEne<<"IntEnergy";
  calibHist=new TH1D(nameCal.str().c_str(),"stat",62.0,0,250);
  Ecorr=new TH1D(ename.str().c_str(),"stat",63.0,0,250);
  integHist=new TH1D(nameInt.str().c_str(),"stat",1875,0,75000);
  intEn=new TH1D(inEne.str().c_str(),"stat",63.0,0,250);
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<15;iModule++){
          std::ostringstream name, name2, name3, name4, name5, name6, tname;
          name<<"IntEne_"; name2<<"Mnfit_"; name3<<"F1fit_"; name4<<"Dhfit_"; name5<<"pHeight", name6<<"Diff";
          tname<<"time";
          name<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          tname<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          h1time[iClock][iFB][iUD][iModule]=new TH1D(tname.str().c_str(),"stat",86.5,0,1300);
          h1cali[iClock][iFB][iUD][iModule]=new TH1D(name.str().c_str(),"stat",1875,0,75000);
        }
      }
    }
  }
  return 0;
}

Long_t covfefe::startup(){
  getBranchObject("treeSing",(TObject **) &csimar); 
  calibcsi=new CATCaliCsI();
  makeBranch("marinCsI",(TObject **) &calibcsi);
  gStyle->SetOptStat(0);

  return 0;
};

Long_t covfefe::process(){
  iclock=csimar->clock-1;
  iModule=csimar->indexCsI-1;
  iUD=csimar->ud; iFB=csimar->fb;
  adcVal=csimar->kmu2;
  intVal=csimar->calInt;
  if(adcVal > 10){
    //std::cout<<" value of clock:  "<<iclock<<std::endl;
    //std::cout<<" value of Module: "<<iModule<<std::endl;
    //std::cout<<" value of iUD:    "<<iUD<<std::endl;
    //std::cout<<" value of iFB:    "<<iFB<<std::endl;
    //std::cout<<" value of adcVal: "<<adcVal<<std::endl;
    h1time[iclock][iFB][iUD][iModule]->Fill(adcVal);
    std::cout<< "  value for integral as seen: "<<intVal<<std::endl;
    if(intVal < 75001)
      h1cali[iclock][iFB][iUD][iModule]->Fill(intVal);
  }

  return 0; // 0 = all ok
};

Long_t covfefe::finalize(){
  TSpectrum *s;
  gStyle->SetOptStat(0);
  double xmax, xx, calpar;
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<15;iModule++){
          xmax=h1time[iClock][iFB][iUD][iModule]->GetMaximumBin();
	  xx=h1time[iClock][iFB][iUD][iModule]->GetXaxis()->GetBinCenter(xmax);
	  nbins=h1time[iClock][iFB][iUD][iModule]->GetXaxis()->GetNbins();
	  lowRange=xx-110;
          upRange=xx+100;
	  TF1* f1=new TF1("f1","gaus",lowRange,upRange);
	  //f1->SetParLimits(0,lowRange,upRange);
	  h1time[iClock][iFB][iUD][iModule]->Fit(f1,"QR");
	  apcsi=f1->GetMaximumX();
	  calpar=dE/apcsi;
	  for(int n=0; n<nbins;n++){
	    double yy=h1time[iClock][iFB][iUD][iModule]->GetBinContent(n);
	    double x=h1time[iClock][iFB][iUD][iModule]->GetBinCenter(n);
	    double xnew=calpar*x;
	    calibHist->Fill(xnew,yy);
	    Ecorr->Fill(xnew+8.9,yy);
	  }
          xmax=h1cali[iClock][iFB][iUD][iModule]->GetMaximumBin();
	  xx=h1cali[iClock][iFB][iUD][iModule]->GetXaxis()->GetBinCenter(xmax);
	  nbins=h1cali[iClock][iFB][iUD][iModule]->GetXaxis()->GetNbins();
	  lowRange=xx-110;
          upRange=xx+100;
	  TF1* f2=new TF1("f2","gaus",lowRange,upRange);
	  //f1->SetParLimits(0,lowRange,upRange);
	  h1cali[iClock][iFB][iUD][iModule]->Fit(f2,"QR");
	  apcsi=f2->GetMaximumX();
	  calpar=dE/apcsi;
	  for(int n=0; n<h1cali[iClock][iFB][iUD][iModule]->GetXaxis()->GetNbins();n++){
	    double yy=h1cali[iClock][iFB][iUD][iModule]->GetBinContent(n);
	    double x=h1cali[iClock][iFB][iUD][iModule]->GetBinCenter(n);
	    double xnew=calpar*x;
	    integHist->Fill(xnew,yy);
	    intEn->Fill(xnew+8.9,yy);
	  }
	  delete f1;
	  //delete f2;
        }
      }
    }
  }
  TCanvas* c1=new TCanvas("MarinateCsI","Marinate",3508,2480);
  TCanvas* c2=new TCanvas("MarinateCsI2","Marinate 2",3508,2480);
  TCanvas* c3=new TCanvas("E_CsI","Energy CsI",808,700);
  c1->Divide(3,4);
  c2->Divide(3,4);
  for(int iClock=0;iClock<1;iClock++){
    for(int iFB=0;iFB<1;iFB++){
      for(int iUD=0;iUD<1;iUD++){
        for(int iModule=0;iModule<12;iModule++){
          c1->cd(iModule+1);
          h1time[iClock][iFB][iUD][iModule]->Draw();
          c2->cd(iModule+1);
          h1time[iClock][iFB][iUD+1][iModule]->Draw();
	}
      }
    }
  }
  c1->Write();
  c2->Write();
  c3->cd();
  Ecorr->SetTitle("CsI Energy for K_{#mu2} ");
  Ecorr->GetXaxis()->SetTitle("T_{#mu} [MeV]");
  Ecorr->SetLineColor(kCyan);
  Ecorr->SetLineWidth(2);
  //double lower=153.4-20, upper=153.4+20;
  //TF1* f2=new TF1("f2","gaus",lower,upper);
  //Ecorr->Fit("f2","QR");
  Ecorr->Draw("hist");
  calibHist->Draw("hist same");
  c3->Write();
  return 0; // 0 = all ok
};

Long_t covfefe::cmdline(char* cmd){
  return 0;
}

extern "C"{
  Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p){
    return (Plugin *) new covfefe(in,out,inf_,outf_,p);
  }
}
