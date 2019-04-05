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
};

covfefe::~covfefe(){
  for(int i=0; i<20; i++)
    delete calibHist[i];
};

Long_t covfefe::histos(){
  for(int iname=0; iname<20; iname++){
    std::ostringstream name;
    name<<"calibH";
    name<<"crysID"<<"_"<<"ang_";
    name<<theta[iname];
    name<<"_3.75";
    calibHist[iname]=new TH1D(name.str().c_str(),"stat",250,0,1250);
  }
  tpeak=new TH1D("tpeak","stat",250,0,250);
  tref=new TH1D("tref","stat",250,0,250);
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
  if(csimar->thSing!=18.75 && csimar->phiSing==3.75){
    if(csimar->thSing==26.25)
      calibHist[1]->Fill(csimar->kmu2);
    if(csimar->thSing==33.75)
      calibHist[2]->Fill(csimar->kmu2);
    if(csimar->thSing==41.25)
      calibHist[3]->Fill(csimar->kmu2);
    if(csimar->thSing==48.75)
      calibHist[4]->Fill(csimar->kmu2);
    if(csimar->thSing==56.25)
      calibHist[5]->Fill(csimar->kmu2);
    if(csimar->thSing==63.25)
      calibHist[6]->Fill(csimar->kmu2);
  }
  tpeak->Fill(csimar->tpeak);
  tref->Fill(csimar->tref);
  if(csimar->thSing!=18.75 && csimar->phiSing==26.25){
    if(csimar->thSing==26.25)
      calibHist[7]->Fill(csimar->kmu2);
    if(csimar->thSing==33.75)
      calibHist[8]->Fill(csimar->kmu2);
    if(csimar->thSing==41.25)
      calibHist[9]->Fill(csimar->kmu2);
    if(csimar->thSing==48.75)
      calibHist[10]->Fill(csimar->kmu2);
    if(csimar->thSing==56.25)
      calibHist[11]->Fill(csimar->kmu2);
    if(csimar->thSing==63.25)
      calibHist[12]->Fill(csimar->kmu2);
  }

  return 0; // 0 = all ok
};

Long_t covfefe::finalize(){
  TSpectrum *s;
  for(int i=1; i<12; i++){
    s=new TSpectrum(4);
    TF1* f1=calibHist[i]->GetFunction("gaus");
    Int_t nfound = s->Search(calibHist[i],2,"",0.10);
    Double_t *xpeaks = s->GetPositionX();
    std::sort(xpeaks,xpeaks+nfound);
    lowRange=xpeaks[0]/1.5;
    upRange=xpeaks[0]*3.0/2.0;
    calibHist[i]->Fit("gaus","Q","",lowRange,upRange);
  }
  TCanvas* c1=new TCanvas("MarinateCsI","Marinate",3508,2480);
  TCanvas* c2=new TCanvas("MarinateCsI2","Marinate 2",3508,2480);
  c1->Divide(3,2);
  c2->Divide(3,2);
  for(int i=1;i<6;i++){
    c1->cd(i);
    calibHist[i]->Draw();
    c2->cd(i);
    calibHist[i+6]->Draw();
  }
  c1->Write();
  c2->Write();
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
