#include <covfefe.h>

#include <iostream>
#include <cmath>
#include "TVector2.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"

#include <stdio.h>
#include <math.h>

covfefe::covfefe(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p){
};

covfefe::~covfefe(){
};

Long_t covfefe::startup(){
  
  return 0;
};

Long_t covfefe::process(){
  int ret=0;

  return 0; // 0 = all ok
};

Long_t covfefe::finalize(){

  return 0; // 0 = all ok
};


extern "C"{
  Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p){
    return (Plugin *) new covfefe(in,out,inf_,outf_,p);
  }
}
