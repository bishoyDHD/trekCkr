#ifndef __COVFEFE__
#define __COVFEFE__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include <iostream>
#include <map>
#include <iterator>

class StrawTubeHits;
class TrackHits;

class TH2D;
class covfefe:public Plugin
{
 private:
  StrawTubeHits *STT;
  TrackHits *Tracks;

  TH2D *hough,*circ,*straws,*line,*both,*lineprime,*chiVangle,*chiVangleprime;
  TH1D *hitresidual,*strawresidual,*hitchisq,*strawchisq,*Angle;




 public:
  covfefe(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~covfefe();

  // add funtions with return value Long_t here:
  
  Long_t startup();
  Long_t process();
  Long_t finalize();

  //static void mychi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

  virtual Long_t cmdline(char * cmd);

  ClassDef(covfefe,1);
    };

#endif
