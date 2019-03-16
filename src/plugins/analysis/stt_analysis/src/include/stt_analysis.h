#ifndef __stt_analysis__
#define __stt_analysis__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <map>
#include <iterator>
#include "StrawTubetree.h"
#include "StrawTubeanalysistree.h"



class stt_analysis:public Plugin
{
 protected:
  TH1D *straw[stt_num];
  
 private:
  StrawTube *rawstt; 
  
 public:
  stt_analysis(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~stt_analysis();

  // add funtions with return value Long_t here:
  
  Long_t defineHistograms();
  Long_t startup();
  Long_t process();
  Long_t finalize();

  virtual Long_t cmdline(char * cmd);

  ClassDef(stt_analysis,1);
};

#endif
