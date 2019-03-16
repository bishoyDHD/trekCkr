#ifndef __DET_CSI__
#define __DET_CSI__
#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "csitree.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <utility>
#include <map>
typedef std::pair<UInt_t,UInt_t> IdCsI;
class Det_CsI:public Plugin
{
 private:

  CRTCaliCsI *treeCali;			/// Output tree for CSI data
  CRTRawCsI *treeRaw;			/// Input tree with CSI raw data
  //
  // Detector parameters set in init file
  //

  // constants 

  std::map<IdCsI,UInt_t> mapCsI;

 public:
  Det_CsI(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_);
  virtual ~Det_CsI();

  // Main CsI analysis methods
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t done();

  // CsI fit
  Long_t histos_fit();
  Long_t startup_fit();
  Long_t process_fit();
  Long_t finalize_fit();
  Long_t setIdCsI(std::map<IdCsI,UInt_t>& mapCsI);
  Long_t set_goodEvents(int, int);
  std::vector<std::pair<int,int> > listGoodEvent;
  //histograms for fit

  TH2D* h2TdcVSCsI;
  TH2D* h2AdcVSCsI;
  TH1D* h1Skipped;
  TH2D* h2BaseVSCsI;  

  virtual Long_t cmdline(char * cmd);

  
  
  ClassDef(Det_CsI,1);
};

#endif
