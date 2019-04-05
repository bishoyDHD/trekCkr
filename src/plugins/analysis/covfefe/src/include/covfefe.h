#ifndef __COVFEFE__
#define __COVFEFE__

#include <marinateCsI.h>
#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TSpectrum.h"
#include <TStyle.h>
#include <TH1D.h>
#include <iostream>
#include <map>
#include <iterator>

class covfefe:public Plugin{
 private:
  CATSingleCsI* csimar;
  CATCaliCsI* calibcsi;

  double theta[20]={18.75, 26.25, 33.75, 41.25, 48.75, 56.25, 63.75, 71.25, 78.75, 86.25,
                    93.75, 101.25, 108.75, 116.25, 123.75, 131.25, 138.75, 146.25, 153.75, 161.25};
  double phi[12][2][2]={{{3.75, 11.25},{18.75, 26.25}},
                        {{33.75,41.25},{48.75, 56.25}},
                      {{63.75,71.25},{78.75, 86.25}},
                      {{93.75,101.25},{108.75,116.25}},
                      {{123.75,131.25},{138.75,146.25}},
                      {{153.75,161.25},{168.75,176.25}},
                      {{183.75,191.25},{198.75,206.25}},
                      {{213.75,221.25},{228.75,236.25}},
                      {{243.75,251.25},{258.75,266.25}},
                      {{273.75,281.25},{288.75,296.25}},
                      {{303.75,311.25},{318.75,326.25}},
                      {{333.75,341.25},{348.75,356.25}}};
  double lowRange, upRange;
 public:
  covfefe(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~covfefe();

  TH1D* calibHist[20];
  TH1D* tpeak;
  TH1D* tref;
  // add funtions with return value Long_t here:
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t finalize();

  virtual Long_t cmdline(char * cmd);

  ClassDef(covfefe,1);
};

#endif
