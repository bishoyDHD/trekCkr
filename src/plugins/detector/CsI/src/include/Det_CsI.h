#ifndef __DET_CSI__
#define __DET_CSI__

#include "csitree.h"
#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TSpectrum.h"
#include <TStyle.h>
#include "TAxis.h"
#include "TAxis.h"
#include "TLine.h"
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <utility>
#include <map>
typedef std::pair<UInt_t,UInt_t> IdCsI;
class Det_CsI:public Plugin{
 private:
  CRTCaliCsI *treeCali;                 /// Output tree for CSI data
  CRTRawCsI *treeRaw;                   /// Input tree with CSI raw data
  BeamInfo* treeBeam;
 public:
  double m1, m2, x1, x2, y1, ymax, xx1, xx2, yy1, yy2;
  Det_CsI(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~Det_CsI();
  // add funtions with return value Long_t here:
  
  char* pName;
  TSpectrum *s;
  Int_t nfound;
  TH1D* h1Mnft[12][2][2][16];
  TH1D* h1Diff[12][2][2][16];
  TH1D* h1rate[12][2][2][16]; // normalized diff b/n fit and histogram
  TH2D* h2Fits[12][2][2][16];
  TH1D* h1Fits[12][2][2][16];
  TH1D* h1Amps[12][2][2][16];
  TH1D* h1time[12][2][2][16];
  TH1D* h1Pamp, *h1ped;
  TH1D* h1kmu2;
  TH1D* h1cali;
  TH2D* h2clus;
  TH1D* h1Intg, *p0, *p10, *p9, *p1, *f1X2;
  TLine *hbox1[22], *hline1[23];
  TLine *hbox2[2], *hline2[26];
  TLine *vbox1[22], *hline3[22];
  TLine *vbox2[4];
  Int_t clock;
  Int_t fb;
  Int_t ud, event, module, multip;
  Double_t ped, kmu2, phei, calInt, dubPed, tpeak, tref, f1chi2;
  Double_t clusE, thClus, phiClus;
  std::vector<double> xpos, val, csThet, csPhi;
  std::vector<int> idCrys, ncrys, nclus;
  std::vector<int> crysID, typeAB, gud, gno, gfb;
  bool clus_csi;
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t finalize();
  // clustering methods
  Long_t histos_clus();
  Long_t startup_clus();
  Long_t process_clus();
  Long_t finalize_clus();
  Long_t setIdCsI(std::map<IdCsI,UInt_t>& mapCsI);
  Long_t set_goodEvents(int, int);
  std::vector<std::pair<int,int> > listGoodEvent;

  virtual Long_t cmdline(char * cmd);

  ClassDef(Det_CsI,1);
};
#endif
