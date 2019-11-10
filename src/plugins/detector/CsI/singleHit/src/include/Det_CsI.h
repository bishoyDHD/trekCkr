#ifndef __DET_CSI__
#define __DET_CSI__

#include "csitree.h"
#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSpectrum.h"
#include <TStyle.h>
#include "TAxis.h"
#include "TLine.h"
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <utility>
#include <map>
#include <algorithm>
#include <numeric>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>

typedef std::pair<UInt_t,UInt_t> IdCsI;
class Det_CsI:public Plugin{
 private:
  CRTCaliCsI *treeCali; // Output branch for CSI data
  CRTRawCsI *treeRaw;   // Input tree with CSI raw data
  CRTRawCsI *treeClus;   // Input tree with CSI raw data
  BeamInfo* treeBeam;
  //CRTClusterCsI *treeClus; // Output branch for CSI cluster var
  CRTSingleCsI *treeSing;  // Output branch for CSI single hit var
  bool loopX,csiT, firedCsI, notfire;
  std::string refT();
  std::ostringstream name, name2, name3, name4, name5, name6, tname;
 public:
  double m1, m2, x1, x2, y1, ymax, xx1, xx2, yy1, yy2, minx;
  double valx1, valx2;
  Det_CsI(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~Det_CsI();
  // add funtions with return value Long_t here:
  
  int clusCrys;
  boost::unordered_map<std::pair<double,double>,double, boost::hash<std::pair<double,double>>> csiph;
  boost::unordered_map<std::pair<double,double>,bool, boost::hash<std::pair<double,double>>> csiClus;
  //arranged theta[fb][crystalNo. 0-15]: f=0,b=1
  double theta[2][16]={{86.25, 78.75, 71.25, 63.75, 56.25, 48.75, 41.25, 33.75, 26.25,
                                            63.75, 56.25, 48.75, 41.25, 33.75, 26.25, 18.75},
	              {93.75, 101.25, 108.75, 116.25, 123.75, 131.25, 138.75, 146.25, 153.75,
                                            116.25, 123.75, 131.25, 138.75, 146.25, 153.75, 161.25}};
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
  Double_t ped, kmu2, phei, calInt, dubPed, tpeak, tref, f1chi2, dummy;
  Double_t clusE, thClus, phiClus, lowRange, upRange, tsigL;
  Double_t T_ref[3], maxfn[3], minfn[3], cf50[3], refgaus[3];
  std::vector<double> xpos, val, csThet, csPhi;
  std::vector<double> param, parm,err,par;
  std::vector<int> idCrys, ncrys, nclus;
  std::vector<int> crysID, typeAB, gud, gno, gfb;
  int xmin=0, xmax=250;
  int bmin,bmax, integral;
  double mnx,mny,max,may, diff, area;
  double rtime, *xpeaks;
  //analysis methods and classes
  bool clus_csi;
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t finalize();
  //additional methods for an easy life
  void initVar(); //initialize storage variables
  void erase(); // erase histos at end of event
  virtual Long_t cmdline(char * cmd);

  ClassDef(Det_CsI,1);
};
#endif
