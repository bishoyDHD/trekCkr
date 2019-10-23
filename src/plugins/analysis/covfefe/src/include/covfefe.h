#ifndef __COVFEFE__
#define __COVFEFE__

#include <marinateCsI.h>
#include "evalClusters.h"
#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TSpectrum.h"
#include "TLegend.h"
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TLine.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TH1D.h>
#include <iostream>
#include <map>
#include <iterator>

class covfefe:public Plugin{
 private:
  CATSingleCsI* csimar;
  CATClusterCsI* clsmar;
  CATCaliCsI* calibcsi;
  const double dE=143.5;
  const double mpi0=0.1349766;
  const double mpip=0.13957018;
  double sigdE=3.25;
  double T_0[3], tcorr[3], trefCorr[3], cdf50[3], refgaus[3];
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
  double lowRange, upRange, apcsi;
  evalClusters* clustEval;
 public:
  covfefe(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~covfefe();

  int iclock, iModule;
  Int_t iUD, iFB, nbins;
  Double_t adcVal, intVal;
  TH1D* calibHist;
  TH1D* integHist;
  TH1D* tpeak;
  TH1D* Ecorr;
  TH1D* intEn;
  TH1D* phdis;
  TH1D* hkmu2;
  TH1D* timing, *phdistr;
  TH1D* h1time[12][2][2][16];
  TH1D* h1cali[12][2][2][16];
  // add funtions with return value Long_t here:
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t finalize();
  // functions for cluster analysis:
  Long_t hist_clust();
  Long_t startup_clust();
  Long_t process_clust();
  Long_t finalize_clust();
  // histograms for cluster analysis
  TH1D* E_pi0[3], *M_pi0[3]; // pi0 total energy
  TH1D* waveID[2], *clustM[3], *id1;
  TH1D* g1px[2], *g1py[2], *g1pz[2];
  TH1D* g2px[2], *g2py[2], *g2pz[2];
  TH1D* pi0px[2], *pi0py[2], *pi0pz[2];
  TH1D* vertpx[2], *vertpy[2], *vertpz[2];
  TH2D* kmass[2], *h2Angle[2];
  // angles
  TH1D *h1theta[2], *h1phi[2];
  TH1D *angPP[3], *piPang[3], *pi0ang[3];
  TH2D *h2corrAng;
  // timining histograms
  TH1D* h1t_0[3], *h1tcorr[3],*h1trefCorr[3],*h1cdf50[3],*h1refgaus[3];

  virtual Long_t cmdline(char * cmd);

  ClassDef(covfefe,1);
};
#endif
