#ifndef __Det_ClusterCsI__
#define __Det_ClusterCsI__

#include "csitree.h"
#include "multiTree.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include <fstream>
#include "TSpectrum.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TCanvas.h"
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
#include <algorithm>
#include <numeric>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#define _USE_MATH_DEFINES

typedef std::pair<UInt_t,UInt_t> IdCsI;
class Det_ClusterCsI:public Plugin{
 private:
  CRTCaliCsI *treeCali; // Output branch for CSI data
  CRTRawCsI *treeRaw;   // Input tree with CSI raw data
  BeamInfo* treeBeam;
  trackingE36* tracktree; // Input from tracking
  //CRTClusterCsI *treeClus; // Output branch for CSI cluster var
  CRTClusterCsI *treeClus;  // Output branch for CSI single hit var
  TCanvas* c1;
  bool resetH, notfire;
  const double M_pi0=0.1349766;
  const double M_piP=0.13957018;
  const double E_kpi2=0.2455612; //total energy for Kpi2
  double px[4], py[4], pz[4];
  double T_pi0, ppip, thetaE, phiE;
  double rtheta, rphi;
  //double masses[3];
 public:
  Det_ClusterCsI(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~Det_ClusterCsI();
  double m1, m2, x1, x2, y1, ymax, xx1, xx2, yy1, yy2, minx;
  double valx1, valx2, calib;
  std::ifstream paramFile;
  // add funtions with return value Long_t here:
  
  int clusCrys;
  boost::unordered_map<std::pair<double,double>,double, boost::hash<std::pair<double,double>>> csiph;
  boost::unordered_map<std::pair<double,double>,bool, boost::hash<std::pair<double,double>>> csiClus;
  //enum tk {TK34,TK32,TK36,TK37,TK38,TK40,TK08,TK50,TK09,TK54,TK31,TK04,TK45,TK33,TK39,TK41};
  int moduleNo;
  // indices start at zero now
  int thetaCsI[16][48];
  int phiCsI[16][48];
  // mapped angles from Osaka Group
  double mapTheta[20]={18.75, 26.25, 33.75, 41.25, 48.75, 56.25,63.75, 71.25, 78.75, 86.25,
	               93.75, 101.25, 108.75, 116.25, 123.75, 131.25, 138.75, 146.25, 153.75,161.25};
  double mapPhi;
  // angles to be used by the clusterFinder table
  double otheta, ophi;
  //arranged theta[fb][crystalNo. 0-15]: f=0,b=1
  double theta[2][16]={86.25, 78.75, 71.25, 63.75, 56.25, 48.75, 41.25, 33.75, 26.25,
                                            63.75, 56.25, 48.75, 41.25, 33.75, 26.25, 18.75,
                       93.75, 101.25, 108.75, 116.25, 123.75, 131.25, 138.75, 146.25, 153.75,
                                            116.25, 123.75, 131.25, 138.75, 146.25, 153.75, 161.25};
  double phi[12][2][2]={{{63.75,71.25},{78.75, 86.25}}, //g1 60
                        {{33.75,41.25},{48.75, 56.25}}, //g2 30
		        {{3.75, 11.25},{18.75, 26.25}},  //g3 -0
                        {{333.75,341.25},{348.75,356.25}}, //g4 330
                        {{303.75,311.25},{318.75,326.25}}, //g5 300
                        {{273.75,281.25},{288.75,296.25}}, //g6 270
                        {{243.75,251.25},{258.75,266.25}}, //g7 240
                        {{213.75,221.25},{228.75,236.25}}, // g8 210
                        {{183.75,191.25},{198.75,206.25}}, // g9 180
                        {{153.75,161.25},{168.75,176.25}}, // g10 150
                        {{123.75,131.25},{138.75,146.25}}, // g11 120
                        {{93.75,101.25},{108.75,116.25}}}; // g12 90
  /*
  double phi[12][2][2]={{{303.75,311.25},{318.75,326.25}},
                      {{333.75,341.25},{348.75,356.25}},
		      {{3.75, 11.25},{18.75, 26.25}},
                      {{33.75,41.25},{48.75, 56.25}},
                      {{63.75,71.25},{78.75, 86.25}},
                      {{93.75,101.25},{108.75,116.25}},
                      {{123.75,131.25},{138.75,146.25}},
                      {{153.75,161.25},{168.75,176.25}},
                      {{183.75,191.25},{198.75,206.25}},
                      {{213.75,221.25},{228.75,236.25}},
                      {{243.75,251.25},{258.75,266.25}},
                      {{273.75,281.25},{288.75,296.25}}};

  // phi back
  double phiB[12][2][2]={{{123.75,131.25},{138.75,146.25}},   //Iclock_1  120
                        {{153.75,161.25},{168.75,176.25}},   //Iclock_2  150
                        {{183.75,191.25},{198.75,206.25}},   //Iclock_3  180
                        {{213.75,221.25},{228.75,236.25}},   //Iclock_4  210
                        {{243.75,251.25},{258.75,266.25}},   //Iclock_5  240
                        {{273.75,281.25},{288.75,296.25}},   //Iclock_6  270
                        {{303.75,311.25},{318.75,326.25}},   //Iclock_7  300
                        {{333.75,341.25},{348.75,356.25}},   //Iclock_8  330
		        {{3.75, 11.25},{18.75, 26.25}},      //Iclock_9  0
                        {{33.75,41.25},{48.75, 56.25}},      //Iclock_10 30
	                {{63.75,71.25},{78.75, 86.25}},      //Iclock_11 60
                        {{93.75,101.25},{108.75,116.25}}};   //Iclock_12 90
  // phi front
  double phiF[12][2][2]={{{93.75,101.25},{108.75,116.25}},    //Iclock_1 90
	                {{63.75,71.25},{78.75, 86.25}},      //Iclock_2 60
                        {{33.75,41.25},{48.75, 56.25}},      //Iclock_3 30
		        {{3.75, 11.25},{18.75, 26.25}},      //Iclock_4 -0
                        {{333.75,341.25},{348.75,356.25}},   //Iclock_5 330
                        {{303.75,311.25},{318.75,326.25}},   //Iclock_6 300
                        {{273.75,281.25},{288.75,296.25}},   //Iclock_7 270
                        {{243.75,251.25},{258.75,266.25}},   //Iclock_8 240
                        {{213.75,221.25},{228.75,236.25}},   //Iclock_9 210
                        {{183.75,191.25},{198.75,206.25}},   //Iclock_10 180
                        {{153.75,161.25},{168.75,176.25}},   //Iclock_11 150
	                {{123.75,131.25},{138.75,146.25}}};  //Iclock_12 120*/
  double phiCh16[12]={60,30,0,330,300,270,240,210,180,150,120,90};


  char* pName;
  TSpectrum *s;
  Int_t nfound;
  std::ofstream outFile;
  std::ifstream parfile;
  TH1D* h1Mnft[12][2][2][16];
  TH1D* h1Diff[12][2][2][16];
  TH1D* h1Fits[12][2][2][16];
  TH1D* h1time[12][2][2][16];
  double ovrpar[10]={1023.*4, 35.76, 26.68, 19.85, 15.83, 0.065, 2.255, 31.21,120,  120.5};
  TH1D* h1Mpi0, *h1Mpi02, *pi0Etot, *E2g;
  TH1D* h1pi0px, *h1pi0py, *h1pi0pz, *h1clust, *h1sclus;
  TH1D* h1vertpx, *h1vertpy, *h1vertpz;
  TH2D* h2clus, *h2Ene, *h2ang, *h2deg;
  TH2D* h2corrAng, *h2theta, *h2phi, *h2Ang;
  TH1D* E_cut, *cosTheta, *vertOp;
  TLine *hbox1[25], *hline1[25];
  TLine *hbox2[2], *hline2[26];
  TLine *vbox1[25], *hline3[25];
  TLine *vbox2[4];
  Int_t clock;
  Int_t fb;
  double calibpar[12][2][2][16];
  Int_t ud, event, module, multip;
  Double_t ped, kmu2, phei, calInt, dubPed, tpeak, tref, f1chi2;
  Double_t clusE, thClus, phiClus, lowRange, upRange, tsigL;
  Double_t T_ref[3], maxfn[3], minfn[3], cf50[3];
  double pi0px, pi0py, pi0pz;
  double g1px, g1py, g1pz;
  double g2px, g2py, g2pz;
  double piPpx, piPpy, piPpz;
  double pi0x, pi0y, pi0z;
  double pipEtot;
  int evtNum; 
  std::vector<double> xpos, val, csThet, csPhi;
  std::vector<int> idCrys, ncrys, nclus;
  std::vector<int> crysID, typeAB, gud, gno, gfb;
  std::vector<double> clusEne,singleEne;
  std::vector<double> clusThetaE, clusPhiE;
  double Eclus;
  bool clus_csi;
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t finalize();
  Long_t angleCsI(int id, int module, int channel, int yy, int zz);
  void empty();
  // clustering methods
  Long_t histos_clus();
  Long_t startup_clus();
  Long_t process_clus();
  Long_t finalize_clus();
  Long_t setIdCsI(std::map<IdCsI,UInt_t>& mapCsI);
  Long_t set_goodEvents(int, int);
  std::vector<std::pair<int,int> > listGoodEvent;
  void initVar();
  void readFiles();

  virtual Long_t cmdline(char * cmd);

  ClassDef(Det_ClusterCsI,1);
};
#endif
