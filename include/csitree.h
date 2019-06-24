/** 
 * This file has the definition of the ToF tree branch objects 
 */

#ifndef __CSITREE_H_
#define __CSITREE_H_
#include "cookerrawtree.h" // for CRTBase
#include <vector>
//static const double E_kpi2=0.10854566040017; // in GeV
class CRTCaliCsI:public CRTBase{
 public:
  UInt_t runNo;
  UInt_t eventNo;
  Int_t isBad;
  UInt_t nChannel;
  std::vector<UInt_t> indexCsI;
  std::vector<float> baseline;
  std::vector<float> tdc;
  std::vector<float> adc;
  std::vector<float> peak;
  CRTCaliCsI();
  virtual ~CRTCaliCsI();
  ClassDef(CRTCaliCsI,1);
};
class CRTRawCsI:public CRTBase{
 public:
  UInt_t runNo;
  UInt_t eventNo;
  Int_t isBad;
  UInt_t nChannel;
  std::vector<UInt_t> nameModule;
  std::vector<UInt_t> indexChannel;
  std::vector<UInt_t> nameCsI;
  std::vector<UInt_t> indexCsI;
  std::vector<ULong64_t> timeStamp;
  std::vector<UInt_t> timeCFD;
  std::vector<ULong64_t> charge;
  std::vector<UInt_t> nSample;
  std::vector<std::vector<UShort_t> > data;  
  CRTRawCsI();
  virtual ~CRTRawCsI();
  ClassDef(CRTRawCsI,1);
};
class CRTSingleCsI:public CRTBase{
 public:
  //Var for all pulse types
  Int_t ud, fb;
  Double_t ped, phei, calInt, tpeak, tref[3];
  Double_t refpk[3], tcorr[3], refmn[3];
  Double_t thSing, phiSing, trise;
  int crysID, typeAB;
  int indexCsI, clock;
  int csiArrange[2];
  int waveID;  // distinguish between 3-different kinds of waves
  double phdstr;
  //single peak
  Double_t sphei; // single peak pulse-height distribution
  Double_t sptime; //timing of single peak
  Double_t sped; // pedestal for single pulse
  //Double peak var
  Double_t kmu2, dubPed, intKmu2;
  Double_t dubphei; //location of second peak
  //Overrange variables
  Double_t ovrpH, ovrpLoc, ovrped;
  CRTSingleCsI();
  virtual ~CRTSingleCsI();
  ClassDef(CRTSingleCsI,1);
};
class CRTClusterCsI:public CRTBase{
 public:
  Int_t evtNo, channel;
  int waveID;  // distinguish between 4-different kinds of waves
  Int_t dubP_1; // pre-pile up with double peak
  Double_t E_pi0;
  double thetaE, phiE;
  Double_t ggCosTheta, piCosTheta;
  Int_t clusterM; // cluster multiplicity
  Int_t Ncrys; // number of fired crystals
  Int_t ClustCrys; // number of crystals within cluster
  Double_t piPpi0, piP2g;
  // state vector information for 2 gammas
  Double_t g1Px, g1Py, g1Pz;
  Double_t g2Px, g2Py, g2Pz;
  // state vector information for pions
  Double_t pi0px, pi0py, pi0pz;
  Double_t piPpx, piPpy, piPpz;/*
  std::vector<double> xpos, val;
  std::vector<double> csThet, csPhi;
  std::vector<std::pair<double,double>> pos;
  std::vector<int> idCrys, ncrys, nclus;
  std::vector<int> crysID, typeAB, gud, gno, gfb;*/
  CRTClusterCsI();
  virtual ~CRTClusterCsI();
  ClassDef(CRTClusterCsI,1);
};

#endif

