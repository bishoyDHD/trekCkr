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
  Double_t refpk[3], tcorr[3],rgaus[3],refmn[3];
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
  // timing determination for CsI(Tl):
  Double_t tpeak, tref[3],rgaus[3];
  Double_t refpk[3], tcorr[3], refmn[3];
  Double_t trise;
  // cluster variables
  Int_t evtNo, channel;
  int waveID;  // distinguish between 4-different kinds of waves
  Int_t dubP_1; // pre-pile up with double peak
  Double_t E_prim1, M_prim1, prim1M2;
  Double_t E_prim2, M_prim2, prim2M2;
  Double_t M_k, kM2;
  double cpid1thetaE, cpid1phiE;
  double cpid2thetaE, cpid2phiE;
  Double_t clCosTheta, prCosTheta;
  Int_t clusterM; // cluster multiplicity
  Int_t Ncrys; // number of fired crystals
  Int_t ClustCrys; // number of crystals within cluster
  Double_t Clus2M, Clus2E, Clus2gAng, Clus2piAng, Clus2;
  Double_t Clus1M, Clus1E, Clus1gAng, Clus1piAng, Clus1;
  // state vector information for 2 gammas
  Double_t cpid1Px, cpid1Py, cpid1Pz;
  Double_t cpid2Px, cpid2Py, cpid2Pz;
  Double_t cpid1E, cpid2E;
  Double_t cpid1theta, cpid2theta, cpid1phi, cpid2phi;
  // position of cluster particles
  Double_t cpid1x, cpid1y, cpid1z;
  Double_t cpid2x, cpid2y, cpid2z;
  Double_t cpid1r, cpid2r;
  // state vector information for primary particles
  // ---> i.e the tracked charged particle & reconstructed particle
  Double_t prim1px, prim1py, prim1pz;
  Double_t prim2px, prim2py, prim2pz;
  /*
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

