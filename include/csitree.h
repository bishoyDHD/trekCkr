/** 
 * This file has the definition of the ToF tree branch objects 
 */

#ifndef __CSITREE_H_
#define __CSITREE_H_
#include "cookerrawtree.h" // for CRTBase
#include <vector>
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
  Int_t ud, fb;
  Double_t ped, kmu2, phei, calInt, dubPed, tpeak, tref;
  Double_t thSing, phiSing;
  int crysID, typeAB;
  int indexCsI, clock;
  int csiArrange[2];
  CRTSingleCsI();
  virtual ~CRTSingleCsI();
  ClassDef(CRTSingleCsI,1);
};
class CRTClusterCsI:public CRTBase{
 public:
  Int_t clock;
  Int_t fb;
  Int_t ud, event, module, multip;
  Double_t ped, kmu2, phei, calInt, dubPed, tpeak, tref;
  Double_t clusE, thClus, phiClus;
  std::vector<double> xpos, val;
  std::vector<double> csThet, csPhi;
  std::vector<std::pair<double,double>> pos;
  std::vector<int> idCrys, ncrys, nclus;
  std::vector<int> crysID, typeAB, gud, gno, gfb;
  CRTClusterCsI();
  virtual ~CRTClusterCsI();
  ClassDef(CRTClusterCsI,1);
};

#endif

