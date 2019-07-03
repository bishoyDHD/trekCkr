/*
 *  This is the definition of low level analysis data structure
 *  CAT = Cooker Analysis Tree
 */
#ifndef __MARINATECSI_H_
#define __MARINATECSI_H_

#include <TObject.h>
#include <TTimeStamp.h>
#include <vector>
#include <map>

/**
 * Base class for all components. We may use this for function prototypes.
 */

class CATBase: public TObject
{
 public:
  CATBase();
  virtual ~CATBase();
  ClassDef(CATBase,1);
};

/* 
 * Out from raw root file -->
 *   Input for to be marinated
 */

class CATSingleCsI:public CATBase{
 public:
  //Var for all pulse types
  //Var for all pulse types
  Int_t ud, fb;
  Double_t ped, phei, calInt, tpeak, tref[3];
  Double_t refpk[3], tcorr[3], refmn[3];
  Double_t thSing, phiSing, trise;
  int crysID, typeAB;
  int indexCsI, clock;
  int csiArrange[2];
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
  /*
  Int_t ud, fb;
  Double_t ped, phei, calInt, tpeak, tref;
  Double_t thSing, phiSing, trise;
  int crysID, typeAB;
  int indexCsI, clock;
  int csiArrange[2];
  double tcsi,phdstr;
  //single peak
  Double_t sphei; // single peak pulse-height distribution
  Double_t sptime; //timing of single peak
  Double_t sped; // pedestal for single pulse
  //Double peak var
  Double_t kmu2, dubPed, intKmu2;
  Double_t dubphei; //location of second peak
  //Overrange variables
  Double_t ovrpH, ovrpLoc, ovrped;*/
  CATSingleCsI();
  virtual ~CATSingleCsI();
  ClassDef(CATSingleCsI,1);
};

class CATClusterCsI:public CATBase{
 public:
  Int_t evtNo, channel;
  int waveID;  // distinguish between 4-different kinds of waves
  Int_t dubP_1; // pre-pile up with double peak
  Double_t E_pi0, M_pi0, pi0M2;
  Double_t M_k, kM2;
  Double_t thetaE, phiE;
  Double_t ggCosTheta, piCosTheta;
  Int_t clusterM; // cluster multiplicity
  Int_t Ncrys; // number of fired crystals
  Int_t ClustCrys; // number of crystals within cluster
  Double_t piPpi0, piP2g;
  // state vector information for 2 gammas
  Double_t g1Px, g1Py, g1Pz;
  Double_t g2Px, g2Py, g2Pz;
  Double_t g1E, g2E;
  // state vector information for pions
  Double_t pi0px, pi0py, pi0pz;
  Double_t pi0vx, pi0vy, pi0vz;
  Double_t piPpx, piPpy, piPpz;
  CATClusterCsI();
  virtual ~CATClusterCsI();
  ClassDef(CATClusterCsI,1);
};

/* 
 * Marinate-er class definition
 */
class CATCaliCsI: public CATBase{
 public:
  Double_t tcmode, edist;
  Double_t thSing, phiSing;
  std::vector<Double_t> icalCrys;
  Double_t ovrpH, ovrpLoc, pHloc, ovrped; //pHloc --> 2nd peak position (no overrange)
  Double_t pkloc;
  CATCaliCsI();
  virtual ~CATCaliCsI();
  ClassDef(CATCaliCsI,1);
};

#endif
