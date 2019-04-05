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
  Int_t ud, fb;
  Double_t ped, kmu2, phei, calInt, dubPed, tpeak, tref;
  Double_t thSing, phiSing;
  Double_t ovrpH, ovrpLoc, pHloc, ovrped; //pHloc --> 2nd peak position (no overrange)
  int crysID, typeAB;
  int indexCsI, clock;
  int csiArrange[2];
  CATSingleCsI();
  virtual ~CATSingleCsI();
  ClassDef(CATSingleCsI,1);
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
