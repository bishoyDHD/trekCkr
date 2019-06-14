/** 
 * This file has the definition of the ToF tree branch objects 
 */

#ifndef __MULTITREE_H_
#define __MULTITREE_H_

#include "cookerrawtree.h" // for CRTBase
#include <vector>
static const double E_kpi2=0.10854566040017; // in GeV
class CRTtrackVar:public CRTBase{
 public:
  UInt_t eventNo;
  CRTtrackVar();
  virtual ~CRTtrackVar();
  ClassDef(CRTtrackVar,1);
};
class trackingE36 : public CRTBase{
public:
  Double_t nxVert, nyVert, nzVert;
  Double_t xVert, yVert, zVert;
  Double_t pVertpi0;
  Int_t tof2Gap, evtNum;
  trackingE36();
  virtual ~trackingE36();
  ClassDef(trackingE36, 2);
};

#endif

