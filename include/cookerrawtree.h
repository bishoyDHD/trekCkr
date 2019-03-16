/*
 * This is the definition of the cooker raw data structures
 * CRT = Cooker Root Tree
 */

#ifndef __COOKERRAWTREE_H__
#define __COOKERRAWTREE_H__


#include <TObject.h>
#include <TTimeStamp.h>
#include <vector>
#include <map>

/**
 * Base class for all components. We may use this for function prototypes.
 */

class CRTBase: public TObject
{
 public:
  CRTBase();
  virtual ~CRTBase();
  ClassDef(CRTBase,1);
};


/**
 * A class to store daq config files
 */

class CRTConfigFile:public CRTBase
{
 public:
  std::string filename;
  std::string content;
  CRTConfigFile();
  virtual ~CRTConfigFile();
  ClassDef(CRTConfigFile,1);
};



/**
 * A class to store information of the current run
 */

class CRTRunInfo:public CRTBase
{
	public:
		UInt_t runNumber;

		TTimeStamp startTime;
		TTimeStamp stopTime;
		UInt_t nrOfEvents;
		UInt_t activeLEVBS;
		// currently unused values
		Double_t targetFlow;
		Double_t magnetCurrent;
		Int_t beamSpecies; // -1=Electrons , 0=Cosmics , 1=Positrons

		CRTRunInfo();
		virtual ~CRTRunInfo();
		ClassDef(CRTRunInfo,2);
};

// Per event info 

class CRTEventInfo:public CRTBase
{
	public:
		UInt_t eventNumber;
		UInt_t eventType;
		TTimeStamp eventTime;
		UInt_t trigFired;
		UInt_t eventCounter;
		UInt_t resetCounter;
		UInt_t liveTimeCounter;
		UInt_t freeClock;
		UInt_t deadTimeCounter;
		double weight;
		CRTEventInfo();
		virtual ~CRTEventInfo();
		ClassDef(CRTEventInfo,2);
};

/**
 * A Class saving a binary blob, with a pointer and a size
 */

class CRTBinaryBlob
{
 public:
  Int_t size; 
  UChar_t *data; //[size]
  CRTBinaryBlob();
  virtual ~CRTBinaryBlob();
  ClassDef(CRTBinaryBlob,1);
};

class EventInfo:public CRTBase{
  /*
  const Int_t nEvTag = 14;
  const Int_t nEvFlag = 40;
  */
 public:
  Int_t run;
  Int_t event;
  Byte_t EvTag[14];//nEvTag
  Bool_t EvFlag[40];//nEvFlag
  EventInfo();
  virtual ~EventInfo();
  ClassDef(EventInfo,2);
};
class BeamInfo:public CRTBase{
 public:
  Int_t run;
  Int_t event;
  std::vector<Int_t> TDC_Trig[2];//nVT48Trig=2
  std::vector<Int_t> TDC_Ck[14];//nVT48Ck=14
  std::vector<Int_t> TDC_Cpi[14];//nVT48Cpi=14
  std::vector<Int_t> TDC_Hodo[24];//nHodo=24
  BeamInfo();
  virtual ~BeamInfo();
  ClassDef(BeamInfo,2);

};

class SftInfo:public CRTBase{
  //  const Int_t nFibrs = 128;
 public:
  Int_t run;
  Int_t event;
  Int_t ADC_High_SFT[128];
  Int_t ADC_Low_SFT[128];
  std::vector<Int_t> TDC_LE_SFT[128];
  std::vector<Int_t> TDC_TE_SFT[128];
  SftInfo();
  virtual ~SftInfo();
  ClassDef(SftInfo,2);


};
class TargetInfo:public CRTBase{
  //  const Int_t nBars = 256;
 public:
  Int_t run;
  Int_t event;
  Int_t ADC_High_TARGET[256];
  Int_t ADC_Low_TARGET[256];
  std::vector<Int_t> TDC_LE_TARGET[256];
  std::vector<Int_t> TDC_TE_TARGET[256];
  TargetInfo();
  virtual ~TargetInfo();
  ClassDef(TargetInfo,2);
  
};

class MwpcInfo:public CRTBase{
  /*
  const Int_t nMWPCADC = 512;
  const Int_t nC2X = 56;
  const Int_t nC2Y = 16;
  const Int_t nC3X = 64;
  const Int_t nC3Y = 16;
  const Int_t nC4X = 72;
  const Int_t nC4Y = 16;
  */
 public:
  Int_t run;
  Int_t event;
  Int_t ADC_C2X_L[56];
  Int_t ADC_C2X_R[56];
  Int_t ADC_C2Y_L[16];
  Int_t ADC_C2Y_R[16];
  Int_t ADC_C3X_L[64];
  Int_t ADC_C3X_R[64];
  Int_t ADC_C3Y_L[16];
  Int_t ADC_C3Y_R[16];
  Int_t ADC_C4X_L[72];
  Int_t ADC_C4X_R[72];
  Int_t ADC_C4Y_L[16];
  Int_t ADC_C4Y_R[16];
  MwpcInfo();
  virtual ~MwpcInfo();
  ClassDef(MwpcInfo,2);

};

class Tof1Info:public CRTBase{
 public:
  Int_t run;
  Int_t event;
  Int_t ADC_TOF1U[12];
  Int_t ADC_TOF1D[12];
  
  Int_t TDC_TOF1U[12];
  Int_t TDC_TOF1D[12];
  Tof1Info();
  virtual ~Tof1Info();
  ClassDef(Tof1Info,2);

};
class Tof2Info:public CRTBase{
 public:
  Int_t run;
  Int_t event;
  Int_t TDC_TOF2AO[12];
  Int_t TDC_TOF2AI[12];
  Int_t TDC_TOF2BO[12];
  Int_t TDC_TOF2BI[12];
  
  Int_t ADC_TOF2AO[12];
  Int_t ADC_TOF2AI[12];
  Int_t ADC_TOF2BO[12];
  Int_t ADC_TOF2BI[12];
  Tof2Info();
  virtual ~Tof2Info();
  ClassDef(Tof2Info,2);

};

class AcInfo:public CRTBase{
 public:
  Int_t run;
  Int_t event;
  Int_t ADC_ACU[12];
  Int_t ADC_ACD[12];
  std::vector<Int_t> TDC_ACU[12];
  std::vector<Int_t> TDC_ACD[12];
  AcInfo();
  virtual ~AcInfo();
  ClassDef(AcInfo,2);

};

class GvInfo:public CRTBase{
 public:
  Int_t run;
  Int_t event;
  Int_t ADC_GV[12];
  std::vector<Int_t> TDC_GV[12];
  GvInfo();
  virtual ~GvInfo();
  ClassDef(GvInfo,2);

};

// TTC
class TtcInfo:public CRTBase{
 public:
  Int_t run;
  Int_t event;
  Int_t ADC_TTC[12];
  std::vector<Int_t> TDC_TTC[12];
  TtcInfo();
  virtual ~TtcInfo();
  ClassDef(TtcInfo,2);

};
// PGC
class PgcInfo:public CRTBase{
 public:
  Int_t run;
  Int_t event;
  Int_t ADC_PGC_Gap[12][8];
  std::vector<Int_t> TDC_PGC_Gap[12][8];
  PgcInfo();
  virtual ~PgcInfo();
  ClassDef(PgcInfo,2);

};


#endif


