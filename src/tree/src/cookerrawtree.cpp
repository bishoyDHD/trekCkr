/** 
 * This file has the implementations of the CRT 
 */

#include "../../../include/cookerrawtree.h"

CRTBase::CRTBase() {};
CRTBase::~CRTBase() {};
ClassImp(CRTBase);

CRTConfigFile::CRTConfigFile() {};
CRTConfigFile::~CRTConfigFile() {};
ClassImp(CRTConfigFile);

CRTRunInfo::CRTRunInfo() {};
CRTRunInfo::~CRTRunInfo() {};
ClassImp(CRTRunInfo);

CRTEventInfo::CRTEventInfo() {};
CRTEventInfo::~CRTEventInfo() {};
ClassImp(CRTEventInfo);


CRTBinaryBlob::CRTBinaryBlob() {
data=NULL;
};
CRTBinaryBlob::~CRTBinaryBlob() {};
ClassImp(CRTBinaryBlob);

EventInfo::EventInfo(){}
EventInfo::~EventInfo(){}
ClassImp(EventInfo);

BeamInfo::BeamInfo(){}
BeamInfo::~BeamInfo(){}
ClassImp(BeamInfo);

SftInfo::SftInfo(){}
SftInfo::~SftInfo(){}
ClassImp(SftInfo);

TargetInfo::TargetInfo(){}
TargetInfo::~TargetInfo(){}
ClassImp(TargetInfo);

MwpcInfo::MwpcInfo(){}
MwpcInfo::~MwpcInfo(){}
ClassImp(MwpcInfo);

Tof1Info::Tof1Info(){}
Tof1Info::~Tof1Info(){}
ClassImp(Tof1Info);

Tof2Info::Tof2Info(){}
Tof2Info::~Tof2Info(){}
ClassImp(Tof2Info);

AcInfo::AcInfo(){}
AcInfo::~AcInfo(){}
ClassImp(AcInfo);

GvInfo::GvInfo(){}
GvInfo::~GvInfo(){}
ClassImp(GvInfo);

TtcInfo::TtcInfo(){}
TtcInfo::~TtcInfo(){}
ClassImp(TtcInfo);

PgcInfo::PgcInfo(){}
PgcInfo::~PgcInfo(){}
ClassImp(PgcInfo);

