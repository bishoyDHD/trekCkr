#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class std::vector<unsigned short>+;
#pragma link C++ class std::vector<UInt_t>+;
#pragma link C++ class std::pair<unsigned short, std::string>+;
#pragma link C++ class std::map<unsigned short, std::string>+;

#pragma link C++ class CRTBase+;
#pragma link C++ class CRTConfigFile+;
#pragma link C++ class CRTRunInfo+;

#pragma link C++ class CRTEventInfo+;
#pragma link C++ class CRTBinaryBlob+;

#pragma link C++ class CATBase+;
#pragma link C++ class CATCaliCsI+;
#pragma link C++ class CATSingleCsI+;

#pragma link C++ class CRTRawCsI+;
#pragma link C++ class CRTCaliCsI+;
#pragma link C++ class CRTClusterCsI+;
#pragma link C++ class CRTSingleCsI+;

#pragma link C++ class std::map<std::string, double>+;

#pragma link C++ class EventInfo+;
#pragma link C++ class BeamInfo+;
#pragma link C++ class SftInfo+;
#pragma link C++ class TargetInfo+;
#pragma link C++ class MwpcInfo+;
#pragma link C++ class Tof1Info+;
#pragma link C++ class Tof2Info+;
#pragma link C++ class AcInfo+;
#pragma link C++ class GvInfo+;
#pragma link C++ class TtcInfo+;
#pragma link C++ class PgcInfo+;


#endif
