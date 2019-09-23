// This is a class that does cluster scoring as well as 
// calculating combination sums of clusters
// Some function methods are superfluous because the 
// orignal author was curious about such structures

#ifndef clusterScore_h
#define clusterScore_h 1

#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <math.h>
#include <numeric>
#include <map>
// container class for cluster variables
struct clusterVar{
  // cluster particle pid are clus1 etc...
  // vectors are need in order to return index
  std::string prmid,mcid,scid;
  double mcpx, mcpy, mcpz, mcE;
  double clpx, clpy, clpz, clE;
  double clx, cly, clz, clt, clr;
  double prpx, prpy, prpz, prE;
  double scTheta, scPhi;
  double mcTheta, mcPhi;
  // 4 and 3-vector definitions
  TLorentzVector cpidlv;
  TVector3 pi0v3, cpidv3;
};
// Inherit cluster vars here for purposes of container storage
class clusterScore{
public:
  clusterScore();
  ~clusterScore();
  // Really the scoring function could be effective with taking only the vector
  // I really wanted to try a pass by map function,
  // hence the slight complication
  void scoring(std::vector<double> &invmass);
  void scoring(UInt_t size,std::vector<double> &invmass);
  void scoring(UInt_t size,std::vector<double> &invmass, std::map<double,double> mass);
  // funtion to merge 2 different vectors
  // will be used for heterogeneous cluster types
  void mergeVect(std::vector<double>& clust1,std::vector<double>& clust2);
  void setKey(double key);
  void setE(const std::vector<double> &totE);
  // calculate the Energy, inv. mass, angles, position, momentum etc.
  //  from given configuration of clusters
  void clusterEval(std::vector<double> &mcCrys,std::vector<double> &scCrys,std::vector<double> r1,std::vector<double> z1,std::vector<double> &theta1,std::vector<double> &phi1,std::vector<double>& r2,std::vector<double>& z2,std::vector<double>& th2,std::vector<double>& phi2);
  // function for many crystal cluster case:
  // function for single crystal cluster case:
  void clusterEval(const std::vector<double> &mCrys,const std::vector<double> r,const std::vector<double> z,const std::vector<double> &mtheta,const std::vector<double> &mphi);
  void setOpangClus(const std::vector<double> &theta,const std::vector<double> &phi);
  void setOpangPrimary(const std::vector<double> &theta,const std::vector<double> &phi);
  // set cluster PID
  void setCpid(int clustPID){cpid=clustPID;};
  // set methods for various cluster variables
  // this if for tracked particle as well as cluster variables
  void setprPx(std::vector<clusterVar> &px);
  void setprPy(std::vector<clusterVar> &py);
  void setprPz(std::vector<clusterVar> &pz);
  void setprE(std::vector<clusterVar> &Energy);
  void setclPx(std::vector<clusterVar> &px);
  void setclPy(std::vector<clusterVar> &py);
  void setclPz(std::vector<clusterVar> &pz);
  void setclE(std::vector<clusterVar> &Energy);
  void setclAngle(const std::map<double,double> angle); 
  void setScoreMass(double);
  // get methods for various variables
  double getInvMass();
  inline double getE(){return Etot;} //return total energy
  double getOpAngleClust();
  double getOpAnglePrimary();
  double getprPx();
  double getprPy();
  double getprPz();
  TLorentzVector getprimLV();
  double getprE();
  double getclPx();
  double getclPy();
  double getclPz();
  double getclR();
  double getclZ();
  double getclTime();
  double getclE();
  void init(); // initialize
  void reset(); // reset container
protected:
  // variables to push items to container
  double px,py,pz,E;
  double x,y,z,t,r;
  double theta,phi;
  TLorentzVector particlelv;
  TVector3 particle3v;
private:
  int clustEvalNo;
  int cpid; // determines which value std::map returns
  const int dummy=-1000;
  double particleM;
  double mass; // will deterine scoring scale
  double diffMass;
  std::vector<clusterVar> cvars;
  clusterVar clustvar;

  std::vector<double> mdiff;
  // merge manyCrys and singleCrys clusters into single vector
  std::vector<double> mergeclust;
  std::map<double,double> clustE;
  std::map<double,double> InvMass;
  std::map<double,double> openAng;
  std::map<double,TLorentzVector> primLV;
  // Scoring variable storage for 1st 2 clusters
  std::map<double,std::pair<double,double>> csipx;
  std::map<double,std::pair<double,double>> csipy;
  std::map<double,std::pair<double,double>> csipz;
  std::map<double,std::pair<double,double>> csix;
  std::map<double,std::pair<double,double>> csiy;
  std::map<double,std::pair<double,double>> csiz;
  std::map<double,std::pair<double,double>> csir;
  std::map<double,std::pair<double,double>> csiE;
  std::map<double,std::pair<double,double>> csitheta;
  std::map<double,std::pair<double,double>> csiphi;
  std::map<double,std::pair<double,double>> csitime;
  // Scoring variable storage for additional clusters
  std::map<double,std::pair<double,double>> csipx_;
  std::map<double,std::pair<double,double>> csipy_;
  std::map<double,std::pair<double,double>> csipz_;
  std::map<double,std::pair<double,double>> csiE_;
  std::map<double,std::pair<double,double>> csitheta_;
  std::map<double,std::pair<double,double>> csiphi_;
  std::map<double,std::pair<double,double>> csitime_;
  std::map<double,std::pair<double,double>> csix_;
  std::map<double,std::pair<double,double>> csiy_;
  std::map<double,std::pair<double,double>> csiz_;
  std::map<double,std::pair<double,double>> csir_;
  double mkey, ival;
  double invMass;
  double clprpx;
  double clprpy;
  double clprpz;
  double clscpx;
  double clscpy;
  double clscpz;
  double clmcpx;
  double clmcpy;
  double clmcpz;
  double clmcE;
  double clscE;
  double clprE;
  double scTheta;
  double mcTheta;
  double scPhi;
  double mcPhi;
  double Etot, energy,opAngle;
};
#endif
