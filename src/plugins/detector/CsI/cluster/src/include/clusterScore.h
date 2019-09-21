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
  double scoring(UInt_t size,std::vector<double> &invmass, std::map<double,double> mass);
  void setIndex(int ind);
  void setE(const std::vector<double> &totE);
  // calculate the Energy of given combination of clusters
  // returns combination sum of cluster configurations
  void clusterEval(const std::vector<double> &mcCrys,const std::vector<double> &scCrys,const std::vector<double> r,const std::vector<double> z,const std::vector<double> &theta,const std::vector<double> &phi);
  // function for many crystal cluster case:
  // function for single crystal cluster case:
  void clusterEval(const std::vector<double> &mCrys,const std::vector<double> r,const std::vector<double> z,const std::vector<double> &mtheta,const std::vector<double> &mphi);
  void setOpangClus(const std::vector<double> &theta,const std::vector<double> &phi);
  void setOpangPrimary(const std::vector<double> &theta,const std::vector<double> &phi);
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
  void setScoreMass(double);
  // get methods for various variables
  double getprPx();
  double getprPy();
  double getprPz();
  double getprE();
  double getclPx();
  double getclPy();
  double getclPz();


  double getclE();
  inline double getE(){return Etot;} //return total energy
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
  const int dummy=-1000;
  double particleM;
  double mass; // will deterine scoring scale
  double diffMass;
  std::vector<clusterVar> cvars;
  clusterVar clustvar;

  std::vector<double> mdiff;
  std::map<double,double> clustE;
  std::map<double,double> InvMass;
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
  int index, ival;
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
  double Etot, energy;
};
#endif
