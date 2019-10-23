#ifndef boostCM_H
#define boostCM_H 1

#include <iostream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <TLorentzVector.h>

class boostCM{
public:
  boostCM();
  ~boostCM();
  void setbxbybz(double betax,double betay,double betaz);
  void setPxPyPzE(double px,double py,double pz,double E);
  double getPx(){return cmpx;};
  double getPy(){return cmpy;};
  double getPz(){return cmpz;};
  double getE(){return E;};
private:
  double cmpx,cmpy,cmpz,E;
  double bx,by,bz;
  double b2,gamma;
};
#endif
