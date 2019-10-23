#include <boostCM.h>

boostCM::boostCM(){

}
boostCM::~boostCM(){

}
void boostCM::setbxbybz(double betax,double betay,double betaz){
  bx=betax;  by=betay;  bz=betaz;
}
void boostCM::setPxPyPzE(double Px,double Py,double Pz, double Energy){
  cmpx=-1*gamma*bx*Energy+(1+(gamma-1)*(bx*bx/b2))*Px+
	  (gamma-1)*(bx*by/b2)*Py+(gamma-1)*(bx*bz/b2)*Pz;
  cmpy=-1*gamma*by*Energy+(gamma-1)*(by*bx/b2)*Px+
	  (1+(gamma-1)*(by*by/b2))*Py+(gamma-1)*(by*bz/b2)*Pz;
  cmpz=-1*gamma*bz*Energy+(gamma-1)*(bz*bx/b2)*Px+
	  (gamma-1)*(bz*by/b2)*Py+(1+(gamma-1)*(bz*bz/b2))*Pz;
}
