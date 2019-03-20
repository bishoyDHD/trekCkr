// Creating function to calculate min Chi2
// and determine the fit parameters for waveform
// fitting function

#ifndef mn2CsIfn_H
#define mn2CsIfn_H 1
#include <stdio.h>
#include "TMath.h"
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnHesse.h>
#include <Minuit2/Minuit2Minimizer.h>
using namespace std;
using namespace ROOT::Minuit2;


double model(double x, vector<double> p){
  //double den= exp(1-(x-p[1])/(pow(p[4]+p[6],2)+p[5]*pow(p[4],2))); //<-- original func.
  //double num= p[0]*2*TMath::Freq((x-p[1]-p[2])/p[3]) *((x-p[1])/p[4] * exp(1-(x-p[1])/p[4])+
  //      p[5]*(x-p[1])/(p[4]+p[6])*exp(1-(x-p[1])/(p[4]+p[6])));
  double den= 1-exp(-(x-p[1])/p[6]);
  double num= (p[0]/den)*TMath::Freq((x-p[1]-p[2])/p[3])*((x-p[1])/p[4]*exp(1-(x-p[1])/p[4])+
        p[5]*(x-p[1])/(p[7])*exp(1-(x-p[1])/(p[7])));
  return num+p[8];
};
double overR(double x, vector<double> p){
  //double den= exp(1-(x-p[1])/(pow(p[4]+p[6],2)+p[5]*pow(p[4],2)));
  double den= 1-exp(-(x-p[1])/p[6]);
  double num= p[0]*2*TMath::Freq((x-p[1]-p[2])/p[3]) *((x-p[1])/p[4] * exp(1-(x-p[1])/p[4])+
        p[5]*(x-p[1])/(p[7])*exp(1-(x-p[1])/(p[7])));
  return (num/(den)+p[8]);//*(x<p[10]||x>=p[9]) + 1023*(x>=p[10] && x<p[9]);
};
double model2(double x, vector<double> p){
  double den1= 1-exp(-1*(x-p[1])/p[6]);
  double den2= 1-exp(-1*(x-p[9])/p[6]);
  double frac=p[0]*2*TMath::Freq((x-p[1]-p[2])/p[3]) *((x-p[1])/p[4] *
          exp(1-(x-p[1])/p[4]) +p[5]*(x-p[1])/(p[7])*exp(1-(x-p[1])/(p[7])))/(den1)+
          p[10]*2*TMath::Freq((x-p[9]-p[2])/p[3])*((x-p[9])/p[4] * exp(1-(x-p[9])/p[12]) +
          p[5]*(x-p[9])/(p[7])*exp(1-(x-p[9])/(p[7])))/(den2);
  return frac+p[8];
};/*
double model2(double x, vector<double> p){
  double den1=exp(1-(x-p[1])/(pow(p[4]+p[6],2)+p[5]*pow(p[4],2)));
  double den2=p[13]*p[12]*(p[13]*p[12]+(p[12]+p[14]))/(pow(p[12]+p[14],2)+p[13]*pow(p[12],2))*
                  exp(1 -p[12]*(p[13]*p[12]+(p[12]+p[14]))/(pow(p[12]+p[14],2)+p[13]*pow(p[12],2))) +
                  (p[12]+p[14])*(p[13]*p[12]+(p[12]+p[14]))/(pow(p[12]+p[14],2) + p[13]*pow(p[12],2))*
                  exp(1 - (p[12]+p[14])*(p[13]*p[12]+(p[12]+p[14])) /(pow(p[12]+p[14],2) + p[13]*pow(p[12],2)));
  double frac=p[0]*2*TMath::Freq((x-p[1]-p[2])/p[3]) *((x-p[1])/p[4] *
          exp(1-(x-p[1])/p[4]) +p[5]*(x-p[1])/(p[4]+p[6])*exp(1-(x-p[1])/(p[4]+p[6])))/(den1)+
          p[8]*2*TMath::Freq((x-p[9]-p[10])/p[11])*((x-p[9])/p[12] * exp(1-(x-p[9])/p[12]) +
          p[13]*(x-p[9])/(p[12]+p[14])*exp(1-(x-p[9])/(p[12]+p[14])))/(den2)+p[7];
  return frac;
};*/
//single pulse fitting
class fitfn: public FCNBase{
public:
  fitfn(const vector<double>& xx1, const vector<double>& meas){
    x1=xx1;  y1=meas;
  }
  fitfn(const vector<double>& xx1, double xl, double xr, const vector<double>& meas, double ymax){
    x1=xx1;  y1=meas; xlw=xl;
    xhg=xr;  y2=ymax;
  }
  double operator()(const vector<double> &par)const override{
    double chisq=0.0, chisq1=0.0, chisq2=0.0;
    //cout<< "  --- Making sure that this thing actually works   ------\n";
    for(int i=0; i<250; i+=1){
      chisq += pow(y1[i]-model(x1[i],par),2);
    }
    return chisq;
  }
  double Up() const override { return 0.5; }
private:
  vector<double> x1,y1;
  double xlw=0.,xhg=0.,y2;
};
//double pulse fitting
class fitfn2: public FCNBase{
public:
  fitfn2(const vector<double>& xx1, const vector<double>& meas){
    x1=xx1;  y1=meas;
  }
  fitfn2(const vector<double>& xx1, double xl, double xr, const vector<double>& meas, double ymax){
    x1=xx1;  y1=meas; xlw=xl;
    xhg=xr;  y2=ymax;
  }
  double operator()(const vector<double> &par)const override{
    double chisq=0.0, chisq1=0.0, chisq2=0.0;
    //cout<< "  --- Making sure that this thing actually works   ------\n";
    for(int i=0; i<250; i+=1){
      chisq += pow(y1[i]-model2(x1[i],par),2);
    }
    return chisq;
  }
  double Up() const override { return 0.5; }
private:
  vector<double> x1,y1;
  double xlw=0.,xhg=0.,y2;
};
//overrange fitting
class ovrfn: public FCNBase{// <---- flat-top cases
public:
  ovrfn(const vector<double>& xx1, const vector<double>& meas){
    x1=xx1;  y1=meas;
  }
  ovrfn(const vector<double>& xx1, double xl, double xr, const vector<double>& meas, double ymax){
    x1=xx1;  y1=meas; xlw=xl;
    xhg=xr;  y2=ymax;
  }
  double operator()(const vector<double> &par)const override{
    double chisq=0.0, chisq1=0.0, chisq2=0.0;
    for(int i=0; i<=xlw; i+=1){
      chisq1 += pow(y1[i]-overR(x1[i],par),2);
    }
    //cout<< "  --- Option2: Making sure that this thing actually works   ------\n";
    for(int i=xhg; i<250; i+=1){
      chisq2 += pow(y1[i]-overR(x1[i],par),2);
    }
    chisq=chisq1+chisq2;
    return chisq;
  }
  double Up() const override { return 0.5; }
private:
  vector<double> x1,y1;
  double xlw=0.,xhg=0.,y2;
};
#endif

// waveform multifitting function
// Fitting function usage depends on certain conditions
std::string singlemodel(){
  char norm[1000],modelname[1000];
/*
  sprintf(norm,"[5]*[4]*([5]*[4]+([4]+[6]))/(([4]+[6])**2+[5]*[4]**2)* \
		  exp(1 - [4]*([5]*[4]+([4]+[6]))/(([4]+[6])**2+[5]*[4]**2)) + \
		  ([4]+[6])*([5]*[4]+([4]+[6]))/(([4]+[6])**2 + [5]*[4]**2)* \
		  exp(1 - ([4]+[6])*([5]*[4]+([4]+[6])) /( ([4]+[6])**2 + [5]*[4]**2))");
*/
  sprintf(norm,"1-exp(-(x-[1])/[6])");

  sprintf(modelname,"[0]*TMath::Freq((x-[1]-[2])/[3])/(%s)*((x-[1])/[4]*exp(1-(x-[1])/[4])+\
        [5]*((x-[1])/([7]))*exp(1-(x-[1])/[7]))+[8]",norm);
  std::string namestring=modelname;
  return namestring;
}

std::string overrangemodel(){
  char norm[1000],modelname[1000];
/*
  sprintf(norm,"[5]*[4]*([5]*[4]+([7]))/(([7])**2+[5]*[4]**2)* \
		  exp(1 - [4]*([5]*[4]+([7]))/(([7])**2+[5]*[4]**2)) + \
		  ([7])*([5]*[4]+([7]))/(([7])**2 + [5]*[4]**2)* \
		  exp(1 - ([7])*([5]*[4]+([7])) /( ([7])**2 + [5]*[4]**2))");
*/
  sprintf(norm,"1-exp(-(x-[1])/[6])");
  sprintf(modelname,"[0]*TMath::Freq((x-[1]-[2])/[3])/(%s)*((x-[1])/[4]*exp(1-(x-[1])/[4])+\
        [5]*((x-[1])/([7]))*exp(1-(x-[1])/[7]))+[8]",norm);/*
  sprintf(modelname,"([0]*2*TMath::Freq((x-[1]-[2])/[3]) *((x-[1])/[4] * exp(1-(x-[1])/[4]) + \
	  [5]*(x-[1])/([7])*exp(1-(x-[1])/([7])))/(%s)+[8])*(x<[10]||x>=[9]) + \
		  1023* (x>=[10] && x<[9])",norm);


  sprintf(modelname,"([0]*2*TMath::Freq((x-[1]-[2])/[3]) *((x-[1])/[4] * exp(1-(x-[1])/[4]) + \
	  [5]*(x-[1])/([7])*exp(1-(x-[1])/([7])))/(%s)+[7])*(x<[10]||x>=[9]) + \
		  1023* (x>=[10] && x<[9])",norm);
*/

  std::string namestring=modelname;
  return namestring;
}

std::string doublemodel(){
  char norm[1000],doublemodel[1000],norm2[1000];/*
  sprintf(norm,"exp(1-(x-[1])/(([4]+[6])**2+[5]*[4]**2))");
  sprintf(norm2,"[13]*[12]*([13]*[12]+([12]+[14]))/(([12]+[14])**2+[13]*[12]**2)* \
		  exp(1 -[12]*([13]*[12]+([12]+[14]))/(([12]+[14])**2+[13]*[12]**2)) + \
		  ([12]+[14])*([13]*[12]+([12]+[14]))/(([12]+[14])**2 + [13]*[12]**2)* \
		  exp(1 - ([12]+[14])*([13]*[12]+([12]+[14])) /( ([12]+[14])**2 + [13]*[12]**2))");
*/
  sprintf(norm,"1-exp(-(x-[1])/[6])");
  sprintf(norm2,"1-exp(-(x-[9])/[6])");
  sprintf(doublemodel,"[0]*2*TMath::Freq((x-[1]-[2])/[3]) *((x-[1])/[4] * \
	  exp(1-(x-[1])/[4]) +[5]*(x-[1])/([7])*exp(1-(x-[1])/([7])))/(%s)+ \
	  [10]*2*TMath::Freq((x-[9]-[2])/[3])*((x-[9])/[4] * exp(1-(x-[9])/[12]) + \
	  [5]*(x-[9])/([7])*exp(1-(x-[9])/([7])))/(%s)+[8]",norm,norm2);

  std::string namestring=doublemodel;
  return namestring;
}

std::string triplemodel(){
  char norm[1000],triplemodel[10000],norm2[1000],norm3[1000];
/*
  sprintf(norm,"[5]*[4]*([5]*[4]+([4]+[6]))/(([4]+[6])**2+[5]*[4]**2)* \
		  exp(1 - [4]*([5]*[4]+([4]+[6]))/(([4]+[6])**2+[5]*[4]**2)) + \
		  ([4]+[6])*([5]*[4]+([4]+[6]))/(([4]+[6])**2 + [5]*[4]**2)* \
		  exp(1 - ([4]+[6])*([5]*[4]+([4]+[6])) /( ([4]+[6])**2 + [5]*[4]**2))");
*/
  sprintf(norm,"exp(1-(x-[1])/(([4]+[6])**2+[5]*[4]**2))");
  sprintf(norm2,"[13]*[12]*([13]*[12]+([12]+[14]))/(([12]+[14])**2+[13]*[12]**2)* \
		  exp(1 -[12]*([13]*[12]+([12]+[14]))/(([12]+[14])**2+[13]*[12]**2)) + \
		  ([12]+[14])*([13]*[12]+([12]+[14]))/(([12]+[14])**2 + [13]*[12]**2)* \
		  exp(1 - ([12]+[14])*([13]*[12]+([12]+[14])) /( ([12]+[14])**2 + [13]*[12]**2))");

  sprintf(norm3,"[20]*[19]*([20]*[19]+([19]+[21]))/(([19]+[21])**2+[20]*[19]**2)* \
		  exp(1 -[19]*([20]*[19]+([19]+[21]))/(([19]+[21])**2+[20]*[19]**2)) + \
		  ([19]+[21])*([20]*[19]+([19]+[21]))/(([19]+[21])**2 + [20]*[19]**2)* \
		  exp(1 - ([19]+[21])*([20]*[19]+([19]+[21])) /( ([19]+[21])**2 + [20]*[19]**2))");

  sprintf(triplemodel,"[0]*2*TMath::Freq((x-[1]-[2])/[3]) *((x-[1])/[4] * \
	  exp(1-(x-[1])/[4]) +[5]*(x-[1])/([4]+[6])*exp(1-(x-[1])/([4]+[6])))/(%s)+ \
	  [8]*2*TMath::Freq((x-[9]-[10])/[11])*((x-[9])/[12] * exp(1-(x-[9])/[12]) + [13]*(x-[9])/([12]+[14])* \
	  exp(1-(x-[9])/([12]+[14])))/(%s) + [15]*2*TMath::Freq((x-[16]-[17])/[18])*((x-[16])/[19] * \
	  exp(1-(x-[16])/[19]) + [20]*(x-[16])/([19]+[21])*exp(1-(x-[16])/([19]+[21])))/(%s)+[7]", \
	  norm,norm2,norm3);

  std::string namestring=triplemodel;
  return namestring;
}

std::string quadruplemodel(){
  char norm[1000],quadruplemodel[10000],norm2[1000],norm3[1000],norm4[1000];
/*
  sprintf(norm,"[5]*[4]*([5]*[4]+([4]+[6]))/(([4]+[6])**2+[5]*[4]**2)* \
  		exp(1 - [4]*([5]*[4]+([4]+[6]))/(([4]+[6])**2+[5]*[4]**2)) + \
  		([4]+[6])*([5]*[4]+([4]+[6]))/(([4]+[6])**2 + [5]*[4]**2)* \
  		exp(1 - ([4]+[6])*([5]*[4]+([4]+[6])) /( ([4]+[6])**2 + [5]*[4]**2))");
*/
  sprintf(norm,"exp(1-(x-[1])/(([4]+[6])**2+[5]*[4]**2))");
  sprintf(norm2,"[13]*[12]*([13]*[12]+([12]+[14]))/(([12]+[14])**2+[13]*[12]**2)* \
		  exp(1 -[12]*([13]*[12]+([12]+[14]))/(([12]+[14])**2+[13]*[12]**2)) + \
		  ([12]+[14])*([13]*[12]+([12]+[14]))/(([12]+[14])**2 + [13]*[12]**2)* \
		  exp(1 - ([12]+[14])*([13]*[12]+([12]+[14])) /( ([12]+[14])**2 + [13]*[12]**2))");

  sprintf(norm3,"[20]*[19]*([20]*[19]+([19]+[21]))/(([19]+[21])**2+[20]*[19]**2)* \
		  exp(1 -[19]*([20]*[19]+([19]+[21]))/(([19]+[21])**2+[20]*[19]**2)) + \
		  ([19]+[21])*([20]*[19]+([19]+[21]))/(([19]+[21])**2 + [20]*[19]**2)* \
		  exp(1 - ([19]+[21])*([20]*[19]+([19]+[21])) /( ([19]+[21])**2 + [20]*[19]**2))");

  sprintf(norm4,"[27]*[26]*([27]*[26]+([26]+[28]))/(([26]+[28])**2+[27]*[26]**2)* \
		  exp(1 -[26]*([27]*[26]+([26]+[28]))/(([26]+[28])**2+[27]*[26]**2)) + \
		  ([26]+[28])*([27]*[26]+([26]+[28]))/(([26]+[28])**2 + [27]*[26]**2)* \
		  exp(1 - ([26]+[28])*([27]*[26]+([26]+[28])) /( ([26]+[28])**2 + [27]*[26]**2))");

  sprintf(quadruplemodel,"[0]*2*TMath::Freq((x-[1]-[2])/[3]) *((x-[1])/[4] * \
	  exp(1-(x-[1])/[4]) +[5]*(x-[1])/([4]+[6])*exp(1-(x-[1])/([4]+[6])))/(%s) + \
	  [8]*2*TMath::Freq((x-[9]-[10])/[11])*((x-[9])/[12] * exp(1-(x-[9])/[12]) + \
	  [13]*(x-[9])/([12]+[14])*exp(1-(x-[9])/([12]+[14])))/(%s)+ \
	  [15]*2*TMath::Freq((x-[16]-[17])/[18])*((x-[16])/[19] * exp(1-(x-[16])/[19]) + \
	  [20]*(x-[16])/([19]+[21])*exp(1-(x-[16])/([19]+[21])))/(%s)+ \
	  [22]*2*TMath::Freq((x-[23]-[24])/[25])*((x-[23])/[26] * exp(1-(x-[23])/[26]) + \
	  [27]*(x-[23])/([26]+[28])*exp(1-(x-[23])/([26]+[28])))/(%s)+[7]",norm,norm2,norm3,norm4);

  std::string namestring=quadruplemodel;
  return namestring;
}

// Parameters for the waveform fitting function
double par(int i){
  double param[] = {1000, 35.76, 26.68, 19.85, 15.83, 0.065, 2.255, 31.21,140,
                    120.5, 800, 28.9, 18.0, 17.1, 0.065, 62.7-17.1,
                    800, 25, 28.9, 18.0, 17.1, 0.065, 62.7-17.1,
                    25, 28.9, 18.0, 17.1, 0.065, 62.7-17.1};
  return param[i];
}
double parlim(int i){
  double param[] = {1e5, 270.1, 60, 45, 20, 1.07, 4, 90, 350,
                    250, 1e5, 80, 70, 50, 1.07, 1000, 250};
  return param[i];
}
double parmin(int i){
  double param[] = {80.0, -900.1, 1, 5, 0., 1e-4, 0., 10,
                    70.1, 15, 10, 10, 10, 1e-4, 10, 20};
  return param[i];
}

double parDouble(int i){
  double param[] = {750, 25, 28.9, 18.0, 17.1, 0.065, 62.7-17.1,119,
	            800, 25, 28.9, 18.0, 17.1, 0.065, 62.7-17.1,
		    800, 25, 28.9, 18.0, 17.1, 0.065, 62.7-17.1,
                    25, 28.9, 18.0, 17.1, 0.065, 62.7-17.1};
  return param[i];
}
double parlimDouble(int i){
  double param[] = {1e5, 30, 80, 70, 50, .07, 1000, 250,
                    1e5, 30, 80, 70, 50, .07, 1000, 250};
  return param[i];
}
double parminDouble(int i){
  double param[] = {0.0, 15, 10, 10, 10, .05, 10, 20,
                    0.0, 15, 10, 10, 10, .05, 10, 20};
  return param[i];
}
std::string nameL(int i){
  std::string parN[] = {"p0", "p1", "p2", "p3", "p4", "p5", "p6",
	  "p7", "p8", "p9", "p10", "p11", "p12", "p13", "p14",
          "p15", "p16"};
  return parN[i];
}

double step(int i){
  double ssize[]{20, 2, 2, 2, 2, 1e-3, 2, 2, 2, 2, 2};
  return ssize[i];
}

double overfit(int i){
  double param[] = {1000, 25, 28.9, 18.0, 17.1, 0.065, 62.7-17.1, 120, 50, 70};
  return param[i];
}
