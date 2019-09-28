#ifndef evalClusters_H
#define evalClusters_H

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <TLine.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include "Plugin.h"

class evalClusters{
public:
  evalClusters();
  ~evalClusters();
  void defHistos();
  // fill histos for sanity checks
  void fillHistos(double InvMass,double OpAng1,double Etot,double primOpAng);
  void fillHistos(double InvMass,double Etot,double primOpAng);
  // fill histos for centroid angles
  //void fillHistos(double theta1,double phi1,double theta2,double phi2);
  void fillHistos(double theta1,double phi1,double theta2,double phi2,double theta3,double phi3);
  void drawHistos();
  void drawCanvas(TH1D* h,int val);
  void drawCanvas(TH2D* h,int val);
  void drawCanvas(TH1D* hist1,TH1D* hist2,TH1D* hist3,TH1D* hist4,int val);
  void setChannel(int val); // select channel no. for eval
private:
  int channelNo; //used for iterating TCanvas and selecting No. of cluster to analyze.
  int Ncrys, Nclust;
  std::string title1, title2, title3, title4, title5;
  std::string prname1, prname2, clname1, clname2, clname3;
  double invMass, cpidOpAng, primOpAng, E2clust, E3clust;
  double cpidtheta1, cpidtheta2, cpidtheta3;
  double cpidphi1, cpidphi2, cpidphi3;
  TH1D* clustAng, *primAng, *Eclust, *invM;
  TH1D* cl1E, *cl2E, *cl3E;
  TH2D* thetaPhi;
  // just in case I need this later down the line
  double cpid1x, cpid1y, cpid1z;
  double cpid2x, cpid2y, cpid2z;
  double cpid3x, cpid3y, cpid3z;
};
#endif
