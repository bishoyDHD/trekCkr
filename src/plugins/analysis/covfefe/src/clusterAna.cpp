#include <covfefe.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "TVector2.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"
#include <TCanvas.h>
#include <stdio.h>
#include <math.h>

Long_t covfefe::hist_clust(){
  gStyle->SetOptStat(0);
  clustEval=new evalClusters();
  clustEval->setChannel(7); // must be called before defHistos (defaul chan7)
  clustEval->defHistos();

  return 0;
}

Long_t covfefe::startup_clust(){
  getBranchObject("treeClus",(TObject **) &clsmar); 
  //calibcsi=new CATCaliCsI();
  //makeBranch("marinCsI",(TObject **) &calibcsi);
  gStyle->SetOptStat(0);

  return 0;
};

Long_t covfefe::process_clust(){
  if(clsmar->E_prim2>0){
      // Construct CM Lorentz vectors for 2gammas
      std::cout<<" making sure this thing works "<<clsmar->clusterM<<" "<<clsmar->E_prim2<<" "<<clsmar->prCosTheta<<"\n";
      double piPpx=-1*clsmar->prim1px;
      double piPpy=-1*clsmar->prim1py;
      double piPpz=-1*clsmar->prim1pz;
      double P=std::sqrt(piPpx*piPpx+piPpy*piPpy+piPpz*piPpz);
      //double E=std::sqrt(P*P+mpip*mpip); //total energy of pi0 CM
      double E=0.493677-std::sqrt(P*P+mpip*mpip); //total energy of pi0 CM
      double bx=piPpx/E;
      double by=piPpy/E;
      double bz=piPpz/E;
      double b2=bx*bx+by*by+bz*bz;
      double gamma=1.0/std::sqrt(1.0-b2);
      clustEval->fillHistos(clsmar->M_prim2,clsmar->clCosTheta,clsmar->E_prim2,clsmar->prCosTheta);
      //clustEval->fillHistos(clsmar->cpid1theta,clsmar->cpid1phi,clsmar->cpid2theta,clsmar->cpid2phi,0,0);
      clustEval->fillMltp(clsmar->clusterM);
    //}
  }

  return 0; // 0 = all ok
};

Long_t covfefe::finalize_clust(){
  clustEval->drawHistos();

  return 0; // 0 = all ok
};
