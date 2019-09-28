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
  //clustEval->setChannel(14); // must be called before defHistos (defaul chan7)
  clustEval->defHistos();
  E_pi0[0]=dH1("Epi0", " E_{total}(2#gamma)",75.5,0.,0.60);
  M_pi0[0]=dH1("Mpi0", " M_{#pi^{0}}",50.5,0.,0.4);
  waveID[0]=dH1("waveID", "Waveform ID", 8,0.,8.);
  g1px[0]=dH1("g1px", "px_{#gamma1}",63.5,-.3,.3);
  g1py[0]=dH1("g1py", "py_{#gamma1}",63.5,-.3,.3);
  g1pz[0]=dH1("g1pz", "pz_{#gamma1}",63.5,-.3,.3);
  g2px[0]=dH1("g2px", "px_{#gamma2}",63.5,-.3,.3);
  g2py[0]=dH1("g2py", "py_{#gamma2}",63.5,-.3,.3);
  g2pz[0]=dH1("g2pz", "pz_{#gamma2}",63.5,-.3,.3);
  pi0px[0]=dH1("pi0px", "px_{#pi^{0}}",63.5,-.3,.3);
  pi0py[0]=dH1("pi0py", "py_{#pi^{0}}",63.5,-.3,.3);
  pi0pz[0]=dH1("pi0pz", "pz_{#pi^{0}}",63.5,-.3,.3);
  vertpx[0]=dH1("vertpx", "px_{#pi^{+}}",63.5,-.3,.3);
  vertpy[0]=dH1("vertpy", "py_{#pi^{+}}",63.5,-.3,.3);
  vertpz[0]=dH1("vertpz", "pz_{#pi^{+}}",63.5,-.3,.3);
  h1theta[0]=dH1("theta", "#theta distribution",24.,0,M_PI);
  h1phi[0]=dH1("phi", "#phi distribution",48.,0,2*M_PI);
  h2Angle[0]=dH2("h2Ang","#phi vs #theta", 24.0,0.,M_PI,48.0,0.,2.*M_PI);
  id1=dH1("id1", "Waveform ID", 8,0.,8.);
  clustM[0]=dH1("clustM", "Cluster multiplicity", 12,0.,12.);
  kmass[0]=dH2("kmass", "E_{tot}(#pi^{+}+#pi^{0}) vs. E_{tot}(2#gamma+#pi^{+})", 62.5,0.,1.,62.5,0.,.7);
  // histograms under various cut conditions
  E_pi0[1]=dH1("Epi0_1", " E_{total}(2#gamma)",75.5,0.,0.60);
  M_pi0[1]=dH1("Mpi0_1", " Invariant Mass of #pi^{0} ",50.5,0.,0.4);
  E_pi0[2]=dH1("Epi0_2", " E_{total}(2#gamma)",75.5,0.,0.60);
  M_pi0[2]=dH1("Mpi0_2", " Invariant Mass of #pi^{0} ",50.5,0.,0.4);
  g1px[1]=dH1("g1px_1", "px_{#gamma1}",63.5,-.3,.3);
  g1py[1]=dH1("g1py_1", "py_{#gamma1}",63.5,-.3,.3);
  g1pz[1]=dH1("g1pz_1", "pz_{#gamma1}",63.5,-.3,.3);
  g2px[1]=dH1("g2px_1", "px_{#gamma2}",63.5,-.3,.3);
  g2py[1]=dH1("g2py_1", "py_{#gamma2}",63.5,-.3,.3);
  g2pz[1]=dH1("g2pz_1", "pz_{#gamma2}",63.5,-.3,.3);
  pi0px[1]=dH1("pi0px_1", "px_{#pi^{0}}",63.5,-.3,.3);
  pi0py[1]=dH1("pi0py_1", "py_{#pi^{0}}",63.5,-.3,.3);
  pi0pz[1]=dH1("pi0pz_1", "pz_{#pi^{0}}",63.5,-.3,.3);
  vertpx[1]=dH1("vertpx_1", "px_{#pi^{+}}",63.5,-.3,.3);
  vertpy[1]=dH1("vertpy_1", "py_{#pi^{+}}",63.5,-.3,.3);
  vertpz[1]=dH1("vertpz_1", "pz_{#pi^{+}}",63.5,-.3,.3);
  h1theta[1]=dH1("theta_1", "#theta distribution",24.,0,M_PI);
  h1phi[1]=dH1("phi_1", "#phi distribution",48.,0,2*M_PI);
  h2Angle[1]=dH2("h2Ang_1","#phi vs #theta", 24.0,0.,M_PI,48.0,0.,2.*M_PI);
  clustM[1]=dH1("clustM_1", "Cluster multiplicity", 12,0.,12.);
  kmass[1]=dH2("kmass_1", "E_{tot}(#pi^{+}+#pi^{0}) vs. E_{tot}(2#gamma+#pi^{+})", 62.5,0.,1.,62.5,0.,.7);
  angPP[0]=dH1("angPP_0","Opening angle of #pi^{+}#pi^{0} cos(#theta_{#pi^{+}#pi^{0}})",50.5,-1.1,1.1);
  angPP[1]=dH1("angPP_1","Opening angle of #pi^{+}#pi^{0} cos(#theta_{#pi^{+}#pi^{0}})",50.5,-1.1,1.1);
  h2corrAng=dH2("h2cAng", "cos(#theta_{#pi^{+}#gamma{1}}) vs. cos(#theta_{#pi^{+}#gamma{2}})",25.,-1.,1.,25,-1.,1.);
  piPang[0]=dH1("piPAng1", "cos(#theta_{#pi^{+}}", 45.,-1.5,1.5);
  pi0ang[0]=dH1("pi0Ang0", "cos(#theta_{#gamma#gamma})", 50.5,-1.1,1.1);
  pi0ang[1]=dH1("pi0Ang1", "cos(#theta_{#gamma#gamma})", 50.5,-1.1,1.1);
  pi0ang[2]=dH1("pi0Ang2", "cos(#theta_{#gamma#gamma})", 50.5,-1.1,1.1);

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
  TLorentzVector piPlv, pi0lv, g1lv, g2lv;
  TLorentzVector cm1, cm2, cmg1, cmg2; // CM lorentz vect for 2gammas
  TVector3 cmg1v, cmg2v; // CM vectors for opening angle calc.
  // 3-vector to calculate the angle between pi+ and pi0
  //TVector3 piPv3, pi0v3, ggv3, g1v3, g2v3;
  //TVector3 
  if(clsmar->E_prim2>0){
    waveID[0]->Fill(clsmar->waveID);
    id1->Fill(clsmar->dubP_1);
    //M_pi0[0]->Fill(clsmar->M_pi0);
    clustM[0]->Fill(clsmar->clusterM);
    g1px[0]->Fill(clsmar->cpid1Px);
    g1py[0]->Fill(clsmar->cpid1Py);
    g1pz[0]->Fill(clsmar->cpid1Pz);
    // gamma2
    g2px[0]->Fill(clsmar->cpid2Px);
    g2py[0]->Fill(clsmar->cpid2Py);
    g2pz[0]->Fill(clsmar->cpid2Pz);
    // pi0 vertex
    pi0px[0]->Fill(clsmar->prim2px);
    pi0py[0]->Fill(clsmar->prim2py);
    pi0pz[0]->Fill(clsmar->prim2pz);
    // pi+ vertex
    vertpx[0]->Fill(clsmar->prim1px);
    vertpy[0]->Fill(clsmar->prim1py);
    vertpz[0]->Fill(clsmar->prim1pz);
    //kmass[0]->Fill(clsmar->piP2g,clsmar->piPpi0);
    //h2Angle[0]->Fill(clsmar->thetaE,clsmar->phiE);
    // angles
    //h1theta[0]->Fill(clsmar->thetaE);
    //h1phi[0]->Fill(clsmar->phiE);
    // apply cut in total energy of pi0
    // replot with new cut condition
    //if(clsmar->E_pi0>=.1 && clsmar->E_pi0<=.3){
      //E_pi0[1]->Fill(clsmar->E_pi0);
      clustM[1]->Fill(clsmar->clusterM);
      g1px[1]->Fill(clsmar->cpid1Py);
      g1py[1]->Fill(clsmar->cpid1Px);
      g1pz[1]->Fill(clsmar->cpid1Pz);
      // gamma2
      g2px[1]->Fill(clsmar->cpid2Py);
      g2py[1]->Fill(clsmar->cpid2Px);
      g2pz[1]->Fill(clsmar->cpid2Pz);
      // pi0 vertex
      pi0px[1]->Fill(clsmar->prim2px);
      pi0py[1]->Fill(clsmar->prim2py);
      pi0pz[1]->Fill(clsmar->prim2pz);
      // pi+ vertex
      vertpx[1]->Fill(clsmar->prim1px);
      vertpy[1]->Fill(clsmar->prim1py);
      vertpz[1]->Fill(clsmar->prim1pz);
      //std::cout<<" checking vert: (px, py, pz) ("<<clsmar->piPpx<<", "<<clsmar->piPpy<<", "<<clsmar->piPpz<<") \n";
      //kmass[1]->Fill(clsmar->piP2g,clsmar->piPpi0);
      //h2Angle[1]->Fill(clsmar->thetaE,clsmar->phiE);
      // angles
      //TVector3 piPv3(clsmar->piPpx,clsmar->piPpy,clsmar->piPpz);
      //TVector3 pi0v3(clsmar->pi0px,clsmar->pi0py,clsmar->pi0pz);
      //TVector3 ggv3(clsmar->pi0vx,clsmar->pi0vy,clsmar->pi0vz);
      //TVector3 g1v3(clsmar->g1Px,clsmar->g1Py,clsmar->g1Pz);
      //TVector3 g2v3(clsmar->g2Px,clsmar->g2Py,clsmar->g2Pz);
      // location
      //TVector3 g1pos(clsmar->g1x,clsmar->g1y,clsmar->g1z);
      //TVector3 g2pos(clsmar->g2x,clsmar->g2y,clsmar->g2z);
      // rotate angles
      //pi0v3.RotateZ(M_PI/2);
      //g1v3.RotateZ(M_PI/6);
      //g2v3.RotateZ(M_PI/6);
      //g2v3.RotateY(M_PI);
      //g2v3.RotateZ(-1*M_PI/2);
      //g1lv.SetVect(g1v3);
      //g2lv.SetVect(g2v3);
      //g1lv.SetE(clsmar->g1E);
      //g2lv.SetE(clsmar->g2E);
      //pi0lv=g1lv + g2lv;
      // ---------------------------------------------------------
      // Construct CM Lorentz vectors for 2gammas
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
      // gamma1
      double gam1px=clsmar->cpid1Px;
      double gam1py=clsmar->cpid1Py;
      double gam1pz=clsmar->cpid1Pz;
      double E1=clsmar->cpid1E;
      double g1prx=-1*gamma*bx*E1+(1+(gamma-1)*(bx*bx/b2))*gam1px+
	      (gamma-1)*(bx*by/b2)*gam1py+(gamma-1)*(bx*bz/b2)*gam1pz;
      double g1pry=-1*gamma*by*E1+(gamma-1)*(by*bx/b2)*gam1px+
	      (1+(gamma-1)*(by*by/b2))*gam1py+(gamma-1)*(by*bz/b2)*gam1pz;
      double g1prz=-1*gamma*bz*E1+(gamma-1)*(bz*bx/b2)*gam1px+
	      (gamma-1)*(bz*by/b2)*gam1py+(1+(gamma-1)*(bz*bz/b2))*gam1pz;
      // gamma2
      double gam2px=clsmar->cpid2Px;
      double gam2py=clsmar->cpid2Py;
      double gam2pz=clsmar->cpid2Pz;
      double E2=clsmar->cpid2E;
      double g2prx=-1*gamma*bx*E2+(1+(gamma-1)*(bx*bx/b2))*gam2px+
	      (gamma-1)*(bx*by/b2)*gam2py+(gamma-1)*(bx*bz/b2)*gam2pz;
      double g2pry=-1*gamma*by*E2+(gamma-1)*(by*bx/b2)*gam2px+
	      (1+(gamma-1)*(by*by/b2))*gam2py+(gamma-1)*(by*bz/b2)*gam2pz;
      double g2prz=-1*gamma*bz*E2+(gamma-1)*(bz*bx/b2)*gam2px+
	      (gamma-1)*(bz*by/b2)*gam2py+(1+(gamma-1)*(bz*bz/b2))*gam2pz;
      cmg1=g1lv;
      cmg1.Boost(bx,by,bz); // perform the L. Boost
      cmg2=g2lv;
      cmg2.Boost(bx,by,bz); // perform the L. Boost
      //cmg1v.SetXYZ(cmg1.Px(),cmg1.Py(),cmg1.Pz()); // construct TVector3
      //cmg2v.SetXYZ(cmg2.Px(),cmg2.Py(),cmg2.Pz()); // construct TVector3
      cmg1v.SetXYZ(g1prx,g1pry,g1prz); // construct TVector3
      cmg2v.SetXYZ(g2prx,g2pry,g2prz); // construct TVector3
      // ---------------------------------------------------------
      //double angpiPg1=piPv3.Angle(g1v3);
      //double angpiPg2=piPv3.Angle(g2v3);
      //double angpiPpi0=piPv3.Angle(pi0v3);
      //double ggAng=g1v3.Angle(g2v3);
      //double ggAng2=g1pos.Angle(g2pos);
      //h2corrAng->Fill(std::cos(angpiPg1),std::cos(angpiPg2));
      angPP[1]->Fill(clsmar->prCosTheta);
      //angPP[1]->Fill(std::cos(angpiPpi0));
      //h1theta[1]->Fill(clsmar->thetaE);
      //h1phi[1]->Fill(clsmar->phiE);
      //piPang[0]->Fill(std::cos(cmg1v.Angle(cmg2v)));
      pi0ang[1]->Fill(clsmar->clCosTheta);
      pi0ang[2]->Fill(clsmar->prCosTheta);
      M_pi0[0]->Fill(clsmar->Clus2M);
      M_pi0[2]->Fill(clsmar->Clus1M);
      M_pi0[1]->Fill(clsmar->M_prim2);
      //E_pi0[0]->Fill(clsmar->Clus2E);
      //E_pi0[2]->Fill(clsmar->Clus1E);
      E_pi0[1]->Fill(clsmar->E_prim2);
      clustEval->fillHistos(clsmar->M_prim2,clsmar->clCosTheta,clsmar->E_prim2,clsmar->prCosTheta);
      clustEval->fillHistos(clsmar->cpid1theta,clsmar->cpid1phi,clsmar->cpid2theta,clsmar->cpid2phi,0,0);
    //}
  }

  return 0; // 0 = all ok
};

Long_t covfefe::finalize_clust(){
  TCanvas* c1=new TCanvas("c1","#theta and #phi distr.",908,980);
  TCanvas* c2=new TCanvas("c2","Cluster multiplicity",908,980);
  TCanvas* c3=new TCanvas("c3","Total energy of 2#gamma",808,700);
  TCanvas* c4=new TCanvas("c4","Total energy of K^{+}",808,700);
  TCanvas* c5=new TCanvas("c5","State vector info.",808,700);
  TCanvas* c6=new TCanvas("c6"," Invariant Mass",808,700);
  TCanvas* c7=new TCanvas("c7"," Angle b/n pi+ and pi0 ",808,700);
  TCanvas* c8=new TCanvas("c8"," Angle b/n pi+ and 2gamma ",808,700);
  TCanvas* c9=new TCanvas("c9"," Parameters to check ",808,700);
  clustEval->drawHistos();
  c1->Divide(2,2);
  c1->cd(3);
  gStyle->SetOptStat(0);
  h2Angle[0]->SetLineColor(kBlue);
  h2Angle[0]->GetXaxis()->SetTitle("#theta [rad]");
  h2Angle[0]->GetYaxis()->SetTitle("#phi [rad]");
  TAxis* a=h2Angle[0]->GetXaxis();
  a->SetNdivisions(-502);
  a->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
  a->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi");
  TAxis* a1=h2Angle[0]->GetYaxis();
  a1->SetNdivisions(-502);
  a1->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
  a1->ChangeLabel(-1,-1,-1,-1,-1,-1,"2#pi");
  h2Angle[0]->Draw("colz");
  gStyle->SetOptStat(0);
  // E-cut phi vs. theta 2D plot
  c1->cd(4);
  gStyle->SetOptStat(0);
  h2Angle[1]->GetXaxis()->SetTitle("#theta [rad]");
  h2Angle[1]->GetYaxis()->SetTitle("#phi [rad]");
  TAxis* a_1=h2Angle[1]->GetXaxis();
  a_1->SetNdivisions(-502);
  a_1->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
  a_1->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi");
  TAxis* a11=h2Angle[1]->GetYaxis();
  a11->SetNdivisions(-502);
  a11->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
  a11->ChangeLabel(-1,-1,-1,-1,-1,-1,"2#pi");
  h2Angle[1]->Draw("colz");
  gStyle->SetOptStat(0);
  // theta
  c1->cd(1);
  h1theta[0]->SetLineColor(kBlue);
  h1theta[0]->GetXaxis()->SetTitle("#theta");
  h1theta[0]->GetYaxis()->SetTitle("counts/bin");
  h1theta[0]->SetLineWidth(2);
  TAxis* ax1=h1theta[0]->GetXaxis();
  ax1->SetNdivisions(-502);
  ax1->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
  ax1->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi");
  h1theta[0]->Draw("hist");
  // E-cut theta
  h1theta[1]->SetLineWidth(2);
  h1theta[1]->SetLineColor(kRed);
  h1theta[1]->Draw("hist same");
  gStyle->SetOptStat(0);
  c1->cd(2);
  h1phi[0]->SetLineColor(kBlue);
  h1phi[0]->GetXaxis()->SetTitle("#phi");
  h1phi[0]->GetYaxis()->SetTitle("counts/bin");
  h1phi[0]->SetLineWidth(2);
  TAxis* ax2=h1phi[0]->GetXaxis();
  ax2->SetNdivisions(-502);
  ax2->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
  ax2->ChangeLabel(-1,-1,-1,-1,-1,-1,"2#pi");
  h1phi[0]->Draw("hist");
  // E-cut phi
  h1phi[1]->SetLineWidth(2);
  h1phi[1]->SetLineColor(kRed);
  h1phi[1]->Draw("hist same");
  gStyle->SetOptStat(0);
  c1->Write();

  c2->cd();
  gStyle->SetOptStat(0);
  piPang[0]->SetLineColor(kBlue);
  piPang[0]->GetXaxis()->SetTitle("cos(#theta_{2#gamma})");
  piPang[0]->GetYaxis()->SetTitle("counts/bin");
  piPang[0]->SetLineWidth(2);
  piPang[0]->Draw("hist");
  c2->Write();

  c3->cd();
  gStyle->SetOptStat(0);
  E_pi0[0]->GetXaxis()->SetTitle("E_{#gamma#gamma} [GeV/c]");
  E_pi0[0]->GetYaxis()->SetTitle("counts/bin");
  E_pi0[0]->SetLineWidth(2);
  E_pi0[0]->Draw("hist");
  // cut E_pi0
  E_pi0[1]->SetLineColor(kRed);
  E_pi0[1]->SetLineWidth(2);
  E_pi0[1]->Draw("hist same");
  c3->Write();

  c4->cd();
  gStyle->SetOptStat(0);
  kmass[0]->GetXaxis()->SetTitle("E_{#gamma#gamma} + E_{#pi^{+}} [GeV/c]");
  kmass[0]->GetYaxis()->SetTitle("E_{#pi^{+}} + E_{#pi^{0}} [GeV/c]");
  kmass[0]->SetLineWidth(2);
  kmass[0]->Draw("colz");
  auto tl3=new TLine(.480,0.0,.480,.695);
  auto tl4=new TLine(.510,0.0,.510,.695);
  tl3->SetLineColor(kRed);
  tl4->SetLineColor(kRed);
  tl3->SetLineWidth(3);
  tl4->SetLineWidth(3);
  tl3->Draw();
  tl4->Draw();
  // along y: pi+ cut
  auto tl5=new TLine(0.00,.490,.955,.490);
  auto tl6=new TLine(0.00,.505,.955,.505);
  tl5->SetLineColor(kRed);
  tl6->SetLineColor(kRed);
  tl5->SetLineWidth(3);
  tl6->SetLineWidth(3);
  // possible shower leakage area
  auto tl7=new TLine(.410,0.0,.410,.695);
  auto tl8=new TLine(.457,0.0,.457,.695);
  tl7->SetLineColor(kCyan);
  tl8->SetLineColor(kCyan);
  tl7->SetLineWidth(3);
  tl8->SetLineWidth(3);
  tl3->Draw();
  tl4->Draw();
  tl5->Draw();
  tl6->Draw();
  tl7->Draw();
  tl8->Draw();
  c4->Write();

  c5->Divide(3,4);
  c5->cd(1);
  gStyle->SetOptStat(0);
  g1px[0]->GetXaxis()->SetTitle("px_{#gamma1} [GeV/c]");
  g1px[0]->GetYaxis()->SetTitle("counts/bin");
  g1px[0]->SetLineWidth(2);
  g1px[0]->Draw("hist");
  // E-cut condition for gamma1
  g1px[1]->SetLineWidth(2);
  g1px[1]->SetLineColor(kRed);
  g1px[1]->Draw("hist same");
  c5->cd(2);
  gStyle->SetOptStat(0);
  g1py[0]->GetXaxis()->SetTitle("py_{#gamma1} [GeV/c]");
  g1py[0]->GetYaxis()->SetTitle("counts/bin");
  g1py[0]->SetLineWidth(2);
  g1py[0]->Draw("hist");
  // E-cut condition for gamma1
  g1py[1]->SetLineWidth(2);
  g1py[1]->SetLineColor(kRed);
  g1py[1]->Draw("hist same");
  c5->cd(3);
  gStyle->SetOptStat(0);
  g1pz[0]->GetXaxis()->SetTitle("pz_{#gamma1} [GeV/c]");
  g1pz[0]->GetYaxis()->SetTitle("counts/bin");
  g1pz[0]->SetLineWidth(2);
  g1pz[0]->Draw("hist");
  // E-cut condition for gamma1
  g1pz[1]->SetLineWidth(2);
  g1pz[1]->SetLineColor(kRed);
  g1pz[1]->Draw("hist same");
  // gamma2
  c5->cd(4);
  gStyle->SetOptStat(0);
  g2px[0]->GetXaxis()->SetTitle("px_{#gamma2} [GeV/c]");
  g2px[0]->GetYaxis()->SetTitle("counts/bin");
  g2px[0]->SetLineWidth(2);
  g2px[0]->Draw("hist");
  // E-cut condition for gamma2
  g2px[1]->SetLineWidth(2);
  g2px[1]->SetLineColor(kRed);
  g2px[1]->Draw("hist same");
  c5->cd(5);
  gStyle->SetOptStat(0);
  g2py[0]->GetXaxis()->SetTitle("py_{#gamma2} [GeV/c]");
  g2py[0]->GetYaxis()->SetTitle("counts/bin");
  g2py[0]->SetLineWidth(2);
  g2py[0]->Draw("hist");
  // E-cut condition for gamma2
  g2py[1]->SetLineWidth(2);
  g2py[1]->SetLineColor(kRed);
  g2py[1]->Draw("hist same");
  c5->cd(6);
  gStyle->SetOptStat(0);
  g2pz[0]->GetXaxis()->SetTitle("pz_{#gamma2} [GeV/c]");
  g2pz[0]->GetYaxis()->SetTitle("counts/bin");
  g2pz[0]->SetLineWidth(2);
  g2pz[0]->Draw("hist");
  // E-cut condition for gamma2
  g2pz[1]->SetLineWidth(2);
  g2pz[1]->SetLineColor(kRed);
  g2pz[1]->Draw("hist same");
  c5->cd(7);
  gStyle->SetOptStat(0);
  pi0px[0]->GetXaxis()->SetTitle("px_{#pi^{0}} [GeV/c]");
  pi0px[0]->GetYaxis()->SetTitle("counts/bin");
  pi0px[0]->SetLineWidth(2);
  pi0px[0]->Draw("hist");
  // E-cut condition for pi0px
  pi0px[1]->SetLineWidth(2);
  pi0px[1]->SetLineColor(kRed);
  pi0px[1]->Draw("hist same");
  c5->cd(8);
  gStyle->SetOptStat(0);
  pi0py[0]->GetXaxis()->SetTitle("py_{#pi^{0}} [GeV/c]");
  pi0py[0]->GetYaxis()->SetTitle("counts/bin");
  pi0py[0]->SetLineWidth(2);
  pi0py[0]->Draw("hist");
  // E-cut condition for pi0py
  pi0py[1]->SetLineWidth(2);
  pi0py[1]->SetLineColor(kRed);
  pi0py[1]->Draw("hist same");
  c5->cd(9);
  gStyle->SetOptStat(0);
  pi0pz[0]->GetXaxis()->SetTitle("pz_{#pi^{0}} [GeV/c]");
  pi0pz[0]->GetYaxis()->SetTitle("counts/bin");
  pi0pz[0]->SetLineWidth(2);
  pi0pz[0]->Draw("hist");
  // E-cut condition for pi0pz
  pi0pz[1]->SetLineWidth(2);
  pi0pz[1]->SetLineColor(kRed);
  pi0pz[1]->Draw("hist same");
  c5->cd(10);
  gStyle->SetOptStat(0);
  vertpx[1]->GetXaxis()->SetTitle("px_{#pi^{+}} [GeV/c]");
  vertpx[1]->GetYaxis()->SetTitle("counts/bin");
  vertpx[1]->SetLineWidth(2);
  vertpx[1]->Draw("hist");
  c5->cd(11);
  gStyle->SetOptStat(0);
  vertpy[1]->GetXaxis()->SetTitle("py_{#pi^{+}} [GeV/c]");
  vertpy[1]->GetYaxis()->SetTitle("counts/bin");
  vertpy[1]->SetLineWidth(2);
  vertpy[1]->Draw("hist");
  c5->cd(12);
  gStyle->SetOptStat(0);
  vertpz[1]->GetXaxis()->SetTitle("pz_{#pi^{+}} [GeV/c]");
  vertpz[1]->GetYaxis()->SetTitle("counts/bin");
  vertpz[1]->SetLineWidth(2);
  vertpz[1]->Draw("hist");
  c5->Write();

  c6->cd();
  gStyle->SetOptStat(0);
  M_pi0[0]->GetXaxis()->SetTitle("M_{#pi^{0}} [GeV/c]");
  M_pi0[0]->GetYaxis()->SetTitle("counts/bin");
  M_pi0[0]->SetLineWidth(2);
  M_pi0[0]->Draw("hist");
  // E-cut condition for M_pi0
  M_pi0[1]->SetLineWidth(2);
  M_pi0[1]->SetLineColor(kRed);
  M_pi0[1]->Draw("hist same");
  c6->Write();

  c7->Divide(2,1);
  c7->cd(1);
  M_pi0[2]->GetXaxis()->SetTitle("M_{#gamma#gamma}");
  M_pi0[2]->GetYaxis()->SetTitle("counts/bin");
  M_pi0[2]->SetLineWidth(2);
  M_pi0[2]->SetLineWidth(2);
  M_pi0[2]->Draw("hist");
  c7->cd(2);
  E_pi0[2]->GetXaxis()->SetTitle("E_{#gamma#gamma}");
  E_pi0[2]->GetYaxis()->SetTitle("counts/bin");
  E_pi0[2]->SetLineWidth(2);
  E_pi0[2]->SetLineWidth(2);
  E_pi0[2]->Draw("hist");
  c7->Write();

  c8->Divide(2,2);
  c8->cd(1);
  pi0ang[0]->GetXaxis()->SetTitle("cos(#theta)");
  pi0ang[0]->GetYaxis()->SetTitle("counts/bin");
  pi0ang[0]->SetLineWidth(2);
  pi0ang[0]->SetLineWidth(2);
  pi0ang[0]->Draw("hist");
  c8->cd(2);
  M_pi0[0]->GetXaxis()->SetTitle("M_{#gamma#gamma}");
  M_pi0[0]->GetYaxis()->SetTitle("counts/bin");
  M_pi0[0]->SetLineWidth(2);
  M_pi0[0]->SetLineWidth(2);
  M_pi0[0]->Draw("hist");
  auto ymax1=M_pi0[0]->GetMaximum();
  auto al1=new TLine(.095,0.0,.095,1.02*ymax1);
  auto al2=new TLine(.145,0.0,.145,1.02*ymax1);
  al1->SetLineColor(kBlue);
  al2->SetLineColor(kBlue);
  al1->SetLineWidth(3);
  al2->SetLineWidth(3);
  al1->Draw();
  al2->Draw();
  c8->cd(3);
  E_pi0[0]->GetXaxis()->SetTitle("E_{#gamma#gamma}");
  E_pi0[0]->GetYaxis()->SetTitle("counts/bin");
  E_pi0[0]->SetLineWidth(2);
  E_pi0[0]->SetLineWidth(2);
  E_pi0[0]->Draw("hist");
  c8->cd(4);
  angPP[0]->GetXaxis()->SetTitle("cos(#theta_{#pi^{+}#pi^{0}})");
  angPP[0]->GetYaxis()->SetTitle("counts/bin");
  angPP[0]->SetLineWidth(2);
  angPP[0]->SetLineWidth(2);
  angPP[0]->Draw("hist");
  c8->Write();

  c9->Divide(2,2);
  c9->cd(1);
  gStyle->SetOptStat(0);
  pi0ang[1]->GetXaxis()->SetTitle("cos(#theta)");
  pi0ang[1]->GetYaxis()->SetTitle("counts/bin");
  pi0ang[1]->SetLineWidth(2);
  pi0ang[1]->SetLineWidth(2);
  pi0ang[1]->Draw("hist");
  c9->cd(2);
  M_pi0[1]->GetXaxis()->SetTitle("M_{#gamma#gamma}");
  M_pi0[1]->GetYaxis()->SetTitle("counts/bin");
  M_pi0[1]->SetLineWidth(2);
  M_pi0[1]->SetLineWidth(2);
  M_pi0[1]->Draw("hist");
  auto ymax=M_pi0[1]->GetMaximum();
  auto ml1=new TLine(.095,0.0,.095,1.02*ymax);
  auto ml2=new TLine(.145,0.0,.145,1.02*ymax);
  ml1->SetLineColor(kBlue);
  ml2->SetLineColor(kBlue);
  ml1->SetLineWidth(3);
  ml2->SetLineWidth(3);
  ml1->Draw();
  ml2->Draw();
  c9->cd(3);
  E_pi0[1]->GetXaxis()->SetTitle("E_{#gamma#gamma}");
  E_pi0[1]->GetYaxis()->SetTitle("counts/bin");
  E_pi0[1]->SetLineWidth(2);
  E_pi0[1]->SetLineWidth(2);
  E_pi0[1]->Draw("hist");
  c9->cd(4);
  angPP[1]->GetXaxis()->SetTitle("cos(#theta_{#pi^{+}#pi^{0}})");
  angPP[1]->GetYaxis()->SetTitle("counts/bin");
  angPP[1]->SetLineWidth(2);
  angPP[1]->SetLineWidth(2);
  angPP[1]->Draw("hist");
  c9->Write();

  return 0; // 0 = all ok
};
