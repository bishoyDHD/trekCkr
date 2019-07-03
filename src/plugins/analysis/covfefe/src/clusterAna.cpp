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
  E_pi0[0]=dH1("Epi0", " E_{total}(2#gamma)",40.5,0.,0.70);
  M_pi0[0]=dH1("Mpi0", " M_{#pi^{0}}",40.5,0.,0.4);
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
  E_pi0[1]=dH1("Epi0_1", " E_{total}(2#gamma)",40.5,0.,0.70);
  M_pi0[1]=dH1("Mpi0_1", " M_{#pi^{0}}",40.5,0.,0.4);
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
  angPP[0]=dH1("angPP_1","Angle between #pi^{+} and #pi^{0}",26.0,-1.5,1.5);
  h2corrAng=dH2("h2cAng", "cos(#theta_{#pi^{+}}) vs. cos(#theta_{#gamma#gamma})",25.,-1.,1.,25,-1.,1.);
  piPang[0]=dH1("piPAng1", "cos(#theta_{#pi^{+}}", 45.,-1.5,1.5);
  pi0ang[0]=dH1("pi0Ang0", "cos(#theta_{#gamma#gamma}", 45.,-1.1,1.1);
  pi0ang[1]=dH1("pi0Ang1", "cos(#theta_{#gamma#gamma}", 45.,-1.1,1.1);

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
  TLorentzVector piPlv, pi0lv;
  // 3-vector to calculate the angle between pi+ and pi0
  TVector3 piPv3, pi0v3, ggv3, g1v3, g2v3;
  if(clsmar->waveID>0){
    waveID[0]->Fill(clsmar->waveID);
    id1->Fill(clsmar->dubP_1);
    E_pi0[0]->Fill(clsmar->E_pi0);
    M_pi0[0]->Fill(clsmar->M_pi0);
    clustM[0]->Fill(clsmar->clusterM);
    g1px[0]->Fill(clsmar->g1Px);
    g1py[0]->Fill(clsmar->g1Py);
    g1pz[0]->Fill(clsmar->g1Pz);
    // gamma2
    g2px[0]->Fill(clsmar->g2Px);
    g2py[0]->Fill(clsmar->g2Py);
    g2pz[0]->Fill(clsmar->g2Pz);
    // pi0 vertex
    pi0px[0]->Fill(clsmar->pi0px);
    pi0py[0]->Fill(clsmar->pi0py);
    pi0pz[0]->Fill(clsmar->pi0pz);
    // pi+ vertex
    vertpx[0]->Fill(clsmar->piPpx);
    vertpy[0]->Fill(clsmar->piPpy);
    vertpz[0]->Fill(clsmar->piPpz);
    kmass[0]->Fill(clsmar->piP2g,clsmar->piPpi0);
    h2Angle[0]->Fill(clsmar->thetaE,clsmar->phiE);
    // angles
    h1theta[0]->Fill(clsmar->thetaE);
    h1phi[0]->Fill(clsmar->phiE);
    // apply cut in total energy of pi0
    // replot with new cut condition
    if(clsmar->E_pi0>=.1 && clsmar->E_pi0<=.3){
      E_pi0[1]->Fill(clsmar->E_pi0);
      M_pi0[1]->Fill(clsmar->M_pi0);
      clustM[1]->Fill(clsmar->clusterM);
      g1px[1]->Fill(clsmar->g1Px);
      g1py[1]->Fill(clsmar->g1Py);
      g1pz[1]->Fill(clsmar->g1Pz);
      // gamma2
      g2px[1]->Fill(clsmar->g2Px);
      g2py[1]->Fill(clsmar->g2Py);
      g2pz[1]->Fill(clsmar->g2Pz);
      // pi0 vertex
      pi0px[1]->Fill(clsmar->pi0vx);
      pi0py[1]->Fill(clsmar->pi0vy);
      pi0pz[1]->Fill(clsmar->pi0vz);
      // pi+ vertex
      vertpx[1]->Fill(clsmar->piPpx);
      vertpy[1]->Fill(clsmar->piPpy);
      vertpz[1]->Fill(clsmar->piPpz);
      kmass[1]->Fill(clsmar->piP2g,clsmar->piPpi0);
      h2Angle[1]->Fill(clsmar->thetaE,clsmar->phiE);
      // angles
      piPv3.SetXYZ(clsmar->piPpx,clsmar->piPpy,clsmar->piPpz);
      pi0v3.SetXYZ(clsmar->pi0px,clsmar->pi0py,clsmar->pi0pz);
      ggv3.SetXYZ(clsmar->pi0vx,clsmar->pi0vy,clsmar->pi0vz);
      g1v3.SetXYZ(clsmar->g1Px,clsmar->g1Py,clsmar->g1Pz);
      g2v3.SetXYZ(clsmar->g2Px,clsmar->g2Py,clsmar->g2Pz);
      // rotate angles
      pi0v3.RotateZ(M_PI/2);
      double ang=piPv3.Angle(ggv3);
      double ggAng=g1v3.Angle(g2v3);
      h2corrAng->Fill(pi0v3.CosTheta(),std::cos(ang));
      angPP[0]->Fill(std::cos(ang));
      h1theta[1]->Fill(clsmar->thetaE);
      h1phi[1]->Fill(clsmar->phiE);
      piPang[0]->Fill(piPv3.CosTheta());
      pi0ang[0]->Fill(pi0v3.CosTheta());
      pi0ang[1]->Fill(ggAng);
    }
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
  TCanvas* c9=new TCanvas("c9"," Angle cosTheta ",808,700);
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
  clustM[0]->SetLineColor(kBlue);
  clustM[0]->GetXaxis()->SetTitle("Cluster multiplicity");
  clustM[0]->GetYaxis()->SetTitle("counts/bin");
  clustM[0]->SetLineWidth(2);
  clustM[0]->Draw("hist");
  clustM[1]->SetLineWidth(2);
  clustM[1]->SetLineColor(kRed);
  clustM[1]->Draw("hist same");
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
  vertpx[0]->GetXaxis()->SetTitle("px_{#pi^{+}} [GeV/c]");
  vertpx[0]->GetYaxis()->SetTitle("counts/bin");
  vertpx[0]->SetLineWidth(2);
  vertpx[0]->Draw("hist");
  c5->cd(11);
  gStyle->SetOptStat(0);
  vertpy[0]->GetXaxis()->SetTitle("py_{#pi^{+}} [GeV/c]");
  vertpy[0]->GetYaxis()->SetTitle("counts/bin");
  vertpy[0]->SetLineWidth(2);
  vertpy[0]->Draw("hist");
  c5->cd(12);
  gStyle->SetOptStat(0);
  vertpz[0]->GetXaxis()->SetTitle("pz_{#pi^{+}} [GeV/c]");
  vertpz[0]->GetYaxis()->SetTitle("counts/bin");
  vertpz[0]->SetLineWidth(2);
  vertpz[0]->Draw("hist");
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

  c7->cd();
  gStyle->SetOptStat(0);
  angPP[0]->GetXaxis()->SetTitle("angle");
  angPP[0]->GetYaxis()->SetTitle("counts/bin");
  angPP[0]->SetLineWidth(2);
  angPP[0]->Draw("hist");
  c7->Write();

  c8->cd();
  gStyle->SetOptStat(0);
  h2corrAng->GetXaxis()->SetTitle("cos(#theta_{#gamma#gamma})");
  h2corrAng->GetYaxis()->SetTitle("cos(#theta_{#pi^{+}})");
  h2corrAng->SetLineWidth(2);
  h2corrAng->Draw("colz");
  c8->Write();

  c9->cd();
  gStyle->SetOptStat(0);
  pi0ang[1]->GetXaxis()->SetTitle("cos(#theta)");
  pi0ang[1]->GetYaxis()->SetTitle("counts/bin");
  pi0ang[1]->SetLineWidth(2);
  //piPang[0]->Draw("hist");
  //pi0ang[0]->SetLineWidth(2);
  //pi0ang[0]->SetLineColor(kRed);
  //pi0ang[0]->Draw("hist same");
  pi0ang[1]->SetLineWidth(2);
  //pi0ang[1]->SetLineColor(kCyan);
  pi0ang[1]->Draw("hist");
  c9->Write();

  return 0; // 0 = all ok
};
