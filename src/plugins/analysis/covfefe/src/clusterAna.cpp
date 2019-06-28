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
  E_pi0=dH1("Epi0", " E_{total}(2#gamma)",110.5,0.,0.9);
  M_pi0=dH1("Mpi0", " M_{#pi^{0}}",110.5,0.,0.9);
  waveID=dH1("waveID", "Waveform ID", 8,0.,8.);
  g1px=dH1("g1px", "px_{#gamma1}",63.5,-.3,.3);
  g1py=dH1("g1py", "py_{#gamma1}",63.5,-.3,.3);
  g1pz=dH1("g1pz", "pz_{#gamma1}",63.5,-.3,.3);
  g2px=dH1("g2px", "px_{#gamma2}",63.5,-.3,.3);
  g2py=dH1("g2py", "py_{#gamma2}",63.5,-.3,.3);
  g2pz=dH1("g2pz", "pz_{#gamma2}",63.5,-.3,.3);
  pi0px=dH1("pi0px", "px_{#pi^{0}}",63.5,-.3,.3);
  pi0py=dH1("pi0py", "py_{#pi^{0}}",63.5,-.3,.3);
  pi0pz=dH1("pi0pz", "pz_{#pi^{0}}",63.5,-.3,.3);
  vertpx=dH1("vertpx", "px_{#pi^{+}}",63.5,-.3,.3);
  vertpy=dH1("vertpy", "py_{#pi^{+}}",63.5,-.3,.3);
  vertpz=dH1("vertpz", "pz_{#pi^{+}}",63.5,-.3,.3);
  h1theta=dH1("theta", "#theta distribution",24.,0,M_PI);
  h1phi=dH1("phi", "#phi distribution",48.,0,2*M_PI);
  h2Angle=dH2("h2Ang","#phi vs #theta", 24.0,0.,M_PI,48.0,0.,2.*M_PI);
  id1=dH1("id1", "Waveform ID", 8,0.,8.);
  clustM=dH1("clustM", "Cluster multiplicity", 12,0.,12.);
  kmass=dH2("kmass", "E_{tot}(#pi^{+}+#pi^{0}) vs. E_{tot}(2#gamma+#pi^{+})", 62.5,0.,1.,62.5,0.,.7);

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
  if(clsmar->waveID>0){
    waveID->Fill(clsmar->waveID);
    id1->Fill(clsmar->dubP_1);
    E_pi0->Fill(clsmar->E_pi0);
    M_pi0->Fill(clsmar->M_pi0);
    clustM->Fill(clsmar->clusterM);
    g1px->Fill(clsmar->g1Px);
    g1py->Fill(clsmar->g1Py);
    g1pz->Fill(clsmar->g1Pz);
    // gamma2
    g2px->Fill(clsmar->g2Px);
    g2py->Fill(clsmar->g2Py);
    g2pz->Fill(clsmar->g2Pz);
    // pi0 vertex
    pi0px->Fill(clsmar->pi0px);
    pi0py->Fill(clsmar->pi0py);
    pi0pz->Fill(clsmar->pi0pz);
    // pi+ vertex
    vertpx->Fill(clsmar->piPpx);
    vertpy->Fill(clsmar->piPpy);
    vertpz->Fill(clsmar->piPpz);
    kmass->Fill(clsmar->piP2g,clsmar->piPpi0);
    h2Angle->Fill(clsmar->thetaE,clsmar->phiE);
  }
  if(clsmar->thetaE>0){
    //std::cout<<" try to get theta and phi "<<clsmar->thetaE<<", "<<clsmar->phiE<<std::endl;
    // angles
    h1theta->Fill(clsmar->thetaE);
    h1phi->Fill(clsmar->phiE);
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
  //TCanvas* c7=new TCanvas("c7"," Peak time ",808,700);
  c1->Divide(1,3);
  c1->cd(3);
  gStyle->SetOptStat(0);
  h2Angle->SetLineColor(kBlue);
  h2Angle->GetXaxis()->SetTitle("#theta [rad]");
  h2Angle->GetYaxis()->SetTitle("#phi [rad]");
  TAxis* a=h2Angle->GetXaxis();
  a->SetNdivisions(-502);
  a->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
  a->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi");
  TAxis* a1=h2Angle->GetYaxis();
  a1->SetNdivisions(-502);
  a1->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
  a1->ChangeLabel(-1,-1,-1,-1,-1,-1,"2#pi");
  h2Angle->Draw("colz");
  gStyle->SetOptStat(0);
  c1->cd(1);
  h1theta->SetLineColor(kBlue);
  h1theta->GetXaxis()->SetTitle("#theta");
  h1theta->GetYaxis()->SetTitle("counts/bin");
  h1theta->SetLineWidth(2);
  TAxis* ax1=h1theta->GetXaxis();
  ax1->SetNdivisions(-502);
  ax1->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
  ax1->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi");
  h1theta->Draw("hist");
  gStyle->SetOptStat(0);
  c1->cd(2);
  h1phi->SetLineColor(kBlue);
  h1phi->GetXaxis()->SetTitle("#phi");
  h1phi->GetYaxis()->SetTitle("counts/bin");
  h1phi->SetLineWidth(2);
  TAxis* ax2=h1phi->GetXaxis();
  ax2->SetNdivisions(-502);
  ax2->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
  ax2->ChangeLabel(-1,-1,-1,-1,-1,-1,"2#pi");
  h1phi->Draw("hist");
  gStyle->SetOptStat(0);
  c1->Write();

  c2->cd();
  gStyle->SetOptStat(0);
  clustM->SetLineColor(kBlue);
  clustM->GetXaxis()->SetTitle("Cluster multiplicity");
  clustM->GetYaxis()->SetTitle("counts/bin");
  clustM->SetLineWidth(2);
  clustM->Draw("hist");
  c2->Write();

  c3->cd();
  gStyle->SetOptStat(0);
  E_pi0->GetXaxis()->SetTitle("E_{#gamma#gamma} [GeV/c]");
  E_pi0->GetYaxis()->SetTitle("counts/bin");
  E_pi0->SetLineWidth(2);
  E_pi0->Draw("hist");
  c3->Write();

  c4->cd();
  gStyle->SetOptStat(0);
  kmass->GetXaxis()->SetTitle("E_{#gamma#gamma} + E_{#pi^{+}} [GeV/c]");
  kmass->GetYaxis()->SetTitle("E_{#pi^{+}} + E_{#pi^{0}} [GeV/c]");
  kmass->SetLineWidth(2);
  kmass->Draw("colz");
  auto tl3=new TLine(.485,0.0,.485,.695);
  auto tl4=new TLine(.515,0.0,.515,.695);
  tl3->SetLineColor(kRed);
  tl4->SetLineColor(kRed);
  tl3->SetLineWidth(3);
  tl4->SetLineWidth(3);
  tl3->Draw();
  tl4->Draw();
  // along y: pi+ cut
  auto tl5=new TLine(0.00,.485,.955,.485);
  auto tl6=new TLine(0.00,.505,.955,.505);
  tl5->SetLineColor(kRed);
  tl6->SetLineColor(kRed);
  tl5->SetLineWidth(3);
  tl6->SetLineWidth(3);
  tl5->Draw();
  tl6->Draw();
  c4->Write();

  c5->Divide(3,4);
  c5->cd(1);
  gStyle->SetOptStat(0);
  g1px->GetXaxis()->SetTitle("px_{#gamma1} [GeV/c]");
  g1px->GetYaxis()->SetTitle("counts/bin");
  g1px->SetLineWidth(2);
  g1px->Draw("hist");
  c5->cd(2);
  gStyle->SetOptStat(0);
  g1py->GetXaxis()->SetTitle("py_{#gamma1} [GeV/c]");
  g1py->GetYaxis()->SetTitle("counts/bin");
  g1py->SetLineWidth(2);
  g1py->Draw("hist");
  c5->cd(3);
  gStyle->SetOptStat(0);
  g1pz->GetXaxis()->SetTitle("pz_{#gamma1} [GeV/c]");
  g1pz->GetYaxis()->SetTitle("counts/bin");
  g1pz->SetLineWidth(2);
  g1pz->Draw("hist");
  // gamma2
  c5->cd(4);
  gStyle->SetOptStat(0);
  g2px->GetXaxis()->SetTitle("px_{#gamma2} [GeV/c]");
  g2px->GetYaxis()->SetTitle("counts/bin");
  g2px->SetLineWidth(2);
  g2px->Draw("hist");
  c5->cd(5);
  gStyle->SetOptStat(0);
  g2py->GetXaxis()->SetTitle("py_{#gamma2} [GeV/c]");
  g2py->GetYaxis()->SetTitle("counts/bin");
  g2py->SetLineWidth(2);
  g2py->Draw("hist");
  c5->cd(6);
  gStyle->SetOptStat(0);
  g2pz->GetXaxis()->SetTitle("pz_{#gamma2} [GeV/c]");
  g2pz->GetYaxis()->SetTitle("counts/bin");
  g2pz->SetLineWidth(2);
  g2pz->Draw("hist");
  c5->cd(7);
  gStyle->SetOptStat(0);
  pi0px->GetXaxis()->SetTitle("px_{#pi^{0}} [GeV/c]");
  pi0px->GetYaxis()->SetTitle("counts/bin");
  pi0px->SetLineWidth(2);
  pi0px->Draw("hist");
  c5->cd(8);
  gStyle->SetOptStat(0);
  pi0py->GetXaxis()->SetTitle("py_{#pi^{0}} [GeV/c]");
  pi0py->GetYaxis()->SetTitle("counts/bin");
  pi0py->SetLineWidth(2);
  pi0py->Draw("hist");
  c5->cd(9);
  gStyle->SetOptStat(0);
  pi0pz->GetXaxis()->SetTitle("pz_{#pi^{0}} [GeV/c]");
  pi0pz->GetYaxis()->SetTitle("counts/bin");
  pi0pz->SetLineWidth(2);
  pi0pz->Draw("hist");
  c5->cd(10);
  gStyle->SetOptStat(0);
  vertpx->GetXaxis()->SetTitle("px_{#pi^{+}} [GeV/c]");
  vertpx->GetYaxis()->SetTitle("counts/bin");
  vertpx->SetLineWidth(2);
  vertpx->Draw("hist");
  c5->cd(11);
  gStyle->SetOptStat(0);
  vertpy->GetXaxis()->SetTitle("py_{#pi^{+}} [GeV/c]");
  vertpy->GetYaxis()->SetTitle("counts/bin");
  vertpy->SetLineWidth(2);
  vertpy->Draw("hist");
  c5->cd(12);
  gStyle->SetOptStat(0);
  vertpz->GetXaxis()->SetTitle("pz_{#pi^{+}} [GeV/c]");
  vertpz->GetYaxis()->SetTitle("counts/bin");
  vertpz->SetLineWidth(2);
  vertpz->Draw("hist");
  c5->Write();

  c6->cd();
  gStyle->SetOptStat(0);
  M_pi0->GetXaxis()->SetTitle("M_{#pi^{0}} [GeV/c]");
  M_pi0->GetYaxis()->SetTitle("counts/bin");
  M_pi0->SetLineWidth(2);
  M_pi0->Draw("hist");
  c6->Write();

  return 0; // 0 = all ok
};
