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
  waveID=dH1("waveID", "Waveform ID", 8,0.,8.);
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
    clustM->Fill(clsmar->clusterM);
    kmass->Fill(clsmar->piP2g,clsmar->piPpi0);
  }

  return 0; // 0 = all ok
};

Long_t covfefe::finalize_clust(){
  TCanvas* c1=new TCanvas("c1","Waveform ID",908,980);
  TCanvas* c2=new TCanvas("c2","Cluster multiplicity",908,980);
  TCanvas* c3=new TCanvas("c3","Total energy of 2#gamma",808,700);
  TCanvas* c4=new TCanvas("c4","Total energy of K^{+}",808,700);
  //TCanvas* c5=new TCanvas("c5","Pulse height Distribution",808,700);
  //TCanvas* c6=new TCanvas("c6"," CsI timing",808,700);
  //TCanvas* c7=new TCanvas("c7"," Peak time ",808,700);
  c1->cd();
  gStyle->SetOptStat(0);
  //waveID->SetTitle(" Opening angle between #pi^{+} and #pi^{0} ");
  waveID->SetLineColor(kBlue);
  waveID->GetXaxis()->SetTitle("waveform ID");
  waveID->SetLineWidth(2);
  waveID->Draw("hist");
  id1->SetLineWidth(2);
  id1->SetFillColor(kMagenta);
  id1->SetFillStyle(3244);
  id1->SetLineWidth(2);
  id1->Draw("hist ][ same");
  auto leg1=new TLegend(0.1,0.7,0.48,0.9);
  leg1->SetHeader("Key:","C");
  leg1->AddEntry(waveID, "Waveform classification");
  leg1->AddEntry(id1, "Pre-pileup with good timing");
  leg1->Draw();
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

  return 0; // 0 = all ok
};
