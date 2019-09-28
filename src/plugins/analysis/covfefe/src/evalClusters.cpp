#include "evalClusters.h"

evalClusters::evalClusters():channelNo(7){

}
evalClusters::~evalClusters(){

}
void evalClusters::setChannel(int val){
  channelNo=val;
}
void evalClusters::defHistos(){
  switch(channelNo){
    case 7:
      title1="Opening angle of 2#gamma";
      title2="Invariant Mass of #pi^{0}";
      title3="Energy of 2#gamma";
      title4="Opening angle of  #pi^{+}#pi^{0}";
      prname1="#pi^{+}";
      prname2="#pi^{0}";
      clname1="#gamma_{1}";
      clname2="#gamma_{2}";
      break;
    case 71: // should really be 7.1 for pi0->e+e-gamma. switch only takes int
      title1="Opening angle of 2#gamma";
      title2="Invariant Mass of #pi^{0}";
      title3="Energy of 2#gamma";
      title4="Opening angle of  #pi^{+}#pi^{0}";
      prname1="#pi^{+}";
      prname2="#pi^{0}";
      clname1="e^{+}";
      clname2="e^{-}";
      clname3="#gamma";
      break;
    case 14: // for A' (channel 14)
      title1="Opening angle of e^{+}e^{-}";
      title2="Invariant Mass of A'";
      title3="Energy of e^{+}e^{-}";
      title4="Opening angle of #mu^{+}A'";
      prname1="#mu^{+}";
      prname2="A'";
      clname1="e^{+}";
      clname2="e^{-}";
      break;
    case 15: // for A' (channel 14)
      title1="Opening angle of e^{+}e^{-}";
      title2="Invariant Mass of A'";
      title3="Energy of e^{+}e^{-}";
      title4="Opening angle of #pi^{+}A'";
      prname1="#pi^{+}";
      prname2="A'";
      clname1="e^{+}";
      clname2="e^{-}";
      break;
    case 16: // for SM (channel 16)
      title1="Opening angle of e^{+}e^{-}";
      //title2="Invariant Mass of ";
      title3="Energy of e^{+}e^{-}";
      //title4="Opening angle of ";
      prname1="#mu^{+}";
      clname1="e^{+}";
      clname2="e^{-}";
      break;
    default:
      std::cout<<" *** Warning!\n --- Channel definition is out of range!!\n";
      std::cout<<" *** Warning!\n";
      break;
  }// end of switch statement
  // Define various histograms
  clustAng=new TH1D("clAng",title1.c_str(),55.0,-1.1,1.1);
  invM=new TH1D("invM",title3.c_str(),37.5,0.0,.150);
  Eclust=new TH1D("Eclust",title3.c_str(),75.0,0.0,.300);
  primAng=new TH1D("primAng",title4.c_str(),55.0,-1.1,1.1);
  // needed regardless of channel no.
  thetaPhi=new TH2D("h2ang", "#theta Vs. #phi",24.0,0.,M_PI,48.0,0.,2.*M_PI);
}
// function to plot actual histos
void evalClusters::drawHistos(){
  switch(channelNo){
    case 7:
      drawCanvas(clustAng,invM,Eclust,primAng,1);
      drawCanvas(thetaPhi,1);
      break;
    case 14:
      std::cout<<" will get to this soon!! \n";
      break;
  }// end of switch statement
}
void evalClusters::fillHistos(double theta1,double phi1,double theta2,double phi2,double theta3,double phi3){
  switch(channelNo){
    case 71:
      thetaPhi->Fill(theta1,phi1);
      thetaPhi->Fill(theta2,phi2);
      thetaPhi->Fill(theta3,phi3);
      break;
    case 7:
      thetaPhi->Fill(theta1,phi1);
      thetaPhi->Fill(theta2,phi2);
      break;
  }
}
void evalClusters::fillHistos(double InvMass,double opAng1,double Etot,double primOpAng){
  invM->Fill(InvMass);
  Eclust->Fill(Etot);
  clustAng->Fill(opAng1);
  primAng->Fill(primOpAng);
}
void evalClusters::drawCanvas(TH2D* h2,int val){
  std::string cname="c_";
  cname+=val;
  TCanvas* c1=new TCanvas(cname.c_str(),"Theta vs Phi",900,800);
  c1->cd();
  std::string xtitle="#theta [rad]";
  std::string ytitle="#phi [rad]";
  h2->GetXaxis()->SetTitle(xtitle.c_str());
  h2->GetYaxis()->SetTitle(ytitle.c_str());
  h2->Draw("colz");
}
void evalClusters::drawCanvas(TH1D* hist1,TH1D* hist2,TH1D* hist3,TH1D* hist4,int val=1){
  std::string cname="c";
  cname+=val;
  TCanvas* c1=new TCanvas(cname.c_str(),"Sanity plots",900,800);
  c1->Divide(2,2);
  std::string xtitle="cos(#theta_{";
  xtitle+=clname1;
  xtitle+=clname2;
  xtitle+="})";
  c1->cd(1);
  hist1->GetXaxis()->SetTitle(xtitle.c_str());
  hist1->GetYaxis()->SetTitle("conts/bin");
  hist1->Draw("hist");
  c1->cd(2);
  xtitle="M_{";
  xtitle+=clname1;
  xtitle+=clname2;
  xtitle+="}";
  hist2->GetXaxis()->SetTitle(xtitle.c_str());
  hist2->GetYaxis()->SetTitle("conts/bin");
  hist2->Draw("hist");
  c1->cd(3);
  xtitle="E_{";
  xtitle+=clname1;
  xtitle+=clname2;
  xtitle+="}";
  //hist3->GetXaxis()->SetLabelFont(0.35);
  hist3->GetXaxis()->SetTitle(xtitle.c_str());
  hist3->GetYaxis()->SetTitle("conts/bin");
  hist3->Draw("hist");
  c1->cd(4);
  xtitle="cos(#theta_{";
  xtitle+=prname1;
  xtitle+=prname2;
  xtitle+="})";
  hist4->GetXaxis()->SetTitle(xtitle.c_str());
  hist4->GetYaxis()->SetTitle("conts/bin");
  hist4->Draw("hist");
  c1->Write();
}
