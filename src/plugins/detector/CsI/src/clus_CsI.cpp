#include <Det_CsI.h>
//#include <mn2CsIfn.h>
#include<iostream>
#include<cmath>

Long_t Det_CsI::histos_clus(){
  for(int iClock=0;iClock<12;iClock++){
    for(int iFB=0;iFB<2;iFB++){
      for(int iUD=0;iUD<2;iUD++){
        for(int iModule=0;iModule<16;iModule++){
          std::ostringstream name, name2, name3, name4, name5, name6, tname;
          name<<"stat_"; name2<<"Mnfit_"; name3<<"F1fit_"; name4<<"Dhfit_"; name5<<"pHeight", name6<<"Diff";
          tname<<"time";
          name<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name2<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name3<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name4<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name5<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          name6<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          tname<<iClock<<"_"<<iFB<<"_"<<iUD<<"_"<<iModule;
          h1time[iClock][iFB][iUD][iModule]=dH1(tname.str().c_str(),"stat",250,0,250);
          h1Fits[iClock][iFB][iUD][iModule]=dH1(name.str().c_str(),"stat",250,0,250);
          h1Amps[iClock][iFB][iUD][iModule]=dH1(name5.str().c_str(),"stat",250,0,250);
          h1Mnft[iClock][iFB][iUD][iModule]=dH1(name2.str().c_str(),"stat",250,0,250);
          h1Diff[iClock][iFB][iUD][iModule]=dH1(name6.str().c_str(),"stat",250,0,250);
        }
      }
    }
  }
  // parameter specific histos
  s = new TSpectrum(4);
  p0=dH1("p0","Parameter p0", 250, 0, 1100);
  p1=dH1("p1","Parameter p1", 250, 0, 250);
  p9=dH1("p9","Parameter p9", 250, 0, 300);
  p10=dH1("p10","Parameter p10", 200, 0, 800);
  h2clus=dH2("hclust","Clusters in the CsI(Tl)  ", 50,-25,25,27,0,54);
  h1Pamp=dH1("hpulse","Pulse height distribution", 250, 0, 1000);
  h1kmu2=dH1("kmu2DP","Pulse height distribution", 250, 0, 1000);
  h1Intg=dH1("Integr","Integrated pulse height distribution", 250, 0, 100000);
  h1cali=dH1("Calibr","Integrated pulse height distribution", 250, 0, 1000);
  h1ped=dH1("Ped","Pedestals for the waveform ", 250, 0, 1000);
  for(int i=0;i<11;i++){
    hbox1[i]=new TLine(-2.0,4*(i+1),4.00,4*(i+1));
    vbox1[i]=new TLine(-2.0,4*(i+1),-2.0,4*(i+1)+2);
    hbox1[i+11]=new TLine(-2.0,4*(i+1)+2,4.00,4*(i+1)+2);
    vbox1[i+11]=new TLine(4.0,4*(i+1),4.00,4*(i+1)+2);
    hline1[i]=new TLine(-9.0,4*(i+1)-1,11.0,4*(i+1)-1);
    hline1[i+11]=new TLine(-9.0,4*(i+1)-1,11.0,4*(i+1)-1);
    hline2[i]=new TLine(-9.0,4*(i+1)+1,-2.0,4*(i+1)+1);
    hline2[i+11]=new TLine(-9.0,4*(i+1)+1,-2.0,4*(i+1)+1);
    hline3[i]=new TLine(4.00,4*(i+1)+1,11.0,4*(i+1)+1);
    hline3[i+11]=new TLine(4.00,4*(i+1)+1,11.0,4*(i+1)+1);
  }
  hline2[22]=new TLine(-9.0,1,-2.0,1);
  hline2[23]=new TLine(4.00,1,11.0,1);
  hline2[24]=new TLine(-9.0,49,-2.0,49);
  hline2[25]=new TLine(4.00,49,11.0,49);
  hbox2[0]=new TLine(-2.0,2.0,4,2.0);
  hbox2[1]=new TLine(-2.0,48.,4,48.);
  vbox2[0]=new TLine(-2.0,1.0,-2,2.0);
  vbox2[1]=new TLine(4.00,1.0,4.,2.0);
  vbox2[2]=new TLine(-2.0,48.,-2,49.);
  vbox2[3]=new TLine(4.00,48.,4.,49.);
  hline1[22]=new TLine(-9.0,47,11.0,47);
  for(int n=0;n<22;n++){
    hbox1[n]->Draw("l");
    vbox1[n]->Draw("l");
    hline1[n]->Draw("l");
    hline2[n]->Draw("l");
    hline3[n]->Draw("l");
  }
  for(int n=0;n<2;n++){
    hbox2[n]->Draw("l");
  }
  for(int n=0;n<4;n++){
    vbox2[n]->Draw("l");
    hline2[n+22]->Draw("l");
  }
  hline1[22]->Draw("l");
  return 0;
}

Long_t Det_CsI::startup_clus(){
  getBranchObject("vf48",(TObject **) &treeRaw);
  getBranchObject("RawBeamInfo",(TObject **) &treeBeam);
  gStyle->SetOptStat(0);

  return 0;
}


Long_t Det_CsI::finalize_clus(){
  return 0;
}
