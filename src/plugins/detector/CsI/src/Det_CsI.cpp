#include <Det_CsI.h>

#include<cmath>

#include <fstream>
#include <iostream>
using namespace std;

Det_CsI::Det_CsI(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_):Plugin(in_,out_,inf_,outf_,p_)
{
  // Set defaults for various options
  treeCali=0;
  setIdCsI(mapCsI);
};


Det_CsI::~Det_CsI()
{
};
Long_t Det_CsI::setIdCsI(map<IdCsI,UInt_t> & map){
  char fb[2]={'f','b'};
  char ud[2]={'u','d'};
  UInt_t index=1;
  for(UInt_t iHole=0;iHole<12;iHole++){
    for(int ifb=0;ifb<2;ifb++){
      for(int iud=0;iud<2;iud++){
	for(UInt_t iCsI=0;iCsI<16;iCsI++){
	  UInt_t name;
	  name=0x30000000+(((iHole+1)/10)<<24);
	  name+=0x00300000+(((iHole+1)%10)<<16);
	  name+=fb[ifb]<<8;
	  name+=ud[iud];
	  IdCsI first(name,iCsI+1);
	  map[first]=index;
	  index++;
	}
      }
    }
  }
}

Long_t Det_CsI::histos()
{

  return 0;

}

Long_t Det_CsI::startup()
{
  /*
   debug(1,"TOF Startup called.\n ToF Config:\nTDC offsets:");
  debug(1,"\nADC offsets:");
  for (int i=0;i<72;i++)
    debug(1,"%i ",adcoffset[i]);
  debug(1,"\nADC gains:");
  for (int i=0;i<72;i++)
    debug(1,"%g ",adcgain[i]);
  for (int i=0;i<36;i++)
    debug(1,"\n Speed Of Light:  %g  mm/count\n",channeldist[i]);


  // Generate our ToF object
  tofo=new ToF();
  
  // Make new branch on output tree
  makeBranch("ToFhits",(TObject **)&tofo);

  // Look up input tof:
  getBranchObject("ToF",(TObject **) &tofi);
  
  getBranchObject("EventInfo", (TObject**)&eventinfo);

  // Check for pedestals from pedestal running
  ped = (ORTPedestal*) getFileObject("TOF_PEDS");
  double avg=0;
  if (ped != NULL) // Pedestals exist in raw file
  {
    debug(1,"############# Found pedestal data in raw file\n");
    for (int a = 0; a < 72; a++) 
    {
      // Indexing used by saved pedestals is different from indexing used everywhere else
      // ped 0-71: L_top, R_top, L_bottom, R_bottom 
      int pindex = (a%18) + ((a %36)>17?36:0)+ (a>35?18:0);
      double peda= ped->pedestals[pindex].mean;

      avg+=peda;
      set_adcoffset(a, peda);
      char buf[1000];
      sprintf(buf,"ToF:ADC:ped:%i",a);
      saveInfo(buf,peda);
    }
    saveInfo("ToF:ADC:ped:all",avg/72);
  }
  debug(1,"ADC pedestals:");
  for (int b = 0; b < 72; b++) 
    debug(1,"%i ",adcoffset[b]);
  debug(1,"\n");

  switch (getRunType())
  {
  case Plugin::PositronPositive:
    debug(1,"ToF: Beam Species:  Positrons\nToF: Toroid Current: positive\n");
    for(int i = 0; i< 72; i++)
      TDCoffset[i] = pTDCoffset[i];
    std::cout << "positive " << std::endl;
    break;
  case Plugin::PositronNegative:
    debug(1,"ToF: Beam Species:  Positrons\nToF: Toroid Current: negative\n");
    for(int i = 0; i< 72; i++)
      TDCoffset[i] = pTDCoffset[i];
    break;
  case Plugin::ElectronPositive:
    debug(1,"ToF: Beam Species:  Electrons\nToF: Toroid Current: positive\n");
    for(int i = 0; i< 72; i++)
      TDCoffset[i] = eTDCoffset[i];
    break;
  case Plugin::ElectronNegative:
    debug(1,"ToF: Beam Species:  Electrons\nToF: Toroid Current: negative\n");
    for(int i = 0; i< 72; i++)
      TDCoffset[i] = eTDCoffset[i];
    break;
  case Plugin::unknown:
    debug(0," *** Slow control does not know run type! Cosmic run? Using Positron/positive.\n");
    for(int i = 0; i< 72; i++)
      TDCoffset[i] = pTDCoffset[i];
    break;
  }


  //setup hv tests
  eventshvoff=0;
  brokenchannels=0;
  scmanager= (slowctrl::manager*)getMemoryObject("SlowCtrl Manager");
  for (int i=1;i<19;i++)
  {
    char channame[1000];
    sprintf(channame,"TOF:LB%02iSTI",i);
    scmanager->addWatch(scmanager->getCurrentByName(channame),0.5,1,&brokenchannels);
    sprintf(channame,"TOF:LT%02iSTI",i);
    scmanager->addWatch(scmanager->getCurrentByName(channame),0.5,1,&brokenchannels);
    sprintf(channame,"TOF:RB%02iSTI",i);
    scmanager->addWatch(scmanager->getCurrentByName(channame),0.5,1,&brokenchannels);
    sprintf(channame,"TOF:RT%02iSTI",i);
    scmanager->addWatch(scmanager->getCurrentByName(channame),0.5,1,&brokenchannels);
  }

  loadGeo(); // Load geometry for phi reconstruction
  */
  return 0;
};


Long_t Det_CsI::process()
{
  /*
  if (brokenchannels>0)
  {
    eventshvoff++;
    debug(1,"TOF HV dropped out\n");
  }

  int tofmulti=0;
  // GetEntry already called by cooker
  tofo->hits.clear();
  for (int i = 0; i < 2; i++)
    tofo->hitpattern[i] = 0;

  // Scan all ToF bars
  tofo->firstHitRight=0;
  for (int b=0;b<36;b++)
  {  
    alltdc_bypmt->Fill(b,tofi->tofs[b].tdc[0]);
    alltdc_bypmt->Fill(b+36,tofi->tofs[b].tdc[1]);
    
    alladc_bypmt->Fill(b,tofi->tofs[b].adc[0]);
    alladc_bypmt->Fill(b+36,tofi->tofs[b].adc[1]);

    // Top PMTs
    alltdc_bypmt->Fill(b,tofi->tofs[b].tdc[0]);
    if (tdc_valid(tofi->tofs[b].tdc[0]))
    {
      validadc_bypmt->Fill(b,tofi->tofs[b].adc[0]);
      validtdc_bypmt->Fill(b,tofi->tofs[b].tdc[0]);
    }
    else
      invalidadc_bypmt->Fill(b, (double)tofi->tofs[b].adc[0]-adcoffset[b]);

    // Bottom PMTs
    alltdc_bypmt->Fill(b+36,tofi->tofs[b].tdc[1]);
    if (tdc_valid(tofi->tofs[b].tdc[1]))
    {
      validadc_bypmt->Fill(b+36,tofi->tofs[b].adc[1]);
      validtdc_bypmt->Fill(b+36,tofi->tofs[b].tdc[1]);
    }
    else
      invalidadc_bypmt->Fill(b+36, (double)tofi->tofs[b].adc[1]-adcoffset[b+36]);

    // We have a coincidence hit
    if (tdc_valid(tofi->tofs[b].tdc[0]) && tdc_valid(tofi->tofs[b].tdc[1]))
    {
      tup_vs_tbottom[b] ->Fill(tofi->tofs[b].tdc[0]+TDCoffset[b]-tofi->reftime,tofi->tofs[b].tdc[1]+TDCoffset[b+36]-tofi->reftime);

      // Fill the hit info
      ToFhit hit;
      hit.bar = b;
      hit.meantime = get_tdc_meantime(b);
      hit.hitpos = get_tdc_position(b);
      hit.qsum = get_energy_deposition(b);
		
		double loc[3] = {0.0, hit.hitpos/10, 0.0};
		double glob[3];
		tofvolumes[hit.bar]->LocalToMaster(loc,glob);
		TVector3 ptg(10*glob[0],10*glob[1],10*glob[2]);
		hit.phi = ptg.Phi();
		if (hit.phi<0) hit.phi = 2*M_PI+hit.phi;

//      hit.altphi = atan((hit.hitpos+V[b].Y())/fabs(V[b].X()));
//      if((b < 18) && (hit.altphi < 0))
//			hit.altphi = 2*M_PI+hit.altphi;
//      if(b > 18)
//			hit.altphi = M_PI - hit.altphi;
      
		tofo->hits.push_back(hit);

      if (b<18) tofo->firstHitRight++;
      tofo->hitpattern[b < 18 ? 0 : 1] |= 0x1 << (b%18);

      // Fill the histograms
  
      hits_bybar->Fill(hit.bar);
      mtime_bybar->Fill(hit.bar,hit.meantime);
      edep_bybar->Fill(hit.bar,hit.qsum);
      e_vs_t_all->Fill(hit.meantime,hit.qsum);
      if (hit.bar %18 >14) e_vs_t_back->Fill(hit.meantime,hit.qsum);
      e_vs_t_bybar[b]->Fill(hit.meantime,hit.qsum);
      RefTime->Fill(tofi->reftime);
      ypos_bybar->Fill(b,hit.hitpos);
     }
  }
  H1(tofmulti, "ToF/TOFmulti", "TOFmultiplicity", 37, -0.5, 36.5);

//	int nhitl = 0;
//	int nhitr = 0;
//   double bnr = -999;
//	double bnl = -999;
//	// Go through the registered hits
//	for (int i=0;i<(int)tofo->hits.size();i++)
//	{
//		// Get which bar was hit
//		int bn = tofo->hits[i].bar;
//		if (bn<18) {nhitl++; bnl=bn;}
//		else {nhitr++; bnr=bn;}
//	}
//	if ((eventinfo->trigFired & 0x2) && nhitr==1) lumitofl->Fill(bnr+0.5);
//	if ((eventinfo->trigFired & 0x4) && nhitl==1) lumitofr->Fill(bnl+0.5);
*/
  return 0;
};


Long_t Det_CsI::done()
{
  /*
  saveInfo("ToF:HVoff",eventshvoff*1.0/in->GetEntries());

  setDetectorFlag("ToF",eventshvoff>0?2:1);

  // Save ADC means for each PMT
  TProfile *adc_mean = validadc_bypmt->ProfileX();
  for (int i=0;i<72;i++)
  {
    char buf[1000];
    sprintf(buf,"ToF:ADC:mean:%i",i);
    saveInfo(buf,adc_mean->GetBinContent(i+1));
  }
  saveInfo("ToF:ADC:mean:all", validadc_bypmt->GetMean());

  // Save energy deposition means for each ToF
  TProfile *edep_mean = edep_bybar->ProfileX();
  for (int i=0;i<36;i++)
  {
    char buf[1000];
    sprintf(buf,"ToF:EnergyDeposited:%i",i);
    saveInfo(buf,edep_mean->GetBinContent(i+1));
  }
  // Done with ToFs
  delete tofo;
  */
  return 0;
};

Long_t Det_CsI::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};

/////////////////////////////////////////////////////////////////////////////////////////////////// 

extern "C"{
  Plugin *factory(TTree *in_, TTree *out_, TFile *inf_, TFile * outf_,TObject *p_)
  {
      return (Plugin *) new Det_CsI(in_,out_,inf_,outf_,p_);
  }
}


ClassImp(Det_CsI);
