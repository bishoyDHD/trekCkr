#include <stt_analysis.h>

#include<iostream>
#include<cstdlib>
#include<cmath>


stt_analysis::stt_analysis(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{ 

};

stt_analysis::~stt_analysis()
{
};


Long_t stt_analysis::defineHistograms()
{
 for (int i=0; i<stt_num; i++) straw[stt_num] = new TH1D("htrb_sbl",Form("%i", i),1000,-100,100);
 return ok;
}



Long_t stt_analysis::startup()
{
  rawstt = new StrawTube();
  getBranchObject("StrawTube",(TObject **) & rawstt);
  if (!rawstt)
    { 
      debug(0,"Could not find the STRAW  branch in the raw tree. Bailing out\n");
      exit(-1);
    }

  return ok;
}


Long_t stt_analysis::process()
{
  

  return ok;
}

Long_t stt_analysis::finalize()
{
  return ok;
}


Long_t stt_analysis::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};

extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new stt_analysis(in,out,inf_,outf_,p);
}
}

ClassImp(stt_analysis);
