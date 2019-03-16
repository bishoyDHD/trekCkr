#include <Det_example.h>

#include<iostream>
#include<cmath>


Det_example::Det_example(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

Det_example::~Det_example()
{
};


Long_t Det_example::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};


extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new Det_example(in,out,inf_,outf_,p);
}
}


ClassImp(Det_example);

