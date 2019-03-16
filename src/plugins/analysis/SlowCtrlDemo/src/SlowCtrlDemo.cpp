#include <SlowCtrlDemo.h>

#include<iostream>
#include<cmath>


SlowCtrlDemo::SlowCtrlDemo(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

SlowCtrlDemo::~SlowCtrlDemo()
{
};




Long_t SlowCtrlDemo::startup()
{

  auto scmanager= (slowctrl::manager*) getMemoryObject("SlowCtrl Manager");

  if (!scmanager)
    {
      debug(0,"Couldn't find SlowCtrl manager. Plugin included?\n");
      exit(-1);
    }

  slctrl_fs11_u=scmanager->getLastValidByName(.cooker:BL:FS11-U:POS:AV");
  return 0;
}

Long_t SlowCtrlDemo::process()
{

  debug(1,"Value of fs11-u: %g status %i at timestamp %lli\n",slctrl_fs11_u->value,slctrl_fs11_u->status,slctrl_fs11_u->timestamp);



  return 0;

}





Long_t SlowCtrlDemo::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};


extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new SlowCtrlDemo(in,out,inf_,outf_,p);
}
}


ClassImp(SlowCtrlDemo);

