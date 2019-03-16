#ifndef __SLOWCTRLDEMO__
#define __SLOWCTRLDEMO__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include <iostream>
#include "slowctrl.h"


class SlowCtrlDemo:public Plugin
{
 private:
  slowctrl::datum * slctrl_fs11_u;

 public:
  SlowCtrlDemo(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~SlowCtrlDemo();
  // add funtions with return value Long_t here:
  
  Long_t startup();
    Long_t process();
  // Long_t finalize()

  virtual Long_t cmdline(char * cmd);

  ClassDef(SlowCtrlDemo,1);
    };

#endif
