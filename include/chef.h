#ifndef __CHEF_H__
#define __CHEF_H__

#include <vector>
#include <map>

#include "Plugin.h"
#include "RecipeReader.h"

#include "TMethodCall.h"
#include "TTree.h"



class TBufferFile;

///
// This class is a helper class for all cooker type applications
//
class callinfo
{
  TMethodCall mc;
  Plugin *plug;
public:
  callinfo(Plugin * p,TClass *cl, std::string method);
  Long_t execute();
  //Long_t IsValid();
};


class Chef:public FrameworkCallbacks
{
 private:
  std::map<std::string,Plugin *> plugins;
  std::map<std::string, TObject *>repo;
  std::vector<callinfo*> cldefineHistograms,clstartup,clexecute,clpostprocess,clexecute2,clfinalize;
  std::vector<callinfo*> compilelist(std::vector<class_method> &list);
  std::map<TString,TH1*> histograms;
  
  double weight;
  int lastresult;
  int genskip,genseed,globseed;

 public:
  RecipeReader recipe;
  TTree *in,*out;
  TFile *infile,*outfile;

  Chef(std::string recipename,unsigned int seed=0,unsigned int gskip=0, unsigned int gseed=0);
  void loadPlugins(int rank=-1);
  void processInit(int debug,std::map<std::string,std::vector<std::pair<std::string,std::string> > > pluginoptions);

  void prepareTrees(std::string input, std::string output,bool empty=false);
  void startup();
  void defineHistograms();
  int processEvent(int i);
  int processEvent2(int i);
  void postprocess();
  int secondpasssize();
  void finalize();
  void addRepo(std::string name,TObject *obj);
  void serializeBranches(TBufferFile &buf);
  void serializeHistogramList(TBufferFile &buf);
  void serializeHistogram(TString path,TBufferFile &buf);
  
  void addHisto(const char *path,TH1D *h);

  };


#endif
