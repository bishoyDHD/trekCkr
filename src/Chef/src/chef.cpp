#include "chef.h"
#include "TTimeStamp.h"
#include "TFile.h"

#include "TSystem.h"
#include "TFriendElement.h"
#include "cookerrawtree.h"
#include "InitReader.h"
#include "TBranchElement.h"
#include <boost/algorithm/string.hpp>
#include "MIDASreader.h"

#include <iostream>
#include <dlfcn.h>

#include <gsl/gsl_rng.h>


#include "TBufferFile.h"
#include "TH1.h"

CRTEventInfo dummy;

callinfo::callinfo(Plugin * p,TClass *cl, std::string method)
{
  mc.Init(cl,method.c_str(),"");
  plug=p;
}

Long_t callinfo::execute()
{
  Long_t rval;
  Long_t kval = 10;
  mc.Execute(plug,rval);
  if (!mc.IsValid()) {kval = -9001;}
  if ((rval<0) || (kval<0))
    {
      std::cout<<"\n********************************************************************************\n";
      std::cerr<<"\nError while calling: "<<plug->IsA()->GetName()<<"::"<<mc.GetMethodName()<<"\n";
      std::cerr<<"The function returned value: "<<rval<<"\n\n";
      std::cerr<<"Check above for errors from the function and/or make sure the function exists!\n\n";
      std::cout<<"********************************************************************************\n\n";
      abort(); // Kill all processes
    }
  return rval;
}

std::vector<callinfo*> Chef::compilelist(std::vector<class_method> &list)
{
  std::vector<callinfo*> erg;
  for (std::vector<class_method>::iterator iter=list.begin();iter!=list.end();iter++)
    {
      std::cout<<"  "<<iter->classname<<"::"<<iter->method<<std::endl; 
      callinfo *mc=new callinfo(plugins[iter->classname],plugins[iter->classname]->IsA(),iter->method);
      erg.push_back(mc);
    }
  return erg;
}


Chef::Chef(std::string recipename,unsigned int seed,unsigned int gskip, unsigned int gseed):recipe(recipename)
{
  // set dynamic path
  gSystem->AddDynamicPath("~/.cooker/" ARCHDIR "/lib");
  
  //set COOKERHOME
  setenv("COOKERHOME",getenv("HOME"),0);
  printf("COOKERHOME: %s\n",getenv("COOKERHOME"));

  weight=1;
  addRepo("Weight",(TObject *) &weight);
  lastresult=Plugin::ok;
  addRepo("LastEventResult",(TObject *) &lastresult);
  
  //pseudo random
  gsl_rng *random=gsl_rng_alloc(gsl_rng_mt19937);
  if (seed==0)
    {
      FILE *f=fopen("/dev/urandom","r");
      fread(&seed,sizeof(unsigned int),1,f);
      fclose(f);
    }
  if (gseed==0)
    {
      FILE *f=fopen("/dev/urandom","r");
      fread(&gseed,sizeof(unsigned int),1,f);
      fclose(f);
    }
  gsl_rng_set(random,seed); 

  genskip=gskip;
  genseed=gseed;
  globseed=seed;
  
  printf("Used seeds: %i %i %i\n\n",seed,gskip,gseed);
  addRepo("Random",(TObject *) random);
  addRepo("GeneratorSkip",(TObject *) &genskip);
  addRepo("GeneratorSeed",(TObject *) &genseed);
  addRepo("GlobalSeed",(TObject *) &globseed);
  addRepo("Callbacks",(TObject *)this);
}

void Chef::prepareTrees(std::string input, std::string output,bool empty)
{

  // split up input string and recipe str by ":"

  std::vector <std::string> inputfilenames,inputtreetypes;
  std::cout<<input<<" "<<output<<" "<<empty<<"\n";
  
  if (empty)
    {
      infile=new TFile();
      in=new TTree;
    }
  else
    {
      boost::split(inputtreetypes,recipe.srctree,boost::is_any_of(":"));
      boost::split(inputfilenames,input,boost::is_any_of(":"));
      
      /* <-- Temporary disabled in case MIDAS libs not needed -Bishoy
      if (boost::algorithm::ends_with(inputfilenames[0],".mid"))
	//infile=new MIDASfile(inputfilenames[0].c_str()); 
      else
	inifile=TFile::Open(inputfilenames[0].c_str(),"READ");*/
      in=(TTree*)infile->Get(inputtreetypes[0].c_str());
      if (!in)
	{
	  std::cerr<<"Tree "<<recipe.srctree<<" not found in "<<input<<std::endl;
	  std::cerr<<"Maybe wrong recipe for this type of tree?"<<std::endl;
	  exit(-1);
	}
      // adding all friends:
      for (unsigned int i=1;i<inputtreetypes.size();i++)
	if (!in->AddFriend(inputtreetypes[i].c_str(),inputfilenames[i].c_str())->GetTree())
	  {
	    std::cerr<<"Could not add friend tree of type "<<inputtreetypes[i]<<" from file "<<inputfilenames[i]<<". Bailing out."<<std::endl;
	    exit(-1);
	      }
      //disabling all branches
      in->SetBranchStatus("*",0);
    }
  
  outfile=0;

  if (output!="")
    {
      std::cout<<"Creating output file:"<<output;
      if (output.find("root://",0)==0)
	outfile=TFile::Open(output.c_str(),"NEW");  // root network protocol. Does not allow recreate (maybe because of dcache backend?
      else
	outfile=TFile::Open(output.c_str(),"RECREATE"); //normal file -> Recreate
      if (outfile) std::cout<<"................success"<<std::endl;
      else {
	std::cout<<"................fail"<<std::endl;
	exit(-5);
      }
      out=new TTree(recipe.dsttree.c_str(),recipe.dsttree.c_str());

      // try to copy runinfo and eventinfo
      
      if (infile->Get("RunInfo"))
	infile->Get("RunInfo")->Write("RunInfo");

      TObject **p = new (TObject *);
      in->SetBranchStatus("EventInfo",1);
      TBranchElement *br=(TBranchElement *)in->GetBranch("EventInfo");
      if (br)
	{
	  *p=(TObject *)br->GetObject();
	  if (*p==0)
	    in->SetBranchAddress("EventInfo",p); 
	  out->Branch("EventInfo",p,16000,0); 
	  addRepo("Weight",(TObject *) & ((CRTEventInfo *) *p)->weight);
	}      
    }
  else
    {
      std::cout<<"No output file given. Will use dummy tree and discard."<<std::endl;
      out=new TTree("dummy","dummy");
    }
}

void Chef::loadPlugins(int rank)
{
  std::cout<<std::endl<<"Plugins:"<<std::endl;
  
  for (std::map<std::string,std::string>::iterator iter=recipe.plugins.begin();iter!=recipe.plugins.end();iter++)
    { 
      std::cout<<"Loading: "<<iter->second<<" as "<<iter->first;
      // DESPITE THE DOCUMENTATION WHICH SAYS IT CAN

      //  void *handle=dlopen(gSystem->ExpandPathName((std::string("~/.cooker/"ARCHDIR"/lib/")+iter->second+std::string(LIBSUFFIX)).c_str()),RTLD_LAZY);

      void *handle=dlopen((iter->second+std::string(LIBSUFFIX)).c_str(),RTLD_LAZY);
      if (handle) std::cout<<"................success"<<std::endl;
      else  {
	std::cout<<"................fail"<<std::endl<<dlerror()<<std::endl;
	exit(-2);
      }
      std::cout<<"Creating factory";
      
      Plugin *(*factory)(TTree *,TTree*,TFile *, TFile*, TObject *)=(Plugin* ( *)(TTree *,TTree*,TFile*, TFile *, TObject *)) dlsym(handle,"factory");
      if (factory) std::cout<<"................success"<<std::endl;
      else {
	std::cout<<"................fail"<<std::endl;
	exit(-3);
      }
      std::cout<<"Creating object";
      Plugin *plug=factory(in,out,infile,outfile,(TObject *) &repo);
      if (plug) std::cout<<"................success"<<std::endl;
      else {
	std::cout<<"................fail"<<std::endl;
	exit(-4);
      }
      plugins[iter->first]=plug;
      plug->rank=rank;
    }

  std::cout<<std::endl<<"defineHistograms:"<<std::endl;
  cldefineHistograms=compilelist(recipe.defineHistograms);
  std::cout<<std::endl<<"Startup:"<<std::endl;
  clstartup=compilelist(recipe.startup);
  std::cout<<std::endl<<"Execute:"<<std::endl;
  clexecute=compilelist(recipe.commands);
  std::cout<<std::endl<<"Post Process:"<<std::endl;
  clpostprocess=compilelist(recipe.postprocess);
  std::cout<<std::endl<<"Execute2:"<<std::endl;
  clexecute2=compilelist(recipe.commands2);
  std::cout<<std::endl<<"Finalize:"<<std::endl;
  clfinalize=compilelist(recipe.finalize); 
}



void Chef::processInit(int debug,std::map<std::string,std::vector<std::pair<std::string,std::string> > > pluginoptions)
{
  // Read init file
  CRTRunInfo *ri=0;
  if (infile)
    ri=(CRTRunInfo*)  infile->Get("RunInfo");

  TTimeStamp ttsstarttime;
  if (ri==0)
    {
      std::cerr<<"Could not read runinfo to get start time. Defaulting to *now*.\n";
      addRepo("Runnumber",(TObject *) 0);
      ri= new CRTRunInfo;
      ri->runNumber=0;
    }
  else
    {
      ttsstarttime=ri->startTime;
      addRepo("Runnumber",(TObject *) ri->runNumber);
    }  
  

  std::string starttimestr =std::string(ttsstarttime.AsString("c"));
  // root is "almost" ISO8601 compliant. Great. Almost. 
  starttimestr[10]='T'; // hopefully it is now.

  std::cout<< "Using run start time:"<<starttimestr<<" Run number:"<< ri->runNumber<<std::endl;
  InitReader init(recipe.InitXML,starttimestr,ri->runNumber);
  std::cout<<std::endl;


  for (std::map<std::string,Plugin *>::iterator piter=plugins.begin();piter!=plugins.end();piter++)       // iterate over all plugins
    { 

      std::cout<<"Configuring "<<piter->first<<std::flush;

      // set debug level from command line
      if (debug)
	piter->second->setDebug(debug);

      //set init options:
      std::map<std::string,std::vector<std::string> > calls=init.getConfig(piter->first);
      std::cout<<"#"<<std::flush;

      for ( std::map<std::string,std::vector<std::string> >::iterator iter=calls.begin();iter!=calls.end();iter++)
	{
	  TMethodCall mc(piter->second->IsA(),iter->first.c_str(),iter->second[0].c_str());
	  for ( std::vector<std::string>::iterator iter2=iter->second.begin() ; iter2!=iter->second.end();iter2++) 	 
	    {
	      Long_t rval;
	      mc.Execute(piter->second,iter2->c_str(),rval);
	      if (rval<0)
		{
		  std::cerr<<std::endl<<"Error with Nr. "<<rval<<" calling "<<piter->first<<"::"<<iter->first<<" with:"<<(*iter2)<<std::endl;
		  exit(rval-1000);
		}
	    }
	}
      // set debug level from command line (overrides init settings)
      if (debug)
	piter->second->setDebug(debug);

      std::cout<<"................done"<<std::endl;
    }
  for (std::map<std::string,std::vector<std::pair<std::string,std::string> > >::iterator iter=pluginoptions.begin();iter!=pluginoptions.end();iter++)
    {
      // try to call commandline specified routines of the plugin, including cmdline
      if (plugins[iter->first])
	for (std::vector<std::pair<std::string,std::string> >::iterator iter2=iter->second.begin();iter2!=iter->second.end();iter2++)
	  {
	  
	    TMethodCall mc(plugins[iter->first]->IsA(),iter2->first.c_str(),iter2->second.c_str());
	    Long_t rval;
	    mc.Execute(plugins[iter->first],rval);
	    if (rval<0)
	      {
		std::cerr<<std::endl<<"Error with Nr. "<<rval<<" calling "<<iter2->first  <<" of "<<iter->first<<std::endl;
		exit( rval-1000);
	      }
	  }
	
      else
	{
	  std::cerr<<"Unknown option "<<iter->first<<" No plugin with that name defined!"<<std::endl;
	}
    }
   
}

void Chef::defineHistograms()
{
  int rval=0;
  for (std::vector<callinfo*>::iterator iter=cldefineHistograms.begin();iter!=cldefineHistograms.end();iter++)
    {
      gDirectory->cd("/");
      if ((rval=(*iter)->execute())<0)
	exit(rval-1000);
    }
}

void Chef::startup()
{
  int rval=0;
  for (std::vector<callinfo*>::iterator iter=clstartup.begin();iter!=clstartup.end();iter++)
    {
      gDirectory->cd("/");
      if ((rval=(*iter)->execute())<0)
	exit(rval-1000);
    }
}

int Chef::processEvent(int i)
{
  int rval=0;
  if (i>=0)
    in->GetEntry(i);
  bool mskip=true;
  int returncode=0;
  weight=1; // set default weight.

  for (std::vector<callinfo*>::iterator iter=clexecute.begin();iter!=clexecute.end();iter++)
    {
      gDirectory->cd("/");
      if ((rval=(*iter)->execute())<0)
	exit(rval-1000);
      if (rval & Plugin::skip)
	{
	  lastresult=Plugin::skip;
	  return Plugin::skip; // We skip, no fill
	}
      if (rval & Plugin::stop)
	returncode|= Plugin::stop; // We stop, no fill, and loop should stop.
      if (!(rval & Plugin::maySkip)) // we skip if ALL plugins give maySkip
	mskip=false;      
      if (rval & Plugin::redo)
	returncode|=Plugin::redo;
    }

  //  if (!mskip && outfile) out->Fill();
  //if (outfile) out->Fill();
  int pass2 = clexecute2.size();
  if (!pass2 && outfile) out->Fill();
  lastresult=returncode;
  return returncode;
}

int Chef::processEvent2(int i)
{
  int rval=0;
  if (i>=0)
    in->GetEntry(i);
  bool mskip=true;
  int returncode=0;
  weight=1; // set default weight.

  for (std::vector<callinfo*>::iterator iter=clexecute2.begin();iter!=clexecute2.end();iter++)
    {
      gDirectory->cd("/");
      if ((rval=(*iter)->execute())<0)
	exit(rval-1000);
      if (rval & Plugin::skip)
	{
	  lastresult=Plugin::skip;
	  return Plugin::skip; // We skip, no fill
	}
      if (rval & Plugin::stop)
	returncode|= Plugin::stop; // We stop, no fill, and loop should stop.
      if (!(rval & Plugin::maySkip)) // we skip if ALL plugins give maySkip
	mskip=false;      
      if (rval & Plugin::redo)
	returncode|=Plugin::redo;
    }

  //  if (!mskip && outfile) out->Fill();
  if (outfile) out->Fill();
  lastresult=returncode;
  return returncode;
}

void Chef::postprocess()
{
  //std::cout << "postprocess size: " << clpostprocess.size() << std::endl;
  int rval=0;
  for (std::vector<callinfo*>::iterator iter=clpostprocess.begin();iter!=clpostprocess.end();iter++)
    {
      gDirectory->cd("/");
      if ((rval=(*iter)->execute()<0))
	exit(rval-1000);
    }
}

int Chef::secondpasssize()
{
  //std::cout << "execute2 size: " << clexecute2.size() << std::endl;
  return clexecute2.size();
}

void Chef::finalize()
{
  int rval=0;
  for (std::vector<callinfo*>::iterator iter=clfinalize.begin();iter!=clfinalize.end();iter++)
    {
      gDirectory->cd("/");
      if ((rval=(*iter)->execute()<0))
	exit(rval-1000);
    }
  std::cout<<"--Saving Tree--"<<std::endl;
  std::cout<<"Source tree:"<<std::endl;
  in->Print();
  std::cout<<"Destination tree:"<<std::endl;
  out->Print();
  if (outfile) outfile->Write();
}

void  Chef::addRepo(std::string name,TObject *obj)
{
  repo[name]=obj;
}


void Chef::serializeBranches(TBufferFile &buf)
{
  auto branches=out->GetListOfBranches();
  for (int i=0;i<branches->GetEntries();i++)
    {
      TBranchElement *br=(TBranchElement*) branches->At(i);
      TObject *obj=(TObject*)br->GetObject();
      //printf("%s %s\n",br->GetClassName(),obj->GetName());
      obj->Streamer(buf);
      //
      //br->Streamer(buf);
    }
}

void Chef::serializeHistogramList(TBufferFile &buf)
{

  buf.WriteInt(histograms.size());
  for (auto&& h :histograms)
    {
      TString name=h.first;
      name.Streamer(buf);
    }
}

void Chef::serializeHistogram(TString name, TBufferFile &buf)
{
  if (buf.IsReading())
    {
      TH1 *h=(TH1 *)buf.ReadObject(NULL);
      if (histograms.find(name)!=histograms.end())
	{
	  histograms[name]->Add(h);
	  //	delete h;
	}
      else
	{
	  gDirectory->Cd("/");
	  histograms[name]=h;
	  TObjArray *dirs=TString(name).Tokenize("/");
	  for (int i=0;i<dirs->GetEntriesFast()-1;i++)
	    {
	      if (!gDirectory->GetDirectory(((TObjString *) dirs->At(i))->GetString().Data()))
		gDirectory->mkdir(((TObjString *) dirs->At(i))->GetString().Data());
	      gDirectory->cd(((TObjString *) dirs->At(i))->GetString().Data());
	    }
	  delete dirs;
	  gDirectory->Add(h);
	}
    }
  else
    {
      buf.WriteObject(histograms[name]);
      //      histograms[name]->Streamer(buf);
    }
}

void Chef::addHisto(const char *path,TH1D *h)
{
  TString tsp(path);
  tsp.Remove(0,tsp.First(":")+1);
  histograms[tsp]=h;
}
