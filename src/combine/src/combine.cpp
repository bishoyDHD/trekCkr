
#include <iostream>
#include <boost/program_options.hpp>
#include<vector>
#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"

#include <fstream>
//#include "sclumi.h"
namespace po = boost::program_options;


class hisinfo
{
public:
  std::string path;
  TH1* his;
};


std::map< TString, hisinfo> histos;

void chdir(TFile *f, std::string path)
{
  f->cd("/"); // go back to root
  std::stringstream ss(path);
  std::string item;
  while(std::getline(ss, item, '/')) 
    if (item!="")
      {
	if (!f->GetDirectory(item.c_str()))
	  f->mkdir(item.c_str());
	f->cd(item.c_str());
      }
}


void recurse(std::string path,double lum)
{
  TList* list =gDirectory->GetListOfKeys() ;
  TIter next(list);
  TKey *key;
  TObject* obj ;
  
  while ( key = (TKey*)next() ) {
    obj = key->ReadObj() ;
    if (obj->InheritsFrom("TDirectory"))
      { 
	gDirectory->cd(obj->GetName());
	recurse(path+"/"+obj->GetName(),lum);
	gDirectory->cd("..");
      }
    if (obj->InheritsFrom("TH1"))
      if ( histos.find(path+"/"+obj->GetName()) !=histos.end())
	{
	  histos[path+"/"+obj->GetName()].his->Add((TH1*) obj,1);
	  histos[path+"/"+obj->GetName()].his->SetNormFactor(histos[path+"/"+obj->GetName()].his->GetNormFactor()+lum);
	    }
      else
	{
	  ((TH1*) obj)->SetDirectory(0);
	  histos[path+"/"+obj->GetName()].path=path;
	  histos[path+"/"+obj->GetName()].his=(TH1*) obj;
	  histos[path+"/"+obj->GetName()].his->SetNormFactor(lum);
	}
    else
      delete obj;
  }
  
}

int main(int argc, char **argv)
{
  //setup options
  po::options_description desc("Combine -- Combine Histograms and export to gnuplot\nSyntax: combine [options] <input files>\n\tAllowed options:");
  desc.add_options()
    ("help,h","print help message")
    ("merge,m",po::value<std::string>(),"Merge trees with given name")
    ("root,R",po::value<std::string> (),"export histograms to root file")
    ("gnuplot,G",po::value<std::string> (),"export histograms to gnuplot file")
    ("only,o",po::value< std::vector<std::string> > (),"only export these histograms")
    ("normalize,N","Normalize to slow control luminosity (needs SCLumiInfo in *all* files)")
    ("inputfile,i",po::value< std::vector<std::string> > (),"input file")
    ;
  
  po::positional_options_description p;
  p.add("inputfile", -1);
  
  po::variables_map vm;
  po::parsed_options parsopt=po::command_line_parser(argc, argv).options(desc).positional(p).run();
  po::store(parsopt, vm); 
  po::notify(vm);
  if (vm.count("help") || ( vm.count("root")==0 && vm.count("gnuplot")==0))
    {
      std::cout<<desc<<std::endl;
      return 1;
    }
  
  bool only=vm.count("only")>0;
  std::vector<std::string> onlys;
  if (only) onlys=vm["only"].as<std::vector<std::string> >();

  std::set<TString> onlyset;

  if (only)
    for (unsigned int i=0;i<onlys.size();i++)
	   onlyset.insert(onlys[i]);
  bool normalize=vm.count("normalize")>0;
  if (normalize) std::cout<<"Will normalize to slow control luminosity\n";

  TList trees;


  std::vector< std::string > filelist=vm["inputfile"].as< std::vector<std::string> >();
  for (unsigned int i=0;i<filelist.size();i++)
    {
      std::cout<<"Reading "<<filelist[i]<<std::endl;
      TFile *f=new TFile(filelist[i].c_str(),"READ");
      if (vm.count("merge")>0)
	{
	  if (TTree *tree =(TTree *)f->Get(vm["merge"].as<std::string>().c_str()))
	    trees.Add(tree);
	  else
	    std::cout<<"File "<<filelist[i]<<" has no tree named "<<vm["merge"].as<std::string>()<<std::endl;
	}

      if (normalize)
	{
	  // doesn't work in.cooker because sclumi.h is unknown...
	  // SCLumiInfo *li=(SCLumiInfo*)f->Get("SCLumiInfo");
	  // if (!li)
	  //   std::cout<<"Skipping file "<<filelist[i]<<": No slowctrl info\n";
	  // else
	  //   {
	  //     std::cout<<"Luminosity:"<<li->dtcLuminosity<<std::endl;
	  //     recurse("",li->dtcLuminosity);
	  //   }
	}
      else
	recurse("",0);
      //  f.Close();
    }



  if (vm.count("root"))
    {
      std::cout<<"Writing histograms to root file: "<<vm["root"].as<std::string>()<<std::endl;
      TFile f(vm["root"].as<std::string>().c_str(),"RECREATE");
       

      for (std::map<TString,hisinfo>::iterator it=histos.begin();it!=histos.end();it++)
	{
	  chdir(&f,it->second.path);
	  it->second.his->Write();
	}
      f.cd();
      f.cd("/");
      if (vm.count("merge")>0)
	{
	  std::cout<<"Merging trees... this might take a while\n";
	  TTree::MergeTrees(&trees)->Write();
	}
      
      f.Close();
    }

  if (vm.count("gnuplot"))
    {
      std::cout<<"Writing histograms to gnuplot file: "<<vm["gnuplot"].as<std::string>()<<std::endl;
      double norm=1;
      std::ofstream out(vm["gnuplot"].as<std::string>().c_str());

      int i=0;
      out<<"#Directory:\n";
      for (std::map<TString,hisinfo>::iterator it=histos.begin();it!=histos.end();it++)
	if ( !only || onlyset.count(it->first)>0)
	  out<<"#Index "<<i++<<": "<<it->first<<" "<<it->second.his->GetName()<<" "<<it->second.his->IsA()->GetName()<<std::endl;

      

      for (std::map<TString,hisinfo>::iterator it=histos.begin();it!=histos.end();it++)
	{
	  if ( only && onlyset.count(it->first)==0) continue;
	  out<<"# "<<it->second.his->GetName()<<std::endl;
	  if (normalize)
	    norm=it->second.his->GetNormFactor();

	  if (it->second.his->InheritsFrom("TH3"))
	    {
	      out<<"# 3d histos not possible with gnuplot"<<std::endl;
	      
	      // 3 dim histo
	      /*      TH3 *h=(TH3*) it->second;
	      for (int z=0;z<h->GetNbinsZ();z++)
		for (int y=0;y<h->GetNbinsY();y++)
		  {
		    for (int x=0;x<h->GetNbinsX();x++)
		      out<<h->GetXaxis()->GetBinCenter(x)<<" "<<h->GetYaxis()->GetBinCenter(y)<<" "<<h->GetZaxis()->GetBinCenter(z)<<" "<<h->GetBinContent(x,y,z)<<std::endl;
		      out<<std::endl;
		}
	      */
	    }
	  else
	  if (it->second.his->InheritsFrom("TH2"))
	    {
	      // 2 dim histo
	      TH2 *h=(TH2*) it->second.his;
	      for (int y=0;y<h->GetNbinsY();y++)
		{
		  for (int x=0;x<h->GetNbinsX();x++)
		    out<<h->GetXaxis()->GetBinCenter(x)<<" "<<h->GetYaxis()->GetBinCenter(y)<<" "<<h->GetBinContent(x,y)<<" "<<norm<<std::endl;
		    out<<std::endl;
		}
	    }
	   else
	     {
	       TH1 *h=it->second.his;
	       for (int x=0;x<h->GetNbinsX();x++)
		 out<<h->GetXaxis()->GetBinCenter(x)<<" "<<h->GetBinContent(x)<<" "<<norm<<std::endl;
	       out<<std::endl;
	     }   
	  out<<std::endl; // next record
	}


      out.close();

    }



  return 0;
}
