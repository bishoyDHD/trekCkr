
#include "Plugin.h"

#include <iostream>
#include <sstream>
#include <boost/program_options.hpp>
#include <sys/time.h>
#include <algorithm>
#include "chef.h"

#include <boost/array.hpp>
#include <boost/asio.hpp>
#include <boost/format.hpp>

using boost::asio::ip::udp;


#ifdef WITHMPI
#include <cstdlib>
#include <boost/mpi.hpp>
#include "TBufferFile.h"
namespace mpi = boost::mpi;
#endif


timeval costarttime;

double getdifftime()
{
  timeval tv;
  gettimeofday(&tv,0);
  return (tv.tv_sec-costarttime.tv_sec)*1000+(tv.tv_usec-costarttime.tv_usec)*1.0/1000;
}


namespace po = boost::program_options;



#ifdef WITHMPI
mpi::environment *env=NULL;
mpi::communicator *world=NULL;

auto mpirecv(int s,int t)
{
  auto status=world->probe(s,t);
  boost::optional<int> size=status.count<char>();
  std::unique_ptr<TBufferFile> buf(new TBufferFile(TBuffer::kRead,*size));
  world->recv(status.source(),status.tag(),buf->Buffer(),*size);
  return std::make_pair(status,std::move(buf));
}
#endif

int main(int argc, char **argv)
{
  //setup options
  po::options_description desc("Cooker -- Raw to cooked with a recipe\nSyntax: cooker [options] <recipe> <source tree> <target tree>\n\tAllowed options:");
  desc.add_options()
    ("help,h","print help message")
    ("init,i",po::value<std::string> (),"override initialization xml file")
    ("start,s",po::value<unsigned int> ()->default_value(0),"Number of first event to process (starts with 0)")
    ("split,S",po::value<unsigned int>(),"Split in this many parts")
    ("part,p",po::value<unsigned int>(),"do this part")
    ("num,n",po::value<unsigned int> (),"Number of events to process")
    ("recipe,r",po::value<std::string>(),"recipe")
    ("input_tree,I",po::value<std::string>(),"input root tree")
    ("output_tree,O",po::value<std::string>(),"output root tree")
    ("verbose,v",po::value<int>()->implicit_value(1),"Verbose mode (optionally specify level)")
    ("every,e",po::value<int>()->default_value(500),"Print out status every x events")
    ("call,c",po::value<std::vector<std::string> >(),"Call a plugin's function. Needs argument with format: <plugin>:<function>:<arguments>") 
    ("monitor,m",po::value<std::string>(),"send monitor udp to this host (port 5555)")
    ("remote","Switches on remote mode -- internal use only")
#ifdef WITHMPI
    ("mpi,M","Start as MPI master program")
    ("mpiworkers,W",po::value<unsigned int>()->default_value(4),"Use N workers (default: 4)")
    ("mpislave","Start as MPI worker program (internal use)")
#endif
    ("fake_input,F","Fake input")
    ;


  po::positional_options_description pod;
  pod.add("recipe",1).add("input_tree",1).add("output_tree",1);
  po::variables_map vm;
  po::parsed_options parsopt=po::command_line_parser(argc, argv).options(desc).positional(pod).allow_unregistered().run();
  po::store(parsopt, vm); 
  po::notify(vm);
 
#ifdef WITHMPI
  if (vm.count("mpi")>vm.count("mpislave"))
    {
      std::cerr<<"Starting as MPI master\n";
      char **args = (char **)calloc(argc + 5, sizeof(char *));
      args[0] = "mpirun";
      args[1] = "-np";
     
      char parnum[100];
      sprintf(parnum,"%u",vm["mpiworkers"].as<unsigned int>());
      args[2] = parnum;
      for (int i=0;i<argc;i++)
	{
	  args[3+i]=argv[i];
	  printf("%s\n",argv[i]);
	}
      args[3+argc]= "--mpislave";
      args[4+argc]=0;
      int i=execvp("mpirun", args);
      std::cerr<<"Got:"<<i<<" "<<strerror(errno)<<"\n";
      return EXIT_SUCCESS;
    }
  if (vm.count("mpislave")>0)
    {
      env=new mpi::environment(argc, argv);
      world=new mpi::communicator; 
      for (int i=0;i<argc;i++)
	{
	  std::cerr<<world->rank()<<" args: "<<i<<" "<<argv[i]<<"\n";
	}
      std::cerr<<"Starting as MPI slave rank:"<<world->rank()<<"\n";
    }
#endif

  if (vm.count("help") ||( vm.count("recipe")!=1 )|| ((vm.count("input_tree")!=1) && (vm.count("fake_input")==0)))
    {
      std::cerr<<desc<<std::endl;
      return 1;
    }

  int every=vm["every"].as<int>();

  if (vm.count("split")!=vm.count("part"))
    {
      std::cerr<<"Need both split and part, or neither\n";
      exit(-100);
    }

  int debug= 0;

  if (vm.count("verbose"))
    {
      debug=vm["verbose"].as<int>();
      std::cerr<<"Setting debug level to:"<<debug<<std::endl;
    }

  Chef chef(argv[1]);

  // init file  
  if (vm.count("init"))
    chef.recipe.InitXML=vm["init"].as<std::string>();
   
  std::cerr<<"Init file:  "<<chef.recipe.InitXML<<std::endl;

  bool remote=false;
  std::string input;
  if (vm.count("fake_input")==0)
    input=vm["input_tree"].as<std::string>();
  if (vm.count("remote"))
    {
      remote=true;
      char *tmpfilename=tempnam(NULL,"cooker.temp");
      std::cout<<"##REMOTE## filename "<<tmpfilename<<" \n";
      chef.prepareTrees(input,tmpfilename,vm.count("fake_input"));
    }
  else
#ifdef WITHMPI
    if (world && (world->rank()>0))
      chef.prepareTrees(input,"",vm.count("fake_input"));
    else
#endif 
      {
	if (vm.count("output_tree"))  
	  chef.prepareTrees(input,vm["output_tree"].as<std::string>(),vm.count("fake_input"));
	else
	  chef.prepareTrees(input,"",vm.count("fake_input"));
      }

  bool sendudp=false;
  boost::asio::io_service io_service;
  boost::asio::ip::udp::socket socket(io_service);
  boost::asio::ip::udp::resolver resolver(io_service);
  boost::asio::ip::udp::endpoint dest_endpoint;
  if (vm.count("monitor"))
    {
    
      dest_endpoint= *resolver.resolve(udp::resolver::query(udp::v4(),vm["monitor"].as<std::string>(), "5555"));
      socket.open(dest_endpoint.protocol());
      sendudp=true;
    }

  if (vm.count("fake_input"))
    {
      if (vm.count("num"))
	chef.in->SetEntries(vm["num"].as<unsigned int>());
      else
	chef.in->SetEntries(1e6);
    }

#ifdef WITHMPI
  if (world)
    chef.loadPlugins(world->rank());
  else
    chef.loadPlugins();
#else
  chef.loadPlugins();
#endif WITHMPI


  // make map from unknown options
  std::map<std::string,std::vector<std::pair<std::string,std::string> > > pluginoptions;
  for (std::vector<po::basic_option<char> >::iterator iter=parsopt.options.begin();iter!=parsopt.options.end();iter++)

    if (iter->unregistered )
      if  (iter->value.size()==1)
	{
	  pluginoptions[iter->string_key].push_back(std::pair<std::string,std::string>("cmdline",iter->value[0]));
	}
      else
	std::cerr<<"Unknown token:"<< iter->string_key<<std::endl;

    
  if (vm.count("call"))
    {
      std::vector<std::string> ol=vm["call"].as<std::vector<std::string> >();
      for (std::vector<std::string>::iterator iter=ol.begin();iter!=ol.end();iter++)
	{
	  int part1, part2;
	  part1=iter->find(":");
	  part2=iter->rfind(":");
	  if (part2>part1)
	    {
	      std::cerr<<iter->substr(0,part1)<<" - "<<iter->substr(part1+1,part2-part1-1)<<" - "<<iter->substr(part2+1,iter->length())<<std::endl;
	      pluginoptions[iter->substr(0,part1)].push_back(std::pair<std::string,std::string>(iter->substr(part1+1,part2-part1-1),iter->substr(part2+1,iter->length())));
	    }
	  else
	    {std::cerr<<"Could not parse option "<<*iter<<std::endl;
	      exit(-20);
	    }
	}	   
    }


#ifdef WITHMPI
  
  if ((world) && (chef.secondpasssize()>0))
    {
      std::cerr<<"Can not have a second pass in MPI mode\n";
      exit(-1);
    }
#endif
	    chef.processInit(debug,pluginoptions);


  std::cerr<<"--Define Histograms--"<<std::endl;
  chef.defineHistograms();
  std::cerr<<"--Startup--"<<std::endl;
  chef.startup();

  std::cerr<<"--Looping--"<<std::endl;


  int count=chef.in->GetEntries();
  int num=count;
  
  if (vm.count("num"))
    num=vm["num"].as<unsigned int>();

  int start=vm["start"].as<unsigned int>();
  if (vm.count("split"))
    {
      num=(int) (0.5+num*1.0/vm["split"].as<unsigned int>());
      start=num*(vm["part"].as<unsigned int>());

    }
  count=std::max(std::min(count-start,num),0);

#define blocksize 500
#ifdef WITHMPI
  TBufferFile *sbuf; //send and receive buffers;
  std::map<int,std::unique_ptr<TBufferFile> > doneblocks;
  if (world && world->rank()==0)
    {
      gettimeofday(&costarttime,0);
      // We don't do main loop, but distribute work and write the output file
      sbuf=new TBufferFile(TBuffer::kWrite,10000);
      int blocks=(count+blocksize-1)/blocksize;
      int finishedblock=0;
      int sendblocks=0;
      int current=start;
      for (;;)
	{
	  //	  std::cerr<<"Main receive\n";
	  auto rec=mpirecv(mpi::any_source, mpi::any_tag);
	  //  std::cerr<<"Main receive done\n";
	    
	  switch(rec.first.tag())
	    {
	    case 0:
	      {
		if (count ==0)
		  break;
		sbuf->Reset();
		sbuf->WriteInt(start+blocksize*sendblocks);
		sbuf->WriteInt(std::min(blocksize,count));
		sbuf->WriteInt(sendblocks);
		count-=std::min(blocksize,count);
		sendblocks++;
		world->send(rec.first.source(),0,sbuf->Buffer(),sbuf->Length());
	      }
	      break;
	    case 1:
	      int start,count,block;
	      finishedblock++;

	      rec.second->ReadInt(start);
	      doneblocks[start]=std::move(rec.second);

	     
	      while(doneblocks.find(current)!=doneblocks.end())
		{
		  std::unique_ptr<TBufferFile> buf=std::move(doneblocks[current]);
  		  doneblocks.erase(current);

		  buf->ReadInt(count);
		  buf->ReadInt(block);
		  for (int i=0;i<count ;i++)
		    {
		      chef.serializeBranches(*buf);
		      chef.out->Fill();
		      current++;
		    }
		}
	      unsigned int tdiff=getdifftime();
	      unsigned int cps=(unsigned int) (finishedblock*blocksize*1.0/tdiff*1000.0);
	      unsigned int eta= (unsigned int) (tdiff*1.0*(blocks-finishedblock)/finishedblock/1000);
	      tdiff/=1000;
	      //	      std::cerr<<"Work result received from rank "<<status.source()<<" Size:"<<(*size)<<"\n";
	      std::cerr<<"\r"<<finishedblock<<"/"<<blocks<<" + "<<doneblocks.size()<<" "<<finishedblock*100/blocks<<"%  "<<cps<<" Hz   ETA: ";
	      if (eta / 3600 >0) std::cerr<<eta/3600<<"h ";
	      if (eta / 60 >0) std::cerr<<eta % 3600 / 60<<"min ";
	      std::cerr<<eta % 60<<"s  Elapsed: ";
	      if (tdiff / 3600 >0) std::cerr<<tdiff/3600<<"h ";
	      if (tdiff / 60 >0) std::cerr<<tdiff % 3600 / 60<<"min ";
	      std::cerr<<tdiff % 60<<"s     "<<std::flush;
	     
	      break;
	    }
	  if (finishedblock==blocks)
	    break;
	}
      
      std::cerr<<"--Getting histogram list--\n";
      for (int i=1;i<world->size();i++)
	{
	  // get list of histograms
	  sbuf->Reset();
	  world->send(i,1,sbuf->Buffer(),sbuf->Length());
	  auto rec=mpirecv(i,2);
	  int hcount;
	  rec.second->ReadInt(hcount);
	  std::cerr<<"Rank "<<i<<" has "<<hcount<<" histograms\n";
	  std::vector<TString>paths;
	  for (int h=0;h<hcount;h++)
	    {
	      TString path;
	      rec.second->ReadTString(path);
	      paths.push_back(path);
	    }
	  std::cerr<<"Starting collection\n";
	  for (auto name:paths)
	    {
	      sbuf->Reset();
	      name.Streamer(*sbuf);
	      std::cerr<<name<<"\n";
	      world->send(i,2,sbuf->Buffer(),sbuf->Length());
	      auto rech=mpirecv(i,2);
	      chef.serializeHistogram(name,*rech.second);
	    }
	  

	}
	 
	  
      std::cerr<<"--Killing workers--\n";
      
      //kill other 
      for (int i=1;i<world->size();i++)
	world->send(i,10,sbuf->Buffer(),0);
      
      std::cerr<<std::endl<<"--Finalize--"<<std::endl;
      
      chef.finalize();
   
      std::cerr<<"--Done--"<<std::endl;
      delete env;
      delete world;
      exit(0);
    }

  if (world)
    {
      sbuf=new TBufferFile(TBuffer::kWrite,10*1024*1024);
      world->send(0,0,sbuf->Buffer(),0);
    }
  for(;;)
    {
      int blocknum=0;
      if (world)
	{
	  auto rec=mpirecv(0,mpi::any_tag);
	  if (rec.first.tag()==1) // get histograms
	    {
	      sbuf->Reset();
	      chef.serializeHistogramList(*sbuf);
	      world->send(0,2,sbuf->Buffer(),sbuf->Length());	    
	      continue;
	    }
	  if (rec.first.tag()==2) //get one histogram
	    {
	      TString name;
	      rec.second->ReadTString(name);
	      sbuf->Reset();
	      chef.serializeHistogram(name,*sbuf);
	      world->send(0,2,sbuf->Buffer(),sbuf->Length());
	      continue;
	    }
	  if (rec.first.tag()==10)
	    break;

	  rec.second->ReadInt(start);
	  rec.second->ReadInt(count);
	  rec.second->ReadInt(blocknum);
	  
	  
	  sbuf->Reset();
	  sbuf->WriteInt(start);
	  sbuf->WriteInt(count);
	  sbuf->WriteInt(blocknum);
	}
#endif
      gettimeofday(&costarttime,0);

      unsigned int ltdiff=0;
      for (int i=0;i<count;i++)
	{ 
	  int code=chef.processEvent(i+start);
#ifdef WITHMPI
	  if (world)
	    chef.serializeBranches(*sbuf);
#endif
	  if (debug>=100 && code!=0)
	    std::cerr<<"Loop gave:"<<code<<std::endl;

	  if ( (i % every) ==0  || sendudp) {
	    unsigned int tdiff=getdifftime();
	    unsigned int cps=(unsigned int) (i*1.0/tdiff*1000.0);
	    unsigned int eta= (unsigned int) (tdiff*1.0*(count-i)/i/1000);
	    tdiff/=1000;
#ifdef WITHMPI
	    if (!world)
#endif
	      if( (i %every)==0)
		{
		  std::cerr<<"\r"<<i<<"/"<<count<<" @ "<<i+start<<" " <<i*100/count<<"%  "<<cps<<" Hz   ETA: ";
		  if (eta / 3600 >0) std::cerr<<eta/3600<<"h ";
		  if (eta / 60 >0) std::cerr<<eta % 3600 / 60<<"min ";
		  std::cerr<<eta % 60<<"s  Elapsed: ";
		  if (tdiff / 3600 >0) std::cerr<<tdiff/3600<<"h ";
		  if (tdiff / 60 >0) std::cerr<<tdiff % 3600 / 60<<"min ";
		  std::cerr<<tdiff % 60<<"s     "<<std::flush;
		}
	
	    if (sendudp && (tdiff>ltdiff))
	      {
		ltdiff=tdiff;
		std::string msg;
#ifdef WITHMPI
		if (world)
		  msg= boost::str(boost::format("Rank %i: Block: %i ETA: %i s,  Elapsed: %i s\n") %world->rank() % blocknum % eta % tdiff);
		else
#endif
		  msg= boost::str(boost::format("%s: ETA: %i s,  Elapsed: %i s\n") % input % eta % tdiff);
		socket.send_to(boost::asio::buffer(msg.c_str(), msg.size()), dest_endpoint);
	      }
	


	    if (remote)
	      std::cout<<"\n##REMOTE## stat "<<i<<" "<<count<<" "<<cps<<" "<<eta<<std::endl;
	  }

	  if (code & Plugin::redo) i--;
	}
#ifdef WITHMPI
      if (!world)
#endif  
	std::cerr<<"\r"<<count<<"/"<<count<<" "<<100<<"%                       "<<std::flush;
#ifdef WITHMPI
      if (world)
	{
	  world->send(0,1,sbuf->Buffer(),sbuf->Length());
  	  world->send(0,0,sbuf->Buffer(),0);
	}
      else
	break;
    }
  if(world)
    {
      delete env;
      delete world;
    }
  else
    {
#endif
      
   
      std::cerr<<std::endl<<"--Post Process--"<<std::endl;
  
      chef.postprocess();
  
      if (chef.secondpasssize()) {
	std::cerr<<std::endl<<"--2nd Pass Looping--"<<std::endl;
	gettimeofday(&costarttime,0);

	int count=chef.in->GetEntries();
	int num=count;
	int start=vm["start"].as<unsigned int>();
	if (vm.count("num"))
	  num=vm["num"].as<unsigned int>();
	if (vm.count("split"))
	  {
	    num=(int) (0.5+num*1.0/vm["split"].as<unsigned int>());
	    start=num*(vm["part"].as<unsigned int>());

	  }

	count=std::max(std::min(count-start,num),0);

	unsigned int ltdiff=0;
	for (int i=0;i<count;i++)
	  { 
	    int code=chef.processEvent2(i+start);
	    if (debug>=100 && code!=0)
	      std::cerr<<"Loop gave:"<<code<<std::endl;

	    if ( (i % every) ==0  || sendudp) {
	      unsigned int tdiff=getdifftime();
	      unsigned int cps=(unsigned int) (i*1.0/tdiff*1000.0);
	      unsigned int eta= (unsigned int) (tdiff*1.0*(count-i)/i/1000);
	      tdiff/=1000;
	      if( (i %every)==0)
		{
		  std::cerr<<"\r"<<i<<"/"<<count<<" @ "<<i+start<<" " <<i*100/count<<"%  "<<cps<<" Hz   ETA: ";
		  if (eta / 3600 >0) std::cerr<<eta/3600<<"h ";
		  if (eta / 60 >0) std::cerr<<eta % 3600 / 60<<"min ";
		  std::cerr<<eta % 60<<"s  Elapsed: ";
		  if (tdiff / 3600 >0) std::cerr<<tdiff/3600<<"h ";
		  if (tdiff / 60 >0) std::cerr<<tdiff % 3600 / 60<<"min ";
		  std::cerr<<tdiff % 60<<"s     "<<std::flush;
		}
	
	      if (sendudp && (tdiff>ltdiff))
		{
		  ltdiff=tdiff;
		  std::string msg = boost::str(boost::format("%s: ETA: %i s,  Elapsed: %i s\n") % vm["input_tree"].as<std::string>() % eta % tdiff);
		  socket.send_to(boost::asio::buffer(msg.c_str(), msg.size()), dest_endpoint);
		}

	      if (remote)
		std::cout<<"\n##REMOTE## stat "<<i<<" "<<count<<" "<<cps<<" "<<eta<<std::endl;
	    }

	    if (code & Plugin::redo) i--;
	  }
  
  
	std::cerr<<"\r"<<count<<"/"<<count<<" "<<100<<"%                       "<<std::flush;
      }

      std::cerr<<std::endl<<"--Finalize--"<<std::endl;

      chef.finalize();

      std::cerr<<"--Done--"<<std::endl;

#ifdef WITHMPI
    }
#endif

  return 0;
}
