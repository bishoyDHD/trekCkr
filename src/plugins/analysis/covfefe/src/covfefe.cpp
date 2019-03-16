#include <covfefe.h>

#include <iostream>
#include <cmath>
#include "TVector2.h"
#include "StrawTubehittree.h"
#include "StrawTubetree.h"
#include "Trackhittree.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"

#include <stdio.h>
#include <math.h>
//#include "~.cooker/external/cminpack/cminpack.h"
#include "../../../../../external/cminpack/cminpack.h"

int fcn(void *p, int m, int n, const double *x, double *fvec, int iflag);

covfefe::covfefe(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

covfefe::~covfefe()
{
};


TVector2 getStrawPos(int id)
{
  // ask gdml or something, for now, make it up.
  
  TVector2 res;
  int side, plane, straw_in_plane;
  STT_internal_to_logic(id,&side,&plane,&straw_in_plane);
  double tube_diameter = 10;
  double distance_between_planes = 8.7;
  double pitch=10.1; //10.1 mm spacing
  double zoffset=0;
  
  if (side==2)
  {
    double pos=straw_in_plane*pitch+(plane % 2)*pitch/2;
    res.Set(pos,zoffset+plane*distance_between_planes);
  }

  if(side == 1)
  {
    if(plane < 5) //these planes are vertical
    {
      //zoffset = 300;
      double pos = straw_in_plane*pitch +(plane % 2)*pitch/2;//original plus
      res.Set(pos,zoffset+plane*distance_between_planes);

    }
    // else if(plane >4 && plane<10)
    // {
    //   //zoffset = 380;
    //   double pos = straw_in_plane*pitch +(plane%2)*pitch/2;
    //   res.Set(pos,zoffset+plane*distance_between_planes);
    // }
  }
  return res;
}
//Used by TMinuit Minimizer
void mychi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  double thechi2 = 0.;
  const TVector2 off(80,20);
  double numhits = 0;

  // for(auto hit:STT->hits)
  // {
  //   int side, plane, straw_in_plane;
  //   STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
  //   auto pos=getStrawPos(hit.id)-off;
  //   if (side==1 && plane <5&& STT->hits.size()>3 && hit.dist>0.6 && plane!=0)
  //   {
  //     double slope = par[0];
  //     double intercept = par[1];
  //     double lineX = (pos.X()/slope+pos.Y()- intercept)/(slope+1/slope);//this is given by solving equation for two perpendicular intersecting lines 
  //     double lineY = par[0]*lineX+par[1];
  //     double res = sqrt(pow(pos.X()-lineX,2)+pow(pos.Y()-lineY,2))-hit.dist;
  //     thechi2 = thechi2 + pow(res,2)/pow(0.15,2);
  //     numhits++;
  //   }
    
  // }
  f = thechi2/(numhits-2);
  //return;
};


Long_t covfefe::startup()
{
  STT = NULL;
  getOutBranchObject("StrawTubeHits",(TObject **) &STT);
  if(!STT) getBranchObject("StrawTubeHits",(TObject **) &STT); 
  if(!STT) debug(0,"Could not find STT hits in file\n");

  Tracks = new TrackHits;
  makeBranch("TrackHits",(TObject **) &Tracks);

  hough=dH2("hough","hough",180,-0.5,179.5,401,-200,200);
  circ=dH2("circ","circ",401,0,400,401,-40,40);
  line=dH2("line","line",151,100,250,401,-40,40);
  lineprime=dH2("line prime","line prime",151,100,250,401,-40,40);

  both=dH2("Tracks","Tracks;Pos X (mm);Pos Z (mm)",151,100,250,401,-40,40);
  hitresidual=dH1("Hit Residual","Hit Residual; Residual (mm);counts",401,-10,10);
  hitchisq=dH1("Hit Chi Squared","Hit Chi Squared",401,-40,40);
  //strawresidual=dH1("Straw Residual","Straw Residual",401,-40,40);
  //strawchisq=dH1("Straw Chi Squared","Straw Chi Squared",401,-40,40);
  straws=dH2("Plane","Plane",51,150,200,41,-20,20);
  Angle=dH1("Angle","Angle",180,-0.5,179.5);
  chiVangle=dH2("Chi sq V Angle","Chi sq V angle;angle (degree);chisq",180,-0.5,179.5,500,0,500);
  chiVangleprime=dH2("Chi sq V Angle prime","Chi sq V angle prime;angle (degree);chisq",180,-0.5,179.5,500,0,500);
  
  return 0;
  
};



//CMin Stuff below

// int fcn(void *p, int m, int n, const double *x, double *fvec, int iflag)
// {
	
// 	//fcn is the name of the user-supplied subroutine which calculates
// 	//the functions. fcn should be written as follows:

// 	//subroutine fcn(m,n,x,fvec,iflag)
// 	//integer m, n, iflag
// 	//double precision x(n), fvec(m)

// 	//calculate the functions at x and return this vector in fvec

// 	//return

// 		double numhits = 0;
// 		double chisq = 0;
// 		double res = 0;
//   	const TVector2 off(80,20);
//   	TF1 *f = new TF1("f", "[0]*x+[1]",100,250);
//     int side, plane, straw_in_plane;
//     STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
// 		//STT = NULL
//   	//getOutBranchObject("StrawTubeHits",(TObject **) &STT);
//   	//if(!STT) getBranchObject("StrawTubeHits",(TObject **) &STT); 
//   	//if(!STT) debug(0,"Could not find STT hits in file\n");

//     for (auto hit:STT->hits)
//     {
//     		if (side==1 && plane <5&& STT->hits.size()>3 && hit.dist>0.6 && hit.dist<4.5){
//         numhits++;
// 				}
//         //ret=Plugin::stop;
//         int side, plane, straw_in_plane;
//         STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
//         if (side==1 && plane <5 && STT->hits.size()>3 && hit.dist>0.6 && hit.dist<4.5)
//         {
//           auto pos=getStrawPos(hit.id)-off;//this is wire position x',y'
//           //now need intercept of line perpendicular to fit that passes through the straw
//           //the x coordinate of intercept is given by
//           double slope = f->GetParameter(0);
//           double intercept = f->GetParameter(1);
//           double lineX = (pos.X()/slope+pos.Y()- intercept)/(slope+1/slope);//this is given by solving equation for two perpendicular intersecting lines 
//           double lineY = f->Eval(lineX);
//           res = sqrt(pow(pos.X()-lineX,2)+pow(pos.Y()-lineY,2))-hit.dist;
//           chisq = chisq + pow(res,2)/pow(0.15,2)/(numhits-2);
//         }
//     }

// }

int minimizer()
{
	int m, n, info, lwa;
	double  fnorm;
	
	m = 0; //This should be number of functions
	n = 0; //This should be number of variables, < m
	lwa = m*n + 5*n + m;
	
	int iwa[n];
	double x[n], fvec[m], wa[lwa];	
	//x should contain initial estimate, so need to initialize that

	double tol = 0; //This should be a precision

	//info = lmdif1(fcn, 0, m, n, fvec, tol, iwa, wa, lwa);

}



/*
Long_t covfefe::process()
{
  int ret=0;
	double tol = 0; //This should be a precision

	//info = lmdif1(fcn, 0, m, n, fvec, tol, iwa, wa, lwa);

	return 0;
}
*/

//end cmin stuff

/*
Long_t covfefe::startup()
{
  STT = NULL;
  getOutBranchObject("StrawTubeHits",(TObject **) &STT);
  if(!STT) getBranchObject("StrawTubeHits",(TObject **) &STT); 
  if(!STT) debug(0,"Could not find STT hits in file\n");

  Tracks = new TrackHits;
  makeBranch("TrackHits",(TObject **) &Tracks);

  hough=dH2("hough","hough",180,-0.5,179.5,401,-200,200);
  circ=dH2("circ","circ",401,0,400,401,-40,40);
  line=dH2("line","line",151,100,250,401,-40,40);
  lineprime=dH2("line prime","line prime",151,100,250,401,-40,40);

  both=dH2("Tracks","Tracks;Pos X (mm);Pos Z (mm)",151,100,250,401,-40,40);
  hitresidual=dH1("Hit Residual","Hit Residual; Residual (mm);counts",401,-10,10);
  hitchisq=dH1("Hit Chi Squared","Hit Chi Squared",401,-40,40);
  //strawresidual=dH1("Straw Residual","Straw Residual",401,-40,40);
  //strawchisq=dH1("Straw Chi Squared","Straw Chi Squared",401,-40,40);
  straws=dH2("Plane","Plane",51,150,200,41,-20,20);
  Angle=dH1("Angle","Angle",180,-0.5,179.5);
  chiVangle=dH2("Chi sq V Angle","Chi sq V angle;angle (degree);chisq",180,-0.5,179.5,500,0,500);
  chiVangleprime=dH2("Chi sq V Angle prime","Chi sq V angle prime;angle (degree);chisq",180,-0.5,179.5,500,0,500);
  
  return 0;
  
};

*/


Long_t covfefe::process()
{
  int ret=0;
  const TVector2 off(80,20);



  hough->Reset();
  //circ->Reset();
  straws->Reset();
  line->Reset();
  both->Reset();
  //Angle->Reset();
  chiVangle->Reset();
  // double plane0[33]={0};
  // double plane1[33]={0};
  // double plane2[33]={0};
  // double plane3[33]={0};
  // double plane4[33]={0};
  if (STT->hits.size()<3) return ok;
  //residual->Reset();
  // for(auto hit:STT->hits)
  // {
  //   int side, plane, straw_in_plane;
  //   STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
  //   if (side==1 && plane <5&& STT->hits.size()>3 && hit.dist>0.6 && plane!=0)
  //     numhits++;

  //   if(plane==0)
  //     plane0[straw_in_plane]=1;
  //   if(plane==1)
  //     plane1[straw_in_plane]=1;
  //   if(plane==2)
  //     plane2[straw_in_plane]=1;
  //   if(plane==3)
  //     plane3[straw_in_plane]=1;
  //   if(plane==4)
  //     plane4[straw_in_plane]=1;

  // }

  // for(int i=0;i<33;i++)
  // {
  //   for(int j=0;j<33;j++)
  //   {
  //     if(plane0[i]&&plane1[j])
  //       H2(i,j,"Plane0-1","Plane0-1",33,-0.5,32.5,33,-0.5,32.5);
  //     if(plane1[i]&&plane2[j])
  //       H2(i,j,"Plane1-2","Plane1-2",33,-0.5,32.5,33,-0.5,32.5);
  //     if(plane2[i]&&plane3[j])
  //       H2(i,j,"Plane2-3","Plane2-3",33,-0.5,32.5,33,-0.5,32.5);
  //     if(plane3[i]&&plane4[j])
  //       H2(i,j,"Plane3-4","Plane3-4",33,-0.5,32.5,33,-0.5,32.5);

  //}
  //}

  double numhits =0;
  for (auto hit:STT->hits)
  {
      //ret=Plugin::stop;
      int side, plane, straw_in_plane;
      STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
      if (side==1 && plane <5&& STT->hits.size()>3 && hit.dist>0.6 && hit.dist<4.5)
	    {
        numhits++;
	      auto pos=getStrawPos(hit.id)-off;
	      straws->Fill(pos.X()+hit.dist,pos.Y());
        straws->Fill(pos.X()-hit.dist,pos.Y());
        straws->Fill(pos.X(),pos.Y()-hit.dist);
        straws->Fill(pos.X(),pos.Y()+hit.dist);
//        straws->Fill(pos.X(),pos.Y());

	      //std::cout<<"Plane: "<<plane << " Straw: " << straw_in_plane;
        //std::cout <<" hit dist: "<<hit.dist<<std::endl;

	      for (int i =0; i<180;i++)// probably should do -90 to 90 but that means things have to be rewritten more intelligently for the negative angles.
	      {

	       TVector2 d=pos+TVector2(hit.dist,0).Rotate(2*i*M_PI/180);
	       circ->Fill(d.X(),d.Y());
         both->Fill(d.X(),d.Y());
	      
	       double x=pos.Rotate(i*M_PI/180).X();
         for(int j =0;j<10;j++)
         {
          double r = gRandom->Gaus(0,0.5);//random smearing. This makes regions with a high density of Hough points stand out in the search for best Hough points
          double r1 = gRandom->Gaus(0,0.7);
	        hough->Fill(i+r,(x+hit.dist)+r1);
	        if (trunc(x+hit.dist+r1+0.5)!=trunc(x-hit.dist-r1+0.5)) 
          {
            hough->Fill(i+r,x-hit.dist-r1);
          }
         }
	      } 
	    }
  }
  if (numhits<3)
    return ok;

  TF1 *f1 = new TF1("f1","[0]*x+[1]",100,200);
  f1->SetParameter(0,1);
  f1->SetParameter(1,1);
  straws->Fit("f1","q","",100,200);

  double strawangle = 90-TMath::ATan(f1->GetParameter(0))*180/M_PI;//90 because we measure from 0 differently then ROOT does, for us 0 is +y for ROOT it is +x
  int x,y,z;
  double res = 100000;
  double chisq = 100000;
  double dx=100000;
  double dy = 100000;
  int counter =0;
  hough->GetXaxis()->SetRange(0,180);
  std::vector<double> *chi = new std::vector<double>;
  std::vector<double> *par0 = new std::vector<double>;
  std::vector<double> *par1 = new std::vector<double>;
  std::vector<double> *maxX = new std::vector<double>;
  std::vector<double> *maxY = new std::vector<double>;
  std::vector<double> *angledx = new std::vector<double>;
  std::vector<double> *rho = new std::vector<double>;
  TF1 *f = new TF1("f", "[0]*x+[1]",100,250);
  while(chisq>2 && counter < 100)
  {
    chisq=0;
    //line->Reset();
    hough->GetXaxis()->SetRange(strawangle-40,strawangle+40);
      hough->GetMaximumBin(x,y,z);
      dx= hough->GetXaxis()->GetBinCenter(x);
      dy= hough->GetYaxis()->GetBinCenter(y);
    hough->SetBinContent(x,y,0);
    maxX->push_back(x);
    maxY->push_back(y);
    //std::cout << x <<" "<<y<<" "<<z<<std::endl;
    for (int i=-200;i<200;i++)
    {
      auto p=TVector2(dy,i).Rotate(-x*M_PI/180);
      line->Fill(p.X(),p.Y());
    }

    f->SetParameter(0,TMath::Cos(x*M_PI/180)/TMath::Sin(x*M_PI/180));//I don't remember Cotangent in ROOT and am too lazy to guess
    f->SetParameter(1,dy/(TMath::Sin(-x*M_PI/180)));
    //std::cout << " Par 0: " << f->GetParameter(0) << " Par 1: " << f->GetParameter(1) << std::endl;
    //std::cout << " Hough theta: " << x << " Hough R: " << dy << std::endl;
    //line->Fit("f","","MUL",100,250);


  //derp don't be an idiot, do it this way.
    //finding residuals and chi sq from fit
    for (auto hit:STT->hits)
    {
        //ret=Plugin::stop;
        int side, plane, straw_in_plane;
        STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
        if (side==1 && plane <5 && STT->hits.size()>3 && hit.dist>0.6 && hit.dist<4.5)
        {
          auto pos=getStrawPos(hit.id)-off;//this is wire position x',y'
          //now need intercept of line perpendicular to fit that passes through the straw
          //the x coordinate of intercept is given by
          double slope = f->GetParameter(0);
          double intercept = f->GetParameter(1);
          double lineX = (pos.X()/slope+pos.Y()- intercept)/(slope+1/slope);//this is given by solving equation for two perpendicular intersecting lines 
          double lineY = f->Eval(lineX);
          res = sqrt(pow(pos.X()-lineX,2)+pow(pos.Y()-lineY,2))-hit.dist;
          chisq = chisq + pow(res,2)/pow(0.15,2)/(numhits-2);
          // std::cout << " lineX: " << lineX << " lineY: " << lineY<< " WireX: " << pos.X() << " WireY: " << pos.Y() << " hit.dist " << hit.dist;
          // std::cout<<" Residual: " << res << std::endl;
        }
    }
    chiVangle->Fill(dx,chisq);
    chi->push_back(chisq);
    par0->push_back(f->GetParameter(0));
    par1->push_back(f->GetParameter(1));
    angledx->push_back(x);
    rho->push_back(dy);
    //std::cout <<" optimizing chi sq: " << chisq << std::endl;
    counter++;
  }


  double min = 1000000;
  int loc;
  for(int k=0;k<chi->size();k++)
  {
    if(min > chi->at(k))
    {
      min=chi->at(k);
      loc=k;
    }
  }




  // for(int i=-200;i<250;i++)
  // {
  //   //auto p=TVector2(dy,i).Rotate(-x*M_PI/180);
  //   //both->Fill(p.X(),p.Y());
  //   both->Fill(i,par0->at(loc)*i+par1->at(loc));
  // }


//Do it again but use the starting point of fit as the optimal value found above.
//Hough Transform isn't perfect for us (only accurate to level of bin size) so have to kind of look around optimal guessed value.

  std::vector<double> *chiprime = new std::vector<double>;
  std::vector<double> *par0prime = new std::vector<double>;
  std::vector<double> *par1prime = new std::vector<double>;
  int counter1=0;
  double chisqprime=100000;
  double resprime=100000;
  TF1 *f2 = new TF1("f2", "[0]*x+[1]",100,250);
  dx= hough->GetXaxis()->GetBinCenter(maxX->at(loc));
  dy= hough->GetYaxis()->GetBinCenter(maxY->at(loc));
  double gran = 0.1;//granularity of search around best hough params
  chiVangleprime->Reset();
  lineprime->Reset();
  double param0 = TMath::Cos(  (angledx->at(loc)) *M_PI/180)/TMath::Sin(  (angledx->at(loc)) *M_PI/180);
  double param1 = (rho->at(loc))/(TMath::Sin(- (angledx->at(loc))*M_PI/180));
  for(int i = -100;i<100;i++)
  {
    for(int j=-100;j<100;j++)
    {
      param0 = TMath::Cos(  (angledx->at(loc) +gran*i) *M_PI/180)/TMath::Sin(  (angledx->at(loc)+gran*i) *M_PI/180);
      param1 = (rho->at(loc)+gran*j)/(TMath::Sin(- (angledx->at(loc)+gran*i)*M_PI/180));
      double slope = param0;
      double intercept = param1;
      f2->SetParameter(0,param0);
      f2->SetParameter(1,param1);
      for (int k=0;k<400;k++)
      {
        lineprime->Fill(k,f2->Eval(k));
      }
      // lineprime->Fit("f2","q","MUL",100,250);
      chisqprime = 0;
      for (auto hit:STT->hits)
      {
          //ret=Plugin::stop;
          int side, plane, straw_in_plane;
          STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
          if (side==1 && plane <5 && STT->hits.size()>3 && hit.dist>0.6 && hit.dist<4.5)
          {
            auto pos=getStrawPos(hit.id)-off;//this is wire position x',y'
            //now need intercept of line perpendicular to fit that passes through the straw
            //the x coordinate of intercept is given by
            double lineX = (pos.X()/slope+pos.Y()- intercept)/(slope+1/slope);//this is given by solving equation for two perpendicular intersecting lines 
            double lineY = slope*lineX+intercept;
            resprime = sqrt(pow(pos.X()-lineX,2)+pow(pos.Y()-lineY,2))-hit.dist;
            chisqprime = chisqprime + pow(resprime,2)/pow(0.15,2)/(numhits-2);
            // std::cout << " lineX: " << lineX << " lineY: " << lineY<< " WireX: " << pos.X() << " WireY: " << pos.Y() << " hit.dist " << hit.dist;
            // std::cout<<" Residual: " << resprime << std::endl;
          }
      }

      //std::cout <<" optimizing chi sq: " << chisqprime << std::endl;
      chiprime->push_back(chisqprime);
      par0prime->push_back(slope);
      par1prime->push_back(intercept);
      chiVangleprime->Fill(angledx->at(loc)+gran*i,chisqprime);

    }
  }

  // //The minimization by TMinuit
  // Int_t iflag =0;
  // TMinuit *gMinuit = new TMinuit(2);//2 param fit
  // gMinuit->SetFCN(mychi2);
  // Double_t arglist[10];
  // arglist[0] = 1;
  // gMinuit->mnexcm("SET ERR",arglist,1,iflag);

  // gMinuit->mnparm(0,"Slope",param0,0.1,-1000,1000,iflag);//par num, par name, star val, step size,min val, max val, errflag
  // gMinuit->mnparm(1,"Intercept",param1,0.1,-1000,1000,iflag);

  // gMinuit->mnexcm("CALL FCN", arglist,1,iflag);

  // gMinuit->mnexcm("MIGRAD",arglist,2,iflag);

      

  double minprime = 1000000000;
  int locprime;
  for(int k=0;k<chiprime->size();k++)
  {
    if(minprime>chiprime->at(k))
    {
      minprime=chiprime->at(k);
      locprime=k;
    }
  }
  double finalangle = 90-TMath::ATan(par0prime->at(locprime))*180/M_PI;
  Angle->Fill(finalangle);
  // if (finalangle < 30)
  //   ret = Plugin::stop;


  //std::cout<<"new min chisq is: " << minprime << " old was: " << min<<std::endl;
  //plot the track
  for (int i=0;i<400;i++)
  {
    both->Fill(i,par0prime->at(locprime)*i+par1prime->at(locprime));
  }
  Tracks->clear();
  if(chiprime->at(locprime) < 10)
  {
    chisqprime=0;
    //do it one last time to actually fill plots with the smallest chisq and its params
			int kk = 0;
      for (auto hit:STT->hits)
      {
				kk++;
				std::cout << "\n\nLOOK HERE:" << kk << std::endl;
        //ret=Plugin::stop;
        int side, plane, straw_in_plane;
        STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
        if (side==1 && plane <5 && STT->hits.size()>3 && hit.dist>0.6 && hit.dist<4.5)
        {
          auto pos=getStrawPos(hit.id)-off;//this is wire position x',y'
          //now need intercept of line perpendicular to fit that passes through the straw
          //the x coordinate of intercept is given by
          double slope = par0prime->at(locprime);
          double intercept = par1prime->at(locprime);
          double lineX = (pos.X()/slope+pos.Y()- intercept)/(slope+1/slope);//this is given by solving equation for two perpendicular intersecting lines 
          double lineY = slope*lineX+intercept;
          resprime = sqrt(pow(pos.X()-lineX,2)+pow(pos.Y()-lineY,2))-hit.dist;
          chisqprime = chisqprime + pow(resprime,2)/pow(0.15,2)/(numhits-2);
          hitresidual->Fill(resprime);  
        }
    
      }
    TrackHit outtrack;
    outtrack.chisq=chisqprime;
    outtrack.slope=par0prime->at(locprime);
    outtrack.intercept=par1prime->at(locprime);
    Tracks->tracks.emplace_back(std::move(outtrack));
    hitchisq->Fill(chisqprime); 
  }
}

Long_t covfefe::cmdline(char *cmd)
{

  return 0; // 0 = all ok
};

Long_t covfefe::finalize()
{

  return 0; // 0 = all ok
};


extern "C"{
  Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
  {
    return (Plugin *) new covfefe(in,out,inf_,outf_,p);
  }
}
