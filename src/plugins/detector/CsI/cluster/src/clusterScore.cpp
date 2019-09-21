#include "clusterScore.h"

//clusterVar::clusterVar(){}

clusterScore::clusterScore():mass(0.1349766),clustEvalNo(2){
 
}
clusterScore::~clusterScore(){

}
// called for each valid clustering event
// must empty bucket before filling 
// with new entries
void clusterScore::init(){
  cvars.clear();
  mdiff.clear();
  clustE.clear();
  InvMass.clear();
  csipx.clear();         csipx_.clear();      csix.clear();      csix_.clear();
  csipy.clear();         csipy_.clear();      csiy.clear();      csiy_.clear();
  csipz.clear();         csipz_.clear();      csiz.clear();      csiz_.clear();
  csiE.clear();          csitheta_.clear();   csir.clear();      csir_.clear();
  csitheta.clear();      csiphi_.clear();
  csiphi.clear();        csitime_.clear();
  csitime.clear();
}
// will be needed for sanity check plots and scoring
// default is pi0
void clusterScore::setScoreMass(double particleMass=0.1349766){
  mass=particleMass;
}
// must empty the bucket after every fill
void clusterScore::reset(){
  cvars.clear();
  mdiff.clear();
  clustE.clear();
  InvMass.clear();
  csipx.clear();         csipx_.clear();      csix.clear();      csix_.clear();
  csipy.clear();         csipy_.clear();      csiy.clear();      csiy_.clear();
  csipz.clear();         csipz_.clear();      csiz.clear();      csiz_.clear();
  csiE.clear();          csitheta_.clear();   csir.clear();      csir_.clear();
  csitheta.clear();      csiphi_.clear();
  csiphi.clear();        csitime_.clear();
  csitime.clear();
}
// combination clusters
// evaluate cluster variables for different channels
void clusterScore::clusterEval(const std::vector<double> &mCrys,const std::vector<double> &sCrys,const std::vector<double> r, const std::vector<double> z,const std::vector<double> &theta,const std::vector<double> &phi){

  
}
// single and many crystal cluster
void clusterScore::clusterEval(const std::vector<double> &eneCrys,const std::vector<double> Wr,const std::vector<double> Wz,const std::vector<double> &theta,const std::vector<double> &phi){
  // this is the case for single and many crysCluster;
  // =======================================
  // Case for either n many crysCluster or
  // n single crysCluster evaluation
  // ========================================
  //energy=(eneCrys[0]+eneCrys[1]);
  //std::cout<<" ... new energy eval.: "<<energy<<std::endl;
  for(UInt_t i=0; i<eneCrys.size(); i++){
    // cluster storage:
    px=eneCrys[i]*std::sin(theta[i])*std::cos(phi[i]);
    py=eneCrys[i]*std::sin(theta[i])*std::sin(phi[i]);
    pz=eneCrys[i]*std::cos(theta[i]);
    x=Wr[i]*std::cos(phi[i]);
    y=Wr[i]*std::sin(phi[i]);
    z=Wz[i];
    r=Wr[i];
    clustvar.clpx=px;      clustvar.clx=x;
    clustvar.clpy=py;      clustvar.cly=y;
    clustvar.clpz=pz;      clustvar.clz=z;
    clustvar.clr=r;
    // evaluate 4-momenta and corresponding vectors
    clustvar.cpidlv.SetPxPyPzE(px,py,pz,eneCrys[i]);
    cvars.push_back(clustvar);
  }
  switch(clustEvalNo){
    case 2:
      // =========================================
      // Calculate the combinatorial sum assuming
      // 2 clusters. 
      // ==========================================
      for(UInt_t i=0; i<cvars.size()-1; i++){
        UInt_t m=i+1;
        for(UInt_t n=m; n<cvars.size(); n++){
          std::cout<<" we have: "<<i<<" + "<<n<<"\n";
          energy=(eneCrys[i]+eneCrys[n]);
          particlelv=cvars[i].cpidlv+cvars[n].cpidlv;
          invMass=particlelv.M();
          diffMass=std::abs(mass-invMass);
          // fill vars for scoring:
          mdiff.push_back(diffMass);
          InvMass[diffMass]=invMass;
          clustE[diffMass]=energy;
          csix[diffMass]=std::make_pair(cvars[i].clx,cvars[n].clx);
          csiy[diffMass]=std::make_pair(cvars[i].cly,cvars[n].cly);
          csiz[diffMass]=std::make_pair(cvars[i].clz,cvars[n].clz);
          csir[diffMass]=std::make_pair(cvars[i].clr,cvars[n].clr);
          csipx[diffMass]=std::make_pair(cvars[i].clpx,cvars[n].clpx);
          csipy[diffMass]=std::make_pair(cvars[i].clpy,cvars[n].clpy);
          csipz[diffMass]=std::make_pair(cvars[i].clpz,cvars[n].clpz);
          csiE[diffMass]=std::make_pair(eneCrys[i],eneCrys[n]);
          csitheta[diffMass]=std::make_pair(theta[i],theta[n]);
          csiphi[diffMass]=std::make_pair(phi[i],phi[n]);
          std::cout<<" ... Inv. Mass of pi0 is: "<<diffMass<<std::endl;
        }
      }
      break;
    case 3:
      // case for evaluation of and scoring of 3 different cluster
      // =========================================
      // Calculate the combinatorial sum assuming
      // 3 clusters. 
      // ==========================================
      for(UInt_t i=0; i<cvars.size()-2; i++){
        UInt_t m=i+1;
        for(UInt_t n=m; n<cvars.size()-1; n++){
          UInt_t s=n+1;
          for(UInt_t a=s; a<cvars.size(); a++){
            std::cout<<" we have: "<<i<<" + "<<n<<" + "<<a<<"\n";
            energy=(eneCrys[i]+eneCrys[n]+eneCrys[a]);
            particlelv=cvars[i].cpidlv+cvars[n].cpidlv+cvars[a].cpidlv;
            invMass=particlelv.M();
            diffMass=std::abs(mass-invMass);
            // fill vars for scoring:
            mdiff.push_back(diffMass);
            InvMass[diffMass]=invMass;
            clustE[diffMass]=energy;
	    // fill scoring vars for first 2 hist
            csix[diffMass]=std::make_pair(cvars[i].clx,cvars[n].clx);
            csiy[diffMass]=std::make_pair(cvars[i].cly,cvars[n].cly);
            csiz[diffMass]=std::make_pair(cvars[i].clz,cvars[n].clz);
            csir[diffMass]=std::make_pair(cvars[i].clr,cvars[n].clr);
            csipx[diffMass]=std::make_pair(cvars[i].clpx,cvars[n].clpx);
            csipy[diffMass]=std::make_pair(cvars[i].clpy,cvars[n].clpy);
            csipz[diffMass]=std::make_pair(cvars[i].clpz,cvars[n].clpz);
            csiE[diffMass]=std::make_pair(eneCrys[i],eneCrys[n]);
            csitheta[diffMass]=std::make_pair(theta[i],theta[n]);
            csiphi[diffMass]=std::make_pair(phi[i],phi[n]);
	    // fill scoring vars for first 2 hist
            csix_[diffMass]=std::make_pair(cvars[a].clx,dummy);
            csiy_[diffMass]=std::make_pair(cvars[a].cly,dummy);
            csiz_[diffMass]=std::make_pair(cvars[a].clz,dummy);
            csir_[diffMass]=std::make_pair(cvars[a].clr,dummy);
            csipx_[diffMass]=std::make_pair(cvars[a].clpx,dummy);
            csipy_[diffMass]=std::make_pair(cvars[a].clpy,dummy);
            csipz_[diffMass]=std::make_pair(cvars[a].clpz,dummy);
            csiE_[diffMass]=std::make_pair(eneCrys[a],dummy);
            csitheta_[diffMass]=std::make_pair(theta[a],dummy);
            csiphi_[diffMass]=std::make_pair(phi[a],dummy);
            std::cout<<" ... Inv. Mass of pi0 is: "<<diffMass<<std::endl;
          }
	}
      }
      break;
  } // end of switch statement

  std::cout<<" ... new energy eval.: "<<clustE[diffMass]<<std::endl;
  //std::cout<<" ..... size of cvars: "<<cvars.size()<<std::endl;
  //std::cout<<" ..... new x eval: "<<csipx[diffMass].first<<"\t"<<csipx[diffMass].second<<std::endl;
  //std::cout<<" ..... new y eval: "<<csipy[diffMass].first<<"\t"<<csipy[diffMass].second<<std::endl;
  //std::cout<<" ..... new z eval: "<<csipz[diffMass].first<<"\t"<<csipz[diffMass].second<<std::endl;
  //std::cout<<" ..... new Inv. mass eval: "<<cvars[0].pi0lv.M()<<std::endl;
  std::cout<<" ..... new Inv. mass scored return: "<<InvMass[diffMass]<<std::endl;
  std::abort();
}
// cluster scoring
double clusterScore::scoring(UInt_t size, std::vector<double> &invmass,std::map<double,double> mass){
  // obtain the index of the lowest entry of the invMass 
  // new size of vector to determine scoring
  std::vector<double>::iterator min=std::min_element(invmass.begin(),invmass.end());
  ival=std::distance(invmass.begin(),min);
  invMass=mass[*min];
  setIndex(ival);
  std::cout<<" ---- Invariant mass is: "<<invMass<<"\n";
  return invMass;
}
void clusterScore::setE(const std::vector<double> &totE){
  Etot=totE[index];
}

void clusterScore::setIndex(int I){
  index=I;
}
void clusterScore::setOpangClus(const std::vector<double> &theta,const std::vector<double> &phi){

}
void clusterScore::setOpangPrimary(const std::vector<double> &theta,const std::vector<double> &phi){

}
void clusterScore::setprPx(std::vector<clusterVar> &px){

}
void clusterScore::setprPy(std::vector<clusterVar> &py){

}
void clusterScore::setprPz(std::vector<clusterVar> &pz){

}
void clusterScore::setprE(std::vector<clusterVar> &Energy){

}
void clusterScore::setclPx(std::vector<clusterVar> &px){

}
void clusterScore::setclPy(std::vector<clusterVar> &py){

}
void clusterScore::setclPz(std::vector<clusterVar> &pz){

}
void clusterScore::setclE(std::vector<clusterVar> &Energy){

}
