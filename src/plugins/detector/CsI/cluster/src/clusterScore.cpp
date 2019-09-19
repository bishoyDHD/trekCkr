#include "clusterScore.h"

//clusterVar::clusterVar(){}

clusterScore::clusterScore(){
 
}
clusterScore::~clusterScore(){

}
// called for each valid clustering event
void clusterScore::init(){
  //clustvar=new clusterVar();
  cvars.clear();
}
// will be needed for sanity check plots and scoring
// default is pi0
void clusterScore::setScoreMass(double particleMass=0.1349766){
  mass=particleMass;
}
void clusterScore::reset(){
  cvars.clear();
}
// combination clusters
// evaluate cluster variables for different channels
void clusterScore::clusterEval(const std::vector<double> &mCrys,const std::vector<double> &sCrys,const std::vector<double> r, const std::vector<double> z,const std::vector<double> &theta,const std::vector<double> &phi){

  
}
// many crystal cluster
void clusterScore::clusterEval(const std::vector<double> &eneCrys,const std::vector<double> Wr,const std::vector<double> Wz,const std::vector<double> &theta,const std::vector<double> &phi){
  // this is the case for single and multiple cluster;
    // =======================================
    // Case for either n many crysCluster or
    // n single crysCluster evaluation
    // ========================================
    energy=(eneCrys[0]+eneCrys[1]);
    std::cout<<" ... new energy eval.: "<<energy<<std::endl;
    for(UInt_t i=0; i<eneCrys.size(); i++){
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
    //std::cout<<" ..... size of cvars: "<<cvars.size()<<std::endl;
    //std::cout<<" ..... new x eval: "<<cvars[0].clx<<"\t"<<cvars[1].clx<<std::endl;
    //std::cout<<" ..... new y eval: "<<cvars[0].cly<<"\t"<<cvars[1].cly<<std::endl;
    //std::cout<<" ..... new z eval: "<<cvars[0].clz<<"\t"<<cvars[1].clz<<std::endl;
    //std::cout<<" ..... new Inv. mass eval: "<<cvars[0].pi0lv.M()<<std::endl;
    //std::abort();
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
