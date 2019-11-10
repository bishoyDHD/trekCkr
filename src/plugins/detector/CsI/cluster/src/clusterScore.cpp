#include "clusterScore.h"

//clusterVar::clusterVar(){}

clusterScore::clusterScore():/*mass(0.1349766),*/cpid(1),clustEvalNo(2){
 
}
clusterScore::~clusterScore(){

}
// called at the begining of each valid clustering event
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
void clusterScore::setScoreMass(double particleMass=0.0){ //0.1349766){
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
//   ------ This function has the sole purpose of appending 2 std::vector's
void clusterScore::mergeVect(std::vector<double>& clust1,std::vector<double>& clust2){
  clust1.insert(
         clust1.end(),
         std::make_move_iterator(clust2.begin()),
         std::make_move_iterator(clust2.end())
         );
};
void clusterScore::clusterEval(std::vector<double> &mCrys,std::vector<double> &sCrys,std::vector<double> Wr,std::vector<double> Wz,std::vector<double> &theta,std::vector<double> &phi,std::vector<double> &r2,std::vector<double> &z2,std::vector<double> &theta2,std::vector<double> &phi2){
  // this is the case for single and many crysCluster;
  // =======================================
  // Case for mixed clusters
  // ========================================
  // manyCrysE and singleCrysE        // manyCrysR and singleCrysR
  mergeVect(mCrys,sCrys);             mergeVect(Wr,r2);
  // manyCrysZ and singleCrysZ        // manyCrysTheta and singleCrysTheta
  mergeVect(Wz,z2);                   mergeVect(theta,theta2);
  // manyCrysPhi and singleCrysPhi
  mergeVect(phi,phi2);
  // set Cluster multiplicity
  setClustM(mCrys.size());
  // push cluster variables to container
  for(UInt_t i=0; i<mCrys.size(); i++){
    // cluster storage:
    px=mCrys[i]*std::sin(theta[i])*std::cos(phi[i]);
    py=mCrys[i]*std::sin(theta[i])*std::sin(phi[i]);
    pz=mCrys[i]*std::cos(theta[i]);
    x=Wr[i]*std::cos(phi[i]);
    y=Wr[i]*std::sin(phi[i]);
    z=Wz[i];
    r=Wr[i];
    clustvar.clpx=px;      clustvar.clx=x;
    clustvar.clpy=py;      clustvar.cly=y;
    clustvar.clpz=pz;      clustvar.clz=z;
    clustvar.clr=r;
    // evaluate 4-momenta and corresponding vectors
    clustvar.cpidlv.SetPxPyPzE(px,py,pz,mCrys[i]);
    clustvar.cpidv3.SetXYZ(px,py,pz);
    cvars.push_back(clustvar);
  }
  // cluster evaluation and scoring
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
          energy=(mCrys[i]+mCrys[n]);
          particlelv=cvars[i].cpidlv+cvars[n].cpidlv;
          opAngle=std::cos(cvars[i].cpidv3.Angle(cvars[n].cpidv3));
          invMass=particlelv.M();
          diffMass=std::abs(mass-invMass);
          // fill vars for scoring:
          mdiff.push_back(diffMass);
          InvMass[diffMass]=invMass;
          clustE[diffMass]=energy;
	  openAng[diffMass]=opAngle;
	  primLV[diffMass]=particlelv;
          csix[diffMass]=std::make_pair(cvars[i].clx,cvars[n].clx);
          csiy[diffMass]=std::make_pair(cvars[i].cly,cvars[n].cly);
          csiz[diffMass]=std::make_pair(cvars[i].clz,cvars[n].clz);
          csir[diffMass]=std::make_pair(cvars[i].clr,cvars[n].clr);
          csipx[diffMass]=std::make_pair(cvars[i].clpx,cvars[n].clpx);
          csipy[diffMass]=std::make_pair(cvars[i].clpy,cvars[n].clpy);
          csipz[diffMass]=std::make_pair(cvars[i].clpz,cvars[n].clpz);
          csiE[diffMass]=std::make_pair(mCrys[i],mCrys[n]);
          csitheta[diffMass]=std::make_pair(theta[i],theta[n]);
          csiphi[diffMass]=std::make_pair(phi[i],phi[n]);
          std::cout<<" ... Inv. Mass of pi0 is: "<<diffMass<<std::endl;
	  std::cout<<" ... primpidLV().Px() "<<primLV[diffMass].Px()<<std::endl;
	  std::cout<<" ... opening Ang: "<<openAng[diffMass]<<" - "<<opAngle<<std::endl;
          std::cout<<" ... new energy eval.: "<<clustE[diffMass]<<std::endl;
          std::cout<<" ..... new x eval: "<<csipx[diffMass].first<<"\t"<<csipx[diffMass].second<<std::endl;
          std::cout<<" ..... new y eval: "<<csipy[diffMass].first<<"\t"<<csipy[diffMass].second<<std::endl;
          std::cout<<" ..... new z eval: "<<csipz[diffMass].first<<"\t"<<csipz[diffMass].second<<std::endl;
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
            energy=(mCrys[i]+mCrys[n]+mCrys[a]);
            particlelv=cvars[i].cpidlv+cvars[n].cpidlv+cvars[a].cpidlv;
            invMass=particlelv.M();
            diffMass=std::abs(mass-invMass);
            // fill vars for scoring:
            mdiff.push_back(diffMass);
            InvMass[diffMass]=invMass;
            clustE[diffMass]=energy;
	    primLV[diffMass]=particlelv;
	    // fill scoring vars for first 2 hits
            csix[diffMass]=std::make_pair(cvars[i].clx,cvars[n].clx);
            csiy[diffMass]=std::make_pair(cvars[i].cly,cvars[n].cly);
            csiz[diffMass]=std::make_pair(cvars[i].clz,cvars[n].clz);
            csir[diffMass]=std::make_pair(cvars[i].clr,cvars[n].clr);
            csipx[diffMass]=std::make_pair(cvars[i].clpx,cvars[n].clpx);
            csipy[diffMass]=std::make_pair(cvars[i].clpy,cvars[n].clpy);
            csipz[diffMass]=std::make_pair(cvars[i].clpz,cvars[n].clpz);
            csiE[diffMass]=std::make_pair(mCrys[i],mCrys[n]);
            csitheta[diffMass]=std::make_pair(theta[i],theta[n]);
            csiphi[diffMass]=std::make_pair(phi[i],phi[n]);
	    // fill scoring vars for 3rd cluster
            csix_[diffMass]=std::make_pair(cvars[a].clx,dummy);
            csiy_[diffMass]=std::make_pair(cvars[a].cly,dummy);
            csiz_[diffMass]=std::make_pair(cvars[a].clz,dummy);
            csir_[diffMass]=std::make_pair(cvars[a].clr,dummy);
            csipx_[diffMass]=std::make_pair(cvars[a].clpx,dummy);
            csipy_[diffMass]=std::make_pair(cvars[a].clpy,dummy);
            csipz_[diffMass]=std::make_pair(cvars[a].clpz,dummy);
            csiE_[diffMass]=std::make_pair(mCrys[a],dummy);
            csitheta_[diffMass]=std::make_pair(theta[a],dummy);
            csiphi_[diffMass]=std::make_pair(phi[a],dummy);
            std::cout<<" ... Inv. Mass of pi0 is: "<<diffMass<<std::endl;
            std::cout<<" ..... new x eval: "<<csipx[diffMass].first<<"\t"<<csipx[diffMass].second<<std::endl;
            std::cout<<" ..... new y eval: "<<csipy[diffMass].first<<"\t"<<csipy[diffMass].second<<std::endl;
            std::cout<<" ..... new z eval: "<<csipz[diffMass].first<<"\t"<<csipz[diffMass].second<<std::endl;
          }
	}
      }
      break;
  } // end of switch statement
  scoring(mdiff);
  std::cout<<" ... Merged vector size: "<<mCrys.size()<<std::endl;
  std::cout<<" ... Merged ang/rz size: "<<theta.size()<<" "<<phi.size()<<" "<<Wr.size()<<Wz.size()<<std::endl;
  std::cout<<" ... new energy eval.: "<<clustE[diffMass]<<std::endl;
  //std::cout<<" ..... size of cvars: "<<cvars.size()<<std::endl;
  //std::cout<<" ..... new x eval: "<<csipx[diffMass].first<<"\t"<<csipx[diffMass].second<<std::endl;
  //std::cout<<" ..... new y eval: "<<csipy[diffMass].first<<"\t"<<csipy[diffMass].second<<std::endl;
  //std::cout<<" ..... new z eval: "<<csipz[diffMass].first<<"\t"<<csipz[diffMass].second<<std::endl;
  //std::cout<<" ..... new Inv. mass eval: "<<cvars[0].pi0lv.M()<<std::endl;
  std::cout<<" ..... new Inv. mass scored return: "<<InvMass[diffMass]<<std::endl;
  //std::abort();
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
  // set Cluster multiplicity
  setClustM(eneCrys.size());
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
    clustvar.cpidv3.SetXYZ(px,py,pz);
    cvars.push_back(clustvar);
  }
  // cluster evaluation prodcedure
  switch(clustEvalNo){
    case 2:
      // =========================================
      // Calculate the combinatorial sum assuming
      // 2 clusters. 
      // ==========================================
      std::cout<<"  checking 2 clusters: 2 clusters: 2 clusters \n";
      for(UInt_t i=0; i<cvars.size()-1; i++){
        UInt_t m=i+1;
        for(UInt_t n=m; n<cvars.size(); n++){
          std::cout<<" we have: "<<i<<" + "<<n<<"\n";
          energy=(eneCrys[i]+eneCrys[n]);
	  std::cout<<"+++ Energy ++ Energy ++ Energy: "<<energy<<std::endl;
          particlelv=cvars[i].cpidlv+cvars[n].cpidlv;
          opAngle=std::cos(cvars[i].cpidv3.Angle(cvars[n].cpidv3));
          invMass=particlelv.M();
          diffMass=std::abs(mass-invMass);
          // fill vars for scoring:
          mdiff.push_back(diffMass);
          InvMass[diffMass]=invMass;
          clustE[diffMass]=energy;
	  openAng[diffMass]=opAngle;
	  primLV[diffMass]=particlelv;
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
	    primLV[diffMass]=particlelv;
	    // fill scoring vars for first 2 hits
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
	    // fill scoring vars for 3rd cluster
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
  scoring(mdiff);
  std::cout<<" ... new energy eval.: "<<clustE[diffMass]<<std::endl;
  std::cout<<" ... new angle eval. : "<<openAng[mkey]<<std::endl;
  //std::cout<<" ..... size of cvars: "<<cvars.size()<<std::endl;
  //std::cout<<" ..... new x eval: "<<csipx[diffMass].first<<"\t"<<csipx[diffMass].second<<std::endl;
  //std::cout<<" ..... new y eval: "<<csipy[diffMass].first<<"\t"<<csipy[diffMass].second<<std::endl;
  //std::cout<<" ..... new z eval: "<<csipz[diffMass].first<<"\t"<<csipz[diffMass].second<<std::endl;
  //std::cout<<" ..... new Inv. mass eval: "<<cvars[0].pi0lv.M()<<std::endl;
  std::cout<<" ..... new Inv. mass scored return: "<<InvMass[diffMass]<<std::endl;
  //std::abort();
}
// cluster scoring
void clusterScore::scoring(std::vector<double> &invmass){
  // obtain the index of the lowest entry of the invMass 
  // new size of vector to determine scoring
  std::vector<double>::iterator min=std::min_element(invmass.begin(),invmass.end());
  ival=std::distance(invmass.begin(),min);
  setKey(*min);
  std::cout<<" ---- best Invariant mass diff: "<<*min<<" corresponding to "<<InvMass[*min]<<" & "<<clustE[*min]<<"\n";
}
void clusterScore::scoring(UInt_t size, std::vector<double> &invmass,std::map<double,double> mass){
  // obtain the index of the lowest entry of the invMass 
  // new size of vector to determine scoring
  auto min=std::min_element(invmass.begin(),invmass.end());
  ival=std::distance(invmass.begin(),min);
  invMass=mass[*min];
  setKey(*min);
  std::cout<<" ---- Invariant mass is: "<<invMass<<"\n";
}
double clusterScore::getE(){
  return clustE[mkey];
}
// obtain map key corresponding to distance of closest approach
void clusterScore::setKey(double key){
  mkey=key;
}
TLorentzVector clusterScore::getprimLV(){
  return primLV[mkey];
}
// get methods for cluster variables
double clusterScore::getprE(){
  // return the total energy of recontructed particle
  return primLV[mkey].E();
}
double clusterScore::getprM(){
  // return the invariant mass of recontructed particle
  return primLV[mkey].M();
}
double clusterScore::getprPx(){
  // return the px of recontructed particle
  return primLV[mkey].Px();
}
double clusterScore::getprPy(){
  // return the py of recontructed particle
  return primLV[mkey].Py();
}
double clusterScore::getprPz(){
  // return the pz of recontructed particle
  return primLV[mkey].Pz();
}
double clusterScore::getInvMass(){
  return InvMass[mkey];
}
double clusterScore::getOpAngleClust(){
  return openAng[mkey];
}
double clusterScore::getclPx(){
  switch(cpid){
    case 1:
      return csipx[mkey].first;
      break;
    case 2:
      return csipx[mkey].second;
      break;
    case 3:
      return csipx_[mkey].first;
      break;
    case 4:
      return csipx_[mkey].second;
      break;
  }
  return 0; // non void function
}
double clusterScore::getclPy(){
  switch(cpid){
    case 1:
      return csipy[mkey].first;
      break;
    case 2:
      return csipy[mkey].second;
      break;
    case 3:
      return csipy_[mkey].first;
      break;
    case 4:
      return csipy_[mkey].second;
      break;
  }
  return 0; // non void function
}
double clusterScore::getclPz(){
  switch(cpid){
    case 1:
      return csipz[mkey].first;
      break;
    case 2:
      return csipz[mkey].second;
      break;
    case 3:
      return csipz_[mkey].first;
      break;
    case 4:
      return csipz_[mkey].second;
      break;
  }
  return 0; // non void function
}
double clusterScore::getclR(){
  switch(cpid){
    case 1:
      return csir[mkey].first;
      break;
    case 2:
      return csir[mkey].second;
      break;
    case 3:
      return csir_[mkey].first;
      break;
    case 4:
      return csir_[mkey].second;
      break;
  }
  return 0; // non void function
}
double clusterScore::getclZ(){
  switch(cpid){
    case 1:
      return csiz[mkey].first;
      break;
    case 2:
      return csiz[mkey].second;
      break;
    case 3:
      return csiz_[mkey].first;
      break;
    case 4:
      return csiz_[mkey].second;
      break;
  }
  return 0; // non void function
}
double clusterScore::getclTime(){
  switch(cpid){
    case 1:
      return csitime[mkey].first;
      break;
    case 2:
      return csitime[mkey].second;
      break;
    case 3:
      return csitime_[mkey].first;
      break;
    case 4:
      return csitime_[mkey].second;
      break;
  }
  return 0; // non void function
}
double clusterScore::getclE(){
  switch(cpid){
    case 1:
      return csiE[mkey].first;
      break;
    case 2:
      return csiE[mkey].second;
      break;
    case 3:
      return csiE_[mkey].first;
      break;
    case 4:
      return csiE_[mkey].second;
      break;
  }
  return 0; // non void function
}
double clusterScore::getclTheta(){
  switch(cpid){
    case 1:
      return csitheta[mkey].first;
      break;
    case 2:
      return csitheta[mkey].second;
      break;
    case 3:
      return csitheta_[mkey].first;
      break;
    case 4:
      return csitheta_[mkey].second;
      break;
  }
  return 0; // non void function
}
double clusterScore::getclPhi(){
  switch(cpid){
    case 1:
      return csiphi[mkey].first;
      break;
    case 2:
      return csiphi[mkey].second;
      break;
    case 3:
      return csiphi_[mkey].first;
      break;
    case 4:
      return csiphi_[mkey].second;
      break;
  }
  return 0; // non void function
}
