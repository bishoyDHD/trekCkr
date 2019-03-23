#ifndef clus_var_h
#define clus_var_h 1
#include <stdio.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>
#include <fstream>
#include <execinfo.h>
#include <map>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
using namespace std;
int clusCrys;
string NameFilelist;
UInt_t RunNo=0;
UInt_t MaxEvent=0;
string NameRoot="vf48_dump.root";
vector<string> ListNameFile;
boost::unordered_map<std::pair<double,double>,double, boost::hash<pair<double,double>>> csiph;
//std::unordered_map<std::pair<double,double>,double, boost::hash<pair<double,double>>> csiph;
boost::unordered_map<std::pair<double,double>,bool, boost::hash<pair<double,double>>> csiClus;
typedef vector<double> ve;
extern ve indexph;
ve indexph;
vector<double> *phval, *clusth, *clusphi;
std::size_t get_nthIndex(ve, std::size_t k){
  std::vector<std::size_t> indexes(indexph.size());
  std::iota(indexes.begin(), indexes.end(), 0);

  std::nth_element(indexes.begin(), indexes.begin() + k, indexes.end(),
    [&](int lhs, int rhs){
      return indexph[lhs] > indexph[rhs];
    }
  );
  return indexes[k];
}
//arranged theta[fb][crystalNo. 0-15]: f=0,b=1
double theta[2][16]={86.25, 78.75, 71.25, 63.75, 56.25, 48.75, 41.25, 33.75, 26.25,
                                          63.75, 56.25, 48.75, 41.25, 33.75, 26.25, 18.75,
                     93.75, 101.25, 108.75, 116.25, 123.75, 131.25, 138.75, 146.25, 153.75,
		                            116.25, 123.75, 131.25, 138.75, 146.25, 153.75, 161.25};
double phi[12][2][2]={{{3.75, 11.25},{18.75, 26.25}},
                      {{33.75,41.25},{48.75, 56.25}},
		      {{63.75,71.25},{78.75, 86.25}},
		      {{93.75,101.25},{108.75,116.25}},
		      {{123.75,131.25},{138.75,146.25}},
		      {{153.75,161.25},{168.75,176.25}},
		      {{183.75,191.25},{198.75,206.25}},
		      {{213.75,221.25},{228.75,236.25}},
		      {{243.75,251.25},{258.75,266.25}},
		      {{273.75,281.25},{288.75,296.25}},
		      {{303.75,311.25},{318.75,326.25}},
		      {{333.75,341.25},{348.75,356.25}}};

#endif
