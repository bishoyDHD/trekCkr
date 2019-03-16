#ifndef __SLOWCTRL_H__
#define __SLOWCTRL_H__
#include <map>
#include <vector>
#include <string>

namespace slowctrl{

class SCwatch
{
 public:
  double low,high;
  bool state;
  int *invalid;
};

 class datum{
 private:
   std::vector<slowctrl::SCwatch> watches;
 public:
   
   unsigned long long timestamp;
   double value;
   unsigned short status;
   friend class manager;
  };


class manager
{
 private:
  std::vector<slowctrl::datum*> current;
  std::vector<slowctrl::datum*> lastValid;
  std::map<std::string,int> namemap;
  std::map<unsigned short ,int> idmap;
 public:
  datum* getCurrentByName(std::string name);
  datum* getCurrentByID(unsigned short id);
  datum* getLastValidByName(std::string name);
  datum* getLastValidByID(unsigned short id);
  bool datumExists(std::string name);
  void add(std::string,unsigned short id);
  void addWatch(datum*,double low, double high, int *inv);
  void updateWatches(datum *);
};

}
#endif
