#include <slowctrl.h>
#include<iostream>


void slowctrl::manager::add(std::string name,unsigned short id)
{
  // check if it's really new:
  if (namemap.find(name)!=namemap.end())
    return;
  current.push_back(new datum());
  lastValid.push_back(new datum());
  namemap[name]=current.size()-1;
  idmap[id]=current.size()-1;

}


slowctrl::datum* slowctrl::manager::getCurrentByName(std::string name)
{
  //std::cout<<"Look up:"<<name<<std::endl;
  return current[namemap[name]];
}

slowctrl::datum* slowctrl::manager::getCurrentByID(unsigned short id)
{
  return current[idmap[id]];
}


slowctrl::datum* slowctrl::manager::getLastValidByName(std::string name)
{
  return lastValid[namemap[name]];
}

slowctrl::datum* slowctrl::manager::getLastValidByID(unsigned short id)
{
  return lastValid[idmap[id]];
}

void slowctrl::manager::addWatch(slowctrl::datum *d,double low, double high, int *inv)
{
  slowctrl::SCwatch w;
  w.low=low;
  w.high=high;
  w.invalid=inv;
  w.state= ((d->value < low) ||(d->value > high));
  if (w.state) (*inv)++;
  d->watches.push_back(w);
}

void slowctrl::manager::updateWatches(slowctrl::datum *d)
{
  for (std::vector<SCwatch >::iterator iter=d->watches.begin();iter!=d->watches.end();++iter)
  {
    bool newstate=(d->value<iter->low || d->value>iter->high);
    if (newstate && !iter->state)  (*iter->invalid)++;
    if (!newstate && iter->state)  (*iter->invalid)--;
    iter->state=newstate;
   }
}
