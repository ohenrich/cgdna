
#ifndef COUNTBOUNT_H
#define COUNTBOUNT_H

#include <set>

using namespace std;

struct simluationsystem;

void count_bound(int,char*[],const double&);

int get_bound(const simluationsystem &,ifstream &,const set<int> &,const double&);


#endif
