
#ifndef COMPAREMAPS_H
#define COMPAREMAPS_H

using namespace std;

struct cmap;
struct randoms;


void compare_maps(int,char*[],randoms&);
void shuffle_maps(int,char*[],randoms&);

double get_chi2(const cmap&, const cmap&,const int&);

cmap get_shuffled_cmap(const int&,cmap*[],randoms&); 

#endif
