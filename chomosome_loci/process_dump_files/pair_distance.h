
#ifndef PDIST_H
#define PDIST_H

#include <vector>

using namespace std;

struct simluationsystem;
struct atom;

void pairwise_dist(int, char *[],const double&);

double get_distance(const simluationsystem&, const atom[], const atom[], const double&);
double get_distance(const simluationsystem&, const atom[], const atom[], const double&, const vector<int> &);

#endif

