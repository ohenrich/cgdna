
#ifndef RG_H
#define RG_H

using namespace std;

struct simluationsystem;

void gyration_radius_t(int, char*[]);

void gyration_radius_r(int, char*[]);

double get_Rg(const simluationsystem&, ifstream&, const int&, const int&);

#endif
