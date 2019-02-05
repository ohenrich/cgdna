// Program to change the linking number of a DNA loop
// Reads a lammps input file which has been generate by 
// restart to data, and outputs a lammps input file

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>

#include "tools.h"
#include "tools_functions.cc"

#define PI 3.14159265358979

using namespace std;

int main(int argv, char* argc[]) {

  if (argv!=4) {cout<<"Usage : a.out input.file output.file Lk"<<endl<<" where Lk is the number of turns to add"<<endl;exit(0);}

  char fn[50];
  ifstream inf;
  ofstream ouf;
  string line;
  stringstream sline;
  string first,second;

  int N,                  // number of atoms
    Ne,                   // number of ellipsoids
    Lk;                   // number of turns to add

  double lx,ly,lz,        // box size
    d_theta,              // the increment in the angle to add
    theta;                // angle

  xyz axis;               // an axis vector

  atom *atoms;              // an array of atoms

  inf.open(argc[1]);

  do {  
    string third;
    double hi,lo;
    getline(inf,line); sline.str(line); sline>>first>>second>>third; sline.clear();
    if (second.compare("atoms")==0) {sline.str(line); sline>>N; sline.clear();}
    if (second.compare("ellipsoids")==0) {sline.str(line); sline>>Ne; sline.clear();}
    if (third.compare("xlo")==0) {sline.str(line); sline>>lo>>hi; lx=hi-lo; sline.clear();}
    if (third.compare("ylo")==0) {sline.str(line); sline>>lo>>hi; ly=hi-lo; sline.clear();}
    if (third.compare("zlo")==0) {sline.str(line); sline>>lo>>hi; lz=hi-lo; sline.clear();}
  } while ( first.compare("Atoms")!=0 );
  getline(inf,line); // a blank line
 
  if (Ne!=N) {
    cout<<"Warning - number of ellipsoids is not the same as number of atoms."<<endl; 
    cout<<"This program assumes that the DNA loop had atoms ids 1-N, and that only DNA atoms are ellipsoids."<<endl;
  }

  atoms = new atom[Ne];

  // now get the atoms lines
  for (int i=0;i<N;i++) {
    int id, type, mol, eflag, density, ix, iy, iz;
    double x,y,z;
    getline(inf,line); sline.str(line);
    sline>>id>>type;
    sline>>x>>y>>z; 
    sline>>mol>>eflag>>density>>ix>>iy>>iz;
    if (eflag) { // ignore the atom if its not an elipsoid
      atoms[id-1].x=x; atoms[id-1].y=y; atoms[id-1].z=z;
      atoms[id-1].x+=ix*lx; atoms[id-1].y+=iy*ly; atoms[id-1].z+=iz*lz;
    }
    sline.clear();
  }

  do {
    getline(inf,line); sline.str(line); sline>>first>>second; sline.clear();
  } while ( first.compare("Ellipsoids")!=0 );
  getline(inf,line); // a blank line
  
  // now get the ellipsoid lines 
  for (int i=0;i<Ne;i++) {
    int id;
    double shape;
    getline(inf,line); sline.str(line);
    sline>>id>>shape>>shape>>shape;
    sline>>atoms[id-1].q0>>atoms[id-1].q1>>atoms[id-1].q2>>atoms[id-1].q3;
    sline.clear();
  }
  inf.close();

  // do the rotations
  Lk = atoi(argc[3]);
  d_theta = Lk * 2 * PI / double(Ne);
  theta=0.0;
  for (int i=0;i<Ne;i++) {
    axis=atoms[i].zaxis();  // rotate about the z axis
    atoms[i].rotate(axis,theta);
    theta+=d_theta;         // by an angle which increases each time
  }

  // now output
  inf.open(argc[1]);
  ouf.open(argc[2]);
  do {
    getline(inf,line); sline.str(line); sline>>first; sline.clear();
    ouf<<line<<endl;
  } while ( first.compare("Ellipsoids")!=0 );
  getline(inf,line); ouf<<line<<endl;; // a blank line
  for (int i=0;i<Ne;i++) {
    getline(inf,line);
    ouf<<i+1<<" 1.0 1.0 1.0 "<<atoms[i].q0<<" "<<atoms[i].q1<<" "<<atoms[i].q2<<" "<<atoms[i].q3<<endl;
  }
  do {
    getline(inf,line);
    ouf<<line<<endl;
  } while ( !inf.eof() );
  ouf.close();
  inf.close();

  // clean up
  delete [] atoms;

}

