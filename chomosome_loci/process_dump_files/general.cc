
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "general.h"

using namespace std;

//### Global functions

string to_string(const int &a) { stringstream ss; ss<<a; return ss.str(); }

void fastforward_frames(const simluationsystem &thesystem, ifstream &inf) {
  // skip forward all but the last frame
  string line;
  for (int fr=0;fr<thesystem.frames_per_file-1;fr++){
    for (int j=0;j<thesystem.Natoms+9;j++) {
      getline(inf,line);
    }
  }
}

void load_conf_unwrapped(atom *bead, const simluationsystem &thesystem, ifstream &inf) {
  // get an unwrapped conformation of polymer beads
  string line;
  stringstream sline;

  int id, type, ix, iy, iz;

  // clear bead
  for (int i=0;i<thesystem.Ndna;i++) {
    bead[i].clear();
  }

  // lead configuration from file
  for (int i=0;i<9;i++) {getline(inf,line);} // lines of junk
  for (int i=0;i<thesystem.Natoms;i++) {
    getline(inf,line);
    sline.str(line);
    sline>>id>>type;
    if (id<=thesystem.Ndna) {
      sline>>bead[id-1].x>>bead[id-1].y>>bead[id-1].z>>ix>>iy>>iz;; 
      bead[id-1].x*=thesystem.lx; bead[id-1].y*=thesystem.ly;bead[id-1].z*=thesystem.lz;
      // unwrap pbcs
      bead[id-1].x+=ix*thesystem.lx; bead[id-1].y+=iy*thesystem.ly; bead[id-1].z+=iz*thesystem.lz;
    }
    sline.clear();
  }

}

void load_conf_wrapped(atom *bead, const simluationsystem &thesystem, ifstream &inf) {
  // get an unwrapped conformation of polymer beads
  string line;
  stringstream sline;

  int id, type, ix, iy, iz;

  // clear bead
  for (int i=0;i<thesystem.Ndna;i++) {
    bead[i].clear();
  }

  // lead configuration from file
  for (int i=0;i<9;i++) {getline(inf,line);} // lines of junk
  for (int i=0;i<thesystem.Natoms;i++) {
    getline(inf,line);
    sline.str(line);
    sline>>id>>type;
    if (id<=thesystem.Ndna) {
      sline>>bead[id-1].x>>bead[id-1].y>>bead[id-1].z>>ix>>iy>>iz;; 
      bead[id-1].x*=thesystem.lx; bead[id-1].y*=thesystem.ly;bead[id-1].z*=thesystem.lz;
      // DO NOT unwrap pbcs
      //bead[id-1].x+=ix*thesystem.lx; bead[id-1].y+=iy*thesystem.ly; bead[id-1].z+=iz*thesystem.lz;
    }
    sline.clear();
  }

}

//### Memeber functions of random
randoms::randoms(void) {
 initialized=0;
}

void randoms::initialize(const int &seed) {
  initialized=1;
  srand(seed);
}

int randoms::integer(const int &start, const int &end) {
  // return a random integer between start and end inclusive
  if ( initialized==0 ) {
    cout<<"Error : using random numbers which have not been initialized."<<endl;
    exit(0);
  }
  return start + ( rand() % (int)(end-start+1) ) ;
}

//### Memeber functions of simluationsystem
simluationsystem::simluationsystem(const string &filename, const int &Ndna) : Ndna(Ndna) {
  // read some system properties from a dump file
  ifstream inf;
  ofstream ouf;
  string line;
  stringstream sline;
    
  // some miscilaneous variables
  double dub1, dub2;

  inf.open( filename.c_str() );
  if ( !inf.good() ) {
    cout<<" ERROR while trying to open file "<<filename<<endl;
    exit(0);
  }
  
  getline(inf,line);
  getline(inf,line); sline.str(line); sline>>starttime; sline.clear();
  getline(inf,line);
  getline(inf,line); sline.str(line); sline>>Natoms; sline.clear();
  getline(inf,line);
  getline(inf,line); sline.str(line); sline>>dub1>>dub2; lx=dub2-dub1; sline.clear();
  getline(inf,line); sline.str(line); sline>>dub1>>dub2; ly=dub2-dub1; sline.clear();
  getline(inf,line); sline.str(line); sline>>dub1>>dub2; lz=dub2-dub1; sline.clear();
  getline(inf,line);
  for (int i=0;i<Natoms;i++) {getline(inf,line);} // read first frame
  getline(inf,line);
  getline(inf,line); sline.str(line); sline>>step; sline.clear();
  step=step-starttime;
  inf.close();
  Nprot=Natoms-Ndna;

  // get number of frames in the file
  FILE *pee;
  char command[200];
  int flen;

  sprintf(command,"wc -l %s | awk '{ print $1 }'",filename.c_str());
  pee=popen(command,"r");
  fscanf(pee,"%u",&flen);
  pclose(pee);
  flen=flen;

  frames_per_file=flen/(Natoms+9);

  cout<<" Reading from dump files with "<<Natoms<<" atoms,"<<endl;
  cout<<" including "<<Ndna<<" DNA and "<<Nprot<<" proteins,"<<endl;
  cout<<" which have "<<frames_per_file<<" frames per file"<<endl<<endl;

  configuration_loaded_flag=0;

}

void simluationsystem::getconfiguration(const string &filename) {
  // read some system properties from a simulation configuration file
  ifstream inf;
  ofstream ouf;
  string line;
  stringstream sline;

  inf.open( filename.c_str() );
  if ( !inf.good() ) {
    cout<<" ERROR while trying to open file "<<filename<<endl;
    exit(0);
  }

  // get chromosome
  getline(inf,line);
  chromosome=line;

  // get region start in bp
  getline(inf,line);
  for (int i=0;i<line.size();i++) {
    if (!isdigit(line[i])) {cout<<"ERROR : non-integer region start in file "<<filename<<endl; exit(0);}
  }
  sline.clear(); sline.str(line);
  sline>>start_bp;

  // get region end in bp
  getline(inf,line);
  for (int i=0;i<line.size();i++) {
    if (!isdigit(line[i])) {cout<<"ERROR : non-integer region end in file "<<filename<<endl; exit(0);}
  }
  sline.clear(); sline.str(line);
  sline>>end_bp;

  // get bp per bead
  getline(inf,line);
  for (int i=0;i<line.size();i++) {
    if (!isdigit(line[i])) {cout<<"ERROR : non-integer bp per bead in file "<<filename<<endl; exit(0);}
  }
  sline.clear(); sline.str(line);
  sline>>bp_per_bead;

  inf.close();

  configuration_loaded_flag=1;

}

double simluationsystem::bp_to_doublebead(const int &bp) const {
  // return the bead number as a double (bead numbers start from 1)
  if ( !configuration_loaded_flag ) {
    cout<<" ERROR : simulation configuration has not been loaded"<<endl;
    exit(0);
  }
  return double(bp-start_bp)/double(bp_per_bead)+1;
}

//### Memeber functions of atom
atom atom::operator-(const atom &b) const { atom c; c.x=x-b.x; c.y=y-b.y; c.z=z-b.z; return c; }
void atom::operator+=(const atom &a) { x+=a.x; y+=a.y; z+=a.z; };
void atom::operator/=(const double &a) { x/=a; y/=a; z/=a; };
double atom::length() { return sqrt(x*x + y*y + z*z); }
void atom::clear() {x=0.0; y=0.0; z=0.0; type=0;}
double atom::distance_wrapped(const simluationsystem &thesystem, const atom &b) const { 
  // find distance to another atom. assumes distances have not been unwrapped
  double sepx,sepy,sepz;

  sepx=x - b.x;  if (sepx>0.5*thesystem.lx) {sepx=thesystem.lx-sepx;}
  sepy=y - b.y;  if (sepy>0.5*thesystem.ly) {sepy=thesystem.ly-sepy;}
  sepz=z - b.z;  if (sepz>0.5*thesystem.lz) {sepz=thesystem.lz-sepz;}

  return sqrt( sepx*sepx + sepy*sepy + sepz*sepz );

}
