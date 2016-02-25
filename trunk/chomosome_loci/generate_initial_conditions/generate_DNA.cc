// Generate a random linear DNA 
// with spherical or dumbell proteins

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#define PI 3.14159265358979

#include "generate_DNA.h"

using namespace std;

void add_DNA(atom,atom,vector<atom>&,vector<bond>&,vector<angle>&,double[3],int);

int main() {

  int N,              // DNA length
    kuhn,
    Nprot,            // number of proteins
    prot_flag,        // proteins
    Ppatch_flag,      // patches on proteins?
    Nppatch,          // number of patches on proteins
    Npatch_types,     // are the patches of the same atom type
    Nfiles,           // number of files to generate
    seed;             // for randoms

  double lx,ly,lz,    // box size
    pp_sep;           // separation of protein and patch

  stringstream fn;    // file name
  ofstream ouf;

  vector<atom> atoms;
  vector<bond> bonds;
  vector<angle> angles;


  // get parameters
  cout<<"Length of DNA : "<<endl; 
  cin>>N;
 
 kuhn:
  cout<<"Kuhn length of DNA : "<<endl;
  cin>>kuhn;
  if ( kuhn>0.5*N ) goto kuhn;
  
 lbox:
  cout<<"size of box,x,y,x"<<endl; 
  cin>>lx>>ly>>lz;

  cout<<"Number of proteins : "<<endl; cin>>Nprot; 
  if (Nprot>0) {
    prot_flag=1;
  lppat:
    cout<<"Add patches to proteins?  Enter 0 for no, 1 for yes : "<<endl; 
    cin>>Ppatch_flag;
    if (Ppatch_flag!=0 && Ppatch_flag!=1) {goto lppat;}
    if (Ppatch_flag==1) {
    nppat:
      cout<<"Number of patches on proteins? Enter 1 or 2 : "<<endl;
      cin>>Nppatch;
      if (Nppatch!=1 && Nppatch!=2) {goto nppat;}
      cout<<"Separation of patch and centre"<<endl;
      cin>>pp_sep;
      if (Nppatch==2) {
      nppattype:
	cout<<"Number of patch types? Enter 1 or 2 : "<<endl;
	cin>>Npatch_types;
	if (Npatch_types!=1 && Npatch_types!=2) {goto nppattype;}
      } else {
	Npatch_types=1;
      }
    }
  } else {prot_flag=0;}

  cout<<"Number of files to generate : "<<endl; cin>>Nfiles;

  cout<<"Enter seed for random numbers : "<<endl; cin>>seed;


  // set parameters
  srand(seed);


  // loop round files
  for (int filenum=1;filenum<=Nfiles;filenum++) {

    // file name
    fn.clear();
    fn.str("");
    if (Nfiles==1) {
      fn<<"lammps.input";
    } else {
      fn<<"lammps.input_"<<filenum;
    }


    // do DNA
    atom last,lastlast;       // previous DNA atoms

    // first bead
    atoms.push_back( atom(0.0,0.0,0.0) );
    atoms.back().type=TYPE.DNA;
    atoms.back().id=1;
    atoms.back().mol=1;
    last=atoms.back();       // last and lastlast are the previous two core beads
    lastlast=atom(0.0,0.0,0.0);
    lastlast.id=0;
    // middle beads
    for (int i=1;i<N;i++) {
      if (i % kuhn == 0) {
	add_DNA(last,lastlast,atoms,bonds,angles,(double [3]){lx,ly,lz},1);
      } else {
	add_DNA(last,lastlast,atoms,bonds,angles,(double [3]){lx,ly,lz},0);
      }
      lastlast=last;
      last=atoms.back();  
    }


    // do proteins
    if (Nprot>0) {
      int  id=atoms.back().id,
	mol=atoms.back().mol;
      double x,y,z;
      for (int i=0;i<Nprot;i++) {
	x=lx*double(rand())/double(RAND_MAX)-lx*0.5;
	y=ly*double(rand())/double(RAND_MAX)-ly*0.5;
	z=lz*double(rand())/double(RAND_MAX)-lz*0.5;
	// make the new atom
	id++;
	mol++;
	atoms.push_back( atom(x,y,z) );
	atoms.back().id=id;
	atoms.back().type=TYPE.PCORE;
	atoms.back().mol=mol;
	if (Ppatch_flag) {
	  id++;
	  atoms.push_back( atom(x+pp_sep,y,z) );
	  atoms.back().id=id;
	  atoms.back().type=TYPE.PPAT;
	  atoms.back().mol=mol;
	  if (Nppatch==2) {
	    id++;
	    atoms.push_back( atom(x-pp_sep,y,z) );
	    atoms.back().id=id;
	    if (Npatch_types>1) {
	      atoms.back().type=TYPE.PPAT2;
	    } else {
	      atoms.back().type=TYPE.PPAT;
	    }
	    atoms.back().mol=mol;
	  }
	}
      }
    }

    // output
    ouf.open(fn.str().c_str());
    ouf<<" LAMMPS data file"<<endl;
    ouf<<endl;
    ouf<<" "<<atoms.size()<<" atoms"<<endl;
    ouf<<" "<<bonds.size()<<" bonds"<<endl;
    int count=0;
    for (int i=0;i<angles.size();i++) {
      if (angles[i].type==TYPE.BEND) count++;
    }
    ouf<<" "<<count<<" angles"<<endl;
    ouf<<endl;

    ouf<<" "<<1+prot_flag*(1+Ppatch_flag*Npatch_types)<<" atom types"<<endl;
    ouf<<" 1 bond types"<<endl;
    ouf<<" 1 angle types"<<endl;
    ouf<<endl;

    ouf<<" "<<-0.5*lx<<" "<<0.5*lx<<" xlo xhi"<<endl;
    ouf<<" "<<-0.5*ly<<" "<<0.5*ly<<" ylo yhi"<<endl;
    ouf<<" "<<-0.5*lz<<" "<<0.5*lz<<" zlo zhi"<<endl;
    ouf<<endl;

    ouf<<endl<<" Masses"<<endl<<endl;
    ouf<<TYPE.DNA<<" "<<1<<endl;
    if (prot_flag) {
      ouf<<TYPE.PCORE<<" "<<1.0-0.5*Ppatch_flag<<endl;
      if (Ppatch_flag) {
	ouf<<TYPE.PPAT<<" "<<0.5/Nppatch<<endl;
      }
      if (Npatch_types==2) {
	ouf<<TYPE.PPAT2<<" "<<0.5/Nppatch<<endl;
      }

    } 
    ouf<<endl;

    ouf<<" Atoms"<<endl<<endl;
    for (int i=0;i<atoms.size();i++) {
      ouf<<" "<<atoms[i].id<<" "<<atoms[i].mol<<" "<<atoms[i].type<<" "<<atoms[i].x<<" "<<atoms[i].y<<" "<<atoms[i].z<<"  0 0 0"<<endl;
    }
    ouf<<endl;

    ouf<<" Velocities"<<endl<<endl;
    for (int i=0;i<atoms.size();i++) {
      ouf<<" "<<atoms[i].id<<" 0 0 0"<<endl;
    }

    ouf<<endl<<" Bonds"<<endl<<endl;
    for (int i=0;i<bonds.size();i++) {
      ouf<<" "<<i+1<<" "<<bonds[i].type<<" "<<bonds[i].a<<" "<<bonds[i].b<<endl;
    }
    ouf<<endl;

    ouf<<" Angles"<<endl<<endl;
    for (int i=0;i<angles.size();i++) {
      if (angles[i].type==TYPE.BEND) {ouf<<" "<<i+1<<" "<<angles[i].type<<" "<<angles[i].a<<" "<<angles[i].b<<" "<<angles[i].c<<endl;}  
    }
    ouf.close();

    // clean up vectors
    atoms.clear(); vector<atom>().swap(atoms);
    bonds.clear(); vector<bond>().swap(bonds);
    angles.clear(); vector<angle>().swap(angles);

  }



}










void add_DNA(atom last,atom lastlast,vector<atom> &atoms,vector<bond> &bonds,vector<angle> &angles, double box[3], int changedir) {
  // Add a DNA bead to the random configuration

  static double theta=0,phi=0;
  double  dx,dy,dz,
    hbox[3],id;

  for (int i=0;i<3;i++) {
    hbox[i]=box[i]*0.5;
  }

  if (changedir==0) {
    dx=last.x+sin(theta)*cos(phi);
    dy=last.y+sin(theta)*sin(phi);
    dz=last.z+cos(theta);
  } else if (changedir==1) {
    theta=double(rand())/double(RAND_MAX)*PI;
    phi=double(rand())/double(RAND_MAX)*2.0*PI;
    dx=last.x+sin(theta)*cos(phi);
    dy=last.y+sin(theta)*sin(phi);
    dz=last.z+cos(theta);
  } else {
    cout<<"Error in add_DNA."<<endl;
    exit(0);
  }

  id=atoms.back().id;

  id++;
  atoms.push_back( atom(dx,dy,dz) );
  atoms.back().type=TYPE.DNA;
  atoms.back().id=id;
  atoms.back().mol=last.mol;

  // bond
  bonds.push_back( bond(last.id,atoms.back().id,TYPE.DNADNA) );

  // angle
  if (lastlast.id==0) { // this is bead number 2

  } else { 
    angles.push_back( angle(lastlast.id,last.id,atoms.back().id,TYPE.BEND) );
  } 
   

}

