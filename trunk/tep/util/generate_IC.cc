// Generate a random linear DNA 
// with spherical proteins or nucleosome cores

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#define PI 3.14159265358979

#include "generate_DNA.h"
#include "generate_functions.cc"

//#include "nucleosome.cc"

using namespace std;

int main() {

  int N,              // DNA length
    S,                // size of side of grid if loop
    Nnuc,             // number of nucleosomes
    Natpernuc,        // number of atoms per nucleosome
    Nprot,            // number of proteins
    nuc_flag,         // nucleosomes
    prot_flag,        // proteins
    Dpatch_flag,      // patches on DNA?
    Ppatch_flag,      // patches on proteins?
    Nppatch,          // number of patches on proteins
    tors_flag,        // flag for if DNA is to have torsional stifness
    loop_flag,        // flag for DNA arrangement
    Nfiles,           // number of files to generate
    seed;             // for randoms

  double lx,ly,lz,    // box size
    density,          // density of DNA beads
    dp_sep,           // separation of DNA bead and patch
    pp_sep,           // separation of protein and patch
    dp_angle;         // DNA patch orientation offset angle

  stringstream fn;    // file name
  ofstream ouf;

  vector<atom> atoms;
  vector<bond> bonds;
  vector<angle> angles;


  // get parameters
 lconfig:
  cout<<"DNA configuration : Enter 1 for loop, 0 for random linear"<<endl; 
  cin>>loop_flag;
  if (loop_flag!=1 && loop_flag!=0) {goto lconfig;}

  cout<<"Length of DNA : "<<endl; 
  cin>>N;
  if (loop_flag==1) {
    S=sqrt(N); if (S%2 == 1) S++;  N=S*S;
    cout<<"Using "<<N<<" DNA beads,"<<endl;
  }

 lbox:
  cout<<"size of box,x,y,x"<<endl; 
  cin>>lx>>ly>>lz;
  if (loop_flag==1 && (lx<1.5*S || ly<1.5*S || lz<1.5*S) ) {cout<<"Box too small."<<endl; goto lbox;}

 ltors:
  cout<<"Add torsional rigidity to DNA? Enter 0 for no, 1 for yes : "<<endl; 
  cin>>tors_flag;
  if (tors_flag!=0 && tors_flag!=1) {goto ltors;}

 ldpat:
  cout<<"Add patches to DNA? Enter 0 for no, 1 for yes : "<<endl; 
  cin>>Dpatch_flag;
  if (Dpatch_flag!=0 && Dpatch_flag!=1) {goto ldpat;}
  if (Dpatch_flag) {
    cout<<"Enter orientational offset angle (degrees) for patch compared to DNA bead orientation : "<<endl;
    cin>>dp_angle;
  }

  //cout<<"Number of nucleosome cores : "<<endl; cin>>Nnuc; 
  //if (Nnuc>0) {nuc_flag=1;} else {nuc_flag=0;}
nuc_flag=0;

  cout<<"Number of other proteins : "<<endl; cin>>Nprot; 
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
      cout<<"Distance to patch (default 0.4)"<<endl;
      cin>>pp_sep;
    }
  } else {prot_flag=0;}

  cout<<"Number of files to generate : "<<endl; cin>>Nfiles;

  cout<<"Enter seed for random numbers : "<<endl; cin>>seed;


  // set up atom types
  if (Dpatch_flag==0) {TYPE.DPAT=0; TYPE.NCORE--; TYPE.NPAT1--; TYPE.NPAT2--; TYPE.PCORE--; TYPE.PPAT--;}
  if (Nnuc==0) {TYPE.NCORE=0; TYPE.NPAT1=0; TYPE.NPAT2=0; TYPE.PCORE-=3; TYPE.PPAT-=3;}

  // set other parameters
  density=(1-0.1*Dpatch_flag)*6.0/PI;
  dp_sep=0.4;  // set dna-patch separation
  //pp_sep=0.4;  // set protein-patch separation
  srand(seed);
  // if (nuc_flag) { Natpernuc=nucleosome().natom(); }


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
    if (loop_flag) {  // DNA is a loop  ******************************************************************************************************************
     
      double x,y,z;
      atom last,lastlast;       // previous DNA atoms

      // 1st core
      x=-S*0.75; y=-S*0.75; z=0.0;
      atoms.push_back( atom(x,y,z) );
      atoms.back().type=TYPE.DNA;
      atoms.back().id=1;
      atoms.back().mol=1;
      atoms.back().ellipse_flag=1;
      atoms.back().density=density;
      atoms.back().q[0]=1; atoms.back().q[1]=0; atoms.back().q[2]=0; atoms.back().q[3]=0;
      last=atoms.back();
      lastlast=atom(0.0,0.0,0.0);
      lastlast.id=0;
      // rest of 1st row
      for (int i=1;i<=S-1;i++) {
	y+=1.5;
	add_DNA(last,lastlast,atoms,bonds,angles,density,x,y,z); 
	lastlast=last; 
	last=atoms.back();
      }
      // next S-2 rows
      for (int j=1;j<=(S-2)/2;j++) {
	// down
	x+=1.5;
	for (int i=1;i<=S-1;i++) {
	  add_DNA(last,lastlast,atoms,bonds,angles,density,x,y,z); 
	  lastlast=last; last=atoms.back();
	  y-=1.5;
	}
	y+=1.5; // we went one too far
	// up
	x+=1.5;
	for (int i=1;i<=S-1;i++) {
	  add_DNA(last,lastlast,atoms,bonds,angles,density,x,y,z); 
	  lastlast=last; last=atoms.back();
	  y+=1.5;
	}
	y-=1.5; // we went one too far
      }
      // final run down
      x+=1.5;
      for (int i=1;i<=S-1;i++) {
	add_DNA(last,lastlast,atoms,bonds,angles,density,x,y,z); 
	lastlast=last; last=atoms.back();
	y-=1.5;
      }
      // and along back to the start
      for (int i=1;i<=S-1;i++) {
	add_DNA(last,lastlast,atoms,bonds,angles,density,x,y,z); 
	lastlast=last; last=atoms.back();
	x-=1.5;
      }
      // complete the loop
      bonds.push_back( bond(last.id,atoms.front().id,TYPE.DNADNA) );
      angles.push_back( angle(lastlast.id,last.id,atoms[0].id,TYPE.BEND) );
      angles.push_back( angle(last.id,atoms[0].id,atoms[1].id,TYPE.BEND) );
      angles.push_back( angle(last.id,atoms[0].id,atoms[0].id,TYPE.TORS) );

    } else {    // DNA is linear ******************************************************************************************************************
      
      atom last,lastlast;       // previous DNA atoms

      // first bead
      atoms.push_back( atom(0.0,0.0,0.0) );
      atoms.back().type=TYPE.DNA;
      atoms.back().id=1;
      atoms.back().mol=1;
      atoms.back().ellipse_flag=1;
      atoms.back().q[0]=1; atoms.back().q[1]=0; atoms.back().q[2]=0; atoms.back().q[3]=0;
      atoms.back().density=density;
      last=atoms.back();       // last and lastlast are the previous two core beads
      lastlast=atom(0.0,0.0,0.0);
      lastlast.id=0;
      // middle beads
      for (int i=1;i<N;i++) {
	add_DNA(last,lastlast,atoms,bonds,angles,(double [3]){lx,ly,lz},dp_sep,density);
	lastlast=last;
	last=atoms.back();  
      }
      // add a final angle for the end of the DNA chain
      angles.push_back( angle(atoms[atoms.size()-2].id,atoms.back().id,atoms.back().id,TYPE.TORSEND) ); // the third id here does nothing 

    } // ******************************************************************************************************************


    // do orientations
    if (loop_flag) { // DNA is a loop 

      evec zax,xax,yax;  // vectors
      quaternion q;      // quaternion

      yax.ei=0; yax.ej=0; yax.ek=-1;
      for (int i=0;i<N-1;i++) {
	zax.ei=atoms[i+1].x-atoms[i].x; zax.ej=atoms[i+1].y-atoms[i].y; zax.ek=atoms[i+1].z-atoms[i].z;
	zax.make_unit();
	xax=yax.cross(zax);
	q.make_quat(xax,yax,zax);
	atoms[i].quat(q);
      }
      zax.ei=atoms[0].x-atoms[N-1].x; zax.ej=atoms[0].y-atoms[N-1].y; zax.ek=atoms[0].z-atoms[N-1].z;
      zax.make_unit();
      xax=yax.cross(zax);
      q.make_quat(xax,yax,zax);
      atoms[N-1].quat(q);

    } else {    // DNA is linear
      // Just leave all oreintations pointing in z direction
    }

    // do DNA patches
    if (Dpatch_flag) {

      int id;        // an intger
      double x,y,z;  // coordinates
      evec xax;      // a vector
      quaternion q;  // a quaternion

      atoms[0].mol=1;
      for (int i=1;i<N;i++) { // each dna bead is a different molecule
	atoms[i].mol=atoms[i-1].mol+1;
      }

      if (dp_angle!=0) {
	cout<<"DNA patch orientation angles not yet implemented."<<endl;
	exit(0);
      }

      id=atoms.back().id;
      for (int i=0;i<N;i++) { // loop round each core
	q.q0=atoms[i].q[0];
	q.q1=atoms[i].q[1];
	q.q2=atoms[i].q[2];
	q.q3=atoms[i].q[3];
	xax=q.xaxis();
	if (i%2==0) {
	  xax.ei=-xax.ei;
	  xax.ej=-xax.ej;
	  xax.ek=-xax.ek;
	}
	x=atoms[i].x+dp_sep*xax.ei;
	y=atoms[i].y+dp_sep*xax.ej;
	z=atoms[i].z+dp_sep*xax.ek;
	// make the new atom
	id++;
	atoms.push_back( atom(x,y,z) );
	atoms.back().ellipse_flag=0;
	atoms.back().id=id;
	atoms.back().type=TYPE.DPAT;
	atoms.back().mol=atoms[i].mol;
	atoms.back().density=1.0;
      }

    }

    // do nucleosomes
    /*
    if (nuc_flag) {
 
      double x,y,z;
      int id=atoms.back().id,
	mol=atoms.back().mol;

      for (int n=0;n<Nnuc;n++) {
	x=lx*double(rand())/double(RAND_MAX)-lx*0.5;
	y=ly*double(rand())/double(RAND_MAX)-ly*0.5;
	z=lz*double(rand())/double(RAND_MAX)-lz*0.5;

	nucleosome nuc(x,y,z); // makes a new nucleosome
	mol++;  // which is a new molecule

	for (int i=0;i<nuc.Ncores;i++) {
	  id++;
	  atoms.push_back( nuc.cores[i] );
	  atoms.back().ellipse_flag=0;
	  atoms.back().id=id;
	  atoms.back().type=TYPE.NCORE;
	  atoms.back().mol=mol;
	  atoms.back().density=1.0;
	}
	for (int i=0;i<nuc.Npatches;i++) {
	  id++;
	  atoms.push_back( nuc.patches[i] );
	  atoms.back().ellipse_flag=0;
	  atoms.back().id=id;
	  // atoms.back().type;  type is already set for patches
	  atoms.back().mol=mol;
	  atoms.back().density=1.0;
	}

      }

    }
*/

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
	atoms.back().ellipse_flag=0;
	atoms.back().id=id;
	atoms.back().type=TYPE.PCORE;
	atoms.back().mol=mol;
	atoms.back().density=1.0;
	if (Ppatch_flag) {
	  id++;
	  atoms.push_back( atom(x+pp_sep,y,z) );
	  atoms.back().ellipse_flag=0;
	  atoms.back().id=id;
	  atoms.back().type=TYPE.PPAT;
	  atoms.back().mol=mol;
	  atoms.back().density=1.0;
	  if (Nppatch==2) {
	    id++;
	    atoms.push_back( atom(x-pp_sep,y,z) );
	    atoms.back().ellipse_flag=0;
	    atoms.back().id=id;
	    atoms.back().type=TYPE.PPAT;
	    atoms.back().mol=mol;
	    atoms.back().density=1.0;
	  }
	}
      }
    }

    // output
    ouf.open(fn.str().c_str());
    ouf<<" LAMMPS data file"<<endl;
    ouf<<endl;
    ouf<<" "<<atoms.size()<<" atoms"<<endl;
    if (tors_flag==1) {ouf<<" "<<N<<" ellipsoids"<<endl;}
    ouf<<" "<<bonds.size()<<" bonds"<<endl;
    if (tors_flag==1) {
      ouf<<" "<<angles.size()<<" angles"<<endl;
    } else {
      int count=0;
      for (int i=0;i<angles.size();i++) {
	if (angles[i].type==TYPE.BEND) count++;
      }
      ouf<<" "<<count<<" angles"<<endl;
    }
    ouf<<endl;
    ouf<<" "<<1+Dpatch_flag+prot_flag*(1+Ppatch_flag)+3*nuc_flag<<" atom types"<<endl;
    ouf<<" 1 bond types"<<endl;
    ouf<<1+2*tors_flag<<" angle types"<<endl;
    ouf<<endl;
    ouf<<" "<<-0.5*lx<<" "<<0.5*lx<<" xlo xhi"<<endl;
    ouf<<" "<<-0.5*ly<<" "<<0.5*ly<<" ylo yhi"<<endl;
    ouf<<" "<<-0.5*lz<<" "<<0.5*lz<<" zlo zhi"<<endl;
    ouf<<endl;
    ouf<<endl<<" Masses"<<endl<<endl;
    ouf<<TYPE.DNA<<" "<<1-0.1*Dpatch_flag<<endl; // mass is ignored for ellipsoids
    if (Dpatch_flag) {
      ouf<<TYPE.DPAT<<" "<<"0.1"<<endl;
    } 
    if (nuc_flag) {
      ouf<<TYPE.NCORE<<" 0.714"<<endl;
      ouf<<TYPE.NPAT1<<" 0.714"<<endl;
      ouf<<TYPE.NPAT2<<" 0.714"<<endl;
    }  
    if (prot_flag) {
      ouf<<TYPE.PCORE<<" "<<1.0-0.5*Ppatch_flag<<endl;
      if (Ppatch_flag) {
	ouf<<TYPE.PPAT<<" "<<0.5/Nppatch<<endl;
      }
    } 
    if (tors_flag==1) { // output atoms, ellipsoids and velocities
      output_tors( ouf,atoms );
    } else {
      output_notors( ouf,atoms );
    }
    ouf<<endl<<" Bonds"<<endl<<endl;
    for (int i=0;i<bonds.size();i++) {
      ouf<<" "<<i+1<<" "<<bonds[i].type<<" "<<bonds[i].a<<" "<<bonds[i].b<<endl;
    }
    ouf<<endl<<" Angles"<<endl<<endl;
    for (int i=0;i<angles.size();i++) {
      if (tors_flag==1) {
	ouf<<" "<<i+1<<" "<<angles[i].type<<" "<<angles[i].a<<" "<<angles[i].b<<" "<<angles[i].c<<endl;
      } else {
	if (angles[i].type==TYPE.BEND) {ouf<<" "<<i+1<<" "<<angles[i].type<<" "<<angles[i].a<<" "<<angles[i].b<<" "<<angles[i].c<<endl;}
      }
    }
    ouf.close();

    // clean up vectors
    atoms.clear(); vector<atom>().swap(atoms);
    bonds.clear(); vector<bond>().swap(bonds);
    angles.clear(); vector<angle>().swap(angles);

  }



}



