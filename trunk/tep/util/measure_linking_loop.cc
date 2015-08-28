// Program to measure the linking number
// twist and write in a loop of DNA

// dump file atom lines need to have the following format:
// ITEM: ATOMS id type xs ys zs ix iy iz c_quat[1] c_quat[2] c_quat[3] c_quat[4]

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "tools.h"
#include "tools_functions.cc"

#define SMALL 1e-10
#define LARGE 1e9

using namespace std;

int main(int argv, char* argc[]) {

  if (!( argv==3 || argv==4 || argv==5)) {
    cout<<"Usage : a.out input.file output.file [-h] [-noWr]"<<endl;
    cout<<"Where -h option generates histograms"<<endl;
    cout<<"and -noWr skips calculation of writhe"<<endl;
    exit(0);
  }

  int noWr=0,doHist=0;
  if (argv>3) {
    if (strcmp(argc[3],"-h")==0) {doHist=1; cout<<"Doing histogram."<<endl;}
    if (strcmp(argc[3],"-noWr")==0) {noWr=1; cout<<"Skipping Wr calculation."<<endl;}
  }
  if (argv>4) {
    if (strcmp(argc[4],"-h")==0) {doHist=1; cout<<"Doing histogram."<<endl;}
    if (strcmp(argc[4],"-noWr")==0) {noWr=1; cout<<"Skipping Wr calculation."<<endl;}
  }

  char fn[50], outfn[50];
  ifstream inf;
  ofstream ouf;
  string line;
  stringstream sline;

  int integer;              // some miscilaneous variables
  double dub1, dub2;

  atom *atoms;              // an array of atoms
  xyz *t, *p, *v, *u;       // arrays of tangent and twist vectors

  int N,                    // number of atoms
    Nframes,                // number of frames in the dump file
    Nbins=50,               // number of bins in histogram
    *Tbins,*Wbins,*Lbins,   // for histograms
    step;                   // the time step

  double lx,ly,lz,          // box size
    Tw,                     // twist
    Wr;                     // writhe

  double *allTw, *allWr, *allLk;     // for storing all twist and writhe

  int id, type, ix, iy, iz;  // for reading files

  // get some information from header of dump file
  sprintf(fn,argc[1]);
  inf.open(fn);
  for (int i=0;i<3;i++) {getline(inf,line);}
  getline(inf,line); sline.str(line); sline>>N; sline.clear();
  getline(inf,line);
  getline(inf,line); sline.str(line); sline>>dub1>>dub2; lx=dub2-dub1; sline.clear();
  getline(inf,line); sline.str(line); sline>>dub1>>dub2; ly=dub2-dub1; sline.clear();
  getline(inf,line); sline.str(line); sline>>dub1>>dub2; lz=dub2-dub1; sline.clear();
  inf.close();
  Nframes=file_length(fn,N);

  // set up arrays
  atoms = new atom[N];
  t = new xyz[N];
  p = new xyz[N];
  v = new xyz[N];
  u = new xyz[N];

  allTw = new double[Nframes];
  allWr = new double[Nframes];
  allLk = new double[Nframes];

  // set up output file
  sprintf(outfn,argc[2]);
  ouf.open(outfn);
  ouf<<"# step, Tw, Wr, Tw+Wr"<<endl;

  inf.open(fn);
  // loop round frames
  for (int fr=0;fr<Nframes;fr++) {
 
    // read atoms
    getline(inf,line);
    getline(inf,line); sline.str(line); sline>>step; sline.clear();
    for (int i=0;i<7;i++) {getline(inf,line);} // lines of junk
    for (int i=0;i<N;i++) {
      getline(inf,line);
      sline.str(line);
      sline>>id>>type;
      sline>>atoms[id-1].x>>atoms[id-1].y>>atoms[id-1].z; 
      atoms[id-1].x*=lx; atoms[id-1].y*=ly;atoms[id-1].z*=lz;
      sline>>ix>>iy>>iz; 
      atoms[id-1].x+=ix*lx; atoms[id-1].y+=iy*ly; atoms[id-1].z+=iz*lz;
      sline>>atoms[id-1].q0>>atoms[id-1].q1>>atoms[id-1].q2>>atoms[id-1].q3;
      sline.clear();
    }


    // get tangent and twist vectors
    for (int i=0;i<N-1;i++) {
      t[i]=atoms[i+1]-atoms[i];
      t[i].make_unit();
    }
    t[N-1]=atoms[0]-atoms[N-1];
    t[N-1].make_unit();
    for (int i=0;i<N;i++) {
      p[i]=atoms[i].xaxis();
      v[i]=atoms[i].yaxis();
      u[i]=atoms[i].zaxis();
    }

    // Calculate twist
    {
      Tw=0.0;
      /*
      int iminus1;
      xyz n, cvec;
      double alpha,gamma;
      for (int i=0;i<N;i++) {

	if (i==0) {iminus1=N-1;} else {iminus1=i-1;}

	n=t[iminus1].cross(t[i]);
     

	if (n.length()<SMALL) { // t[i-1] and t[i] are parallel
	  cvec = p[iminus1].cross(n);
	  alpha = acos( p[iminus1].dot(n) ) * sign( cvec.dot(t[iminus1]) );
	  gamma = 0.0;
	  cout<<"WARNING: parallel t[i-1] and t[i]"<<endl;
	} else {
	  cvec = p[iminus1].cross(n);
	  alpha = acos( p[iminus1].dot(n) ) * sign( cvec.dot(t[iminus1]) );
	  cvec = n.cross(p[i]);
	  gamma = acos( n.dot(p[i]) ) * sign( cvec.dot(t[i]) );
	}

	if (alpha+gamma>PI) {
	  Tw+=-(2*PI-abs(alpha+gamma));
	} else if(alpha+gamma<-PI) {
	  Tw+=(2*PI-abs(alpha+gamma));
	} else {
	  Tw+=alpha+gamma;
	}
      }
      Tw/=(2.0*PI); // want twist in number of turns, not angle
      */
      double cosag,ag;
      xyz n;
      for (int i=0;i<N-1;i++) {
	cosag= p[i+1].dot(p[i]) + v[i+1].dot(v[i]);
	cosag/= 1.0 + u[i+1].dot(u[i]);
	if (cosag>1.0 && cosag<1.1) {cosag=1.0;}
	if (cosag<-1.0 && cosag>-1.1) {cosag=-1.0;}
	if (cosag>1.0) {cout<<i<<":  "<<cosag<<" "<<acos(cosag)<<" "<<(p[i+1].dot(p[i]) + v[i+1].dot(v[i]))<<" "<<(1+u[i+1].dot(u[i]))<<endl;}
	n=p[i].cross(p[i+1]);
	ag=acos(cosag);
	ag*=sign(n.dot(t[i+1]));
	Tw+=ag;
      }
      cosag= p[0].dot(p[N-1]) + v[0].dot(v[N-1]);
      cosag/= 1.0 + u[0].dot(u[N-1]);
      if (cosag>1.0 && cosag<1.1) {cosag=1.0;}
      if (cosag<-1.0 && cosag>-1.1) {cosag=-1.0;}
      if (cosag>1.0) {cout<<N<<":  "<<cosag<<" "<<acos(cosag)<<" "<<(p[0].dot(p[N-1]) + v[0].dot(v[N-1]))<<" "<<(1+u[0].dot(u[N-1]))<<endl;}
      n=p[N-1].cross(p[0]);
      ag=acos(cosag);
      ag*=sign(n.dot(t[0]));
      Tw+=ag;
      Tw=Tw/(2.0*PI);
   
    }

    // calculate writhe
    if (noWr==1) { Wr=0.0; } else
    {
      Wr=0.0;
      xyz one,two, three, four,
	r12,r34,r23,r13,r14,r24,
	n1,n2,n3,n4,cvec;
      double omega, n1n2, n2n3, n3n4, n4n1;
      for (int i=1;i<N;i++) {
	for (int j=0;j<i;j++) {

	  if (i==j||j==i-1||(i==N-1&&j==0)) {
	    omega=0.0;
	  } else {
	    if (j==0) {
	      three=atoms[N-1];//D[L-1];
	    } else {
	      three=atoms[j-1];//D[j-1];
	    }
	    one=atoms[i-1]; two=atoms[i]; four=atoms[j];
	    r12=two-one;
	    r34=four-three;
	    r23=three-two;
	    r13=three-one;
	    r14=four-one;
	    r24=four-two;

	    n1=r13.cross(r14); if (n1.length()<SMALL) {cout<<"error in writhe 1 "<<n1.length()<<" ("<<fr<<")"<<endl;}
	    n1.make_unit();

	    n2=r14.cross(r24); if (n2.length()<SMALL) {cout<<"error in writhe 2 "<<n1.length()<<" ("<<fr<<")"<<endl;}
	    n2.make_unit();

	    n3=r24.cross(r23); if (n3.length()<SMALL) {cout<<"error in writhe 3 "<<n1.length()<<" ("<<fr<<")"<<endl;}
	    n3.make_unit();

	    n4=r23.cross(r13); if (n4.length()<SMALL) {cout<<"error in writhe 4 "<<n1.length()<<" ("<<fr<<")"<<endl;}
	    n4.make_unit();

	    n1n2=n1.dot(n2); fix_roundoff11(n1n2); 
	    n2n3=n2.dot(n3); fix_roundoff11(n2n3); 
	    n3n4=n3.dot(n4); fix_roundoff11(n3n4); 
	    n4n1=n4.dot(n1); fix_roundoff11(n4n1); 

	    cvec=r34.cross(r12);

	    omega = ( asin( n1n2 ) + asin( n2n3 ) + asin( n3n4 ) + asin( n4n1 ) )*sign( cvec.dot(r13) ); 
	
	  }	  
	  Wr+= omega/(4.0*PI);
	}      
      }
      Wr*=2.0;
 
    }

   allTw[fr]=Tw;
   allWr[fr]=Wr;
   allLk[fr]=Tw+Wr;

   // Output
   ouf<<step<<" "<<Tw<<" "<<Wr<<" "<<Tw+Wr<<endl;

  } // close loop round frames

  // clean up
  inf.close();
  ouf.close();
  delete [] p;
  delete [] t;
  delete [] atoms;

  if (doHist) {
    double Tmax=-LARGE,Tmin=LARGE,
      Wmax=-LARGE,Wmin=LARGE,
      Lmax=-LARGE,Lmin=LARGE,
      Tbin_width,Wbin_width,Lbin_width;

    Tbins = new int[Nbins];
    Wbins = new int[Nbins];
    Lbins = new int[Nbins];
    for (int i=0;i<Nbins;i++) {
      Tbins[i]=0; Wbins[i]=0; Lbins[i]=0;
    }

    // find max and min
    for (int i=0;i<Nframes;i++) { // ignore nans
      if (allLk[i]==allLk[i] && allTw[i]>Tmax) Tmax=allTw[i];
      if (allLk[i]==allLk[i] && allTw[i]<Tmin) Tmin=allTw[i];
      if (allLk[i]==allLk[i] && allWr[i]>Wmax) Wmax=allWr[i];
      if (allLk[i]==allLk[i] && allWr[i]<Wmin) Wmin=allWr[i];
      if (allLk[i]==allLk[i] && allLk[i]>Lmax) Lmax=allLk[i];
      if (allLk[i]==allLk[i] && allLk[i]<Lmin) Lmin=allLk[i];
    }
    Tmax=Tmax+0.05;
    Tmin=Tmin-0.05;
    Wmax=Wmax+0.05;
    Wmin=Wmin-0.05;
    Lmax=Lmax+0.05;
    Lmin=Lmin-0.05;

    Tbin_width=(Tmax-Tmin)/Nbins;
    Wbin_width=(Wmax-Wmin)/Nbins;
    Lbin_width=(Lmax-Lmin)/Nbins;

    for (int i=0;i<Nframes;i++) { 
      if (allLk[i]==allLk[i]) {// ignore nans
	Tbins[int( floor( (allTw[i]-Tmin) / Tbin_width) )]++;
	Wbins[int( floor( (allWr[i]-Wmin) / Wbin_width) )]++;
	Lbins[int( floor( (allLk[i]-Lmin) / Lbin_width) )]++;
      }
    }

    sprintf(outfn,"hist_%s",argc[2]);
    ouf.open(outfn);
    ouf<<"# top edge of T bin (width="<<Tbin_width<<"), count for Tw, top edge of W bin (width="<<
      Wbin_width<<"), count for Wr, top edge of L bin (width="<<Lbin_width<<"), count for Lk"<<endl;
    for (int i=0;i<Nbins;i++) {
      ouf<<(i+1)*Tbin_width+Tmin<<" "<<Tbins[i]<<" "<<(i+1)*Wbin_width+Wmin<<" "<<
	Wbins[i]<<" "<<(i+1)*Lbin_width+Lmin<<" "<<Lbins[i]<<endl; 
    }
    ouf.close();


    delete [] Tbins;
    delete [] Wbins;

  }

  delete [] allTw;
  delete [] allWr;

}
