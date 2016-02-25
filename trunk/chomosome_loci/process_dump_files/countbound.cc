//############################################################################
//### Count the number of proteins bound to the polymer
//############################################################################

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <set>
//#include <map>

#include "countbound.h"
#include "general.h"

using namespace std;

void count_bound(int argc, char *argv[],const double &threshold) {

  if (argc<6) {
    cout<<"### Command COUNTBOUND ######################################"<<endl;
    cout<<"Usage : "<<endl;
    cout<<"      ./process_dumps.out COUNTBOUND N Nfiles infile outfile [A1 A2 A3...]"<<endl;
    cout<<"where N is the number of chromatin beads,"<<endl;
    cout<<"      Nfiles is the number of files to read,"<<endl;
    cout<<"      infile is the dumpfile name in the format dumpfile.999.suffix"<<endl;
    cout<<"      and outfile is the output filename"<<endl;
    cout<<"      A1,A2,A3 etc. are protein bead types to include."<<endl;
    cout<<"If no list of bead types is given, all proteins will be counted."<<endl;
    exit(0);
  }
  

  string fullfilename=string( argv[4] ),
    filename1=fullfilename.substr(0,(fullfilename.substr(0,fullfilename.find_last_of("."))).find_last_of(".")+1),
    filename2=fullfilename.substr(fullfilename.find_last_of(".")),
    outfile=string( argv[5] );

  int Nfiles=atoi(argv[3]),
    Ndna=atoi(argv[2]);
  set<int> prot_types;
  set<int>::iterator it;
  pair<set<int>::iterator,bool> insertreturn;

  ifstream inf;
  ofstream ouf;

  int **Nbound;

  if ( argc>6 ) {
    for (int i=6;i<argc;i++) {
      insertreturn = prot_types.insert( atoi(argv[i]) );
      if ( insertreturn.second==false ) {
	cout<<"ERROR : repeated atom type in list of protein types."<<endl;
	exit(0);
      }
    }
  }

  // test output file does not exist
  inf.open( outfile.c_str() );
  if ( inf.good() ) {
    cout<<" ERROR : output file "<<outfile<<" aleady exists"<<endl;
    exit(0);
  }
  inf.close();

  if ( prot_types.size()>0 ) {
    cout<<"Counting number of proteins of type ";
    for (it=prot_types.begin(); it!=prot_types.end(); ++it) {
      cout<<*it<<" ";
    }
    cout<<"which are bound to the polymer"<<endl;
  } else {
    cout<<"Counting number of proteins which are bound to the polymer"<<endl;
  }

  // initialize the simulated system
  simluationsystem thesystem(filename1+to_string(1)+filename2,Ndna);

  // set up variables
  Nbound = new int*[ Nfiles ];
  for (int i=0;i<Nfiles;i++) {
    Nbound[i] = new int[ thesystem.frames_per_file ];
    for (int j=0;j<thesystem.frames_per_file;j++) {
      Nbound[i][j]=0;
    }
  }

  // loop round files
  for (int file=0;file<Nfiles;file++) {
    inf.open( (filename1+to_string(file+1)+filename2).c_str() );
    if ( !inf.good() ) {
      cout<<" ERROR while trying to open file "<<filename1+to_string(file+1)+filename2<<endl;
      exit(0);
    } else {
      cout<<"Reading file "<<filename1+to_string(file+1)+filename2<<"\r"<<flush;
    }

    // loop round frames
    for (int fr=0;fr<thesystem.frames_per_file;fr++) {
      Nbound[file][fr]=get_bound(thesystem,inf,prot_types,threshold);
    }

    inf.close();
  }

  // output
  ouf.open( outfile.c_str() );
  cout<<endl<<"Writting file "<<outfile<<endl;
  ouf<<"# timestep, number of bound proteins in each of "<<Nfiles<<" dump files"<<endl;
  for (int fr=0;fr<thesystem.frames_per_file;fr++) {
    ouf<<fr*thesystem.step+thesystem.starttime;
    for (int file=0;file<Nfiles;file++) {
      ouf<<" "<<Nbound[file][fr];
    }
    ouf<<endl;
  }
 ouf.close();


  // clean up memory
  for (int i=0;i<Nfiles;i++) {
    delete [] Nbound[i];
  }
  delete [] Nbound;


}




//####### global functions

int get_bound(const simluationsystem &thesystem,ifstream &inf,const set<int> &prot_types,const double &thresh) {

  string line;
  stringstream sline;

  atom poly[ thesystem.Ndna ],
    prot[ thesystem.Nprot ];  

  int counter, id, type, ix, iy, iz;
  double sepx, sepy, sepz;

  int isbound[ thesystem.Nprot ],
    bonds;

  counter=0;
  bonds=0;
  for (int i=0;i<thesystem.Nprot;i++) {
    isbound[i]=0;
  }

  for (int i=0;i<9;i++) {getline(inf,line);} // lines of junk
  for (int i=0;i<thesystem.Natoms;i++) {
    getline(inf,line);
    sline.str(line);
    sline>>id>>type;
    if (id<=thesystem.Ndna) {
      sline>>poly[id-1].x>>poly[id-1].y>>poly[id-1].z; 
      poly[id-1].x*=thesystem.lx; poly[id-1].y*=thesystem.ly;poly[id-1].z*=thesystem.lz;
    } else {
      sline>>prot[counter].x>>prot[counter].y>>prot[counter].z; 
      prot[counter].x*=thesystem.lx; prot[counter].y*=thesystem.ly;prot[counter].z*=thesystem.lz;
      prot[counter].type=type;
      counter++;
    }
    sline.clear();
  }

  // count bound proteins
  for (int i=0;i<thesystem.Ndna;i++) {
    for (int j=0;j<thesystem.Nprot;j++) {
      sepx=poly[i].x - prot[j].x;  if (sepx>0.5*thesystem.lx) {sepx=thesystem.lx-sepx;}
      sepy=poly[i].y - prot[j].y;  if (sepy>0.5*thesystem.ly) {sepy=thesystem.ly-sepy;}
      sepz=poly[i].z - prot[j].z;  if (sepz>0.5*thesystem.lz) {sepz=thesystem.lz-sepz;}
      if ( sqrt( sepx*sepx + sepy*sepy + sepz*sepz )<thresh && 
	   ( prot_types.size()==0 || prot_types.find(prot[j].type)!=prot_types.end() )  ) {
	isbound[j]++;
      }
    }
  }

  for (int i=0;i<thesystem.Nprot;i++) {
    if ( isbound[i]>0 ) bonds++;
  }

  return bonds;

}
