//############################################################################
//### Make a contact map from the final frame of a series of dump files
//############################################################################

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <map>

#include "contactmap.h"
#include "general.h"

using namespace std;

void generate_contactmap(int argc, char *argv[], double threshold) {

  if (argc!=6) {
    cout<<"### Command CONTACTMAP ######################################"<<endl;
    cout<<"Usage : "<<endl;
    cout<<"      ./process_dumps.out CONTACTMAP N Nfiles infile outfile"<<endl;
    cout<<"where N is the number of chromatin beads,"<<endl;
    cout<<"      Nfiles is the number of files to read,"<<endl;
    cout<<"      infile is the dumpfile name in the format dumpfile.999.suffix"<<endl;
    cout<<"      and outfile is the output filename"<<endl;
    exit(0);
  }
  cout<<"### Generating contact map from final frame of dump files."<<endl;

  string fullfilename=string( argv[4] ),
    filename1=fullfilename.substr(0,(fullfilename.substr(0,fullfilename.find_last_of("."))).find_last_of(".")+1),
    filename2=fullfilename.substr(fullfilename.find_last_of(".")),
    outfile=string( argv[5] );

  int Nfiles=atoi(argv[3]),
    Ndna=atoi(argv[2]);

  ifstream inf;
  ofstream ouf;

  // initialize the simulated system
  simluationsystem thesystem(filename1+to_string(1)+filename2,Ndna);
 
  // set up contact map
  cmap average_contacts(thesystem.Ndna);
  int counter=0;

  // loop round files
  for (int file=0;file<Nfiles;file++) {
    inf.open( (filename1+to_string(file+1)+filename2).c_str() );
    if ( !inf.good() ) {
      cout<<" ERROR while trying to open file "<<filename1+to_string(file+1)+filename2<<endl;
      exit(0);
    } else {
      cout<<"Reading file "<<filename1+to_string(file+1)+filename2<<"\r"<<flush;
    }

    fastforward_frames(thesystem,inf);
  
    average_contacts+=load_contact(thesystem,inf,threshold);
    counter++;

    inf.close();
  }

  // finsh average
  average_contacts/=double(counter);

  // output
  ouf.open( outfile.c_str() );
  cout<<endl<<"Writting file "<<outfile<<endl;
  for (int i=0;i<thesystem.Ndna;i++) {
    for (int j=0;j<thesystem.Ndna;j++) {
      if (j>=i) {
	ouf<<(i+1)<<" "<<(j+1)<<" "<<average_contacts.themap[i][j]<<endl;
      } else {
	ouf<<(i+1)<<" "<<(j+1)<<" "<<average_contacts.themap[j][i]<<endl;
      }
    }
    ouf<<endl;
  }
  ouf.close();

  cout<<"Done."<<endl;

}


//####################################
//####################################
//####################################


void generate_bincontactmap(int argc, char *argv[], double threshold) {

  if (argc!=7) {
    cout<<"### Command BINCONTACTMAP ######################################"<<endl;
    cout<<"Usage : "<<endl;
    cout<<"      ./process_dumps.out CONTACTMAP N Nfiles BPB infile outfile"<<endl;
    cout<<"where N is the number of chromatin beads,"<<endl;
    cout<<"      BPB is the number of beads per bin,"<<endl;
    cout<<"      Nfiles is the number of files to read,"<<endl;
    cout<<"      infile is the dumpfile name in the format dumpfile.999.suffix"<<endl;
    cout<<"      and outfile is the output filename"<<endl;
    exit(0);
  }
  cout<<"### Generating contact map from final frame of dump files."<<endl;

  string fullfilename=string( argv[5] ),
    filename1=fullfilename.substr(0,fullfilename.find_first_of(".")+1),
    filename2=fullfilename.substr(fullfilename.find_last_of(".")),
    outfile=string( argv[6] );

  int Nfiles=atoi(argv[3]),
    bead_per_bin=atoi(argv[4]),
    Ndna=atoi(argv[2]);

  ifstream inf;
  ofstream ouf;

  int Nbins,
    bini,
    binj;


  // calculate number of bins
  Nbins=floor( double(Ndna)/double(bead_per_bin) );
  if (double(Nbins)!=double(Ndna)/double(bead_per_bin) ) {
    cout<<"Warning :: DNA beads do not divide evenly into bins. Last bin will have fewer beads."<<endl;
  }
  cout<<"Bining data into "<<Nbins<<" bins, each representing "<<bead_per_bin<<" beads"<<endl;


  // initialize the simulated system
  simluationsystem thesystem(filename1+to_string(1)+filename2,Ndna);
 
  // set up contact map
  cmap current_contacts(thesystem.Ndna);
  cmap binned_average_contacts(Nbins);
  int counter=0;

  // loop round files
  for (int file=0;file<Nfiles;file++) {
    inf.open( (filename1+to_string(file+1)+filename2).c_str() );
    if ( !inf.good() ) {
      cout<<" ERROR while trying to open file "<<filename1+to_string(file+1)+filename2<<endl;
      exit(0);
    } else {
      cout<<"Reading file "<<filename1+to_string(file+1)+filename2<<"\r"<<flush;
    }

    fastforward_frames(thesystem,inf);
  
    current_contacts=load_contact(thesystem,inf,threshold);

    // do the bining
    for (int i=0;i<thesystem.Ndna;i++) {
      for (int j=0;j<thesystem.Ndna;j++) {
	bini=floor( double(i)/double(bead_per_bin) );
	binj=floor( double(j)/double(bead_per_bin) );
	binned_average_contacts.themap[bini][binj]+=( current_contacts.themap[i][j] );
      }
    }

    counter++;

    inf.close();
  }

  // finsh average
  binned_average_contacts/=double(bead_per_bin*bead_per_bin);
  binned_average_contacts/=counter;

  // output
  ouf.open( outfile.c_str() );
  cout<<endl<<"Writting file "<<outfile<<endl;
  for (int i=0;i<Nbins;i++) {
    for (int j=0;j<Nbins;j++) {
      if (j>=i) {
	ouf<<(i+1)<<" "<<(j+1)<<" "<<binned_average_contacts.themap[i][j]<<endl;
      } else {
	ouf<<(i+1)<<" "<<(j+1)<<" "<<binned_average_contacts.themap[j][i]<<endl;
      }
    }
    ouf<<endl;
  }
  ouf.close();

  cout<<"Done."<<endl;

}







//####### global functions

cmap load_contact(const simluationsystem &thesystem, ifstream &inf, const double &thresh) {

  string line;
  stringstream sline;
  cmap loaded_map(thesystem.Ndna);
  int id, type, ix, iy, iz;
  double sepx, sepy, sepz;

  atom *bead;  
  bead = new atom[thesystem.Ndna];

  // load the beads
  for (int i=0;i<9;i++) {getline(inf,line);} // lines of junk
  for (int i=0;i<thesystem.Natoms;i++) {
    getline(inf,line);
    sline.str(line);
    sline>>id>>type;
    if (id<=thesystem.Ndna) {
      sline>>bead[id-1].x>>bead[id-1].y>>bead[id-1].z; 
      bead[id-1].x*=thesystem.lx; bead[id-1].y*=thesystem.ly;bead[id-1].z*=thesystem.lz;
    }
    sline.clear();
  }

  // generate the map
  for (int i=0;i<thesystem.Ndna;i++) {
    for (int j=i;j<thesystem.Ndna;j++) { // only do each pair once
      sepx=bead[i].x - bead[j].x;  if (sepx>0.5*thesystem.lx) {sepx=thesystem.lx-sepx;}
      sepy=bead[i].y - bead[j].y;  if (sepy>0.5*thesystem.ly) {sepy=thesystem.ly-sepy;}
      sepz=bead[i].z - bead[j].z;  if (sepz>0.5*thesystem.lz) {sepz=thesystem.lz-sepz;}
      if ( sqrt( sepx*sepx + sepy*sepy + sepz*sepz )<thresh ) {
	loaded_map.themap[i][j]++;
	if (i!=j) loaded_map.themap[j][i]++;
      }

    }
  }

  delete [] bead;

  return loaded_map;

}


//####### member functions of cmap

cmap::cmap(const int &N) : N(N) {
  themap = new double*[N];
  for (int i=0;i<N;i++) {
    themap[i] = new double[N];
    for (int j=0;j<N;j++) {
      themap[i][j]=0.0;
    }
  }
}

cmap::~cmap() {
  for (int i=0;i<N;i++) {
    delete [] themap[i];
  }
  delete [] themap;
}

void cmap::operator+=(const cmap &B) {
  // add two maps 
  if (N!=B.N) {
    cout<<"Fatal Error : tried to add maps of unequal size"<<endl;
  }
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++) {
      themap[i][j]+=B.themap[i][j];
    }
  }
}

void cmap::operator=(const cmap &B) {
  // add two maps 
  if (N!=B.N) {
    cout<<"Fatal Error : tried to equate maps of unequal size"<<endl;
  }
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++) {
      themap[i][j]=B.themap[i][j];
    }
  }
}

void cmap::operator/=(const int &a) {
  // devide by an integer
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++) {
      themap[i][j]/=double(a);
    }
  }
}


void cmap::load_cmap(const string &file) {
  // reset contents of map to zero, and load from file.
  ifstream inf;
  string line;
  stringstream sline;
  
  int i,j;
  double val;
  
  // open file, throw error if not there
  cout<<"Reading file "<<file<<" ... ";
  inf.open( file.c_str() );
  if ( !inf.good() ) {
    cout<<" ERROR while trying to open file "<<file<<endl;
    exit(0);
  }

  // set value of cmap to zero
  clear();

  // read map from file, throw error if number of bins incorrect
  getline(inf,line);
  do {

    if ( line!="" && line.substr(0,1)!="#" ) {
      sline.clear();sline.str(line);
      sline>>i>>j>>val;
      if ( i>N || j>N ) {
	cout<<endl<<"Error : contact map in file "<<file<<" does not have size "<<N<<endl;
	exit(0);
      }
      if ( themap[i-1][j-1]!=0 ) {
	cout<<endl<<"Error : contact map in file "<<file<<" has duplicate entries."<<endl;
      }
      themap[i-1][j-1]=val;
    }

  } while ( getline(inf,line) );

  inf.close();

  cout<<"Done."<<endl;

}

void cmap::save_cmap(const string &outfile) {

  ifstream inf;
  ofstream ouf;

  // check file does not exist
  inf.open( outfile.c_str() );
  if ( inf.good() ) {
    cout<<" ERROR : trying to save to file "<<outfile<<" which already exists."<<endl;
    exit(0);
  }
  inf.close();

  ouf.open( outfile.c_str() );
  cout<<endl<<"Writting file "<<outfile<<endl;
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++) {
      if (j>=i) {
	ouf<<(i+1)<<" "<<(j+1)<<" "<<themap[i][j]<<endl;
      } else {
	ouf<<(i+1)<<" "<<(j+1)<<" "<<themap[j][i]<<endl;
      }
    }
    ouf<<endl;
  }
  ouf.close();

}

void cmap::clear() {
  // reset all values in themap to 0
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++) {
      themap[i][j]=0.0;
    }
  }
}
