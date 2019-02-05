//############################################################################
//### Calculate the radius of gyration from a dump file. 
//############################################################################

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>

#include "rg.h"
#include "general.h"

using namespace std;

void gyration_radius_t(int argc, char *argv[]) {
  // calculate the Rg as a function of time, averaged over all dump files

 if (argc!=6) {
    cout<<"### Command RGT ######################################"<<endl;
    cout<<"Usage : "<<endl;
    cout<<"      ./process_dumps.out RG N Nfiles infile outfile"<<endl;
    cout<<"where N is the number of chromatin beads,"<<endl;
    cout<<"      Nfiles is the number of files to read,"<<endl;
    cout<<"      infile is the dumpfile name in the format dumpfile.999.suffix"<<endl;
    cout<<"      and outfile is the output filename"<<endl;
    exit(0);
  }
  cout<<"### Calculating radius of gyration as a function of time."<<endl;

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
 
  // set up Rg
  double *rg,
    *rg2,
    current_rg;
  int counter=0;

  rg = new double[thesystem.frames_per_file];
  rg2 = new double[thesystem.frames_per_file];
  for (int i=0;i<thesystem.frames_per_file;i++) {
    rg[i]=0.0;
    rg2[i]=0.0;
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
      current_rg=get_Rg(thesystem,inf,1,thesystem.Ndna);
      rg[fr]+=current_rg;
      rg2[fr]+=current_rg*current_rg;
    }
    counter++;

    inf.close();
  }

  // finish means
  for (int fr=0;fr<thesystem.frames_per_file;fr++) {
    rg[fr]/=double(counter);
    rg2[fr]/=double(counter);
  }

  // ouput
  ouf.open( outfile.c_str() );
  ouf<<"# frame number (time), mean R_g, mean squared, number of values in mean"<<endl;
  for (int fr=0;fr<thesystem.frames_per_file;fr++) {
    ouf<<fr+1<<" "<<rg[fr]<<" "<<rg2[fr]<<" "<<counter<<endl;
  }
  ouf.close();

}


void gyration_radius_r(int argc, char *argv[]) {
  // calculate Rg of a certain region of the polymer, for each dump file.

 if (argc!=8) {
    cout<<"### Command RGR ######################################"<<endl;
    cout<<"Usage : "<<endl;
    cout<<"      ./process_dumps.out RGR N Nfiles start end infile outfile"<<endl;
    cout<<"where N is the number of chromatin beads,"<<endl;
    cout<<"      Nfiles is the number of files to read,"<<endl;
    cout<<"      start and end at the (1 based) bead numbers of the region for which to caluclate Rg"<<endl;
    cout<<"      infile is the dumpfile name in the format dumpfile.999.suffix"<<endl;
    cout<<"      and outfile is the output filename"<<endl;
    exit(0);
  }

  string fullfilename=string( argv[6] ),
    filename1=fullfilename.substr(0,fullfilename.find_first_of(".")+1),
    filename2=fullfilename.substr(fullfilename.find_last_of(".")),
    outfile=string( argv[7] );

  int Nfiles=atoi(argv[3]),
    Ndna=atoi(argv[2]),
    start=atoi(argv[4]),
    end=atoi(argv[5]);

  cout<<"### Calculating radius of gyration of the region between beads "<<start<<"-"<<end<<" inclusive."<<endl;

  ifstream inf;
  ofstream ouf;

  // initialize the simulated system
  simluationsystem thesystem(filename1+to_string(1)+filename2,Ndna);
 
  // set up Rg
  double *rg;

  rg = new double[Nfiles];
  for (int i=0;i<Nfiles;i++) {
    rg[i]=0.0;
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

    fastforward_frames(thesystem,inf);
    rg[file]=get_Rg(thesystem,inf,start,end);

    inf.close();
  }

  // output
  ouf.open( outfile.c_str() );
  ouf<<"# file number, radius of regions (beads "<<start<<"-"<<end<<" inclusive)"<<endl;
  for (int i=0;i<Nfiles;i++) {
    ouf<<i+1<<" "<<rg[i]<<endl;
  }
  ouf.close();

}

//####### global functions

double get_Rg(const simluationsystem &thesystem, ifstream &inf, const int &start, const int &end) {
  // calulate the Rg of a single frame, for the region from startbead to endbead
  // startbead and endbead are in "start at 1" numbering, so need to convert to "start at 0"

  string line;
  stringstream sline;

  int id, type, ix, iy, iz;
  double comx,comy,comz,    
    rg2x,rg2y,rg2z,
    rg2;
  int counter,
    startbead=start-1,   // change to 0 start numbering
    endbead=end-1;

  atom *bead;  
  bead = new atom[thesystem.Ndna];

  // load the beads
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

  // first find com
  comx=0; comy=0; comz=0;
  counter=0;
  for (int i=startbead;i<=endbead;i++) {
    comx+=bead[i].x;
    comy+=bead[i].y;
    comz+=bead[i].z;
    counter++;
  }
  comx/=double(counter);
  comy/=double(counter);
  comz/=double(counter);

  // calculate Rg
  rg2x=0; rg2y=0; rg2z=0;
  for (int i=startbead;i<=endbead;i++) {
    rg2x+=pow(bead[i].x-comx,2);
    rg2y+=pow(bead[i].y-comy,2);
    rg2z+=pow(bead[i].z-comz,2);
  }
  rg2 = rg2x + rg2y + rg2z;
  rg2/=double(endbead-startbead+1);
 
  delete [] bead;

  return sqrt(rg2);

}
