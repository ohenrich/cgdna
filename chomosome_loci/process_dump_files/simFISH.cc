//############################################################################
//### Read a bed file of FISH probes, and extract distance from a set of dumps
//############################################################################

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cmath>

#include "simFISH.h"
#include "general.h"

using namespace std;

void get_fish(int argc, char *argv[]) {

  if (argc!=8) {
    cout<<"### Command SIMULATEDFISH ######################################"<<endl;
    cout<<"Usage : "<<endl;
    cout<<"      ./process_dumps.out SIMULATEDFISH N Nfiles infile outfile configfile probesfile.bed"<<endl;
    cout<<"where N is the number of chromatin beads,"<<endl;
    cout<<"      Nfiles is the number of files to read,"<<endl;
    cout<<"      infile is the dumpfile name in the format dumpfile.999.suffix"<<endl;
    cout<<"      outfile is the output filename"<<endl;
    cout<<"      configfile is the file used to set up the simulations, detailing how bead map to bp"<<endl;
    cout<<"      and probesfile.bed is a standard bed file listing probe locations"<<endl;
    exit(0);
  }
  cout<<"### Extracting separations of FISH probes from dump files"<<endl;

  string fullfilename=string( argv[4] ),
    filename1=fullfilename.substr(0,(fullfilename.substr(0,fullfilename.find_last_of("."))).find_last_of(".")+1),
    filename2=fullfilename.substr(fullfilename.find_last_of(".")),
    configfile=string( argv[6] ),
    probesfile=string( argv[7] ),
    outfile=string( argv[5] );

  ifstream inf;
  ofstream ouf;
  string line;

  // Check inputs are integers
  for (int i=0;i<strlen(argv[2]);i++) {
    if (!isdigit(argv[2][i])) {cout<<"ERROR : number of DNA beads must be an integer."<<endl; exit(0);}
  }
  for (int i=0;i<strlen(argv[3]);i++) {
    if (!isdigit(argv[3][i])) {cout<<"ERROR : number of files must be an integer."<<endl; exit(0);}
  }

  int Nfiles=atoi(argv[3]),
    Ndna=atoi(argv[2]);

  vector<FISHprobe> probes;

  // set up simulation system
  simluationsystem thesystem(filename1+to_string(1)+filename2,Ndna);
  atom *bead;

  bead = new atom[thesystem.Ndna];

  // load configuration file
  thesystem.getconfiguration(configfile);

  // read probes file
  inf.open( probesfile.c_str() );
  if ( inf.good() ) {
    cout<<"Reading probes file "<<probesfile<<endl;
    FISHprobe aprobe;
    getline(inf,line);
    do {
      aprobe.read_probe(thesystem,line);
      if ( aprobe.inregion ) { probes.push_back( aprobe ); }
    } while( getline(inf,line) );

  } else {
    cout<<" ERROR : cannot open file "<<probesfile<<endl;
    exit(0);
  }
  inf.close();
  if ( probes.size()<2 ) {
    cout<<" ERROR : Not enough valid probe locations in probes file."<<endl;
    exit(0);
  }

  // set up output
  inf.open( outfile.c_str() );
  if ( inf.good() ) {
    cout<<" ERROR : output file "<<outfile<<" aleady exists"<<endl;
    exit(0);
  }
  inf.close();
  ouf.open( outfile.c_str() );
  ouf<<"# file number, distances:";
  for (int i=0;i<probes.size();i++) {
    for (int j=i+1;j<probes.size();j++) {
      ouf<<", "<<probes[i].name<<"-"<<probes[j].name;
    }
  }
  ouf<<endl;

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

    load_conf_unwrapped(bead,thesystem, inf);

    // get separations and output
    ouf<<file+1;
    for (int i=0;i<probes.size();i++) {
      for (int j=i+1;j<probes.size();j++) {
	ouf<<" "<<probes[i].get_sep(bead,probes[j]);
      }
    }
    ouf<<endl;

    inf.close();
  }

  cout<<endl;

  // tidy up
  ouf.close();
  delete [] bead;

}



//### Member functions of FISHprobe

//FISHprobe 
void FISHprobe::read_probe(const simluationsystem &thesystem, const string &line) {

  stringstream sline;

  // get probe
  sline.str(line);
  sline>>chromosome;
  sline>>start_bp;
  sline>>end_bp;
  sline>>name;

  // test it
  if ( !thesystem.configuration_loaded_flag ) {
    cout<<" ERROR : simulation configuration has not been loaded"<<endl;
    exit(0);
  } 
  inregion=1;
  if ( (thesystem.chromosome != chromosome) ||
       !(end_bp >= thesystem.start_bp && start_bp <= thesystem.end_bp)
       ) {
    inregion=0;
    cout<<" WARNING : bedfile entry not in simulated region will be ignored. "<<line<<endl;
  }

  // convert to beads
  if ( inregion ) {
    double dstart,dend;
    dstart=thesystem.bp_to_doublebead(start_bp);
    dend=thesystem.bp_to_doublebead(end_bp);
    // only include beads if more than 20% of it is covered by the probe
    if (dend-dstart>2) {
      if ( dstart-floor(dstart) < 0.8 ) { start_bead=floor(dstart); } else { start_bead=ceil(dstart); }
      if ( dend-floor(dend) > 0.2 ) { end_bead=ceil(dend); } else { end_bead=floor(dend); }
    } else {
      start_bead=floor(dstart);
      end_bead=ceil(dend);
    }
    length_bead = end_bead - start_bead +1;
    cout<<" Probe "<<name<<" beads "<<start_bead<<" to "<<end_bead<<" inclusive ("<<length_bead<<" beads)."<<endl;
  }

}


double FISHprobe::get_sep(const atom beads[], const FISHprobe &B ) const {
  // find sepatation between this and another probe
  // rememer start_bead and end_bead count from 1

  atom centreA, centreB;

  // get centers
  centreA.clear();
  for (int i=this->start_bead;i<=this->end_bead;i++) {
    centreA += beads[i-1];
  }
  centreA /= double(this->length_bead);

  centreB.clear();
  for (int i=B.start_bead;i<=B.end_bead;i++) {
    centreB += beads[i-1];
  }
  centreB /= double(B.length_bead);

  return (centreA-centreB).length();

}
