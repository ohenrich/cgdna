//############################################################################
//### Read a bed file of capc targets, and extract contact profiles 
//### from a set of dumps
//############################################################################

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <vector>
//#include <cmath>

#include "capturec.h"
#include "simFISH.h"
//#include "contactmap.h"
#include "general.h"

using namespace std;

void capture_c(int argc,char *argv[],double threshold) {


  if (argc!=8) {
    cout<<"### Command SIMULATEDCAPC ######################################"<<endl;
    cout<<"Usage : "<<endl;
    cout<<"      ./process_dumps.out SIMULATEDCAPC N Nfiles infile outfileprefix configfile targets.bed"<<endl;
    cout<<"where N is the number of chromatin beads,"<<endl;
    cout<<"      Nfiles is the number of files to read,"<<endl;
    cout<<"      infile is the dumpfile name in the format dumpfile.999.suffix"<<endl;
    cout<<"      outfileprefix is the first part of the output filenames"<<endl;
    cout<<"      configfile is the file used to set up the simulations, detailing how bead map to bp"<<endl;
    cout<<"      targets.bed is a bed file of targeted regions for the capture c"<<endl;
    exit(0);
  }
  cout<<"### Generating contact profiles from final frame of dump files."<<endl;


  string fullfilename=string( argv[4] ),
    filename1=fullfilename.substr(0,(fullfilename.substr(0,fullfilename.find_last_of("."))).find_last_of(".")+1),
    filename2=fullfilename.substr(fullfilename.find_last_of(".")),
    configfile=string( argv[6] ),
    targetfile=string( argv[7] ),
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

  // set up simulation system
  simluationsystem thesystem(filename1+to_string(1)+filename2,Ndna);

  vector<FISHprobe> targets;
  vector< vector<double> > contact_profile;
  //cmap average_contacts(thesystem.Ndna);
  atom *bead;
  //int counter;

  // load configuration file
  thesystem.getconfiguration(configfile);

  // set up variables
  bead = new atom[thesystem.Ndna];

  // read targets file -- looks the same as a FISH targets files, so use those functions
  inf.open( targetfile.c_str() );
  if ( inf.good() ) {
    cout<<"Reading targets file "<<targetfile<<endl;
    FISHprobe atarget;
    getline(inf,line);
    do {
      atarget.read_probe(thesystem,line);
      if ( atarget.inregion ) { targets.push_back( atarget ); }
    } while( getline(inf,line) );

  } else {
    cout<<" ERROR : cannot open file "<<targetfile<<endl;
    exit(0);
  }
  inf.close();


  // set up interaction profiles
  contact_profile = vector< vector<double> >( targets.size(), vector<double>( thesystem.Ndna, 0.0 ) );


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

    // load conformation
    load_conf_wrapped(bead,thesystem, inf);

    // loop round target beads
    for (int t=0;t<targets.size();t++) {
      for (int i=0;i<thesystem.Ndna;i++) {
	bool incontact=0;
	for (int b=targets[t].start_bead;b<=targets[t].end_bead;b++) { 
	  // remember that beads b start at 1, so subtract 1; i starts at 0
	  if ( bead[i].distance_wrapped(thesystem, bead[b-1])<threshold) { incontact=1; }
	}
	if (incontact) { contact_profile[t][i]++; }
      }
    }

    inf.close();
  }

  // finish average
  for (int t=0;t<targets.size();t++) {
    for (int i=0;i<thesystem.Ndna;i++) {
      contact_profile[t][i]/=double(Nfiles);
    }
  }

  // output
  for (int t=0;t<targets.size();t++) {
    // check output file doesn't exist
    inf.open ( (outfile+"_"+targets[t].name+".dat").c_str() );
    if ( inf.good() ) { 
      cout<<" ERROR : file "<<outfile+"_"+targets[t].name+".dat"<<" already exists. Exiting."<<endl;
      exit(0);
    }
    inf.close();
    cout<<"Writting file "<<outfile+"_"+targets[t].name+".dat"<<endl;
    ouf.open( (outfile+"_"+targets[t].name+".dat").c_str() );
    ouf<<"# bead number, contact probability"<<endl;
    for (int i=0;i<thesystem.Ndna;i++) {
      ouf<<i+1<<" "<<contact_profile[t][i]<<endl;
    }
    ouf.close();
  }

}
