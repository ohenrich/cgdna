//############################################################################
//###
//###    Program which processes dump files in a number of different ways.
//###    
//###    Usage:
//###          ./process_dumps.out COMMAND options
//###
//###    Where command is one of the following :
//###               CONTACTMAP                 generate a contact map
//###               BINCONTACTMAP              generate a binned contact map
//###               COMP2MAPS                  compare two contact maps with same number of beads
//###               SHUFFMAP                   reads contact maps and generates a shuffled map for comparisons
//###               RGT                        calculate radius of gyration as a function of time
//###               RGR                        calculate radius of gyration for a specified region
//###               SIMULATEDFISH              read a bed file of probes, and extract simulated FISH data
//###               SIMULATEDCAPC              read a bed file of capc targets, and extract contact profiles
//###
//###
//############################################################################

#include <iostream>
#include <cstdlib>

#include "general.h"
#include "process_dumps.h"
#include "contactmap.h"
#include "comparemaps.h"
#include "rg.h"
#include "pair_distance.h"
#include "countbound.h"
#include "simFISH.h"
#include "capturec.h"

using namespace std;

const double THRESHOLD=2.75;
const double BOUNDTHRESH=1.45;

void help_message() {
  cout<<"//############################################################################"<<endl;
  cout<<"//###"<<endl;
  cout<<"//###    Program which processes dump files in a number of different ways."<<endl;
  cout<<"//###    "<<endl;
  cout<<"//###    Usage:"<<endl;
  cout<<"//###          ./process_dumps.out COMMAND options"<<endl;
  cout<<"//###"<<endl;
  cout<<"//###    Where command is one of the following :"<<endl;
  cout<<"//###               CONTACTMAP                 generate a contact map"<<endl;
  cout<<"//###               BINCONTACTMAP              generate a binned contact map"<<endl;
  cout<<"//###               COMP2MAPS                  compare two contact maps with same number of beads"<<endl;
  cout<<"//###               SHUFFMAP                   reads contact maps and generates a shuffled map for comparisons"<<endl;
  cout<<"//###               COUNTBOUND                 counts the number of proteins of specified type bound to the polymer"<<endl;
  cout<<"//###               PAIRDISTANCE               calculate the distance between all possible pairs of configurations"<<endl;
  cout<<"//###               RGT                        calculate radius of gyration as a function of time"<<endl;
  cout<<"//###               RGR                        calculate radius of gyration for a specified region"<<endl;
  cout<<"//###               RS                         calculate R(s) the separation of beads separated by s"<<endl;
  cout<<"//###               SIMULATEDFISH              read a bed file of probes, and extract simulated FISH data"<<endl;
  cout<<"//###               SIMULATEDCAPC              read a bed file of capc targets, and extract contact profiles"<<endl;
  cout<<"//###"<<endl;
  cout<<"//###"<<endl;
  cout<<"//###    To get help for a specific command, run with no options"<<endl;
  cout<<"//###"<<endl;
  cout<<"//###    Dump files must have same number of atoms in each frame, and polymer beads must be atom numbers 1-N."<<endl;
  cout<<"//###    Some commands require a configuration file, as used as input for the accompanying data_to_beads tool."<<endl;
  cout<<"//###    In most cases will not throw an error if input is in incorrect format, but will give nonsence output."<<endl;
  cout<<"//###"<<endl;
  cout<<"//############################################################################" <<endl;
}

int main(int argc, char *argv[]) {

  // get options from command line
  if (argc<2) {
    help_message();
    exit(0);
  }

  string command=string( argv[1] );
  randoms rnd;
  int seed;

  if ( command=="CONTACTMAP" ) {
    generate_contactmap(argc,argv,THRESHOLD);

  } else if ( command=="BINCONTACTMAP" ) {
    generate_bincontactmap(argc,argv,THRESHOLD);

  } else if ( command=="COMP2MAPS" ) {
    compare_maps(argc,argv,rnd); 

  } else if ( command=="SHUFFMAP" ) {
    shuffle_maps(argc,argv,rnd); 

  } else if ( command=="RGT" ) {
    gyration_radius_t(argc,argv);

  } else if ( command=="RGR" ) {
    gyration_radius_r(argc,argv);

  } else if ( command=="COUNTBOUND" ) {
    count_bound(argc,argv,BOUNDTHRESH);

 } else if ( command=="PAIRDISTANCE" ) {
    pairwise_dist(argc,argv,THRESHOLD); 

 } else if ( command=="SIMULATEDFISH" ) {
    get_fish(argc,argv); 

 } else if ( command=="SIMULATEDCAPC" ) {
    capture_c(argc,argv,THRESHOLD); 

  } else if ( command=="RS" ) {
    cout<<"Oops. Haven't implemented this yet."<<endl;
    exit(0);

  } else {
    cout<<"Command "<<argv[1]<<" not recognized"<<endl;
    help_message();
    exit(0);
  }


}
