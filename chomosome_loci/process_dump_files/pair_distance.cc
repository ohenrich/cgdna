//############################################################################
//### Calculate the distance between pairs of conformations
//############################################################################

///////////////////////////////////////////////////////////
//
//   diff = 2/(n(n-1)) * 
//      sum_{i>j} [ 1 - kroneker_delta( s_ij^c , s_ij^c' ) ] * ( d_ij^c - d_ij^c )**2
//
//   where s_ij is 1 if bead i and j of configuration P are in contact
//                 0 otherwise
//
///////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
//#include <map>

#include "pair_distance.h"
#include "general.h"

void pairwise_dist(int argc, char *argv[],const double &threshold) {

  if ( !(argc==6||argc==7) ) {
    cout<<"### Command PAIRDISTANCE ######################################"<<endl;
    cout<<"Usage : "<<endl;
    cout<<"      ./process_dumps.out PAIRDISTANCE N Nfiles infile outfile [bead_list.dat]"<<endl;
    cout<<"where N is the number of chromatin beads,"<<endl;
    cout<<"      Nfiles is the number of files to read,"<<endl;
    cout<<"      infile is the dumpfile name in the format dumpfile.999.suffix"<<endl;
    cout<<"      and outfile is the output filename"<<endl;
    cout<<"      and the optional file bead_list.dat has a list of beads to include in the distance calculation"<<endl;
    exit(0);
  }

  cout<<"******** WARNING : this command is not fully tested ***********"<<endl;

  string fullfilename=string( argv[4] ),
    filename1=fullfilename.substr(0,(fullfilename.substr(0,fullfilename.find_last_of("."))).find_last_of(".")+1),
    filename2=fullfilename.substr(fullfilename.find_last_of(".")),
    beadfile,
    outfile=string( argv[5] );

  bool beadlist_flag=0;
  // check if bead_list file is present
  if (argc==7) {
    beadlist_flag=1;
    beadfile=string( argv[6] );
  }

  // Check inputs are integers
  for (int i=0;i<strlen(argv[2]);i++) {
    if (!isdigit(argv[2][i])) {cout<<"ERROR : number of DNA beads must be an integer."<<endl; exit(0);}
  }
  for (int i=0;i<strlen(argv[3]);i++) {
    if (!isdigit(argv[3][i])) {cout<<"ERROR : number of files must be an integer."<<endl; exit(0);}
  }

  int Nfiles=atoi(argv[3]),
    Ndna=atoi(argv[2]);

  ifstream inf1, inf2;
  ofstream ouf;

  atom *conf1,
    *conf2;  

  vector<int> beadlist;
  double dist;

  // check output file does not exist
  inf1.open( outfile.c_str() );
  if ( inf1.good() ) {
    cout<<" ERROR : output file "<<outfile<<" aleady exists"<<endl;
    exit(0);
  }
  inf1.close();

  // check if bead_list file is present
  if ( beadlist_flag ) {
    inf1.open( beadfile.c_str() );
    if ( !inf1.good() ) {
      cout<<" ERROR : cannot open file "<<beadfile<<endl;
      exit(0);
    }
    string line;
    stringstream sline;
    int bead;
    getline(inf1,line);
    do {

      for (int i=0;i<line.size();i++) {
	if (!isdigit(line[i])) {cout<<"ERROR : non-integer entry in file "<<beadfile<<endl; exit(0);}
      }
      sline.clear();sline.str(line);
      sline>>bead;
      beadlist.push_back( bead );

    } while ( getline(inf1,line) );
    inf1.close();

  }

  // initialize the simulated system
  simluationsystem thesystem(filename1+to_string(1)+filename2,Ndna);
  conf1 = new atom[thesystem.Ndna];
  conf2 = new atom[thesystem.Ndna];

  // set up output files
  ouf.open( outfile.c_str() );
  ouf<<"# i, j, distance"<<endl;

  // loop round pairs of files
  for (int file1=0;file1<Nfiles;file1++) {

    // read last frame of file1
    inf1.open( (filename1+to_string(file1+1)+filename2).c_str() );
    if ( !inf1.good() ) {
      cout<<" ERROR while trying to open file "<<filename1+to_string(file1+1)+filename2<<endl;
      exit(0);
    }
    fastforward_frames(thesystem,inf1);
    load_conf_unwrapped(conf1,thesystem,inf1);
    inf1.close();

    for (int file2=file1;file2<Nfiles;file2++) {
      cout<<"Comparing files "<<file1+1<<" and "<<file2+1<<"\r"<<flush;
 
      // read last frame of file2
      inf2.open( (filename1+to_string(file2+1)+filename2).c_str() );
      if ( !inf2.good() ) {
	cout<<" ERROR while trying to open file "<<filename1+to_string(file2+1)+filename2<<endl;
	exit(0);
      }
      fastforward_frames(thesystem,inf2);
      load_conf_unwrapped(conf2,thesystem,inf2);
      inf2.close();
 
      // get distance between the two
      if (!beadlist_flag) {
	dist = get_distance(thesystem,conf1,conf2,threshold);
      } else {
	dist = get_distance(thesystem,conf1,conf2,threshold,beadlist);
      }

      ouf<<file1+1<<" "<<file2+1<<" "<<dist<<endl;

    } 

  }

  ouf.close();
  delete [] conf1;
  delete [] conf2;

}





//### Functions

double get_distance(const simluationsystem &thesystem, const atom conf1[], const atom conf2[], const double &thresh) {

  double diff;

  vector< vector<double> > dist1( thesystem.Ndna, vector<double>( thesystem.Ndna,0.0 ) ),
    dist2( thesystem.Ndna, vector<double>( thesystem.Ndna,0.0 ) );
  vector< vector<int> > contact1( thesystem.Ndna, vector<int>( thesystem.Ndna,0 ) ), 
    contact2( thesystem.Ndna, vector<int>( thesystem.Ndna,0 ) );
 
  for (int i=0;i<thesystem.Ndna;i++) {
    for (int j=i;j<thesystem.Ndna;j++) {

      dist1[i][j] = sqrt( 
      			 pow(conf1[i].x-conf1[j].x,2) + 
      			 pow(conf1[i].y-conf1[j].y,2) +  
      			 pow(conf1[i].z-conf1[j].z,2)  );
      dist1[j][i] = dist1[i][j];
      if (dist1[i][j]<thresh) {contact1[i][j]=1;} else {contact1[i][j]=0;}
      contact1[j][i] = contact1[i][j];

      dist2[i][j] = sqrt( 
      			 pow(conf2[i].x-conf2[j].x,2) + 
      			 pow(conf2[i].y-conf2[j].y,2) +  
      			 pow(conf2[i].z-conf2[j].z,2)  );
      dist2[j][i] = dist2[i][j];
      if (dist2[i][j]<thresh) {contact2[i][j]=1;} else {contact2[i][j]=0;}
      contact2[j][i] = contact2[i][j];

    }
  }
  
  // calculate diff
  diff = 0.0;
  for (int i=0;i<thesystem.Ndna;i++) {
    for (int j=i;j<thesystem.Ndna;j++) {
      if (contact1[i][j]!=contact2[i][j]) { diff+= pow( dist1[i][j] - dist2[i][j], 2); }
    }
  }
  diff /= thesystem.Ndna*(thesystem.Ndna-1)/2.0;
  
  return diff;

}







double get_distance(const simluationsystem &thesystem, const atom conf1[], const atom conf2[], const double &thresh,const vector<int> &beadlist) {

  double diff;

  vector< vector<double> > dist1( thesystem.Ndna, vector<double>( thesystem.Ndna,0.0 ) ),
    dist2( thesystem.Ndna, vector<double>( thesystem.Ndna,0.0 ) );
  vector< vector<int> > contact1( thesystem.Ndna, vector<int>( thesystem.Ndna,0 ) ), 
    contact2( thesystem.Ndna, vector<int>( thesystem.Ndna,0 ) );

  int i,j;

  for (int ii=0;ii<beadlist.size();ii++) {
    for (int jj=ii;jj<beadlist.size();jj++) {

      i = beadlist[ii]-1;
      j = beadlist[jj]-1;

      dist1[i][j] = sqrt( 
			 pow(conf1[i].x-conf1[j].x,2) + 
			 pow(conf1[i].y-conf1[j].y,2) +  
			 pow(conf1[i].z-conf1[j].z,2)  );
      dist1[j][i] = dist1[i][j];
      if (dist1[i][j]<thresh) {contact1[i][j]=1;} else {contact1[i][j]=0;}
      contact1[j][i] = contact1[i][j];

      dist2[i][j] = sqrt( 
			 pow(conf2[i].x-conf2[j].x,2) + 
			 pow(conf2[i].y-conf2[j].y,2) +  
			 pow(conf2[i].z-conf2[j].z,2)  );
      dist2[j][i] = dist2[i][j];
      if (dist2[i][j]<thresh) {contact2[i][j]=1;} else {contact2[i][j]=0;}
      contact2[j][i] = contact2[i][j];

    }
  }

  // calculate diff
  diff = 0.0;
  for (int ii=0;ii<beadlist.size();ii++) {
    for (int jj=ii;jj<beadlist.size();jj++) {
      i = beadlist[ii];
      j = beadlist[jj];
      if (contact1[i][j]!=contact2[i][j]) { diff+= pow( dist1[i][j] - dist2[i][j], 2); }
    }
  }
  diff /= beadlist.size()*(beadlist.size()-1)/2.0;
  
  return diff;


}
