//############################################################################
//### Find chi^2 for two contact maps with the same number of bins/beads
//############################################################################

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

#include "general.h"
#include "comparemaps.h"
#include "contactmap.h"

using namespace std;

void compare_maps(int argc, char *argv[],randoms &rnd)  { 

  if (argc!=7) {
    cout<<"### Command COMP2MAPS ######################################"<<endl;
    cout<<"Usage : "<<endl;
    cout<<"      ./process_dumps.out COMP2MAPS N skip infile1 infile2 seed"<<endl;
    cout<<"where N is the number of chromatin beads/bins in the maps (must be the same for both)"<<endl;
    cout<<"      skip is the number of beads/bins along the diagonal to skip,"<<endl;
    cout<<"      infile1 and infile2 are the file names of the two maps to be compared"<<endl;
    cout<<"      seed is used for the random number generator to give shuffled map comparison"<<endl;
    exit(0);
  }

  string filename1=string( argv[4] ),
    filename2=string( argv[5] );

  int N=atoi(argv[2]),
    skip=atoi(argv[3]),
    seed=atoi(argv[6]);

  cmap map1(N),
    map2(N),
    shuffmap1(N),
    shuffmap2(N);

  double chi2,
    chi2_s1,
    chi2_s2,
    chi2_ss;

  // initialize randoms
  rnd.initialize(seed);

  // Load files
  map1.load_cmap(filename1);
  map2.load_cmap(filename2);

  //calculate chi^2 for map1 and map2
  chi2=get_chi2(map1,map2,skip);

  //create shuffled maps
  cmap *map_ptr[2];
  map_ptr[0] = &map1;
  map_ptr[1] = &map2;
  shuffmap1=get_shuffled_cmap(2,map_ptr,rnd);
  shuffmap2=get_shuffled_cmap(2,map_ptr,rnd);

  //calculate chi^2 between map1/2 and a shuffled map
  chi2_s1=get_chi2(map1,shuffmap1,skip);
  chi2_s2=get_chi2(map2,shuffmap1,skip);
  chi2_ss=get_chi2(shuffmap1,shuffmap2,skip);

  // output
  cout<<endl<<"*******************************************************"<<endl<<endl;
  cout<<"chi^2 between map1 and map2 = "<<chi2<<endl;
  cout<<"chi^2 between map1 and a shuffled map = "<<chi2_s1<<endl;
  cout<<"chi^2 between map2 and a shuffled map = "<<chi2_s2<<endl;
  cout<<"chi^2 between two shuffled maps = "<<chi2_ss<<endl;
  cout<<endl<<"ratio between chi^2 for the two maps, and chi^2 for the shuffled case = "<<chi2*3.0/(chi2_s1+chi2_s2+chi2_ss)<<endl;
  cout<<"(a value close to 1 means the maps are as different as two shuhffled maps)"<<endl;
  cout<<endl<<"*******************************************************"<<endl<<endl;

}


void shuffle_maps(int argc, char *argv[],randoms &rnd) {


  if (argc<6) {
    cout<<"### Command COMP2MAPS ######################################"<<endl;
    cout<<"Usage : "<<endl;
    cout<<"      ./process_dumps.out SHUFFMAP N outfile seed list of infiles"<<endl;
    cout<<"where N is the number of chromatin beads/bins in the maps (must be the same for both)"<<endl;
    cout<<"      outfile is the output file for the shuffled map,"<<endl;
    cout<<"      seed is used for the random number generator to give shuffled map comparison"<<endl;
    cout<<"      then include a list of at least 1 contact map"<<endl;
    cout<<"All of the listed maps are read in, and shuffled randomly to generate a single output map."<<endl;
    exit(0);
  }

  string outfile=string( argv[3] ),
    *infilenames;

  int N=atoi(argv[2]),
    seed=atoi(argv[4]),
    nfiles=argc-5;

  cmap *map_ptr[nfiles],
    shuffled(N);

  // set up randoms
  rnd.initialize(seed);

  // get file names
  infilenames = new string[nfiles];
  for (int i=0;i<nfiles;i++) {
    infilenames[i]=string( argv[i+5] );
  }

  // set up pointers and load the maps
  for (int i=0;i<nfiles;i++) {
    map_ptr[i] = new cmap(N);
    map_ptr[i]->load_cmap( infilenames[i] );
  }

  // generate a shuffled map
  shuffled=get_shuffled_cmap(nfiles,map_ptr,rnd);

  // output the map
  shuffled.save_cmap(outfile);

  for(int i=0;i<nfiles;++i) {
    delete map_ptr[i];
  }

}




double get_chi2(const cmap &A,const cmap &B,const int &skip) {

  double chi2;
  int count;

  if (A.N != B.N ) {
    cout<<"Error : cannot compare maps of different size"<<endl;
  }

  chi2=0.0;
  count=0;

  for (int i=0;i<A.N;i++) {
    for (int j=i;j<A.N;j++) {
      if (j-i>skip) { // skip a region near the diagonal
	chi2 += pow( A.themap[i][j] - B.themap[i][j] , 2.0);
      }
    }
  }

  chi2 /= double(A.N*(A.N - skip - 1));

  return chi2;

}


cmap get_shuffled_cmap(const int &n_maps, cmap *map[],randoms &rnd) {

  // check all maps have same size
  for (int i=0;i<n_maps;i++) {
    for (int j=i+1;j<n_maps;j++) {
      if ( map[i]->N != map[j]->N ) {
	cout<<"Error : cannot make shuffled map from maps of different size."<<endl;
	exit(0);
      }
    }
  }
 
  if ( n_maps<=0 ) {
    cout<<"Error : cannot make shuffled map from 0 maps."<<endl;
    exit(0);
  }

  int N=map[0]->N;
  cmap shuff(N);

  int oi,oj,
    pick_i,pick_j,
    pick_map,pick_k;

  for (int s=0;s<N;s++) { // loop through separations

    for (int k=0;k<N-s;k++) { 
      // swap each entry for a random one of 
      // the same separation from a random map
      oi=k;
      oj=k+s;
      pick_map=rnd.integer(0,n_maps-1);
      if (N-s-1>0) {
	pick_k=rnd.integer(0,N-s-1);
      } else {
	pick_k=0;
      }
      pick_i=pick_k;
      pick_j=pick_k+s;
      shuff.themap[oi][oj]=map[pick_map]->themap[pick_i][pick_j];
    }

  }

  // make symetric
  for (int j=0;j<N;j++) {
    for (int i=j+1;i<N;i++) {
      shuff.themap[i][j]=shuff.themap[j][i];
    }
  }

  return shuff;

}
