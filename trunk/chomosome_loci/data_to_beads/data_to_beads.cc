//############################################################################
//###
//###    Program which reads in data sets from bed files, and generates a
//###    list of beads with different types, depending on the overlaps from
//###    the bed files.
//###
//###    Usage :
//###          ./data_to_beads.out config_file o_dir
//###    where o_dir is the directory where output files are written, and 
//###    config_file has the following lines (without comments) :
//###
//###              chr5                           # chromosome of region
//###              28750000                       # start of region (bp)
//###              29750000                       # end of region (bp)
//###              400                            # number of bp per bead
//###              name path/to/bedfile.bed       # input data
//###
//###    where multiple input data sets with unique names can be listed.
//###    
//###    Notes :
//###            * bed files start counting from 0 and end is not inclusive
//###            * internally beads are counted from 0, but from 1 when output
//###
//############################################################################

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>

#include "data_to_beads.h"

#define COVERFRACTION 0.25

using namespace std;

void help_message() {
  cout<<" ############################################################################"<<endl;
cout<<" ###"<<endl;
cout<<" ###    Program which reads in data sets from bed files, and generates a"<<endl;
cout<<" ###    list of beads with different types, depending on the overlaps from"<<endl;
cout<<" ###    the bed files."<<endl;
cout<<" ###"<<endl;
cout<<" ###    Usage :"<<endl;
cout<<" ###          ./data_to_beads.out config_file o_dir"<<endl;
cout<<" ###    where o_dir is the directory where output files are written, and "<<endl;
cout<<" ###    config_file has the following lines (without comments) :"<<endl;
cout<<" ###"<<endl;
cout<<" ###              chr5                           # chromosome of region"<<endl;
cout<<" ###              28750000                       # start of region (bp)"<<endl;
cout<<" ###              29750000                       # end of region (bp)"<<endl;
cout<<" ###              400                            # number of bp per bead"<<endl;
cout<<" ###              name path/to/bedfile.bed       # input data"<<endl;
cout<<" ###"<<endl;
cout<<" ###    where multiple input data sets with unique names can be listed."<<endl;
cout<<" ###    "<<endl;
cout<<" ###    Notes :"<<endl;
cout<<" ###            * bed files start counting from 0 and end is not inclusive"<<endl;
cout<<" ###            * internally beads are counted from 0, but from 1 when output"<<endl;
cout<<" ###"<<endl;
cout<<" ############################################################################"<<endl;
}



int main(int argc, char *argv[]) {

  region theregion(COVERFRACTION);
  map<string,bed> thebeds;
  map<string,vector<int> > beads;
  vector<int> bead_types;

  // get options from command line
  if (argc!=3) {
    help_message();
    exit(0);
  }

  // read config file
  if (read_config(argv[1],theregion,thebeds)==0) {
    cout<<"Exiting due to read_config error."<<endl;
    exit(0);
  }

  // parse bed files, and set bead types
  for (map<string,bed>::iterator it=thebeds.begin(); it!=thebeds.end(); ++it) {
    beads[it->second.name]=vector<int>(theregion.L_bead,0);
    load_bed(theregion,it->second,beads[it->second.name]);
  }

  // output bead list
  if (theregion.output_beadlist(beads)==0) {
    cout<<"Exiting due to output_beadlist error."<<endl;
    exit(0);
  }

  // identify bead types
  theregion.get_types(bead_types,beads);

  // output bead types
  if (theregion.output_typeslist(bead_types)==0) {
    cout<<"Exiting due to output_beadtypes error."<<endl;
    exit(0);
  }

  cout<<"Bead lists generated sucsessfully."<<endl;

}

//############################################################################
//############################################################################
//############################################################################

//###########################---- FUNCTIONS ----##############################

//####### GLOBAL FUNCTIONS ###################################################

bool read_config(const string &filename, region &reg, map<string,bed> &beds) {

  ifstream inf,
    testinf;
  string line;
  stringstream sline;

  string bedname,
    bedfile;

  inf.open( filename.c_str() );
  if ( !inf.good( )) {
    cout<<"Error opening input file "<<filename<<endl;
    return 0;
  } 

  // first lines are chromosome, start, end, resolution
  getline(inf,line); sline.str(line); sline>>reg.chrom; sline.clear();
  getline(inf,line); sline.str(line); sline>>reg.start; sline.clear();
  getline(inf,line); sline.str(line); sline>>reg.end; sline.clear();
  getline(inf,line); sline.str(line); sline>>reg.bp_per_bead; sline.clear();

  reg.L_bp=reg.end-reg.start;
  reg.L_bead=ceil( double(reg.L_bp) / double(reg.bp_per_bead) );

  if (reg.L_bp<1) {
    cout<<"Error : incorrect region specification."<<endl;
    inf.close();
    return 0;
  }

  cout<<"Setting up region..."<<endl;
  cout<<"                    "<<reg.chrom<<":"<<reg.start<<"-"<<reg.end<<endl;
  if (reg.L_bead*reg.bp_per_bead+reg.start>reg.end) {
    cout<<".... adjusting region to keep whole number of beads :"<<endl;
    reg.end=reg.L_bead*reg.bp_per_bead+reg.start;
    reg.L_bp=reg.end-reg.start;
    cout<<" new region :       "<<reg.chrom<<":"<<reg.start<<"-"<<reg.end<<endl;
  }
  cout<<"Length = "<<reg.L_bp<<" bp. Giving "<<reg.L_bead<<" beads at "<<reg.bp_per_bead<<" bp per bead"<<endl;

  // following lines are names and bed files
  cout<<endl<<"Checking bed files :"<<endl;
  while ( getline(inf, line) ) {
    sline.str(line); sline>>bedname; sline>>bedfile; sline.clear();
    cout<<"    dataset "<<bedname<<" in file "<<bedfile<<" ....";
    testinf.open( bedfile.c_str() );
    if (testinf.good()) {
      testinf.close();
      if ( beds.insert( pair<string,bed>( bedname,bed(bedname,bedfile) ) ).second==0 ) {
	cout<<" FAIL."<<endl;
	cout<<"Error : attempting to load multiple bed files with same name"<<endl;
	inf.close();
	return 0;
      }
      cout<<" OK."<<endl;
    } else {
      testinf.close();
      cout<<" FAIL."<<endl;
      inf.close();
      return 0;
    }
  }

  if ( beds.size()>0 ) {
    cout<<endl;
    inf.close();
    return 1;
  } else {
    cout<<"Error : no bed files specified"<<endl<<endl;
    inf.close();
    return 0;
  }

}

//#######
//#######

int overlap(long int s1, long int e1, long int s2, long int e2) {
  // returns the number of bp by which the two regions overlap
  if (s1>=s2 && e1<=e2) { return e1-s1; }
  else if (s2>=s1 && e2<=e1) {return e2-s2; }
  else if (e1>s2 && e1<e2) {return e1-s2; }
  else if (s1>s2 && s1<e2) {return e2-s1; }
  else {return 0; }
}

//#######
//#######

bool load_bed(const region &reg, const bed &abed, vector<int> &beads) {

  ifstream inf,
    testinf;
  string line;
  bedline bline;

  vector<int> beadlist;

  // open the bed file
  inf.open( abed.filename.c_str() );

  while ( getline(inf,line) ) {

    bline=bedline(line);
    if ( reg.in_region(bline) ) {
      if (bline.width<(1+reg.coverfrac)*reg.bp_per_bead) {
	// site only covers one bead
	beads[ reg.bp2bead(bline.mid) ]=1;
      } else {
	beadlist=reg.overlap_beads(bline);
	for (int i=0;i<beadlist.size();i++) {
	  beads[ beadlist[i] ]=1;
	}
      }
    }


  }

  inf.close();
  return 1;

}

//####### MEMBER FUNCTIONS OF struct bedline #################################

bedline::bedline() {};
bedline::bedline(string line) {
  stringstream sline;
  sline.str(line); 
  sline>>chrom>>start>>end;
  sline.clear();

  width=end-start;
  mid=floor( 0.5*double(start+end) );
}

//####### MEMBER FUNCTIONS OF struct region #################################

bool region::in_region(const bedline &bl) const {
  // returns true if bedline is within the region
  return (bl.chrom==chrom && bl.end>start && bl.start<end);
}
int region::bp2bead(const long int &bp) const {
  // convert a bp to a bead number
  return floor( double(bp-start)/double(bp_per_bead) );
}
long int region::pos(const int &bead) const {
  // returns the bp location representing the middle of the bead
  return start+bp_per_bead*bead+0.5*bp_per_bead;
}
long int region::bead_start(const int &bead) const {
  // returns the bp location representing the start of the bead
  return start+bp_per_bead*bead;
}
long int region::bead_end(const int &bead) const {
  // returns the bp location representing the end of the bead
  return start+bp_per_bead*(bead+1)-1;
}
vector<int> region::overlap_beads(const bedline &bl) const {
  // returns a list of overlapping beads
  int first=bp2bead(bl.start),
    last=bp2bead(bl.end);
  vector<int> list;
  for (int i=0;i<L_bead;i++) {
    if ( overlap(bead_start(i),bead_end(i),bl.start,bl.end)>coverfrac*bp_per_bead ) {
      list.push_back(i);
    }
  }
  return list;
}

//#######
//#######

bool region::output_beadlist(const map<string,vector<int> > &beads) {
  // output a list of properties for each bead
  ifstream testinf;
  ofstream ouf;

  map<string,vector<int> >::const_iterator it;

  cout<<"Writing file bead_list.dat"<<endl;
  // test file doesn't exist
  testinf.open( "bead_list.dat" );
  if (testinf.good()) {
    cout<<"Error : bead_list.dat already exists. Exiting."<<endl;
    testinf.close();  
    return 0;
  } else {
    testinf.close();  
  }

  // write file
  ouf.open( "bead_list.dat" );
  ouf<<"# list of bead properties for region "<<chrom<<":"<<start<<"-"<<end<<endl;
  ouf<<"# bead id, bp midpoint";
  for (it=beads.begin();it!=beads.end();++it) {
    ouf<<", "<<it->first;
  }
  ouf<<endl;
  for (int i=0;i<L_bead;i++) {
    ouf<<i+1<<" "<<pos(i);
    for (it=beads.begin();it!=beads.end(); ++it) {
      ouf<<" "<<it->second[i]; 
    }
    ouf<<endl;
  }
  ouf.close();

  return 1;

}

//#######
//#######

bool region::output_typeslist(const vector<int> &bead_types) {
  // output a list of each bead type
  ifstream testinf;
  ofstream ouf;

  cout<<"Writing file bead_types.dat"<<endl;
  // test file doesn't exist
  testinf.open( "bead_types.dat" );
  if (testinf.good()) {
    cout<<"Error : bead_types.dat already exists. Exiting."<<endl;
    testinf.close();  
    return 0;
  } else {
    testinf.close();  
  }

  // write file
  ouf.open( "bead_types.dat" );
  ouf<<"# List of beads with "<<type_list.size()<<" different types"<<endl;
  for (map<string,int>::const_iterator it=type_list.begin(); it!=type_list.end(); ++it) {
    ouf<<"# type "<<it->second<<" : "<<type_description[ it->first ]<<endl;
  }
  ouf<<"#"<<endl<<"# bead id, type"<<endl;

  for (int i=0;i<L_bead;i++) {
    ouf<<i+1<<" "<<bead_types[i]<<endl;
  }
  ouf.close();

  return 1;
}

//#######
//#######

void region::get_types(vector<int> &bead_types, const map<string,vector<int> > &beads) {
  // identify types of bead, and make a list of what type each bead is

  stringstream typestring,
    descript;
  int typecount=1;

  type_list.clear();
  typestring.clear();
  descript.clear();

  // set up bead_types
  bead_types=vector<int>(L_bead,0);

  // set 0 type
  for (map<string,vector<int> >::const_iterator it=beads.begin(); it!=beads.end(); ++it) {
    typestring<<"0";
    descript<<"!"<<it->first<<" ";
  }
    

  type_list[ typestring.str() ]=typecount;
  type_description[ typestring.str() ]=descript.str();
  typecount++;

  // go through beads finding types
  for (int i=0;i<L_bead;i++) {
    typestring.str(std::string()); typestring.clear();
    descript.str(std::string()); descript.clear();
    for (map<string,vector<int> >::const_iterator it=beads.begin(); it!=beads.end(); ++it) {
      typestring<<it->second[i];
      if (it->second[i]==0) {descript<<"!";}
      descript<<it->first<<" ";
    }
    if ( type_list.find( typestring.str() )==type_list.end() ) {
      type_list[ typestring.str() ]=typecount;
      typecount++;
      type_description[ typestring.str() ]=descript.str();
    }
    bead_types[i]=type_list[ typestring.str() ];     
  }

}
