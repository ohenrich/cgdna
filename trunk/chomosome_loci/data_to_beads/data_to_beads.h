#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>

using namespace std;

void help_message();
int overlap(long int, long int, long int, long int);

//#######

struct bed {
  string name,
    filename;
  bed(string n,string fn) : name(n), filename(fn) {}
};

//#######

struct bedline {
  string chrom;
  long int start,
    end,
    mid;
  int width;

  bedline();
  bedline(string);

};

//#######

struct region {

  string chrom;
  long int start,
    end;
  int bp_per_bead,
    L_bp,
    L_bead;
  double coverfrac;

  map<string,int> type_list;
  map<string,string> type_description;

  region(double cf) : coverfrac(cf) {};

  bool in_region(const bedline&) const;
  int bp2bead(const long int&) const;
  long int pos(const int&) const;
  long int bead_start(const int&) const;
  long int bead_end(const int&) const;
  vector<int> overlap_beads(const bedline&) const;

  bool output_beadlist(const map<string,vector<int> >&);
  bool output_typeslist(const vector<int>&);

  void get_types(vector<int>&, const map<string,vector<int> >&);

};

//#######

bool read_config(const string&,region&,map<string,bed>&); 
bool load_bed(const region&,const bed&, vector<int>&);

