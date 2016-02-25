
#ifndef SIMFISH_H
#define SIMFISH_H

using namespace std;

struct simluationsystem;
struct atom;

struct FISHprobe {

  string chromosome,
    name;
  long int start_bp,
    end_bp;
  int start_bead,   // remember that beads start couning from 1
    end_bead,
    length_bead;
  bool inregion;

  //  FISHprobe 
  void read_probe(const simluationsystem&, const string&);
  double get_sep(const atom[], const FISHprobe& ) const;
 
};

void get_fish(int, char*[]);

#endif
