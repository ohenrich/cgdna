
#ifndef GENERAL_H
#define GENERAL_H

using namespace std;

struct simluationsystem {
  int Ndna,
    Natoms,
    Nprot,
    frames_per_file,
    starttime,
    step;

  double lx,ly,lz;

  string chromosome;
  long int start_bp,
    end_bp;
  int bp_per_bead;
  bool configuration_loaded_flag;

  simluationsystem(const string&, const int&);
  void getconfiguration(const string&);
  double bp_to_doublebead(const int&) const;

};

struct atom {
  double x,y,z;
  int type;
  atom operator-(const atom&) const;
  void operator+=(const atom&);
  void operator/=(const double&);
  double length();
  void clear();
  double distance_wrapped(const simluationsystem&, const atom&) const;
};


struct randoms {
  bool initialized;

  randoms(void);
  void initialize(const int&);
  int integer(const int&, const int&);
};


void fastforward_frames(const simluationsystem&, ifstream&);

void load_conf_wrapped(atom*, const simluationsystem&, ifstream&);
void load_conf_unwrapped(atom*, const simluationsystem&, ifstream&);

string to_string(const int &a);





#endif
