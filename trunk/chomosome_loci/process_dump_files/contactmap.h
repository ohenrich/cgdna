
#ifndef CONTACTMAP_H
#define CONTACTMAP_H

using namespace std;

struct simluationsystem;

struct cmap {

  double **themap;
  int N;

  cmap(const int&);
  ~cmap();
  void operator+=(const cmap&);
  void operator/=(const int&);
  void operator=(const cmap &B);
  void clear();

  void load_cmap(const string &);
  void save_cmap(const string &);

};

void generate_contactmap(int, char*[], double);

void generate_bincontactmap(int, char *[], double);

cmap load_contact(const simluationsystem &, ifstream &, const double &);

#endif
