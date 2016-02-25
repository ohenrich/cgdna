
using namespace std;

struct myenums {
  int  BEND,
    DNADNA,
    DNA, DPAT, PCORE, PPAT, PPAT2;
  myenums() {
    BEND=1;
    DNADNA=1;
    DNA=1; PCORE=2; PPAT=3; PPAT2=4;
  }
} TYPE;

struct bond {
bond(int aa, int bb, int t) : a(aa), b(bb), type(t) {};
  int a,b,
    type;
};

struct angle {
angle(int aa, int bb, int cc, int t) : a(aa), b(bb), c(cc), type(t) {};
  int a,b,c,
    type;
};

struct evec {
  double ei,ej,ek;
  evec cross(evec a) {
    evec c;
    c.ei=ej*a.ek - ek*a.ej;
    c.ej=ek*a.ei - ei*a.ek;
    c.ek=ei*a.ej - ej*a.ei;
    return c;
  }
  double length() {
    return sqrt(ei*ei + ej*ej + ek*ek);
  }
  void make_unit() { // makes it a unit vector
    double l;
    l=length();
    ei/=l;
    ej/=l;
    ek/=l;
  }
};

struct atom {
  atom(double a,double b, double c) : x(a), y(b), z(c) {};
  atom() {}
  double x,y,z;
  int id,
    type, 
    mol;
};
