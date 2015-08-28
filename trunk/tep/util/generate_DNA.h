
using namespace std;

struct myenums {
  int  BEND, TORS, TORSEND,
    DNADNA,
    DNA, DPAT, NCORE, NPAT1, NPAT2, PCORE, PPAT;
  myenums() {
    BEND=1; TORS=2; TORSEND=3;
    DNADNA=1;
    DNA=1; DPAT=2; NCORE=3; NPAT1=4; NPAT2=5; PCORE=6; PPAT=7;
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

struct quaternion {
  double q0,q1,q2,q3;
  double norm();
  void make_quat(evec, evec, evec);
  evec xaxis();
private : 
  double sign(double);
};

struct atom {
  atom(double a,double b, double c) : x(a), y(b), z(c) {};
  atom() {}
  double x,y,z,density;
  double q[4];
  int id,
    type, 
    mol,
    ellipse_flag;
  void quat(quaternion quat) {
    q[0]=quat.q0;
    q[1]=quat.q1;
    q[2]=quat.q2;
    q[3]=quat.q3;
  }
};

void output_tors(ostream&,vector<atom>&);
void output_notors(ostream&,vector<atom>&);

void do_dna_randlin(vector<atom>&,vector<bond>&,vector<angle>&);

double quaternion::norm() {

  return sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);

}

double quaternion::sign(double x) {return (x >= 0.0) ? +1.0 : -1.0;}

void quaternion::make_quat(evec xax, evec yax, evec zax) {

  double r11,r12,r13,
    r21,r22,r23,
    r31,r32,r33,
    r;

  r11=xax.ei; r21=xax.ej; r31=xax.ek;
  r12=yax.ei; r22=yax.ej; r32=yax.ek;
  r13=zax.ei; r23=zax.ej; r33=zax.ek;

  q0=( r11 + r22 + r33 + 1.0f) / 4.0;
  q1=( r11 - r22 - r33 + 1.0f) / 4.0;
  q2=(-r11 + r22 - r33 + 1.0f) / 4.0;
  q3=(-r11 - r22 + r33 + 1.0f) / 4.0;

  if(q0 < 0.0) q0 = 0.0;
  if(q1 < 0.0) q1 = 0.0;
  if(q2 < 0.0) q2 = 0.0;
  if(q3 < 0.0) q3 = 0.0;

  q0 = sqrt(q0);
  q1 = sqrt(q1);
  q2 = sqrt(q2);
  q3 = sqrt(q3);

  if(q0 >= q1 && q0 >= q2 && q0 >= q3) {
    q0 *= +1.0f;
    q1 *= sign(r32 - r23);
    q2 *= sign(r13 - r31);
    q3 *= sign(r21 - r12);
  } else if(q1 >= q0 && q1 >= q2 && q1 >= q3) {
    q0 *= sign(r32 - r23);
    q1 *= 1.0;
    q2 *= sign(r21 + r12);
    q3 *= sign(r13 + r31);
  } else if(q2 >= q0 && q2 >= q1 && q2 >= q3) {
    q0 *= sign(r13 - r31);
    q1 *= sign(r21 + r12);
    q2 *= 1.0;
    q3 *= sign(r32 + r23);
  } else if(q3 >= q0 && q3 >= q1 && q3 >= q2) {
    q0 *= sign(r21 - r12);
    q1 *= sign(r31 + r13);
    q2 *= sign(r32 + r23);
    q3 *= 1.0;
  } else {
    cout<<"quaternion error"<<endl;;
  }
  r = norm();
  q0 /= r;
  q1 /= r;
  q2 /= r;
  q3 /= r;

}


evec quaternion::xaxis() {
  // gives the unit vector corresponding to the x-axis of the bead
  evec p;
  p.ei = q0*q0 + q1*q1 - q2*q2 - q3*q3;
  p.ej = 2.0*q1*q2 + 2.0*q0*q3;
  p.ek = 2.0*q1*q3 - 2.0*q0*q2;
  return p;
}
