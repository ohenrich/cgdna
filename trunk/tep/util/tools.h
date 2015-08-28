// header file containing function and structures for 
// analising dump files with ellipsoids

#define PI 3.14159265358979

struct xyz {
  // a data structure which is a 3-vector
  double x,y,z;
  xyz operator- (xyz);
  double dot(xyz);
  xyz cross(xyz);
  double length();
  void make_unit();
};

struct quaternion {
  // data structure for quaternions
  double q0,q1,q2,q3;
  xyz xaxis();
  xyz yaxis();
  xyz zaxis();
  double alpha();
  double gamma();
  double norm();
  void make_unit();
  quaternion mult(quaternion);
  void rotate(xyz,double);
};

struct atom: public xyz, public quaternion {
  // data structure for atoms which have both a 
  // position vector and an orientaion quaternion
  int type;
};






int file_length( char *, int );
double sign(double );
void fix_roundoff11(double &);
