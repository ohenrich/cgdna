// c++ code for functions which are used for data analasis of 


int file_length( char *fn, int N ) {
  // returns the number of frames in a file

  FILE *pee;
  char command[200];
  int flen,nframes;

  sprintf(command,"wc -l %s | awk '{ print $1 }'",fn);
  pee=popen(command,"r");
  fscanf(pee,"%u",&flen);
  pclose(pee);
  flen=flen;
  nframes=flen/(N+9);

  return nframes;

}

double sign(double a) { 
  // returns the sign of the argument, or -1 if a=0
  if (a>0) {return 1.0;} else {return -1.0;} 
}

void fix_roundoff11(double &a) {
  // checks for round off in a number in the interval [-1,1]
  // before using asin or acos

  if (abs(a)>1.001) {
    std::cout<<"Error - number should be in interval [-1,1]"<<std::endl;
    std::exit(0);
  }
  if (a>1) {a=1.0;}
  if (a<-1) {a=-1.0;}
}


// ****************************************************
// member functions of quaternion struct
// ****************************************************

xyz quaternion::xaxis() {
  // gives the unit vector corresponding to the x-axis of the bead
  xyz p;
  p.x = q0*q0 + q1*q1 - q2*q2 - q3*q3;
  p.y = 2.0*q1*q2 + 2.0*q0*q3;
  p.z = 2.0*q1*q3 - 2.0*q0*q2;
  return p;
}

xyz quaternion::yaxis() {
  // gives the unit vector corresponding to the y-axis of the bead
  xyz p;
  p.x = 2.0*q1*q2 - 2.0*q0*q3;
  p.y = q0*q0 - q1*q1 + q2*q2 - q3*q3;  
  p.z = 2.0*q2*q3 + 2.0*q0*q1;
  return p;
}

xyz quaternion::zaxis() {
  // gives the unit vector corresponding to the z-axis of the bead
  xyz p;
  p.x = 2.0*q1*q3 + 2.0*q0*q2;
  p.y = 2.0*q2*q3 - 2.0*q0*q1;
  p.z = q0*q0 - q1*q1 - q2*q2 + q3*q3;
  return p;
}

double quaternion::alpha() {
  // gives the alpha Euler angle from the quaternion
  return atan2( q0*q2+q1*q3 , q0*q3-q1*q2 );
}

double quaternion::gamma() {
  // gives the gamma Euler angle from the quaternion
  return atan2( q0*q2-q1*q3 , q0*q3+q1*q2 );
}

quaternion quaternion::mult(quaternion b) {
  // multiply quaterion a by b to get c
  quaternion c;
  c.q0 = q0*b.q0 - q1*b.q1 - q2*b.q2 - q3*b.q3;
  c.q1 = q0*b.q1 + q1*b.q0 + q2*b.q3 - q3*b.q2;
  c.q2 = q0*b.q2 - q1*b.q3 + q2*b.q0 + q3*b.q1;
  c.q3 = q0*b.q3 + q1*b.q2 - q2*b.q1 + q3*b.q0;
  return c;
}

double quaternion::norm() {
  // returns the norm of the quaternion
  return sqrt( q0*q0 + q1*q1 + q2*q2 + q3*q3 );
}

void quaternion::make_unit() {
  // make it a unit quaternion
  double l=norm();
  q0/=l;
  q1/=l;
  q2/=l;
  q3/=l;
}

void quaternion::rotate(xyz axis, double angle) {
  // rotate the quaternion by angle about axis

  quaternion b,c;
  double sinhangle;
  
  // generate a quaternion for the rotation
  b.q0 = cos( 0.5*angle );
  sinhangle = sin( 0.5*angle );
  b.q1 = sinhangle * axis.x;
  b.q2 = sinhangle * axis.y;
  b.q3 = sinhangle * axis.z;

  // multiply the two quaternions
  c = b.mult(*this);

  // make sure it is still a unit quaternion
  c.make_unit();
  *this = c;

}

// ****************************************************
// member functions of xyz struct
// ****************************************************

xyz xyz::operator- (xyz b) {  // subtract vectors
  xyz c;
  c.x=x-b.x;
  c.y=y-b.y;
  c.z=z-b.z;
  return c;
}

double xyz::dot(xyz b) { // dot product of vectors
  return (x*b.x + y*b.y + z*b.z);
}

xyz xyz::cross(xyz b) { // cross product of vectors
  xyz c;
  c.x = y*b.z - z*b.y;
  c.y = z*b.x - x*b.z;
  c.z = x*b.y - y*b.x;
  return c;
}

double xyz::length() { // return length of vector
  return sqrt(x*x + y*y + z*z);
}

void xyz::make_unit() { // makes it a unit vector
  double l;
  l=length();
  x/=l;
  y/=l;
  z/=l;
}
