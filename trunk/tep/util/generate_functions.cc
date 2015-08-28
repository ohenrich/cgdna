void output_tors(ostream& ouf,vector<atom>& atoms) {

    ouf<<endl<<" Atoms"<<endl<<endl;
    for (int i=0;i<atoms.size();i++) {
      ouf<<" "<<atoms[i].id<<" "<<atoms[i].type<<" "<<atoms[i].x<<" "<<atoms[i].y<<" "<<atoms[i].z<<" "<<atoms[i].mol<<" "<<atoms[i].ellipse_flag<<" "<<atoms[i].density<<endl;
    }
    ouf<<endl<<" Ellipsoids"<<endl<<endl;
    for (int i=0;i<atoms.size();i++) {
      if (atoms[i].ellipse_flag==1) {
	ouf<<" "<<atoms[i].id<<" 1 1 1 "<<atoms[i].q[0]<<" "<<atoms[i].q[1]<<" "<<atoms[i].q[2]<<" "<<atoms[i].q[3]<<endl; 
      }
    }
    ouf<<endl<<" Velocities"<<endl<<endl;
    for (int i=0;i<atoms.size();i++) {
      ouf<<" "<<atoms[i].id<<" 0 0 0 0 0 0"<<endl;
    }

}


void output_notors(ostream& ouf,vector<atom>& atoms) {


    ouf<<endl<<" Atoms"<<endl<<endl;
    for (int i=0;i<atoms.size();i++) {
      ouf<<" "<<atoms[i].id<<" "<<atoms[i].mol<<" "<<atoms[i].type<<" "<<atoms[i].x<<" "<<atoms[i].y<<" "<<atoms[i].z<<"  0 0 0"<<endl;
    }

    ouf<<endl<<" Velocities"<<endl<<endl;
    for (int i=0;i<atoms.size();i++) {
      ouf<<" "<<atoms[i].id<<" 0 0 0"<<endl;
    }

}

 
void add_DNA(atom last,atom lastlast,vector<atom> &atoms,vector<bond> &bonds,vector<angle> &angles, double box[3], double de, double density) {
  // Add a DNA bead to the random configuration

  double theta,phi,
    dx,dy,dz,
    hbox[3],id;

  for (int i=0;i<3;i++) {
    hbox[i]=box[i]*0.5;
  }

  // position
  do {
    theta=double(rand())/double(RAND_MAX)*PI;
    phi=double(rand())/double(RAND_MAX)*2.0*PI;
    dx=last.x+sin(theta)*cos(phi);
    dy=last.y+sin(theta)*sin(phi);
    dz=last.z+cos(theta);
  } while (abs(dx)>hbox[0]||abs(dy)>hbox[1]||abs(dz)>hbox[2]); // reject if outside box

  id=atoms.back().id;

  id++;
  atoms.push_back( atom(dx,dy,dz) );
  atoms.back().type=TYPE.DNA;
  atoms.back().id=id;
  atoms.back().mol=last.mol;
  atoms.back().ellipse_flag=1;
  atoms.back().density=density;   

  // oreintation
  atoms.back().q[0]=1; atoms.back().q[1]=0; atoms.back().q[2]=0; atoms.back().q[3]=0;

  // bond
  bonds.push_back( bond(last.id,atoms.back().id,TYPE.DNADNA) );

  // angle
  if (lastlast.id==0) { // this is bead number 2
    angles.push_back( angle(last.id,atoms.back().id,atoms.back().id,TYPE.TORS) ); // the third id here does nothing
  } else { 
    angles.push_back( angle(lastlast.id,last.id,atoms.back().id,TYPE.BEND) );
    angles.push_back( angle(last.id,atoms.back().id,atoms.back().id,TYPE.TORS) ); // the third id here does nothing   
  } 
   

}



void add_DNA(atom last,atom lastlast,vector<atom> &atoms,vector<bond> &bonds,vector<angle> &angles, double density, double x,double y,double z) {
  // Add a DNA bead to the loop configuration

  atoms.push_back( atom(x,y,z) );
  atoms.back().type=TYPE.DNA;
  atoms.back().id=last.id+1;
  atoms.back().mol=last.mol;
  atoms.back().ellipse_flag=1;
  atoms.back().density=density;  

  // oreintation
  atoms.back().q[0]=1; atoms.back().q[1]=0; atoms.back().q[2]=0; atoms.back().q[3]=0;

  // bond
  bonds.push_back( bond(last.id,atoms.back().id,TYPE.DNADNA) );

  // angle
  if (lastlast.id==0) { // this is bead number 2
    angles.push_back( angle(last.id,atoms.back().id,atoms.back().id,TYPE.TORS) ); // the third id here does nothing
  } else { 
    angles.push_back( angle(lastlast.id,last.id,atoms.back().id,TYPE.BEND) );
    angles.push_back( angle(last.id,atoms.back().id,atoms.back().id,TYPE.TORS) ); // the third id here does nothing   
  } 

}
