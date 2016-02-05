/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing author: Oliver Henrich (EPCC, University of Edinburgh)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "bond_oxdna.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "atom_vec_ellipsoid.h"
#include "math_extra.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondOxdna::BondOxdna(LAMMPS *lmp) : Bond(lmp)
{

}

/* ---------------------------------------------------------------------- */

BondOxdna::~BondOxdna()
{
  if (allocated) {

    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(Delta);
    memory->destroy(r0);

  }
}

/* ----------------------------------------------------------------------
   compute function for oxDNA FENE-bond interaction
   s=sugar-phosphate backbone site, b=base site, st=stacking site
------------------------------------------------------------------------- */
void BondOxdna::compute(int eflag, int vflag)
{
  int i,j,n,type;
  double delf[3],delt[3]; // force, torque increment;;
  double delr[3],ebond,fbond;
  double rsq,Deltasq,rlogarg;
  double r,rr0,rr0sq;
  // distances COM-backbone
  double d_cs=-0.24; 
  // vectors COM-backbone in lab frame
  double r_cs_i[3],r_cs_j[3];

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;

  int *ellipsoid = atom->ellipsoid;
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  double *quat_i,ex_i[3],ey_i[3],ez_i[3];
  double *quat_j,ex_j[3],ey_j[3],ez_j[3];

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

 // loop over FENE bonds 

  for (n = 0; n < nbondlist; n++) {

    i = bondlist[n][0];
    j = bondlist[n][1];
    type = bondlist[n][2];

    quat_i=bonus[i].quat;
    MathExtra::q_to_exyz(quat_i,ex_i,ey_i,ez_i);
    quat_j=bonus[j].quat;
    MathExtra::q_to_exyz(quat_j,ex_j,ey_j,ez_j);

    r_cs_i[0] = d_cs*ex_i[0];
    r_cs_i[1] = d_cs*ex_i[1];
    r_cs_i[2] = d_cs*ex_i[2];
    r_cs_j[0] = d_cs*ex_j[0];
    r_cs_j[1] = d_cs*ex_j[1];
    r_cs_j[2] = d_cs*ex_j[2];

    delr[0] = x[i][0] + r_cs_i[0] - x[j][0] - r_cs_j[0];
    delr[1] = x[i][1] + r_cs_i[1] - x[j][1] - r_cs_j[1];
    delr[2] = x[i][2] + r_cs_i[2] - x[j][2] - r_cs_j[2];
    rsq = delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2];

    r = sqrt(rsq);
    rr0 = r - r0[type];
    rr0sq = rr0*rr0;
    Deltasq = Delta[type] * Delta[type];
    rlogarg = 1.0 - rr0sq/Deltasq;

    // if r -> Delta, then rlogarg < 0.0 which is an error
    // issue a warning and reset rlogarg = epsilon
    // if r > 2*Delta something serious is wrong, abort

    if (rlogarg < 0.1) {
      char str[128];
      sprintf(str,"FENE bond too long: " BIGINT_FORMAT " " 
              TAGINT_FORMAT " " TAGINT_FORMAT " %g",
              update->ntimestep,atom->tag[i],atom->tag[j],sqrt(rsq));
      error->warning(FLERR,str,0);
      if (rlogarg <= -3.0) error->one(FLERR,"Bad FENE bond");
    }

    fbond = -k[type]*rr0/rlogarg/Deltasq/r;
    delf[0] = delr[0]*fbond;
    delf[1] = delr[1]*fbond;
    delf[2] = delr[2]*fbond;

    // energy

    if (eflag) {
      ebond = -0.5 * k[type]*log(rlogarg);
    }

    // apply force and torque to each of 2 atoms

    if (newton_bond || i < nlocal) {

      f[i][0] += delf[0];
      f[i][1] += delf[1];
      f[i][2] += delf[2];

      MathExtra::cross3(r_cs_i,delf,delt);

      torque[i][0] += delt[0];
      torque[i][1] += delt[1];
      torque[i][2] += delt[2];

    }

    if (newton_bond || j < nlocal) {

      f[j][0] -= delf[0];
      f[j][1] -= delf[1];
      f[j][2] -= delf[2];

      MathExtra::cross3(r_cs_j,delf,delt);

      torque[j][0] -= delt[0];
      torque[j][1] -= delt[1];
      torque[j][2] -= delt[2];

    }

    if (evflag) ev_tally(i,j,nlocal,newton_bond,ebond,fbond,delr[0],delr[1],delr[2]);

  }

}

/* ---------------------------------------------------------------------- */

void BondOxdna::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(Delta,n+1,"bond:Delta");
  memory->create(r0,n+1,"bond:r0");
  memory->create(setflag,n+1,"bond:setflag");

  for (int i = 1; i <= n; i++) setflag[i] = 0;

}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void BondOxdna::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double Delta_one = force->numeric(FLERR,arg[2]);
  double r0_one = force->numeric(FLERR,arg[3]);

  int count = 0;

  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    Delta[i] = Delta_one;
    r0[i] = r0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");

}

/* ----------------------------------------------------------------------
   set special_bond settings and check if valid
------------------------------------------------------------------------- */

void BondOxdna::init_style()
{
  /* special bonds have to be lj = 0 1 1 and coul = 1 1 1 to exclude 
     the ss excluded volume interaction between nearest neighbours   */ 

  force->special_lj[1] = 0.0;
  force->special_lj[2] = 1.0;
  force->special_lj[3] = 1.0;
  force->special_coul[1] = 1.0;
  force->special_coul[2] = 1.0;
  force->special_coul[3] = 1.0;

  fprintf(screen,"Finding 1-2 1-3 1-4 neighbors ...\n"
		 " Special bond factors lj:   %-10g %-10g %-10g\n"
		 " Special bond factors coul: %-10g %-10g %-10g\n",
		 force->special_lj[1],force->special_lj[2],force->special_lj[3],
		 force->special_coul[1],force->special_coul[2],force->special_coul[3]);

  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 ||
      force->special_lj[3] != 1.0) {
    if (comm->me == 0)
      error->warning(FLERR,"Use special bonds lj = 0,1,1 and coul = 1,1,1 with bond style oxdna");
  }

}

/* ---------------------------------------------------------------------- */

double BondOxdna::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondOxdna::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&Delta[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondOxdna::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&Delta[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&Delta[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondOxdna::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,k[i],r0[i],Delta[i]);
}

/* ---------------------------------------------------------------------- */

double BondOxdna::single(int type, double rsq, int i, int j,
                        double &fforce)
{
  double r = sqrt(rsq);
  double rr0 = r - r0[type];
  double rr0sq = rr0*rr0;
  double Deltasq = Delta[type] * Delta[type];
  double rlogarg = 1.0 - rr0sq/Deltasq;

  // if r -> Delta, then rlogarg < 0.0 which is an error
  // issue a warning and reset rlogarg = epsilon
  // if r > 2*Delta something serious is wrong, abort

  if (rlogarg < 0.1) {
    char str[128];
    sprintf(str,"FENE bond too long: " BIGINT_FORMAT " %g",
            update->ntimestep,sqrt(rsq));
    error->warning(FLERR,str,0);
    if (rlogarg <= -3.0) error->one(FLERR,"Bad FENE bond");
  }

  double eng = -0.5 * k[type]*log(rlogarg);
  fforce = -k[type]*rr0/rlogarg/Deltasq/r;

  return eng;
}
