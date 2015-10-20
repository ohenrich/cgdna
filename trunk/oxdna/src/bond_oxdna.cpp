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
    memory->destroy(r0);
    memory->destroy(shift);
  }
}

/* ---------------------------------------------------------------------- */

void BondOxdna::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond, delf[3];
  double rsq,r0sq,rlogarg;//,sr2,sr6;
  double r,rshift,rshiftsq;
  double d_bb = -0.24; // distance backbone - COM

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;

  double **torque = atom->torque;
  int *ellipsoid = atom->ellipsoid;
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  double *quat1,ex1[3],ey1[3],ez1[3],e_cmbb1[3];
  double *quat2,ex2[3],ey2[3],ez2[3],e_cmbb2[3];

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    quat1=bonus[i1].quat;
    MathExtra::q_to_exyz(quat1,ex1,ey1,ez1);
    quat2=bonus[i2].quat;
    MathExtra::q_to_exyz(quat2,ex2,ey2,ez2);

    e_cmbb1[0] = d_bb*ex1[0];
    e_cmbb1[1] = d_bb*ex1[1];
    e_cmbb1[2] = d_bb*ex1[2];
    e_cmbb2[0] = d_bb*ex2[0];
    e_cmbb2[1] = d_bb*ex2[1];
    e_cmbb2[2] = d_bb*ex2[2];

    delx = (x[i1][0] + e_cmbb1[0]) - (x[i2][0] + e_cmbb2[0]);
    dely = (x[i1][1] + e_cmbb1[1]) - (x[i2][1] + e_cmbb2[1]);
    delz = (x[i1][2] + e_cmbb1[2]) - (x[i2][2] + e_cmbb2[2]);

    // force from log term

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    rshift = r - shift[type];
    rshiftsq = rshift*rshift;
    r0sq = r0[type] * r0[type];
    rlogarg = 1.0 - rshiftsq/r0sq;

    // if r -> r0, then rlogarg < 0.0 which is an error
    // issue a warning and reset rlogarg = epsilon
    // if r > 2*r0 something serious is wrong, abort

    if (rlogarg < 0.1) {
      char str[128];
      sprintf(str,"FENE bond too long: " BIGINT_FORMAT " " 
              TAGINT_FORMAT " " TAGINT_FORMAT " %g",
              update->ntimestep,atom->tag[i1],atom->tag[i2],sqrt(rsq));
      error->warning(FLERR,str,0);
      if (rlogarg <= -3.0) error->one(FLERR,"Bad FENE bond");
      rlogarg = 0.1;
    }

    fbond = -k[type]*rshift/rlogarg/r0sq/r;
    delf[0] = delx*fbond;
    delf[1] = dely*fbond;
    delf[2] = delz*fbond;

    // energy

    if (eflag) {
      ebond = -0.5 * k[type]*log(rlogarg);
    }

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delf[0];
      f[i1][1] += delf[1];
      f[i1][2] += delf[2];
      torque[i1][0] += e_cmbb1[1]*delf[2] - e_cmbb1[2]*delf[1];
      torque[i1][1] += e_cmbb1[2]*delf[0] - e_cmbb1[0]*delf[2];
      torque[i1][2] += e_cmbb1[0]*delf[1] - e_cmbb1[1]*delf[0];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delf[0];
      f[i2][1] -= delf[1];
      f[i2][2] -= delf[2];
      torque[i2][0] -= e_cmbb2[1]*delf[2] - e_cmbb2[2]*delf[1];
      torque[i2][1] -= e_cmbb2[2]*delf[0] - e_cmbb2[0]*delf[2];
      torque[i2][2] -= e_cmbb2[0]*delf[1] - e_cmbb2[1]*delf[0];
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }

}

/* ---------------------------------------------------------------------- */

void BondOxdna::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(r0,n+1,"bond:r0");
  memory->create(shift,n+1,"bond:shift");
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
  double shift_one = force->numeric(FLERR,arg[2]);
  double r0_one = force->numeric(FLERR,arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    shift[i] = shift_one;
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
  // special bonds should be 0 1 1

  force->special_lj[1] = 0.0;
  force->special_lj[2] = 1.0;
  force->special_lj[3] = 1.0;
  force->special_coul[1] = 0.0;
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
      error->warning(FLERR,"Use special bonds = 0,1,1 with bond style oxdna");
  }


}

/* ---------------------------------------------------------------------- */

double BondOxdna::equilibrium_distance(int i)
{
  return shift[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondOxdna::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&shift[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondOxdna::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&shift[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&shift[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondOxdna::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,k[i],shift[i],r0[i]);
}

/* ---------------------------------------------------------------------- */

double BondOxdna::single(int type, double rsq, int i, int j,
                        double &fforce)
{
  double r = sqrt(rsq);
  double rshift = r - shift[type];
  double rshiftsq = rshift*rshift;
  double r0sq = r0[type] * r0[type];
  double rlogarg = 1.0 - rshiftsq/r0sq;

  // if r -> r0, then rlogarg < 0.0 which is an error
  // issue a warning and reset rlogarg = epsilon
  // if r > 2*r0 something serious is wrong, abort

  if (rlogarg < 0.1) {
    char str[128];
    sprintf(str,"FENE bond too long: " BIGINT_FORMAT " %g",
            update->ntimestep,sqrt(rsq));
    error->warning(FLERR,str,0);
    if (rlogarg <= -3.0) error->one(FLERR,"Bad FENE bond");
    rlogarg = 0.1;
  }

  double eng = -0.5 * k[type]*log(rlogarg);
  fforce = -k[type]*rshift/rlogarg/r0sq/r;

  return eng;
}
