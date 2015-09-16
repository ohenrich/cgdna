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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_oxdna.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "atom_vec_ellipsoid.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairOxdna::PairOxdna(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairOxdna::~PairOxdna()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(epsilon_bb);
    memory->destroy(r0_bb);
    memory->destroy(Delta_bb);
    memory->destroy(epsilon_exc);
    memory->destroy(sigma_exc);
    memory->destroy(cut_exc);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairOxdna::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  int *type = atom->type;
  int *ellipsoid = atom->ellipsoid;
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        fpair = factor_lj*forcelj*r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
            offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairOxdna::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(epsilon_bb,n+1,n+1,"pair:epsilon_bb");
  memory->create(r0_bb,n+1,n+1,"pair:r0_bb");
  memory->create(Delta_bb,n+1,n+1,"pair:Delta_bb");
  memory->create(epsilon_exc,n+1,n+1,"pair:epsilon_exc");
  memory->create(sigma_exc,n+1,n+1,"pair:sigma_exc");
  memory->create(cut_exc,n+1,n+1,"pair:cut_exc");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairOxdna::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairOxdna::coeff(int narg, char **arg)
{
  if (narg != 6)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int count = 0;

  if (strncmp(arg[0],"exc",3) == 0) {

    int ilo,ihi,jlo,jhi;
    double epsilon_exc_one, sigma_exc_one, cut_exc_one;

    force->bounds(arg[1],atom->ntypes,ilo,ihi);
    force->bounds(arg[2],atom->ntypes,jlo,jhi);

    epsilon_exc_one = force->numeric(FLERR,arg[3]);
    sigma_exc_one = force->numeric(FLERR,arg[4]);
    cut_exc_one = force->numeric(FLERR,arg[5]);

    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
	epsilon_exc[i][j] = epsilon_exc_one;
	sigma_exc[i][j] = sigma_exc_one;
	cut_exc[i][j] = cut_exc_one;
	setflag[i][j] = 1;
	count++;
      }
    }

  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairOxdna::init_style()
{
  int irequest;

  // request regular neighbor lists 

  irequest = neighbor->request(this,instance_me);

}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use regular
------------------------------------------------------------------------- */

void PairOxdna::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  if (id  > 0) error->all(FLERR,"Respa not supported");
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairOxdna::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon_exc[i][j] = mix_energy(epsilon_exc[i][i],epsilon_exc[j][j],
                               sigma_exc[i][i],sigma_exc[j][j]);
    sigma_exc[i][j] = mix_distance(sigma_exc[i][i],sigma_exc[j][j]);
    cut_exc[i][j] = mix_distance(cut_exc[i][i],cut_exc[j][j]);
  }

  lj1[i][j] = 48.0 * epsilon_exc[i][j] * pow(sigma_exc[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon_exc[i][j] * pow(sigma_exc[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon_exc[i][j] * pow(sigma_exc[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon_exc[i][j] * pow(sigma_exc[i][j],6.0);

  if (offset_flag) {
    double ratio = sigma_exc[i][j] / cut_exc[i][j];
    offset[i][j] = 4.0 * epsilon_exc[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double sig2 = sigma_exc[i][j]*sigma_exc[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut_exc[i][j]*cut_exc[i][j]*cut_exc[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon_exc[i][j] *
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon_exc[i][j] *
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
  }

  return cut_exc[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdna::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon_exc[i][j],sizeof(double),1,fp);
        fwrite(&sigma_exc[i][j],sizeof(double),1,fp);
        fwrite(&cut_exc[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdna::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon_exc[i][j],sizeof(double),1,fp);
          fread(&sigma_exc[i][j],sizeof(double),1,fp);
          fread(&cut_exc[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon_exc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_exc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_exc[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdna::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdna::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairOxdna::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,epsilon_exc[i][i],sigma_exc[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairOxdna::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,epsilon_exc[i][j],sigma_exc[i][j],cut_exc[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairOxdna::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double r2inv,r6inv,forcelj,philj;

  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
  fforce = factor_lj*forcelj*r2inv;

  philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
    offset[itype][jtype];
  return factor_lj*philj;
}

/* ---------------------------------------------------------------------- */

void *PairOxdna::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon_exc") == 0) return (void *) epsilon_exc;
  if (strcmp(str,"sigma_exc") == 0) return (void *) sigma_exc;
  return NULL;
}
