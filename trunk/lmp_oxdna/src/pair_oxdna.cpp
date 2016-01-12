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
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairOxdna::PairOxdna(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairOxdna::~PairOxdna()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(epsilon_ss);
    memory->destroy(sigma_ss);
    memory->destroy(cut_ss_lj);
    memory->destroy(b_ss);
    memory->destroy(cut_ss_sm);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
    memory->destroy(cutsq);
    memory->destroy(cutsq_ss_lj);
    memory->destroy(cutsq_ss_sm);
  }
}

/* ---------------------------------------------------------------------- */

void PairOxdna::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp_s,ytmp_s,ztmp_s,xtmp_b,ytmp_b,ztmp_b;
  double delx_ss,dely_ss,delz_ss,rsq_ss;
  double delx_sb,dely_sb,delz_sb,rsq_sb;
  double delx_bs,dely_bs,delz_bs,rsq_bs;
  double delx_bb,dely_bb,delz_bb,rsq_bb;
  double evdwl,fpair,delf[3];
  double r,rinv,r2inv,r6inv,forcelj,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double d_coms=-0.24, d_comb=+0.56; // distance COM-backbone and COM-base

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
  double *quat1,ex1[3],ey1[3],ez1[3],e_coms1[3],e_comb1[3];
  double *quat2,ex2[3],ey2[3],ez2[3],e_coms2[3],e_comb2[3];

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

    quat1=bonus[i].quat;
    MathExtra::q_to_exyz(quat1,ex1,ey1,ez1);

    // position of backbone site i
    e_coms1[0] = d_coms*ex1[0];
    e_coms1[1] = d_coms*ex1[1];
    e_coms1[2] = d_coms*ex1[2];

    xtmp_s = x[i][0] + e_coms1[0];
    ytmp_s = x[i][1] + e_coms1[1];
    ztmp_s = x[i][2] + e_coms1[2];

    // position of base site i
    e_comb1[0] = d_comb*ex1[0];
    e_comb1[1] = d_comb*ex1[1];
    e_comb1[2] = d_comb*ex1[2];

    xtmp_b = x[i][0] + e_comb1[0];
    ytmp_b = x[i][1] + e_comb1[1];
    ztmp_b = x[i][2] + e_comb1[2];

    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)]; // = 0 for nearest neighbours
      j &= NEIGHMASK;

      quat2=bonus[j].quat;
      MathExtra::q_to_exyz(quat2,ex2,ey2,ez2);

      e_coms2[0] = d_coms*ex2[0];
      e_coms2[1] = d_coms*ex2[1];
      e_coms2[2] = d_coms*ex2[2];

      e_comb2[0] = d_comb*ex2[0];
      e_comb2[1] = d_comb*ex2[1];
      e_comb2[2] = d_comb*ex2[2];

      // rel. distance backbone i - backbone j
      delx_ss = xtmp_s - x[j][0] - e_coms2[0];
      dely_ss = ytmp_s - x[j][1] - e_coms2[1];
      delz_ss = ztmp_s - x[j][2] - e_coms2[2];
      rsq_ss = delx_ss*delx_ss + dely_ss*dely_ss + delz_ss*delz_ss;

      // rel. distance backbone i - base j
      delx_sb = xtmp_s - x[j][0] - e_comb2[0];
      dely_sb = ytmp_s - x[j][1] - e_comb2[1];
      delz_sb = ztmp_s - x[j][2] - e_comb2[2];
      rsq_sb = delx_sb*delx_sb + dely_sb*dely_sb + delz_sb*delz_sb;

      // rel. distance base i - backbone j
      delx_bs = xtmp_b - x[j][0] - e_coms2[0];
      dely_bs = ytmp_b - x[j][1] - e_coms2[1];
      delz_bs = ztmp_b - x[j][2] - e_coms2[2];
      rsq_bs = delx_bs*delx_bs + dely_bs*dely_bs + delz_bs*delz_bs;

      // rel. distance base i - base j
      delx_bb = xtmp_b - x[j][0] - e_comb2[0];
      dely_bb = ytmp_b - x[j][1] - e_comb2[1];
      delz_bb = ztmp_b - x[j][2] - e_comb2[2];
      rsq_bb = delx_bb*delx_bb + dely_bb*dely_bb + delz_bb*delz_bb;

      jtype = type[j];

      // excluded volume interaction 
      // backbone-backbone LJ part

      if (rsq_ss < cutsq_ss_lj[itype][jtype]) {

        r2inv = 1.0/rsq_ss;
        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        fpair = factor_lj*forcelj*r2inv;

	delf[0] = delx_ss*fpair;
	delf[1] = dely_ss*fpair;
	delf[2] = delz_ss*fpair;

        f[i][0] += delf[0];
        f[i][1] += delf[1];
        f[i][2] += delf[2];
	torque[i][0] += e_coms1[1]*delf[2] - e_coms1[2]*delf[1]; 
	torque[i][1] += e_coms1[2]*delf[0] - e_coms1[0]*delf[2];
	torque[i][2] += e_coms1[0]*delf[1] - e_coms1[1]*delf[0];

        if (newton_pair || j < nlocal) {
          f[j][0] -= delf[0];
          f[j][1] -= delf[1];
          f[j][2] -= delf[2];
	  torque[j][0] -= e_coms2[1]*delf[2] - e_coms2[2]*delf[1]; 
	  torque[j][1] -= e_coms2[2]*delf[0] - e_coms2[0]*delf[2];
	  torque[j][2] -= e_coms2[0]*delf[1] - e_coms2[1]*delf[0];
        }

        if (eflag) {
          evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
            offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx_ss,dely_ss,delz_ss);
      }

       // backbone-backbone smoothing part

      if (cutsq_ss_lj[itype][jtype] <= rsq_ss && rsq_ss < cutsq_ss_sm[itype][jtype]) {

	r = sqrt(rsq_ss);
	rinv = 1.0/r;

	fpair = factor_lj*epsilon_ss[itype][jtype]*2.0*b_ss[itype][jtype]*(cut_ss_sm[itype][jtype]*rinv - 1.0);

	delf[0] = delx_ss*fpair; 
	delf[1] = dely_ss*fpair; 
	delf[2] = delz_ss*fpair; 

        f[i][0] += delf[0];
        f[i][1] += delf[1];
        f[i][2] += delf[2];
        torque[i][0] += e_coms1[1]*delf[2] - e_coms1[2]*delf[1];
        torque[i][1] += e_coms1[2]*delf[0] - e_coms1[0]*delf[2];
        torque[i][2] += e_coms1[0]*delf[1] - e_coms1[1]*delf[0];

        if (newton_pair || j < nlocal) {
          f[j][0] -= delf[0];
          f[j][1] -= delf[1];
          f[j][2] -= delf[2];
          torque[j][0] -= e_coms2[1]*delf[2] - e_coms2[2]*delf[1];
          torque[j][1] -= e_coms2[2]*delf[0] - e_coms2[0]*delf[2];
          torque[j][2] -= e_coms2[0]*delf[1] - e_coms2[1]*delf[0];
        }

        if (eflag) {
          evdwl = b_ss[itype][jtype]*
		(cut_ss_sm[itype][jtype]-r)*(cut_ss_sm[itype][jtype]-r);
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx_ss,dely_ss,delz_ss);
      }

      // excluded volume interaction done

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

  memory->create(epsilon_ss,n+1,n+1,"pair:epsilon_ss");
  memory->create(sigma_ss,n+1,n+1,"pair:sigma_ss");
  memory->create(cut_ss_lj,n+1,n+1,"pair:cut_ss_lj");
  memory->create(b_ss,n+1,n+1,"pair:b_ss");
  memory->create(cut_ss_sm,n+1,n+1,"pair:cut_ss_sm");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cutsq_ss_lj,n+1,n+1,"pair:cutsq_ss_lj");
  memory->create(cutsq_ss_sm,n+1,n+1,"pair:cutsq_ss_sm");
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

  if (strncmp(arg[0],"excv_ss",7) == 0) {

    int ilo,ihi,jlo,jhi;
    double epsilon_ss_one, sigma_ss_one;
    double cut_ss_lj_one, cut_ss_sm_one, b_ss_one;

    force->bounds(arg[1],atom->ntypes,ilo,ihi);
    force->bounds(arg[2],atom->ntypes,jlo,jhi);

    /* LJ parameters */
    epsilon_ss_one = force->numeric(FLERR,arg[3]);
    sigma_ss_one = force->numeric(FLERR,arg[4]);
    cut_ss_lj_one = force->numeric(FLERR,arg[5]);

    /* Smoothing - determined through continuity and differentiability */
    b_ss_one = 4.0/sigma_ss_one
	*(6.0*pow(sigma_ss_one/cut_ss_lj_one,7)-12.0*pow(sigma_ss_one/cut_ss_lj_one,13))
	*4.0/sigma_ss_one*(6.0*pow(sigma_ss_one/cut_ss_lj_one,7)-12.0*pow(sigma_ss_one/cut_ss_lj_one,13))
	/4.0/(4.0*(pow(sigma_ss_one/cut_ss_lj_one,12)-pow(sigma_ss_one/cut_ss_lj_one,6)));

    cut_ss_sm_one = cut_ss_lj_one 
	- 2.0*4.0*(pow(sigma_ss_one/cut_ss_lj_one,12)-pow(sigma_ss_one/cut_ss_lj_one,6))
	/(4.0/sigma_ss_one*(6.0*pow(sigma_ss_one/cut_ss_lj_one,7)-12.0*pow(sigma_ss_one/cut_ss_lj_one,13)));

    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
	epsilon_ss[i][j] = epsilon_ss_one;
	sigma_ss[i][j] = sigma_ss_one;
	cut_ss_lj[i][j] = cut_ss_lj_one;
	b_ss[i][j] = b_ss_one;
	cut_ss_sm[i][j] = cut_ss_sm_one;
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
    epsilon_ss[i][j] = mix_energy(epsilon_ss[i][i],epsilon_ss[j][j],
                               sigma_ss[i][i],sigma_ss[j][j]);
    sigma_ss[i][j]  = mix_distance(sigma_ss[i][i],sigma_ss[j][j]);
    cut_ss_lj[i][j] = mix_distance(cut_ss_lj[i][i],cut_ss_lj[j][j]);
    cut_ss_sm[i][j]  = mix_distance(cut_ss_sm[i][i],cut_ss_sm[j][j]);
  }

  cutsq_ss_lj[i][j] = cut_ss_lj[i][j]*cut_ss_lj[i][j];
  cutsq_ss_sm[i][j]  = cut_ss_sm[i][j]*cut_ss_sm[i][j];

  lj1[i][j] = 48.0 * epsilon_ss[i][j] * pow(sigma_ss[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon_ss[i][j] * pow(sigma_ss[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon_ss[i][j] * pow(sigma_ss[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon_ss[i][j] * pow(sigma_ss[i][j],6.0);

  if (offset_flag) {
    error->all(FLERR,"Offset not supported");
  } 
  else offset[i][j] = 0.0;

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

    /* TODO: Energy and pressure tail corrections */
    double sig2 = sigma_ss[i][j]*sigma_ss[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut_ss_lj[i][j]*cut_ss_lj[i][j]*cut_ss_lj[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon_ss[i][j] *
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon_ss[i][j] *
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
  }

  return cut_ss_lj[i][j];
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
        fwrite(&epsilon_ss[i][j],sizeof(double),1,fp);
        fwrite(&sigma_ss[i][j],sizeof(double),1,fp);
        fwrite(&cut_ss_lj[i][j],sizeof(double),1,fp);
        fwrite(&b_ss[i][j],sizeof(double),1,fp);
        fwrite(&cut_ss_sm[i][j],sizeof(double),1,fp);
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
          fread(&epsilon_ss[i][j],sizeof(double),1,fp);
          fread(&sigma_ss[i][j],sizeof(double),1,fp);
          fread(&cut_ss_lj[i][j],sizeof(double),1,fp);
          fread(&b_ss[i][j],sizeof(double),1,fp);
          fread(&cut_ss_sm[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon_ss[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_ss[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_ss_lj[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&b_ss[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_ss_sm[i][j],1,MPI_DOUBLE,0,world);
   
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
    fprintf(fp,"%d %g %g %g %g %g\n",i,epsilon_ss[i][i],sigma_ss[i][i],
	cut_ss_lj[i][i],b_ss[i][i],cut_ss_sm[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairOxdna::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g\n",i,j,epsilon_ss[i][j],sigma_ss[i][j],
	cut_ss_lj[i][j],b_ss[i][j],cut_ss_sm[i][j]);
}

/* ---------------------------------------------------------------------- */
/* TODO: Check if single can be used and defined */
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
  if (strcmp(str,"epsilon_ss") == 0) return (void *) epsilon_ss;
  if (strcmp(str,"sigma_ss") == 0) return (void *) sigma_ss;
  if (strcmp(str,"cut_ss_lj") == 0) return (void *) cut_ss_lj;
  if (strcmp(str,"b_ss") == 0) return (void *) b_ss;
  if (strcmp(str,"cut_ss_sm") == 0) return (void *) cut_ss_sm;
  return NULL;
}
