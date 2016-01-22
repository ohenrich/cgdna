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
    memory->destroy(cutsq);

    memory->destroy(epsilon_ss);
    memory->destroy(sigma_ss);
    memory->destroy(cut_ss_lj);
    memory->destroy(b_ss);
    memory->destroy(cut_ss_sm);
    memory->destroy(lj1_ss);
    memory->destroy(lj2_ss);
    memory->destroy(lj3_ss);
    memory->destroy(lj4_ss);
    memory->destroy(offset_ss);
    memory->destroy(cutsq_ss_lj);
    memory->destroy(cutsq_ss_sm);

    memory->destroy(epsilon_sb);
    memory->destroy(sigma_sb);
    memory->destroy(cut_sb_lj);
    memory->destroy(b_sb);
    memory->destroy(cut_sb_sm);
    memory->destroy(lj1_sb);
    memory->destroy(lj2_sb);
    memory->destroy(lj3_sb);
    memory->destroy(lj4_sb);
    memory->destroy(offset_sb);
    memory->destroy(cutsq_sb_lj);
    memory->destroy(cutsq_sb_sm);

    memory->destroy(epsilon_bb);
    memory->destroy(sigma_bb);
    memory->destroy(cut_bb_lj);
    memory->destroy(b_bb);
    memory->destroy(cut_bb_sm);
    memory->destroy(lj1_bb);
    memory->destroy(lj2_bb);
    memory->destroy(lj3_bb);
    memory->destroy(lj4_bb);
    memory->destroy(offset_bb);
    memory->destroy(cutsq_bb_lj);
    memory->destroy(cutsq_bb_sm);

  }
}

/* ---------------------------------------------------------------------- */

void PairOxdna::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,ia;
  double rtmp_s[3],rtmp_b[3],delf[3],delt[3];
  double delr_ss[3],rsq_ss,delr_sb[3],rsq_sb;
  double delr_bs[3],rsq_bs,delr_bb[3],rsq_bb;
  double evdwl,fpair;
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

    for (ia = 0; ia < 3; ia++) {
	// position of backbone site i
	e_coms1[ia] = d_coms*ex1[ia];
	rtmp_s[ia] = x[i][ia] + e_coms1[ia];

	// position of base site i
	e_comb1[ia] = d_comb*ex1[ia];
	rtmp_b[ia] = x[i][ia] + e_comb1[ia];
    }

    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)]; // = 0 for nearest neighbours
      j &= NEIGHMASK;

      quat2=bonus[j].quat;
      MathExtra::q_to_exyz(quat2,ex2,ey2,ez2);

      for (ia = 0; ia < 3; ia++) {

	    e_coms2[ia] = d_coms*ex2[ia];
	    e_comb2[ia] = d_comb*ex2[ia];

	    // rel. distance backbone i - backbone j
	    delr_ss[ia] = rtmp_s[ia] - x[j][ia] - e_coms2[ia];
	    // rel. distance backbone i - base j
	    delr_sb[ia] = rtmp_s[ia] - x[j][ia] - e_comb2[ia];
	    // rel. distance base i - backbone j
	    delr_bs[ia] = rtmp_b[ia] - x[j][ia] - e_coms2[ia];
	    // rel. distance base i - base j
	    delr_bb[ia] = rtmp_b[ia] - x[j][ia] - e_comb2[ia];

      }

      rsq_ss = delr_ss[0]*delr_ss[0] + delr_ss[1]*delr_ss[1] + delr_ss[2]*delr_ss[2];
      rsq_sb = delr_sb[0]*delr_sb[0] + delr_sb[1]*delr_sb[1] + delr_sb[2]*delr_sb[2];
      rsq_bs = delr_bs[0]*delr_bs[0] + delr_bs[1]*delr_bs[1] + delr_bs[2]*delr_bs[2];
      rsq_bb = delr_bb[0]*delr_bb[0] + delr_bb[1]*delr_bb[1] + delr_bb[2]*delr_bb[2];

      jtype = type[j];

      // excluded volume interaction

      // backbone-backbone
      if (rsq_ss < cutsq_ss_sm[itype][jtype]) {

	if (rsq_ss < cutsq_ss_lj[itype][jtype]) {

	  calculate_excv_lj(factor_lj,rsq_ss,lj1_ss[itype][jtype],
				  lj2_ss[itype][jtype],r6inv,fpair);

	  if (eflag) {
	    evdwl = r6inv*(lj3_ss[itype][jtype]*r6inv-lj4_ss[itype][jtype]) -
	      offset_ss[itype][jtype];
	    evdwl *= factor_lj;
	  }

	  if (evflag) ev_tally(i,j,nlocal,newton_pair,
			       evdwl,0.0,fpair,delr_ss[0],delr_ss[1],delr_ss[2]);
	}
	else {

	  calculate_excv_sm(factor_lj,rsq_ss,epsilon_ss[itype][jtype],
				  b_ss[itype][jtype],cut_ss_sm[itype][jtype],r,fpair);

	  if (eflag) {
	    evdwl = b_ss[itype][jtype]*
		  (cut_ss_sm[itype][jtype]-r)*(cut_ss_sm[itype][jtype]-r);
	    evdwl *= factor_lj;
	  }

	  if (evflag) ev_tally(i,j,nlocal,newton_pair,
			       evdwl,0.0,fpair,delr_ss[0],delr_ss[1],delr_ss[2]);
	}

	delf[0] = delr_ss[0]*fpair; 
	delf[1] = delr_ss[1]*fpair; 
	delf[2] = delr_ss[2]*fpair; 

        f[i][0] += delf[0];
        f[i][1] += delf[1];
        f[i][2] += delf[2];

	MathExtra::cross3(e_coms1,delf,delt);

	torque[i][0] += delt[0];
	torque[i][1] += delt[1];
	torque[i][2] += delt[2];

        if (newton_pair || j < nlocal) {
          f[j][0] -= delf[0];
          f[j][1] -= delf[1];
          f[j][2] -= delf[2];

	  MathExtra::cross3(e_coms2,delf,delt);

	  torque[j][0] -= delt[0];
	  torque[j][1] -= delt[1];
	  torque[j][2] -= delt[2];
        }

      }

      // backbone-base
      if (rsq_sb < cutsq_sb_sm[itype][jtype]) {

	if (rsq_sb < cutsq_sb_lj[itype][jtype]) {

	  calculate_excv_lj(1.0,rsq_sb,lj1_sb[itype][jtype],lj2_sb[itype][jtype],r6inv,fpair);

	  if (eflag) {
	    evdwl = r6inv*(lj3_sb[itype][jtype]*r6inv-lj4_sb[itype][jtype]) -
	      offset_sb[itype][jtype];
	    evdwl *= factor_lj;
	  }

	  if (evflag) ev_tally(i,j,nlocal,newton_pair,
			       evdwl,0.0,fpair,delr_sb[0],delr_sb[1],delr_sb[2]);
	}
	else {

	  calculate_excv_sm(1.0,rsq_sb,epsilon_sb[itype][jtype],
				  b_sb[itype][jtype],cut_sb_sm[itype][jtype],r,fpair);

	  if (eflag) {
	    evdwl = b_sb[itype][jtype]*
		  (cut_sb_sm[itype][jtype]-r)*(cut_sb_sm[itype][jtype]-r);
	    evdwl *= factor_lj;
	  }

	  if (evflag) ev_tally(i,j,nlocal,newton_pair,
			       evdwl,0.0,fpair,delr_sb[0],delr_sb[1],delr_sb[2]);
	}

	delf[0] = delr_sb[0]*fpair; 
	delf[1] = delr_sb[1]*fpair; 
	delf[2] = delr_sb[2]*fpair; 

	f[i][0] += delf[0];
	f[i][1] += delf[1];
	f[i][2] += delf[2];

	MathExtra::cross3(e_coms1,delf,delt);

	torque[i][0] += delt[0];
	torque[i][1] += delt[1];
	torque[i][2] += delt[2];

	if (newton_pair || j < nlocal) {
	  f[j][0] -= delf[0];
	  f[j][1] -= delf[1];
	  f[j][2] -= delf[2];

	  MathExtra::cross3(e_comb2,delf,delt);

	  torque[j][0] -= delt[0];
	  torque[j][1] -= delt[1];
	  torque[j][2] -= delt[2];
	}

      }

      // base-backbone
      if (rsq_bs < cutsq_sb_sm[itype][jtype]) {

	if (rsq_bs < cutsq_sb_lj[itype][jtype]) {

	  calculate_excv_lj(1.0,rsq_bs,lj1_sb[itype][jtype],lj2_sb[itype][jtype],r6inv,fpair);

	  if (eflag) {
	    evdwl = r6inv*(lj3_sb[itype][jtype]*r6inv-lj4_sb[itype][jtype]) -
	      offset_sb[itype][jtype];
	  }

	  if (evflag) ev_tally(i,j,nlocal,newton_pair,
			       evdwl,0.0,fpair,delr_bs[0],delr_bs[1],delr_bs[2]);
	}
	else {

	  calculate_excv_sm(1.0,rsq_bs,epsilon_sb[itype][jtype],
				  b_sb[itype][jtype],cut_sb_sm[itype][jtype],r,fpair);


	  if (eflag) {
	    evdwl = b_sb[itype][jtype]*
		  (cut_sb_sm[itype][jtype]-r)*(cut_sb_sm[itype][jtype]-r);
	  }

	  if (evflag) ev_tally(i,j,nlocal,newton_pair,
			       evdwl,0.0,fpair,delr_bs[0],delr_bs[1],delr_bs[2]);
	}

	delf[0] = delr_bs[0]*fpair; 
	delf[1] = delr_bs[1]*fpair; 
	delf[2] = delr_bs[2]*fpair; 

	f[i][0] += delf[0];
	f[i][1] += delf[1];
	f[i][2] += delf[2];
 
	MathExtra::cross3(e_comb1,delf,delt);

	torque[i][0] += delt[0];
	torque[i][1] += delt[1];
	torque[i][2] += delt[2];

	if (newton_pair || j < nlocal) {
	  f[j][0] -= delf[0];
	  f[j][1] -= delf[1];
	  f[j][2] -= delf[2];

	  MathExtra::cross3(e_coms2,delf,delt);

	  torque[j][0] -= delt[0];
	  torque[j][1] -= delt[1];
	  torque[j][2] -= delt[2];
	}

      }

      // base-base
      if (rsq_bb < cutsq_bb_sm[itype][jtype]) {

	if (rsq_bb < cutsq_bb_lj[itype][jtype]) {

	  calculate_excv_lj(1.0,rsq_bb,lj1_bb[itype][jtype],lj2_bb[itype][jtype],r6inv,fpair);

	  if (eflag) {
	    evdwl = r6inv*(lj3_bb[itype][jtype]*r6inv-lj4_bb[itype][jtype]) -
	      offset_bb[itype][jtype];
	  }

	  if (evflag) ev_tally(i,j,nlocal,newton_pair,
			       evdwl,0.0,fpair,delr_bb[0],delr_bb[1],delr_bb[2]);
	}
	else {

	  calculate_excv_sm(1.0,rsq_bb,epsilon_bb[itype][jtype],
				  b_bb[itype][jtype],cut_bb_sm[itype][jtype],r,fpair);

	  if (eflag) {
	    evdwl = b_bb[itype][jtype]*
		  (cut_bb_sm[itype][jtype]-r)*(cut_bb_sm[itype][jtype]-r);
	  }

	  if (evflag) ev_tally(i,j,nlocal,newton_pair,
			       evdwl,0.0,fpair,delr_bb[0],delr_bb[1],delr_bb[2]);
	}

	delf[0] = delr_bb[0]*fpair; 
	delf[1] = delr_bb[1]*fpair; 
	delf[2] = delr_bb[2]*fpair; 

	f[i][0] += delf[0];
	f[i][1] += delf[1];
	f[i][2] += delf[2];
 
	MathExtra::cross3(e_comb1,delf,delt);

	torque[i][0] += delt[0];
	torque[i][1] += delt[1];
	torque[i][2] += delt[2];

	if (newton_pair || j < nlocal) {
	  f[j][0] -= delf[0];
	  f[j][1] -= delf[1];
	  f[j][2] -= delf[2];

	  MathExtra::cross3(e_comb2,delf,delt);

	  torque[j][0] -= delt[0];
	  torque[j][1] -= delt[1];
	  torque[j][2] -= delt[2];
	}

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

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(epsilon_ss,n+1,n+1,"pair:epsilon_ss");
  memory->create(sigma_ss,n+1,n+1,"pair:sigma_ss");
  memory->create(cut_ss_lj,n+1,n+1,"pair:cut_ss_lj");
  memory->create(b_ss,n+1,n+1,"pair:b_ss");
  memory->create(cut_ss_sm,n+1,n+1,"pair:cut_ss_sm");
  memory->create(lj1_ss,n+1,n+1,"pair:lj1_ss");
  memory->create(lj2_ss,n+1,n+1,"pair:lj2_ss");
  memory->create(lj3_ss,n+1,n+1,"pair:lj3_ss");
  memory->create(lj4_ss,n+1,n+1,"pair:lj4_ss");
  memory->create(offset_ss,n+1,n+1,"pair:offset_ss");
  memory->create(cutsq_ss_lj,n+1,n+1,"pair:cutsq_ss_lj");
  memory->create(cutsq_ss_sm,n+1,n+1,"pair:cutsq_ss_sm");

  memory->create(epsilon_sb,n+1,n+1,"pair:epsilon_sb");
  memory->create(sigma_sb,n+1,n+1,"pair:sigma_sb");
  memory->create(cut_sb_lj,n+1,n+1,"pair:cut_sb_lj");
  memory->create(b_sb,n+1,n+1,"pair:b_sb");
  memory->create(cut_sb_sm,n+1,n+1,"pair:cut_sb_sm");
  memory->create(lj1_sb,n+1,n+1,"pair:lj1_sb");
  memory->create(lj2_sb,n+1,n+1,"pair:lj2_sb");
  memory->create(lj3_sb,n+1,n+1,"pair:lj3_sb");
  memory->create(lj4_sb,n+1,n+1,"pair:lj4_sb");
  memory->create(offset_sb,n+1,n+1,"pair:offset_sb");
  memory->create(cutsq_sb_lj,n+1,n+1,"pair:cutsq_sb_lj");
  memory->create(cutsq_sb_sm,n+1,n+1,"pair:cutsq_sb_sm");

  memory->create(epsilon_bb,n+1,n+1,"pair:epsilon_bb");
  memory->create(sigma_bb,n+1,n+1,"pair:sigma_bb");
  memory->create(cut_bb_lj,n+1,n+1,"pair:cut_bb_lj");
  memory->create(b_bb,n+1,n+1,"pair:b_bb");
  memory->create(cut_bb_sm,n+1,n+1,"pair:cut_bb_sm");
  memory->create(lj1_bb,n+1,n+1,"pair:lj1_bb");
  memory->create(lj2_bb,n+1,n+1,"pair:lj2_bb");
  memory->create(lj3_bb,n+1,n+1,"pair:lj3_bb");
  memory->create(lj4_bb,n+1,n+1,"pair:lj4_bb");
  memory->create(offset_bb,n+1,n+1,"pair:offset_bb");
  memory->create(cutsq_bb_lj,n+1,n+1,"pair:cutsq_bb_lj");
  memory->create(cutsq_bb_sm,n+1,n+1,"pair:cutsq_bb_sm");

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
  int count;

  if (narg != 11) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  count = 0;

  double epsilon_ss_one, sigma_ss_one;
  double cut_ss_lj_one, cut_ss_sm_one, b_ss_one;

  // LJ parameters
  epsilon_ss_one = force->numeric(FLERR,arg[2]);
  sigma_ss_one = force->numeric(FLERR,arg[3]);
  cut_ss_lj_one = force->numeric(FLERR,arg[4]);

  // Smoothing - determined through continuity and differentiability
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

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

  count = 0;
  
  double epsilon_sb_one, sigma_sb_one;
  double cut_sb_lj_one, cut_sb_sm_one, b_sb_one;

  // LJ parameters
  epsilon_sb_one = force->numeric(FLERR,arg[5]);
  sigma_sb_one = force->numeric(FLERR,arg[6]);
  cut_sb_lj_one = force->numeric(FLERR,arg[7]);

  // Smoothing - determined through continuity and differentiability 
  b_sb_one = 4.0/sigma_sb_one
      *(6.0*pow(sigma_sb_one/cut_sb_lj_one,7)-12.0*pow(sigma_sb_one/cut_sb_lj_one,13))
      *4.0/sigma_sb_one*(6.0*pow(sigma_sb_one/cut_sb_lj_one,7)-12.0*pow(sigma_sb_one/cut_sb_lj_one,13))
      /4.0/(4.0*(pow(sigma_sb_one/cut_sb_lj_one,12)-pow(sigma_sb_one/cut_sb_lj_one,6)));

  cut_sb_sm_one = cut_sb_lj_one 
      - 2.0*4.0*(pow(sigma_sb_one/cut_sb_lj_one,12)-pow(sigma_sb_one/cut_sb_lj_one,6))
      /(4.0/sigma_sb_one*(6.0*pow(sigma_sb_one/cut_sb_lj_one,7)-12.0*pow(sigma_sb_one/cut_sb_lj_one,13)));

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon_sb[i][j] = epsilon_sb_one;
      sigma_sb[i][j] = sigma_sb_one;
      cut_sb_lj[i][j] = cut_sb_lj_one;
      b_sb[i][j] = b_sb_one;
      cut_sb_sm[i][j] = cut_sb_sm_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

  count = 0;
  
  double epsilon_bb_one, sigma_bb_one;
  double cut_bb_lj_one, cut_bb_sm_one, b_bb_one;

  // LJ parameters
  epsilon_bb_one = force->numeric(FLERR,arg[8]);
  sigma_bb_one = force->numeric(FLERR,arg[9]);
  cut_bb_lj_one = force->numeric(FLERR,arg[10]);

  // Smoothing - determined through continuity and differentiability 
  b_bb_one = 4.0/sigma_bb_one
      *(6.0*pow(sigma_bb_one/cut_bb_lj_one,7)-12.0*pow(sigma_bb_one/cut_bb_lj_one,13))
      *4.0/sigma_bb_one*(6.0*pow(sigma_bb_one/cut_bb_lj_one,7)-12.0*pow(sigma_bb_one/cut_bb_lj_one,13))
      /4.0/(4.0*(pow(sigma_bb_one/cut_bb_lj_one,12)-pow(sigma_bb_one/cut_bb_lj_one,6)));

  cut_bb_sm_one = cut_bb_lj_one 
      - 2.0*4.0*(pow(sigma_bb_one/cut_bb_lj_one,12)-pow(sigma_bb_one/cut_bb_lj_one,6))
      /(4.0/sigma_bb_one*(6.0*pow(sigma_bb_one/cut_bb_lj_one,7)-12.0*pow(sigma_bb_one/cut_bb_lj_one,13)));

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon_bb[i][j] = epsilon_bb_one;
      sigma_bb[i][j] = sigma_bb_one;
      cut_bb_lj[i][j] = cut_bb_lj_one;
      b_bb[i][j] = b_bb_one;
      cut_bb_sm[i][j] = cut_bb_sm_one;
      setflag[i][j] = 1;
      count++;
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

    epsilon_sb[i][j] = mix_energy(epsilon_sb[i][i],epsilon_sb[j][j],
                               sigma_sb[i][i],sigma_sb[j][j]);
    sigma_sb[i][j]  = mix_distance(sigma_sb[i][i],sigma_sb[j][j]);
    cut_sb_lj[i][j] = mix_distance(cut_sb_lj[i][i],cut_sb_lj[j][j]);
    cut_sb_sm[i][j]  = mix_distance(cut_sb_sm[i][i],cut_sb_sm[j][j]);

    epsilon_bb[i][j] = mix_energy(epsilon_bb[i][i],epsilon_bb[j][j],
                               sigma_bb[i][i],sigma_bb[j][j]);
    sigma_bb[i][j]  = mix_distance(sigma_bb[i][i],sigma_bb[j][j]);
    cut_bb_lj[i][j] = mix_distance(cut_bb_lj[i][i],cut_bb_lj[j][j]);
    cut_bb_sm[i][j]  = mix_distance(cut_bb_sm[i][i],cut_bb_sm[j][j]);

  }

  cutsq_ss_lj[i][j] = cut_ss_lj[i][j]*cut_ss_lj[i][j];
  cutsq_ss_sm[i][j]  = cut_ss_sm[i][j]*cut_ss_sm[i][j];

  cutsq_sb_lj[i][j] = cut_sb_lj[i][j]*cut_sb_lj[i][j];
  cutsq_sb_sm[i][j]  = cut_sb_sm[i][j]*cut_sb_sm[i][j];

  cutsq_bb_lj[i][j] = cut_bb_lj[i][j]*cut_bb_lj[i][j];
  cutsq_bb_sm[i][j]  = cut_bb_sm[i][j]*cut_bb_sm[i][j];

  lj1_ss[i][j] = 48.0 * epsilon_ss[i][j] * pow(sigma_ss[i][j],12.0);
  lj2_ss[i][j] = 24.0 * epsilon_ss[i][j] * pow(sigma_ss[i][j],6.0);
  lj3_ss[i][j] = 4.0 * epsilon_ss[i][j] * pow(sigma_ss[i][j],12.0);
  lj4_ss[i][j] = 4.0 * epsilon_ss[i][j] * pow(sigma_ss[i][j],6.0);

  lj1_sb[i][j] = 48.0 * epsilon_sb[i][j] * pow(sigma_sb[i][j],12.0);
  lj2_sb[i][j] = 24.0 * epsilon_sb[i][j] * pow(sigma_sb[i][j],6.0);
  lj3_sb[i][j] = 4.0 * epsilon_sb[i][j] * pow(sigma_sb[i][j],12.0);
  lj4_sb[i][j] = 4.0 * epsilon_sb[i][j] * pow(sigma_sb[i][j],6.0);

  lj1_bb[i][j] = 48.0 * epsilon_bb[i][j] * pow(sigma_bb[i][j],12.0);
  lj2_bb[i][j] = 24.0 * epsilon_bb[i][j] * pow(sigma_bb[i][j],6.0);
  lj3_bb[i][j] = 4.0 * epsilon_bb[i][j] * pow(sigma_bb[i][j],12.0);
  lj4_bb[i][j] = 4.0 * epsilon_bb[i][j] * pow(sigma_bb[i][j],6.0);

  if (offset_flag) {
    error->all(FLERR,"Offset not supported");
  } 
  else {
    offset_ss[i][j] = 0.0;
    offset_sb[i][j] = 0.0;
    offset_bb[i][j] = 0.0;
  }

  lj1_ss[j][i] = lj1_ss[i][j];
  lj2_ss[j][i] = lj2_ss[i][j];
  lj3_ss[j][i] = lj3_ss[i][j];
  lj4_ss[j][i] = lj4_ss[i][j];
  offset_ss[j][i] = offset_ss[i][j];

  lj1_sb[j][i] = lj1_sb[i][j];
  lj2_sb[j][i] = lj2_sb[i][j];
  lj3_sb[j][i] = lj3_sb[i][j];
  lj4_sb[j][i] = lj4_sb[i][j];
  offset_sb[j][i] = offset_sb[i][j];

  lj1_bb[j][i] = lj1_bb[i][j];
  lj2_bb[j][i] = lj2_bb[i][j];
  lj3_bb[j][i] = lj3_bb[i][j];
  lj4_bb[j][i] = lj4_bb[i][j];
  offset_bb[j][i] = offset_bb[i][j];

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
        fwrite(&epsilon_sb[i][j],sizeof(double),1,fp);
        fwrite(&sigma_sb[i][j],sizeof(double),1,fp);
        fwrite(&cut_sb_lj[i][j],sizeof(double),1,fp);
        fwrite(&b_sb[i][j],sizeof(double),1,fp);
        fwrite(&cut_sb_sm[i][j],sizeof(double),1,fp);
        fwrite(&epsilon_bb[i][j],sizeof(double),1,fp);
        fwrite(&sigma_bb[i][j],sizeof(double),1,fp);
        fwrite(&cut_bb_lj[i][j],sizeof(double),1,fp);
        fwrite(&b_bb[i][j],sizeof(double),1,fp);
        fwrite(&cut_bb_sm[i][j],sizeof(double),1,fp);
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
          fread(&epsilon_sb[i][j],sizeof(double),1,fp);
          fread(&sigma_sb[i][j],sizeof(double),1,fp);
          fread(&cut_sb_lj[i][j],sizeof(double),1,fp);
          fread(&b_sb[i][j],sizeof(double),1,fp);
          fread(&cut_sb_sm[i][j],sizeof(double),1,fp);
          fread(&epsilon_bb[i][j],sizeof(double),1,fp);
          fread(&sigma_bb[i][j],sizeof(double),1,fp);
          fread(&cut_bb_lj[i][j],sizeof(double),1,fp);
          fread(&b_bb[i][j],sizeof(double),1,fp);
          fread(&cut_bb_sm[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon_ss[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_ss[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_ss_lj[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&b_ss[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_ss_sm[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&epsilon_sb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_sb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_sb_lj[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&b_sb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_sb_sm[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&epsilon_bb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_bb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_bb_lj[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&b_bb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_bb_sm[i][j],1,MPI_DOUBLE,0,world);
   
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
    fprintf(fp,"%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n",i,
	epsilon_ss[i][i],sigma_ss[i][i],cut_ss_lj[i][i],b_ss[i][i],cut_ss_sm[i][i],
	epsilon_sb[i][i],sigma_sb[i][i],cut_sb_lj[i][i],b_sb[i][i],cut_sb_sm[i][i],
	epsilon_bb[i][i],sigma_bb[i][i],cut_bb_lj[i][i],b_bb[i][i],cut_bb_sm[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairOxdna::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",i,j,
	epsilon_ss[i][j],sigma_ss[i][j],cut_ss_lj[i][j],b_ss[i][j],cut_ss_sm[i][j],
	epsilon_sb[i][j],sigma_sb[i][j],cut_sb_lj[i][j],b_sb[i][j],cut_sb_sm[i][j],
	epsilon_bb[i][j],sigma_bb[i][j],cut_bb_lj[i][j],b_bb[i][j],cut_bb_sm[i][j]);
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
  if (strcmp(str,"epsilon_sb") == 0) return (void *) epsilon_sb;
  if (strcmp(str,"sigma_sb") == 0) return (void *) sigma_sb;
  if (strcmp(str,"cut_sb_lj") == 0) return (void *) cut_sb_lj;
  if (strcmp(str,"b_sb") == 0) return (void *) b_sb;
  if (strcmp(str,"cut_sb_sm") == 0) return (void *) cut_sb_sm;
  if (strcmp(str,"epsilon_bb") == 0) return (void *) epsilon_bb;
  if (strcmp(str,"sigma_bb") == 0) return (void *) sigma_bb;
  if (strcmp(str,"cut_bb_lj") == 0) return (void *) cut_bb_lj;
  if (strcmp(str,"b_bb") == 0) return (void *) b_bb;
  if (strcmp(str,"cut_bb_sm") == 0) return (void *) cut_bb_sm;
  return NULL;
}

/* ----------------------------------------------------------------------
   calculates the lj part of the excluded volume interaction 
------------------------------------------------------------------------- */
inline void PairOxdna::calculate_excv_lj(double factor_lj, double rsq, 
	double lj1, double lj2, double & r6inv, double & fpair) 
{
  double r2inv,forcelj;

  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  forcelj = r6inv * (lj1* r6inv - lj2);
  fpair = factor_lj*forcelj*r2inv;
}

/* ----------------------------------------------------------------------
   calculates the smooting part of the excluded volume interaction 
------------------------------------------------------------------------- */
inline void PairOxdna::calculate_excv_sm(double factor_lj, double rsq, 
	double eps, double b, double cut, double & r, double & fpair) 
{
  double rinv;

  r = sqrt(rsq);
  rinv = 1.0/r;
  fpair = factor_lj*eps*2.0*b*(cut*rinv - 1.0);
}
	
