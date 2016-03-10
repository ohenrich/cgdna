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
#include "molecule.h"
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
    memory->destroy(cut_ss_ast);
    memory->destroy(b_ss);
    memory->destroy(cut_ss_c);
    memory->destroy(lj1_ss);
    memory->destroy(lj2_ss);
    memory->destroy(cutsq_ss_ast);
    memory->destroy(cutsq_ss_c);

    memory->destroy(epsilon_sb);
    memory->destroy(sigma_sb);
    memory->destroy(cut_sb_ast);
    memory->destroy(b_sb);
    memory->destroy(cut_sb_c);
    memory->destroy(lj1_sb);
    memory->destroy(lj2_sb);
    memory->destroy(cutsq_sb_ast);
    memory->destroy(cutsq_sb_c);

    memory->destroy(epsilon_bb);
    memory->destroy(sigma_bb);
    memory->destroy(cut_bb_ast);
    memory->destroy(b_bb);
    memory->destroy(cut_bb_c);
    memory->destroy(lj1_bb);
    memory->destroy(lj2_bb);
    memory->destroy(cutsq_bb_ast);
    memory->destroy(cutsq_bb_c);

    memory->destroy(epsilon_st);
    memory->destroy(a_st);
    memory->destroy(cut_st_0);
    memory->destroy(cut_st_c);
    memory->destroy(cut_st_lo);
    memory->destroy(cut_st_hi);
    memory->destroy(cut_st_lc);
    memory->destroy(cut_st_hc);
    memory->destroy(b_st1_lo);
    memory->destroy(b_st1_hi);
    memory->destroy(shift_st);
    memory->destroy(cutsq_st_hc);

    memory->destroy(a_st4);
    memory->destroy(theta_st4_0);
    memory->destroy(dtheta_st4_ast);
    memory->destroy(b_st4);
    memory->destroy(dtheta_st4_c);

    memory->destroy(a_st5);
    memory->destroy(theta_st5_0);
    memory->destroy(dtheta_st5_ast);
    memory->destroy(b_st5);
    memory->destroy(dtheta_st5_c);

    memory->destroy(a_st6);
    memory->destroy(theta_st6_0);
    memory->destroy(dtheta_st6_ast);
    memory->destroy(b_st6);
    memory->destroy(dtheta_st6_c);

  }
}

/* ----------------------------------------------------------------------
   compute function for oxDNA pair interactions
   s=sugar-phosphate backbone site, b=base site, st=stacking site
------------------------------------------------------------------------- */

void PairOxdna::compute(int eflag, int vflag)
{

  double delf[3],delt[3],delta[3],deltb[3]; // force, torque increment;
  double evdwl,fpair,finc,tpair,factor_lj;
  double rtmp_s[3],rtmp_b[3],rtmp_st[3];
  double delr_ss[3],rsq_ss,delr_sb[3],rsq_sb;
  double delr_bs[3],rsq_bs,delr_bb[3],rsq_bb;
  double delr_st[3],delr_st_norm[3],rsq_st,r_st,rinv_st;
  double theta4,t4dir[3],cost4;
  double theta5p,t5pdir[3],cost5p;
  double theta6p,t6pdir[3],cost6p;

  // distances COM-backbone, COM-base, COM-stack
  double d_cs=-0.24, d_cb=0.56, d_cst=0.5; 
  // vectors COM-backbone, -base, -stack in lab frame
  double ra_cs[3],ra_cb[3],ra_cst[3];
  double rb_cs[3],rb_cb[3],rb_cst[3];

  // quaternions and Cartesian unit vectors in lab frame
  double *qa,ax[3],ay[3],az[3];
  double *qb,bx[3],by[3],bz[3];
  double *special_lj = force->special_lj;

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;

  double f1,f4t4,f4t5,f4t6;
  double df1,df4t4,df4t5,df4t6;

  int a,b,ia,ib,in,anum,bnum,atype,btype;

  int *alist,*blist,*numneigh,**firstneigh;
  int *type = atom->type;
  int *molecule = atom->molecule;
  int *ellipsoid = atom->ellipsoid;

  int nlocal = atom->nlocal;

  int newton_pair = force->newton_pair;
  int newton_bond = force->newton_bond;

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;


  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;

  anum = list->inum;
  alist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over pair interaction neighbours of my atoms

  for (ia = 0; ia < anum; ia++) {

    a = alist[ia];

    qa=bonus[a].quat;
    MathExtra::q_to_exyz(qa,ax,ay,az);

    // position of backbone site a
    ra_cs[0] = d_cs*ax[0];
    ra_cs[1] = d_cs*ax[1];
    ra_cs[2] = d_cs*ax[2];
    rtmp_s[0] = x[a][0] + ra_cs[0];
    rtmp_s[1] = x[a][1] + ra_cs[1];
    rtmp_s[2] = x[a][2] + ra_cs[2];

    // position of base site i
    ra_cb[0] = d_cb*ax[0];
    ra_cb[1] = d_cb*ax[1];
    ra_cb[2] = d_cb*ax[2];
    rtmp_b[0] = x[a][0] + ra_cb[0];
    rtmp_b[1] = x[a][1] + ra_cb[1];
    rtmp_b[2] = x[a][2] + ra_cb[2];

    atype = type[a];
    blist = firstneigh[a];
    bnum = numneigh[a];

    for (ib = 0; ib < bnum; ib++) {

      b = blist[ib];
      factor_lj = special_lj[sbmask(b)]; // = 0 for nearest neighbours
      b &= NEIGHMASK;

      btype = type[b];

      qb=bonus[b].quat;
      MathExtra::q_to_exyz(qb,bx,by,bz);

      rb_cs[0] = d_cs*bx[0];
      rb_cs[1] = d_cs*bx[1];
      rb_cs[2] = d_cs*bx[2];
      rb_cb[0] = d_cb*bx[0];
      rb_cb[1] = d_cb*bx[1];
      rb_cb[2] = d_cb*bx[2];
      rb_cst[0] = d_cst*bx[0];
      rb_cst[1] = d_cst*bx[1];
      rb_cst[2] = d_cst*bx[2];

      // vector backbone site b to a
      delr_ss[0] = rtmp_s[0] - (x[b][0] + rb_cs[0]);
      delr_ss[1] = rtmp_s[1] - (x[b][1] + rb_cs[1]);
      delr_ss[2] = rtmp_s[2] - (x[b][2] + rb_cs[2]);
      rsq_ss = delr_ss[0]*delr_ss[0] + delr_ss[1]*delr_ss[1] + delr_ss[2]*delr_ss[2];

      // vector base site b to backbone site a
      delr_sb[0] =  rtmp_s[0] - (x[b][0] + rb_cb[0]);
      delr_sb[1] =  rtmp_s[1] - (x[b][1] + rb_cb[1]);
      delr_sb[2] =  rtmp_s[2] - (x[b][2] + rb_cb[2]);
      rsq_sb = delr_sb[0]*delr_sb[0] + delr_sb[1]*delr_sb[1] + delr_sb[2]*delr_sb[2];

      // vector backbone site b to base site a
      delr_bs[0] = rtmp_b[0] - (x[b][0] + rb_cs[0]);
      delr_bs[1] = rtmp_b[1] - (x[b][1] + rb_cs[1]);
      delr_bs[2] = rtmp_b[2] - (x[b][2] + rb_cs[2]);
      rsq_bs = delr_bs[0]*delr_bs[0] + delr_bs[1]*delr_bs[1] + delr_bs[2]*delr_bs[2];

      // vector base site b to a
      delr_bb[0] = rtmp_b[0] - (x[b][0] + rb_cb[0]);
      delr_bb[1] = rtmp_b[1] - (x[b][1] + rb_cb[1]);
      delr_bb[2] = rtmp_b[2] - (x[b][2] + rb_cb[2]);
      rsq_bb = delr_bb[0]*delr_bb[0] + delr_bb[1]*delr_bb[1] + delr_bb[2]*delr_bb[2];


      // excluded volume interaction

      // backbone-backbone
      if (rsq_ss < cutsq_ss_c[atype][btype]) {

	evdwl = F3(rsq_ss,cutsq_ss_ast[atype][btype],cut_ss_c[atype][btype],lj1_ss[atype][btype],
			lj2_ss[atype][btype],epsilon_ss[atype][btype],b_ss[atype][btype],fpair);

	// knock out nearest-neighbour interaction between ss
	fpair *= factor_lj;
	evdwl *= factor_lj;

	if (evflag) ev_tally(a,b,nlocal,newton_pair,
		evdwl,0.0,fpair,delr_ss[0],delr_ss[1],delr_ss[2]);

	delf[0] = delr_ss[0]*fpair; 
	delf[1] = delr_ss[1]*fpair; 
	delf[2] = delr_ss[2]*fpair; 

        f[a][0] += delf[0];
        f[a][1] += delf[1];
        f[a][2] += delf[2];

	MathExtra::cross3(ra_cs,delf,delta);

	torque[a][0] += delta[0];
	torque[a][1] += delta[1];
	torque[a][2] += delta[2];

        if (newton_pair || b < nlocal) {

          f[b][0] -= delf[0];
          f[b][1] -= delf[1];
          f[b][2] -= delf[2];

	  MathExtra::cross3(rb_cs,delf,deltb);

	  torque[b][0] -= deltb[0];
	  torque[b][1] -= deltb[1];
	  torque[b][2] -= deltb[2];

        }

      }


      // backbone-base
      if (rsq_sb < cutsq_sb_c[atype][btype]) {

	evdwl = F3(rsq_sb,cutsq_sb_ast[atype][btype],cut_sb_c[atype][btype],lj1_sb[atype][btype],
			lj2_sb[atype][btype],epsilon_sb[atype][btype],b_sb[atype][btype],fpair);

	if (evflag) ev_tally(a,b,nlocal,newton_pair,
		evdwl,0.0,fpair,delr_sb[0],delr_sb[1],delr_sb[2]);

	delf[0] = delr_sb[0]*fpair; 
	delf[1] = delr_sb[1]*fpair; 
	delf[2] = delr_sb[2]*fpair; 

	f[a][0] += delf[0];
	f[a][1] += delf[1];
	f[a][2] += delf[2];

	MathExtra::cross3(ra_cs,delf,delta);

	torque[a][0] += delta[0];
	torque[a][1] += delta[1];
	torque[a][2] += delta[2];

	if (newton_pair || b < nlocal) {

	  f[b][0] -= delf[0];
	  f[b][1] -= delf[1];
	  f[b][2] -= delf[2];

	  MathExtra::cross3(rb_cb,delf,deltb);

	  torque[b][0] -= deltb[0];
	  torque[b][1] -= deltb[1];
	  torque[b][2] -= deltb[2];

	}

      }

      // base-backbone
      if (rsq_bs < cutsq_sb_c[atype][btype]) {

	evdwl = F3(rsq_bs,cutsq_sb_ast[atype][btype],cut_sb_c[atype][btype],lj1_sb[atype][btype],
			lj2_sb[atype][btype],epsilon_sb[atype][btype],b_sb[atype][btype],fpair);

	if (evflag) ev_tally(a,b,nlocal,newton_pair,
		evdwl,0.0,fpair,delr_bs[0],delr_bs[1],delr_bs[2]);

	delf[0] = delr_bs[0]*fpair; 
	delf[1] = delr_bs[1]*fpair; 
	delf[2] = delr_bs[2]*fpair; 

	f[a][0] += delf[0];
	f[a][1] += delf[1];
	f[a][2] += delf[2];
 
	MathExtra::cross3(ra_cb,delf,delta);

	torque[a][0] += delta[0];
	torque[a][1] += delta[1];
	torque[a][2] += delta[2];

	if (newton_pair || b < nlocal) {

	  f[b][0] -= delf[0];
	  f[b][1] -= delf[1];
	  f[b][2] -= delf[2];

	  MathExtra::cross3(rb_cs,delf,deltb);

	  torque[b][0] -= deltb[0];
	  torque[b][1] -= deltb[1];
	  torque[b][2] -= deltb[2];

	}

      }

      // base-base
      if (rsq_bb < cutsq_bb_c[atype][btype]) {

	evdwl = F3(rsq_bb,cutsq_bb_ast[atype][btype],cut_bb_c[atype][btype],lj1_bb[atype][btype],
			lj2_bb[atype][btype],epsilon_bb[atype][btype],b_bb[atype][btype],fpair);

	if (evflag) ev_tally(a,b,nlocal,newton_pair, evdwl,0.0,fpair,
			delr_bb[0],delr_bb[1],delr_bb[2]);

	delf[0] = delr_bb[0]*fpair; 
	delf[1] = delr_bb[1]*fpair; 
	delf[2] = delr_bb[2]*fpair; 

	f[a][0] += delf[0];
	f[a][1] += delf[1];
	f[a][2] += delf[2];

	MathExtra::cross3(ra_cb,delf,delta);

	torque[a][0] += delta[0];
	torque[a][1] += delta[1];
	torque[a][2] += delta[2];

	if (newton_pair || b < nlocal) {

	  f[b][0] -= delf[0];
	  f[b][1] -= delf[1];
	  f[b][2] -= delf[2];

	  MathExtra::cross3(rb_cb,delf,deltb);

	  torque[b][0] -= deltb[0];
	  torque[b][1] -= deltb[1];
	  torque[b][2] -= deltb[2];

	}

      }
      // end excluded volume interaction

    }
  }

  // loop over stacking interaction neighours using bond topology 

  for (in = 0; in < nbondlist; in++) {

    a = bondlist[in][1];
    b = bondlist[in][0];

    atype = type[a];
    btype = type[b];

    qa=bonus[a].quat;
    MathExtra::q_to_exyz(qa,ax,ay,az);
    qb=bonus[b].quat;
    MathExtra::q_to_exyz(qb,bx,by,bz);

    // vector COM a - stacking site a
    ra_cst[0] = d_cst*ax[0];
    ra_cst[1] = d_cst*ax[1];
    ra_cst[2] = d_cst*ax[2];

    // vector COM b - stacking site b
    rb_cst[0] = d_cst*bx[0];
    rb_cst[1] = d_cst*bx[1];
    rb_cst[2] = d_cst*bx[2];

    // vector stacking site b to a
    delr_st[0] = x[a][0] + ra_cst[0] - x[b][0] - rb_cst[0];
    delr_st[1] = x[a][1] + ra_cst[1] - x[b][1] - rb_cst[1] ;
    delr_st[2] = x[a][2] + ra_cst[2] - x[b][2] - rb_cst[2] ;

    rsq_st = delr_st[0]*delr_st[0] + delr_st[1]*delr_st[1] + delr_st[2]*delr_st[2];
    r_st = sqrt(rsq_st);
    rinv_st = 1.0/r_st;

    delr_st_norm[0] = delr_st[0] * rinv_st;
    delr_st_norm[1] = delr_st[1] * rinv_st;
    delr_st_norm[2] = delr_st[2] * rinv_st;

    // cosines and corrections

    cost4 = MathExtra::dot3(az,bz);
    if (cost4 >  1.0) cost4 =  1.0;
    if (cost4 < -1.0) cost4 = -1.0;
    theta4 = acos(cost4);

    cost5p  = MathExtra::dot3(delr_st_norm,az);
    if (cost5p >  1.0) cost5p =  1.0;
    if (cost5p < -1.0) cost5p = -1.0;
    theta5p = acos(cost5p);

    cost6p = MathExtra::dot3(delr_st_norm,bz);
    if (cost6p >  1.0) cost6p =  1.0;
    if (cost6p < -1.0) cost6p = -1.0;
    theta6p = acos(cost6p);



    f1 = F1(r_st, epsilon_st[atype][btype], a_st[atype][btype], cut_st_0[atype][btype], 
	cut_st_lc[atype][btype], cut_st_hc[atype][btype], cut_st_lo[atype][btype], cut_st_hi[atype][btype], 
	b_st1_lo[atype][btype], b_st1_hi[atype][btype], shift_st[atype][btype]);

    df1 = DF1(r_st, epsilon_st[atype][btype], a_st[atype][btype], cut_st_0[atype][btype], 
	cut_st_lc[atype][btype], cut_st_hc[atype][btype], cut_st_lo[atype][btype], cut_st_hi[atype][btype], 
	b_st1_lo[atype][btype], b_st1_hi[atype][btype]);

    f4t4 = F4(theta4, a_st4[atype][btype], theta_st4_0[atype][btype], dtheta_st4_ast[atype][btype], 
	b_st4[atype][btype], dtheta_st4_c[atype][btype]);  

    df4t4 = DF4(theta4, a_st4[atype][btype], theta_st4_0[atype][btype], dtheta_st4_ast[atype][btype], 
	b_st4[atype][btype], dtheta_st4_c[atype][btype]);  

    f4t5 = F4(theta5p, a_st5[atype][btype], theta_st5_0[atype][btype], dtheta_st5_ast[atype][btype], 
	b_st5[atype][btype], dtheta_st5_c[atype][btype]);  

    df4t5 = DF4(theta5p, a_st5[atype][btype], theta_st5_0[atype][btype], dtheta_st5_ast[atype][btype], 
	b_st5[atype][btype], dtheta_st5_c[atype][btype]);  

    f4t6 = F4(theta6p, a_st6[atype][btype], theta_st6_0[atype][btype], dtheta_st6_ast[atype][btype], 
	b_st6[atype][btype], dtheta_st6_c[atype][btype]);  

    df4t6 = DF4(theta6p, a_st6[atype][btype], theta_st6_0[atype][btype], dtheta_st6_ast[atype][btype], 
	b_st6[atype][btype], dtheta_st6_c[atype][btype]);  


    /* TODO: Early rejection criteria */

    // radial force
    finc  = -df1 * f4t4 * f4t5 * f4t6;
    fpair = finc;

    delf[0] = delr_st[0]*finc;
    delf[1] = delr_st[1]*finc;
    delf[2] = delr_st[2]*finc;

    // theta5p force
    if (theta5p) {

      finc   = f1 * f4t4 * df4t5 * f4t6 * rinv_st;
      fpair += finc;

      delf[0] += (az[0] - delr_st_norm[0]*cost5p) * finc;
      delf[1] += (az[1] - delr_st_norm[1]*cost5p) * finc;
      delf[2] += (az[2] - delr_st_norm[2]*cost5p) * finc;

    }

    // theta6p force
    if (theta6p) {

      finc   = f1 * f4t4 * f4t5 * df4t6 * rinv_st;
      fpair += finc;

      delf[0] += (bz[0] - delr_st_norm[0]*cost6p) * finc;
      delf[1] += (bz[1] - delr_st_norm[1]*cost6p) * finc;
      delf[2] += (bz[2] - delr_st_norm[2]*cost6p) * finc;

    }

    // increment forces and calculate according torque increments
    if (newton_bond || a < nlocal) {

      f[a][0] += delf[0];
      f[a][1] += delf[1];
      f[a][2] += delf[2];

      MathExtra::cross3(ra_cst,delf,delta);

    }

    if (newton_bond || b < nlocal) {

      f[b][0] -= delf[0];
      f[b][1] -= delf[1];
      f[b][2] -= delf[2];

      MathExtra::cross3(rb_cst,delf,deltb);

    }

    // theta4 pure torque
    if (theta4) {

      tpair = f1 * df4t4 * f4t5 * f4t6;
      MathExtra::cross3(bz,az,t4dir);

      delta[0] -= t4dir[0]*tpair;
      delta[1] -= t4dir[1]*tpair;
      delta[2] -= t4dir[2]*tpair;

      deltb[0] -= t4dir[0]*tpair;
      deltb[1] -= t4dir[1]*tpair;
      deltb[2] -= t4dir[2]*tpair;

    }

    // theta5p pure torque
    if (theta5p) {

      tpair = f1 * f4t4 * df4t5 * f4t6;
      MathExtra::cross3(delr_st_norm,az,t5pdir);

      delta[0] -= t5pdir[0] * tpair;
      delta[1] -= t5pdir[1] * tpair;
      delta[2] -= t5pdir[2] * tpair;

    }

    // theta6p pure torque
    if (theta6p) {

      tpair = f1 * f4t4 * f4t5 * df4t6;
      MathExtra::cross3(delr_st_norm,bz,t6pdir);

      deltb[0] += t6pdir[0] * tpair;
      deltb[1] += t6pdir[1] * tpair;
      deltb[2] += t6pdir[2] * tpair;

    }

    // increment torques
    if (newton_bond || a < nlocal) {

      torque[a][0] += delta[0];
      torque[a][1] += delta[1];
      torque[a][2] += delta[2];

    }
    if (newton_bond || b < nlocal) {

      torque[b][0] -= deltb[0];
      torque[b][1] -= deltb[1];
      torque[b][2] -= deltb[2];

    }


/*
double tau[3], delr0[3], sum[3];

//delf[0] = 0;
//delf[1] = 0;
//delf[2] = 0;

delr0[0] = x[b][0] - x[a][0];
delr0[1] = x[b][1] - x[a][1];
delr0[2] = x[b][2] - x[a][2];

MathExtra::cross3(delr0,delf,tau);

sum[0] = -delta[0] + deltb[0] - tau[0];
sum[1] = -delta[1] + deltb[1] - tau[1];
sum[2] = -delta[2] + deltb[2] - tau[2];

printf("Timestep %d\n", update->ntimestep);
printf("%d %d  %le %le %le\n",a,b,-delta[0],-delta[1],-delta[2]);
printf("%d %d  %le %le %le\n",a,b,deltb[0],deltb[1],deltb[2]);
printf("%d %d  %le %le %le\n",a,b,tau[0],tau[1],tau[2]);
printf("%d %d  %le %le %le\n",a,b,sum[0],sum[1],sum[2]);
printf("\n");
printf("%d %d  %le %le %le\n",a,b,-delf[0],-delf[1],-delf[2]);
printf("\n");
printf("\n");


*/


    /* TODO Energy calculation for stacking pair interaction */
    if (eflag)  evdwl = f1 * f4t4 * f4t5 * f4t6;
    if (evflag) ev_tally(a,b,nlocal,newton_bond,evdwl,0.0,fpair,delr_st[0],delr_st[1],delr_st[2]);

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
  memory->create(cut_ss_ast,n+1,n+1,"pair:cut_ss_ast");
  memory->create(b_ss,n+1,n+1,"pair:b_ss");
  memory->create(cut_ss_c,n+1,n+1,"pair:cut_ss_c");
  memory->create(lj1_ss,n+1,n+1,"pair:lj1_ss");
  memory->create(lj2_ss,n+1,n+1,"pair:lj2_ss");
  memory->create(cutsq_ss_ast,n+1,n+1,"pair:cutsq_ss_ast");
  memory->create(cutsq_ss_c,n+1,n+1,"pair:cutsq_ss_c");

  memory->create(epsilon_sb,n+1,n+1,"pair:epsilon_sb");
  memory->create(sigma_sb,n+1,n+1,"pair:sigma_sb");
  memory->create(cut_sb_ast,n+1,n+1,"pair:cut_sb_ast");
  memory->create(b_sb,n+1,n+1,"pair:b_sb");
  memory->create(cut_sb_c,n+1,n+1,"pair:cut_sb_c");
  memory->create(lj1_sb,n+1,n+1,"pair:lj1_sb");
  memory->create(lj2_sb,n+1,n+1,"pair:lj2_sb");
  memory->create(cutsq_sb_ast,n+1,n+1,"pair:cutsq_sb_ast");
  memory->create(cutsq_sb_c,n+1,n+1,"pair:cutsq_sb_c");

  memory->create(epsilon_bb,n+1,n+1,"pair:epsilon_bb");
  memory->create(sigma_bb,n+1,n+1,"pair:sigma_bb");
  memory->create(cut_bb_ast,n+1,n+1,"pair:cut_bb_ast");
  memory->create(b_bb,n+1,n+1,"pair:b_bb");
  memory->create(cut_bb_c,n+1,n+1,"pair:cut_bb_c");
  memory->create(lj1_bb,n+1,n+1,"pair:lj1_bb");
  memory->create(lj2_bb,n+1,n+1,"pair:lj2_bb");
  memory->create(cutsq_bb_ast,n+1,n+1,"pair:cutsq_bb_ast");
  memory->create(cutsq_bb_c,n+1,n+1,"pair:cutsq_bb_c");

  memory->create(epsilon_st,n+1,n+1,"pair:epsilon_st");
  memory->create(a_st,n+1,n+1,"pair:a_st");
  memory->create(cut_st_0,n+1,n+1,"pair:cut_st_0");
  memory->create(cut_st_c,n+1,n+1,"pair:cut_st_c");
  memory->create(cut_st_lo,n+1,n+1,"pair:cut_st_lo");
  memory->create(cut_st_hi,n+1,n+1,"pair:cut_st_hi");
  memory->create(cut_st_lc,n+1,n+1,"pair:cut_st_lc");
  memory->create(cut_st_hc,n+1,n+1,"pair:cut_st_hc");
  memory->create(b_st1_lo,n+1,n+1,"pair:b_st1_lo");
  memory->create(b_st1_hi,n+1,n+1,"pair:b_st1_hi");
  memory->create(shift_st,n+1,n+1,"pair:shift_st");
  memory->create(cutsq_st_hc,n+1,n+1,"pair:cutsq_st_hc");

  memory->create(a_st4,n+1,n+1,"pair:a_st4");
  memory->create(theta_st4_0,n+1,n+1,"pair:theta_st4_0");
  memory->create(dtheta_st4_ast,n+1,n+1,"pair:dtheta_st4_ast");
  memory->create(b_st4,n+1,n+1,"pair:b_st4");
  memory->create(dtheta_st4_c,n+1,n+1,"pair:dtheta_st4_c");

  memory->create(a_st5,n+1,n+1,"pair:a_st5");
  memory->create(theta_st5_0,n+1,n+1,"pair:theta_st5_0");
  memory->create(dtheta_st5_ast,n+1,n+1,"pair:dtheta_st5_ast");
  memory->create(b_st5,n+1,n+1,"pair:b_st5");
  memory->create(dtheta_st5_c,n+1,n+1,"pair:dtheta_st5_c");

  memory->create(a_st6,n+1,n+1,"pair:a_st6");
  memory->create(theta_st6_0,n+1,n+1,"pair:theta_st6_0");
  memory->create(dtheta_st6_ast,n+1,n+1,"pair:dtheta_st6_ast");
  memory->create(b_st6,n+1,n+1,"pair:b_st6");
  memory->create(dtheta_st6_c,n+1,n+1,"pair:dtheta_st6_c");

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

  if (narg != 26) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  count = 0;

  double epsilon_ss_one, sigma_ss_one;
  double cut_ss_ast_one, cut_ss_c_one, b_ss_one;

  // Excluded volume interaction
  // LJ parameters
  epsilon_ss_one = force->numeric(FLERR,arg[2]);
  sigma_ss_one = force->numeric(FLERR,arg[3]);
  cut_ss_ast_one = force->numeric(FLERR,arg[4]);

  // smoothing - determined through continuity and differentiability
  b_ss_one = 4.0/sigma_ss_one
      *(6.0*pow(sigma_ss_one/cut_ss_ast_one,7)-12.0*pow(sigma_ss_one/cut_ss_ast_one,13))
      *4.0/sigma_ss_one*(6.0*pow(sigma_ss_one/cut_ss_ast_one,7)-12.0*pow(sigma_ss_one/cut_ss_ast_one,13))
      /4.0/(4.0*(pow(sigma_ss_one/cut_ss_ast_one,12)-pow(sigma_ss_one/cut_ss_ast_one,6)));

  cut_ss_c_one = cut_ss_ast_one 
      - 2.0*4.0*(pow(sigma_ss_one/cut_ss_ast_one,12)-pow(sigma_ss_one/cut_ss_ast_one,6))
      /(4.0/sigma_ss_one*(6.0*pow(sigma_ss_one/cut_ss_ast_one,7)-12.0*pow(sigma_ss_one/cut_ss_ast_one,13)));

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon_ss[i][j] = epsilon_ss_one;
      sigma_ss[i][j] = sigma_ss_one;
      cut_ss_ast[i][j] = cut_ss_ast_one;
      b_ss[i][j] = b_ss_one;
      cut_ss_c[i][j] = cut_ss_c_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

  count = 0;
  
  double epsilon_sb_one, sigma_sb_one;
  double cut_sb_ast_one, cut_sb_c_one, b_sb_one;

  // LJ parameters
  epsilon_sb_one = force->numeric(FLERR,arg[5]);
  sigma_sb_one = force->numeric(FLERR,arg[6]);
  cut_sb_ast_one = force->numeric(FLERR,arg[7]);

  // smoothing - determined through continuity and differentiability 
  b_sb_one = 4.0/sigma_sb_one
      *(6.0*pow(sigma_sb_one/cut_sb_ast_one,7)-12.0*pow(sigma_sb_one/cut_sb_ast_one,13))
      *4.0/sigma_sb_one*(6.0*pow(sigma_sb_one/cut_sb_ast_one,7)-12.0*pow(sigma_sb_one/cut_sb_ast_one,13))
      /4.0/(4.0*(pow(sigma_sb_one/cut_sb_ast_one,12)-pow(sigma_sb_one/cut_sb_ast_one,6)));

  cut_sb_c_one = cut_sb_ast_one 
      - 2.0*4.0*(pow(sigma_sb_one/cut_sb_ast_one,12)-pow(sigma_sb_one/cut_sb_ast_one,6))
      /(4.0/sigma_sb_one*(6.0*pow(sigma_sb_one/cut_sb_ast_one,7)-12.0*pow(sigma_sb_one/cut_sb_ast_one,13)));

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon_sb[i][j] = epsilon_sb_one;
      sigma_sb[i][j] = sigma_sb_one;
      cut_sb_ast[i][j] = cut_sb_ast_one;
      b_sb[i][j] = b_sb_one;
      cut_sb_c[i][j] = cut_sb_c_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

  count = 0;
  
  double epsilon_bb_one, sigma_bb_one;
  double cut_bb_ast_one, cut_bb_c_one, b_bb_one;

  // LJ parameters
  epsilon_bb_one = force->numeric(FLERR,arg[8]);
  sigma_bb_one = force->numeric(FLERR,arg[9]);
  cut_bb_ast_one = force->numeric(FLERR,arg[10]);

  // smoothing - determined through continuity and differentiability 
  b_bb_one = 4.0/sigma_bb_one
      *(6.0*pow(sigma_bb_one/cut_bb_ast_one,7)-12.0*pow(sigma_bb_one/cut_bb_ast_one,13))
      *4.0/sigma_bb_one*(6.0*pow(sigma_bb_one/cut_bb_ast_one,7)-12.0*pow(sigma_bb_one/cut_bb_ast_one,13))
      /4.0/(4.0*(pow(sigma_bb_one/cut_bb_ast_one,12)-pow(sigma_bb_one/cut_bb_ast_one,6)));

  cut_bb_c_one = cut_bb_ast_one 
      - 2.0*4.0*(pow(sigma_bb_one/cut_bb_ast_one,12)-pow(sigma_bb_one/cut_bb_ast_one,6))
      /(4.0/sigma_bb_one*(6.0*pow(sigma_bb_one/cut_bb_ast_one,7)-12.0*pow(sigma_bb_one/cut_bb_ast_one,13)));

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon_bb[i][j] = epsilon_bb_one;
      sigma_bb[i][j] = sigma_bb_one;
      cut_bb_ast[i][j] = cut_bb_ast_one;
      b_bb[i][j] = b_bb_one;
      cut_bb_c[i][j] = cut_bb_c_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
  // excluded volume interaction done

  // stacking interaction
  count =0;

  double epsilon_st_one, a_st_one, b_st1_lo_one, b_st1_hi_one;
  double cut_st_0_one, cut_st_c_one, cut_st_lo_one, cut_st_hi_one;
  double cut_st_lc_one, cut_st_hc_one, tmp, shift_st_one;
  double a_st4_one, theta_st4_0_one, dtheta_st4_ast_one;
  double b_st4_one, dtheta_st4_c_one;
  double a_st5_one, theta_st5_0_one, dtheta_st5_ast_one;
  double b_st5_one, dtheta_st5_c_one;
  double a_st6_one, theta_st6_0_one, dtheta_st6_ast_one;
  double b_st6_one, dtheta_st6_c_one;

  epsilon_st_one = force->numeric(FLERR,arg[11]);
  a_st_one = force->numeric(FLERR,arg[12]);
  cut_st_0_one = force->numeric(FLERR,arg[13]);
  cut_st_c_one = force->numeric(FLERR,arg[14]);
  cut_st_lo_one = force->numeric(FLERR,arg[15]);
  cut_st_hi_one = force->numeric(FLERR,arg[16]);

  a_st4_one = force->numeric(FLERR,arg[17]);
  theta_st4_0_one = force->numeric(FLERR,arg[18]);
  dtheta_st4_ast_one = force->numeric(FLERR,arg[19]);
  a_st5_one = force->numeric(FLERR,arg[20]);
  theta_st5_0_one = force->numeric(FLERR,arg[21]);
  dtheta_st5_ast_one = force->numeric(FLERR,arg[22]);
  a_st6_one = force->numeric(FLERR,arg[23]);
  theta_st6_0_one = force->numeric(FLERR,arg[24]);
  dtheta_st6_ast_one = force->numeric(FLERR,arg[25]);

  b_st1_lo_one = 2*a_st_one*exp(-a_st_one*(cut_st_lo_one-cut_st_0_one))*
	2*a_st_one*exp(-a_st_one*(cut_st_lo_one-cut_st_0_one))*
	(1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))*
	(1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))/
	(4*((1-exp(-a_st_one*(cut_st_lo_one -cut_st_0_one)))*
	(1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))-
	(1-exp(-a_st_one*(cut_st_c_one -cut_st_0_one)))*
	(1-exp(-a_st_one*(cut_st_c_one-cut_st_0_one)))));

  cut_st_lc_one = cut_st_lo_one - a_st_one*exp(-a_st_one*(cut_st_lo_one-cut_st_0_one))*
	(1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))/b_st1_lo_one;

  b_st1_hi_one = 2*a_st_one*exp(-a_st_one*(cut_st_hi_one-cut_st_0_one))*
	2*a_st_one*exp(-a_st_one*(cut_st_hi_one-cut_st_0_one))*
	(1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))*
	(1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))/
	(4*((1-exp(-a_st_one*(cut_st_hi_one -cut_st_0_one)))*
	(1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))-
	(1-exp(-a_st_one*(cut_st_c_one -cut_st_0_one)))*
	(1-exp(-a_st_one*(cut_st_c_one-cut_st_0_one)))));

  cut_st_hc_one = cut_st_hi_one - a_st_one*exp(-a_st_one*(cut_st_hi_one-cut_st_0_one))*
	(1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))/b_st1_hi_one;

  tmp = 1 - exp(-(cut_st_c_one-cut_st_0_one) * a_st_one);
  shift_st_one = epsilon_st_one * tmp * tmp;

  b_st4_one = a_st4_one*a_st4_one*dtheta_st4_ast_one*dtheta_st4_ast_one/(1-a_st4_one*dtheta_st4_ast_one*dtheta_st4_ast_one);
  dtheta_st4_c_one = 1/(a_st4_one*dtheta_st4_ast_one);

  b_st5_one = a_st5_one*a_st5_one*dtheta_st5_ast_one*dtheta_st5_ast_one/(1-a_st5_one*dtheta_st5_ast_one*dtheta_st5_ast_one);
  dtheta_st5_c_one = 1/(a_st5_one*dtheta_st5_ast_one);

  b_st6_one = a_st6_one*a_st6_one*dtheta_st6_ast_one*dtheta_st6_ast_one/(1-a_st6_one*dtheta_st6_ast_one*dtheta_st6_ast_one);
  dtheta_st6_c_one = 1/(a_st6_one*dtheta_st6_ast_one);

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {

      epsilon_st[i][j] = epsilon_st_one;
      a_st[i][j] = a_st_one;
      cut_st_0[i][j] = cut_st_0_one;
      cut_st_c[i][j] = cut_st_c_one;
      cut_st_lo[i][j] = cut_st_lo_one;
      cut_st_hi[i][j] = cut_st_hi_one;
      cut_st_lc[i][j] = cut_st_lc_one;
      cut_st_hc[i][j] = cut_st_hc_one;
      b_st1_lo[i][j] = b_st1_lo_one;
      b_st1_hi[i][j] = b_st1_hi_one;
      shift_st[i][j] = shift_st_one;

      a_st4[i][j] = a_st4_one;
      theta_st4_0[i][j] = theta_st4_0_one;
      dtheta_st4_ast[i][j] = dtheta_st4_ast_one;
      b_st4[i][j] = b_st4_one;
      dtheta_st4_c[i][j] = dtheta_st4_c_one;

      a_st5[i][j] = a_st5_one;
      theta_st5_0[i][j] = theta_st5_0_one;
      dtheta_st5_ast[i][j] = dtheta_st5_ast_one;
      b_st5[i][j] = b_st5_one;
      dtheta_st5_c[i][j] = dtheta_st5_c_one;

      a_st6[i][j] = a_st6_one;
      theta_st6_0[i][j] = theta_st6_0_one;
      dtheta_st6_ast[i][j] = dtheta_st6_ast_one;
      b_st6[i][j] = b_st6_one;
      dtheta_st6_c[i][j] = dtheta_st6_c_one;

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
    error->all(FLERR,"Coefficient mixing not defined in oxDNA");
  }
  if (offset_flag) {
    error->all(FLERR,"Offset not supported in oxDNA");
  } 

  // excluded volume auxiliary parameters

  lj1_ss[i][j] = 4.0 * epsilon_ss[i][j] * pow(sigma_ss[i][j],12.0);
  lj2_ss[i][j] = 4.0 * epsilon_ss[i][j] * pow(sigma_ss[i][j],6.0);

  lj1_sb[i][j] = 4.0 * epsilon_sb[i][j] * pow(sigma_sb[i][j],12.0);
  lj2_sb[i][j] = 4.0 * epsilon_sb[i][j] * pow(sigma_sb[i][j],6.0);

  lj1_bb[i][j] = 4.0 * epsilon_bb[i][j] * pow(sigma_bb[i][j],12.0);
  lj2_bb[i][j] = 4.0 * epsilon_bb[i][j] * pow(sigma_bb[i][j],6.0);

  lj1_ss[j][i] = lj1_ss[i][j];
  lj2_ss[j][i] = lj2_ss[i][j];

  lj1_sb[j][i] = lj1_sb[i][j];
  lj2_sb[j][i] = lj2_sb[i][j];

  lj1_bb[j][i] = lj1_bb[i][j];
  lj2_bb[j][i] = lj2_bb[i][j];

  cutsq_ss_ast[i][j] = cut_ss_ast[i][j]*cut_ss_ast[i][j];
  cutsq_ss_c[i][j]  = cut_ss_c[i][j]*cut_ss_c[i][j];

  cutsq_sb_ast[i][j] = cut_sb_ast[i][j]*cut_sb_ast[i][j];
  cutsq_sb_c[i][j]  = cut_sb_c[i][j]*cut_sb_c[i][j];

  cutsq_bb_ast[i][j] = cut_bb_ast[i][j]*cut_bb_ast[i][j];
  cutsq_bb_c[i][j]  = cut_bb_c[i][j]*cut_bb_c[i][j];

  cutsq_ss_ast[j][i] = cutsq_ss_ast[i][j];
  cutsq_ss_c[j][i]  = cutsq_ss_c[i][j];

  cutsq_sb_ast[j][i] = cutsq_sb_ast[i][j];
  cutsq_sb_c[j][i]  = cutsq_sb_c[i][j];

  cutsq_bb_ast[j][i] = cutsq_bb_ast[i][j];
  cutsq_bb_c[j][i]  = cutsq_bb_c[i][j];

  // stacking auxiliary and derived parameters

  shift_st[j][i] = shift_st[i][j];
  b_st1_lo[j][i] = b_st1_lo[i][j];
  b_st1_hi[j][i] = b_st1_hi[i][j];

  cutsq_st_hc[i][j] = cut_st_hc[i][j]*cut_st_hc[i][j];
  cutsq_st_hc[j][i] = cutsq_st_hc[i][j];

  b_st4[j][i] = b_st4[i][j];
  dtheta_st4_c[j][i] = dtheta_st4_c[i][j];

  b_st5[j][i] = b_st5[i][j];
  dtheta_st5_c[j][i] = dtheta_st5_c[i][j];

  b_st6[j][i] = b_st6[i][j];
  dtheta_st6_c[j][i] = dtheta_st6_c[i][j];

  // set the master list distance cutoff
  return cut_ss_ast[i][j];

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
        fwrite(&cut_ss_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_ss[i][j],sizeof(double),1,fp);
        fwrite(&cut_ss_c[i][j],sizeof(double),1,fp);
        fwrite(&epsilon_sb[i][j],sizeof(double),1,fp);
        fwrite(&sigma_sb[i][j],sizeof(double),1,fp);
        fwrite(&cut_sb_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_sb[i][j],sizeof(double),1,fp);
        fwrite(&cut_sb_c[i][j],sizeof(double),1,fp);
        fwrite(&epsilon_bb[i][j],sizeof(double),1,fp);
        fwrite(&sigma_bb[i][j],sizeof(double),1,fp);
        fwrite(&cut_bb_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_bb[i][j],sizeof(double),1,fp);
        fwrite(&cut_bb_c[i][j],sizeof(double),1,fp);

        fwrite(&epsilon_st[i][j],sizeof(double),1,fp);
        fwrite(&a_st[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_0[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_c[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_lo[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_hi[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_lc[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_hc[i][j],sizeof(double),1,fp);
        fwrite(&b_st1_lo[i][j],sizeof(double),1,fp);
        fwrite(&b_st1_hi[i][j],sizeof(double),1,fp);
        fwrite(&shift_st[i][j],sizeof(double),1,fp);

        fwrite(&a_st4[i][j],sizeof(double),1,fp);
        fwrite(&theta_st4_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st4_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_st4[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st4_c[i][j],sizeof(double),1,fp);

        fwrite(&a_st5[i][j],sizeof(double),1,fp);
        fwrite(&theta_st5_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st5_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_st5[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st5_c[i][j],sizeof(double),1,fp);

        fwrite(&a_st6[i][j],sizeof(double),1,fp);
        fwrite(&theta_st6_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st6_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_st6[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st6_c[i][j],sizeof(double),1,fp);

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
          fread(&cut_ss_ast[i][j],sizeof(double),1,fp);
          fread(&b_ss[i][j],sizeof(double),1,fp);
          fread(&cut_ss_c[i][j],sizeof(double),1,fp);
          fread(&epsilon_sb[i][j],sizeof(double),1,fp);
          fread(&sigma_sb[i][j],sizeof(double),1,fp);
          fread(&cut_sb_ast[i][j],sizeof(double),1,fp);
          fread(&b_sb[i][j],sizeof(double),1,fp);
          fread(&cut_sb_c[i][j],sizeof(double),1,fp);
          fread(&epsilon_bb[i][j],sizeof(double),1,fp);
          fread(&sigma_bb[i][j],sizeof(double),1,fp);
          fread(&cut_bb_ast[i][j],sizeof(double),1,fp);
          fread(&b_bb[i][j],sizeof(double),1,fp);
          fread(&cut_bb_c[i][j],sizeof(double),1,fp);

          fread(&epsilon_st[i][j],sizeof(double),1,fp);
          fread(&a_st[i][j],sizeof(double),1,fp);
          fread(&cut_st_0[i][j],sizeof(double),1,fp);
          fread(&cut_st_c[i][j],sizeof(double),1,fp);
          fread(&cut_st_lo[i][j],sizeof(double),1,fp);
          fread(&cut_st_hi[i][j],sizeof(double),1,fp);
          fread(&cut_st_lc[i][j],sizeof(double),1,fp);
          fread(&cut_st_hc[i][j],sizeof(double),1,fp);
          fread(&b_st1_lo[i][j],sizeof(double),1,fp);
          fread(&b_st1_hi[i][j],sizeof(double),1,fp);
          fread(&shift_st[i][j],sizeof(double),1,fp);

	  fread(&a_st4[i][j],sizeof(double),1,fp);
	  fread(&theta_st4_0[i][j],sizeof(double),1,fp);
	  fread(&dtheta_st4_ast[i][j],sizeof(double),1,fp);
	  fread(&b_st4[i][j],sizeof(double),1,fp);
	  fread(&dtheta_st4_c[i][j],sizeof(double),1,fp);

	  fread(&a_st5[i][j],sizeof(double),1,fp);
	  fread(&theta_st5_0[i][j],sizeof(double),1,fp);
	  fread(&dtheta_st5_ast[i][j],sizeof(double),1,fp);
	  fread(&b_st5[i][j],sizeof(double),1,fp);
	  fread(&dtheta_st5_c[i][j],sizeof(double),1,fp);

	  fread(&a_st6[i][j],sizeof(double),1,fp);
	  fread(&theta_st6_0[i][j],sizeof(double),1,fp);
	  fread(&dtheta_st6_ast[i][j],sizeof(double),1,fp);
	  fread(&b_st6[i][j],sizeof(double),1,fp);
	  fread(&dtheta_st6_c[i][j],sizeof(double),1,fp);

        }

        MPI_Bcast(&epsilon_ss[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_ss[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_ss_ast[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&b_ss[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_ss_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&epsilon_sb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_sb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_sb_ast[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&b_sb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_sb_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&epsilon_bb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma_bb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_bb_ast[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&b_bb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_bb_c[i][j],1,MPI_DOUBLE,0,world);
 
        MPI_Bcast(&epsilon_st[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a_st[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_lo[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_hi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_lc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_hc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st1_lo[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st1_hi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&shift_st[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_st4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_st4_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st4_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st4_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_st5[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_st5_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st5_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st5[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st5_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_st6[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_st6_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st6_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st6[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st6_c[i][j],1,MPI_DOUBLE,0,world);

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
    fprintf(fp,"%d\
	 %g %g %g %g %g\
	 %g %g %g %g %g\
	 %g %g %g %g %g\
	 %g %g %g %g %g %g\
	 %g %g %g %g %g\
	 %g %g %g %g %g\
	 %g %g %g %g %g\
	 %g %g %g %g %g\
	 \n",i,
	epsilon_ss[i][i],sigma_ss[i][i],cut_ss_ast[i][i],b_ss[i][i],cut_ss_c[i][i],
	epsilon_sb[i][i],sigma_sb[i][i],cut_sb_ast[i][i],b_sb[i][i],cut_sb_c[i][i],
	epsilon_bb[i][i],sigma_bb[i][i],cut_bb_ast[i][i],b_bb[i][i],cut_bb_c[i][i],
	epsilon_st[i][i],a_st[i][i],cut_st_0[i][i],cut_st_c[i][i],cut_st_lo[i][i],cut_st_hi[i][i],
	cut_st_lc[i][i],cut_st_hc[i][i],b_st1_lo[i][i],b_st1_hi[i][i],shift_st[i][i],
	a_st4[i][i],theta_st4_0[i][i],dtheta_st4_ast[i][i],b_st4[i][i],dtheta_st4_c[i][i],
	a_st5[i][i],theta_st5_0[i][i],dtheta_st5_ast[i][i],b_st5[i][i],dtheta_st5_c[i][i],
	a_st6[i][i],theta_st6_0[i][i],dtheta_st6_ast[i][i],b_st6[i][i],dtheta_st6_c[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairOxdna::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d\
	 %g %g %g %g %g\
	 %g %g %g %g %g\
	 %g %g %g %g %g\
	 %g %g %g %g %g %g\
	 %g %g %g %g %g\
	 %g %g %g %g %g\
	 %g %g %g %g %g\
	 %g %g %g %g %g\
	 \n",i,j,
	epsilon_ss[i][j],sigma_ss[i][j],cut_ss_ast[i][j],b_ss[i][j],cut_ss_c[i][j],
	epsilon_sb[i][j],sigma_sb[i][j],cut_sb_ast[i][j],b_sb[i][j],cut_sb_c[i][j],
	epsilon_bb[i][j],sigma_bb[i][j],cut_bb_ast[i][j],b_bb[i][j],cut_bb_c[i][j],
	epsilon_st[i][j],a_st[i][j],cut_st_0[i][j],cut_st_c[i][j],cut_st_lo[i][j],cut_st_hi[i][j],
	cut_st_lc[i][j],cut_st_hc[i][j],b_st1_lo[i][j],b_st1_hi[i][j],shift_st[i][j],
	a_st4[i][j],theta_st4_0[i][j],dtheta_st4_ast[i][j],b_st4[i][j],dtheta_st4_c[i][j],
	a_st5[i][j],theta_st5_0[i][j],dtheta_st5_ast[i][j],b_st5[i][j],dtheta_st5_c[i][j],
	a_st6[i][j],theta_st6_0[i][j],dtheta_st6_ast[i][j],b_st6[i][j],dtheta_st6_c[i][j]);
}

/* ---------------------------------------------------------------------- */

void *PairOxdna::extract(const char *str, int &dim)
{
  dim = 2;

  if (strcmp(str,"epsilon_ss") == 0) return (void *) epsilon_ss;
  if (strcmp(str,"sigma_ss") == 0) return (void *) sigma_ss;
  if (strcmp(str,"cut_ss_ast") == 0) return (void *) cut_ss_ast;
  if (strcmp(str,"b_ss") == 0) return (void *) b_ss;
  if (strcmp(str,"cut_ss_c") == 0) return (void *) cut_ss_c;
  if (strcmp(str,"epsilon_sb") == 0) return (void *) epsilon_sb;
  if (strcmp(str,"sigma_sb") == 0) return (void *) sigma_sb;
  if (strcmp(str,"cut_sb_ast") == 0) return (void *) cut_sb_ast;
  if (strcmp(str,"b_sb") == 0) return (void *) b_sb;
  if (strcmp(str,"cut_sb_c") == 0) return (void *) cut_sb_c;
  if (strcmp(str,"epsilon_bb") == 0) return (void *) epsilon_bb;
  if (strcmp(str,"sigma_bb") == 0) return (void *) sigma_bb;
  if (strcmp(str,"cut_bb_ast") == 0) return (void *) cut_bb_ast;
  if (strcmp(str,"b_bb") == 0) return (void *) b_bb;
  if (strcmp(str,"cut_bb_c") == 0) return (void *) cut_bb_c;

  
  if (strcmp(str,"epsilon_st") == 0) return (void *) epsilon_st;
  if (strcmp(str,"a_st") == 0) return (void *) a_st;
  if (strcmp(str,"cut_st_0") == 0) return (void *) cut_st_0;
  if (strcmp(str,"cut_st_c") == 0) return (void *) cut_st_c;
  if (strcmp(str,"cut_st_lo") == 0) return (void *) cut_st_lo;
  if (strcmp(str,"cut_st_hi") == 0) return (void *) cut_st_hi;
  if (strcmp(str,"cut_st_lc") == 0) return (void *) cut_st_lc;
  if (strcmp(str,"cut_st_hc") == 0) return (void *) cut_st_hc;
  if (strcmp(str,"b_st1_lo") == 0) return (void *) b_st1_lo;
  if (strcmp(str,"b_st1_hi") == 0) return (void *) b_st1_hi;
  if (strcmp(str,"shift_st") == 0) return (void *) shift_st;

  if (strcmp(str,"a_st4") == 0) return (void *) a_st4;
  if (strcmp(str,"theta_st4_0") == 0) return (void *) theta_st4_0;
  if (strcmp(str,"dtheta_st4_ast") == 0) return (void *) dtheta_st4_ast;
  if (strcmp(str,"b_st4") == 0) return (void *) b_st4;
  if (strcmp(str,"dtheta_st4_c") == 0) return (void *) dtheta_st4_c;

  if (strcmp(str,"a_st5") == 0) return (void *) a_st5;
  if (strcmp(str,"theta_st5_0") == 0) return (void *) theta_st5_0;
  if (strcmp(str,"dtheta_st5_ast") == 0) return (void *) dtheta_st5_ast;
  if (strcmp(str,"b_st5") == 0) return (void *) b_st5;
  if (strcmp(str,"dtheta_st5_c") == 0) return (void *) dtheta_st5_c;

  if (strcmp(str,"a_st6") == 0) return (void *) a_st6;
  if (strcmp(str,"theta_st6_0") == 0) return (void *) theta_st6_0;
  if (strcmp(str,"dtheta_st6_ast") == 0) return (void *) dtheta_st6_ast;
  if (strcmp(str,"b_st6") == 0) return (void *) b_st6;
  if (strcmp(str,"dtheta_st6_c") == 0) return (void *) dtheta_st6_c;

  return NULL;
}

/* ----------------------------------------------------------------------
   f3 modulation factor, force and energy calculation
------------------------------------------------------------------------- */
inline double PairOxdna::F3(double rsq, double cutsq_ast, double cut_c, 
	double lj1, double lj2, double eps, double b, double & fpair) 
{
  double evdwl = 0.0;

  if (rsq < cutsq_ast) {
    double r2inv = 1.0/rsq;
    double r6inv = r2inv*r2inv*r2inv;
    fpair = r2inv*r6inv*(12*lj1*r6inv - 6*lj2);
    evdwl = r6inv*(lj1*r6inv-lj2);
  }
  else {
    double r = sqrt(rsq);
    double rinv = 1.0/r;
    fpair = 2*eps*b*(cut_c*rinv - 1);
    evdwl = eps*b*(cut_c-r)*(cut_c-r);
  }

  return evdwl;
}

/* ----------------------------------------------------------------------
   f1 modulation factor 
------------------------------------------------------------------------- */
inline double PairOxdna::F1(double r, double eps, double a, double cut_0,
	double cut_lc, double cut_hc, double cut_lo, double cut_hi, 
	double b_lo, double b_hi, double shift) 
{

  if (r > cut_hc) {
    return 0.0;
  }
  else if (r > cut_hi) {
    return eps * b_hi * (r-cut_hc) * (r-cut_hc);
  }
  else if (r > cut_lo) {
    double tmp = 1 - exp(-(r-cut_0) * a);
    return eps * tmp * tmp - shift;
  }
  else if (r > cut_lc) {
    return eps * b_lo * (r-cut_lc) * (r-cut_lc);
  }
  else {
    return 0.0;
  }

}

/* ----------------------------------------------------------------------
   derivative of f1 modulation factor 
------------------------------------------------------------------------- */
inline double PairOxdna::DF1(double r, double eps, double a, double cut_0,
	double cut_lc, double cut_hc, double cut_lo, double cut_hi, 
	double b_lo, double b_hi) 
{

  if (r > cut_hc) {
    return 0.0;
  }
  else if (r > cut_hi) {
    return 2 * eps * b_hi * (1 - cut_hc / r);
  }
  else if (r > cut_lo) {
    double tmp = exp(-(r-cut_0) * a);
    return 2 * eps * (1 - tmp) * tmp * a / r;
  }
  else if (r > cut_lc) {
    return 2 * eps * b_lo * (1 - cut_lc / r);
  }
  else {
    return 0.0;
  }

}

/* ----------------------------------------------------------------------
   f4 modulation factor 
------------------------------------------------------------------------- */
inline double PairOxdna::F4(double theta, double a_st, double theta_st_0, 
	double dtheta_st_ast, double b_st, double dtheta_st_c) 
{
  double dtheta = theta-theta_st_0;

  if (fabs(dtheta) > dtheta_st_c) {
    return 0.0;
  }
  else if (dtheta > dtheta_st_ast) {
    return b_st * (dtheta-dtheta_st_c)*(dtheta-dtheta_st_c);
  }
  else if(dtheta > -dtheta_st_ast) {
    return 1 - a_st * dtheta*dtheta;
  }
  else {
    return b_st * (dtheta+dtheta_st_c)*(dtheta+dtheta_st_c);
  }

}

/* ----------------------------------------------------------------------
   derivative of f4 modulation factor 
------------------------------------------------------------------------- */
inline double PairOxdna::DF4(double theta, double a_st, double theta_st_0,
        double dtheta_st_ast, double b_st, double dtheta_st_c) 
{
  double dtheta = theta-theta_st_0;

  if(fabs(dtheta) > dtheta_st_c) {
    return 0.0;
  }
  else if(dtheta > dtheta_st_ast) {
    return 2*b_st * (dtheta-dtheta_st_c) / sin(theta);
  }
  else if(dtheta > -dtheta_st_ast) {
    return -2*a_st * dtheta / sin(theta);
  }
  else {
    return 2*b_st * (dtheta+dtheta_st_c) / sin(theta);
  } 

}
