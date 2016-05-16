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
#include "pair_oxdna_hbond.h"
#include "mf_oxdna.h"
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
using namespace MFOxdna;

/* ---------------------------------------------------------------------- */

PairOxdnaHbond::PairOxdnaHbond(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairOxdnaHbond::~PairOxdnaHbond()
{
  if (allocated) {

    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(epsilon_hb);
    memory->destroy(a_hb);
    memory->destroy(cut_hb_0);
    memory->destroy(cut_hb_c);
    memory->destroy(cut_hb_lo);
    memory->destroy(cut_hb_hi);
    memory->destroy(cut_hb_lc);
    memory->destroy(cut_hb_hc);
    memory->destroy(b_hb_lo);
    memory->destroy(b_hb_hi);
    memory->destroy(shift_hb);

    memory->destroy(a_hb1);
    memory->destroy(theta_hb1_0);
    memory->destroy(dtheta_hb1_ast);
    memory->destroy(b_hb1);
    memory->destroy(dtheta_hb1_c);

    memory->destroy(a_hb4);
    memory->destroy(theta_hb4_0);
    memory->destroy(dtheta_hb4_ast);
    memory->destroy(b_hb4);
    memory->destroy(dtheta_hb4_c);



  }
}

/* ----------------------------------------------------------------------
   compute function for oxDNA pair interactions
   hb=hydrogen bonding site
------------------------------------------------------------------------- */

void PairOxdnaHbond::compute(int eflag, int vflag)
{

  double delf[3],delta[3],deltb[3]; // force, torque increment;
  double evdwl,fpair,finc,tpair,factor_lj;
  double delr_hb[3],delr_hb_norm[3],rsq_hb,r_hb,rinv_hb;
  double theta1,t1dir[3],cost1;
  double theta4,t4dir[3],cost4;

  // distance COM-hbonding site
  double d_chb=0.56; 
  // vectors COM-h-bonding site in lab frame
  double ra_chb[3],rb_chb[3];

  // quaternions and Cartesian unit vectors in lab frame
  double *qa,ax[3],ay[3],az[3];
  double *qb,bx[3],by[3],bz[3];

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  int *type = atom->type;
  int *molecule = atom->molecule;
  int *ellipsoid = atom->ellipsoid;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int *alist,*blist,*numneigh,**firstneigh;
  double *special_lj = force->special_lj;

  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;

  int a,b,ia,ib,anum,bnum,atype,btype;

  double f1,f4t1,f4t4;
  double df1,df4t1,df4t4;
  double tptofp;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  anum = list->inum;
  alist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over pair interaction neighbours of my atoms

  for (ia = 0; ia < anum; ia++) {

    a = alist[ia];
    atype = type[a];

    qa=bonus[a].quat;
    MathExtra::q_to_exyz(qa,ax,ay,az);

    ra_chb[0] = d_chb*ax[0];
    ra_chb[1] = d_chb*ax[1];
    ra_chb[2] = d_chb*ax[2];
  
    blist = firstneigh[a];
    bnum = numneigh[a];

    for (ib = 0; ib < bnum; ib++) {

      b = blist[ib];
      factor_lj = special_lj[sbmask(b)]; // = 0 for nearest neighbours
      b &= NEIGHMASK;

      btype = type[b];

      qb=bonus[b].quat;
      MathExtra::q_to_exyz(qb,bx,by,bz);

      rb_chb[0] = d_chb*bx[0];
      rb_chb[1] = d_chb*bx[1];
      rb_chb[2] = d_chb*bx[2];

      // vector h-bonding site b to a
      delr_hb[0] = x[a][0] + ra_chb[0] - x[b][0] - rb_chb[0];
      delr_hb[1] = x[a][1] + ra_chb[1] - x[b][1] - rb_chb[1];
      delr_hb[2] = x[a][2] + ra_chb[2] - x[b][2] - rb_chb[2];

      rsq_hb = delr_hb[0]*delr_hb[0] + delr_hb[1]*delr_hb[1] + delr_hb[2]*delr_hb[2];
      r_hb = sqrt(rsq_hb);
      rinv_hb = 1.0/r_hb;

      delr_hb_norm[0] = delr_hb[0] * rinv_hb;
      delr_hb_norm[1] = delr_hb[1] * rinv_hb;
      delr_hb_norm[2] = delr_hb[2] * rinv_hb;

      // angles and corrections

      cost1 = MathExtra::dot3(ax,bx);
      if (cost1 >  1.0) cost1 =  1.0;
      if (cost1 < -1.0) cost1 = -1.0;
      theta1 = acos(-1.0*cost1);

      cost4 = MathExtra::dot3(az,bz);
      if (cost4 >  1.0) cost4 =  1.0;
      if (cost4 < -1.0) cost4 = -1.0;
      theta4 = acos(cost4);

      f1 = F1(r_hb, epsilon_hb[atype][btype], a_hb[atype][btype], cut_hb_0[atype][btype],
	    cut_hb_lc[atype][btype], cut_hb_hc[atype][btype], cut_hb_lo[atype][btype], cut_hb_hi[atype][btype],
	    b_hb_lo[atype][btype], b_hb_hi[atype][btype], shift_hb[atype][btype]);

      f4t1 = F4(theta1, a_hb1[atype][btype], theta_hb1_0[atype][btype], dtheta_hb1_ast[atype][btype],
	    b_hb1[atype][btype], dtheta_hb1_c[atype][btype]);

      f4t4 = F4(theta4, a_hb4[atype][btype], theta_hb4_0[atype][btype], dtheta_hb4_ast[atype][btype],
	    b_hb4[atype][btype], dtheta_hb4_c[atype][btype]);

      df1 = DF1(r_hb, epsilon_hb[atype][btype], a_hb[atype][btype], cut_hb_0[atype][btype],
	    cut_hb_lc[atype][btype], cut_hb_hc[atype][btype], cut_hb_lo[atype][btype], cut_hb_hi[atype][btype],
	    b_hb_lo[atype][btype], b_hb_hi[atype][btype]);

      df4t1 = DF4(theta1, a_hb1[atype][btype], theta_hb1_0[atype][btype], dtheta_hb1_ast[atype][btype],
	    b_hb1[atype][btype], dtheta_hb1_c[atype][btype]);

      df4t4 = DF4(theta4, a_hb4[atype][btype], theta_hb4_0[atype][btype], dtheta_hb4_ast[atype][btype],
	    b_hb4[atype][btype], dtheta_hb4_c[atype][btype]);

      evdwl = f1 * f4t1 * f4t4 * factor_lj;

      // early rejection criterium
      if (evdwl) {

      // force, torque and virial contribution for forces between h-bonding sites

      fpair = 0.0;

      delf[0] = 0.0;
      delf[1] = 0.0;
      delf[2] = 0.0;

      delta[0] = 0.0;
      delta[1] = 0.0;
      delta[2] = 0.0;

      deltb[0] = 0.0;
      deltb[1] = 0.0;
      deltb[2] = 0.0;

      // radial force
      finc  = -df1 * f4t1 * f4t4 * factor_lj;
      fpair += finc;

      delf[0] += delr_hb[0] * finc;
      delf[1] += delr_hb[1] * finc;
      delf[2] += delr_hb[2] * finc;

      // increment forces, torques

      f[a][0] += delf[0];
      f[a][1] += delf[1];
      f[a][2] += delf[2];

      MathExtra::cross3(ra_chb,delf,delta);

      torque[a][0] += delta[0];
      torque[a][1] += delta[1];
      torque[a][2] += delta[2];

      if (newton_pair || b < nlocal) {

	f[b][0] -= delf[0];
	f[b][1] -= delf[1];
	f[b][2] -= delf[2];


	MathExtra::cross3(rb_chb,delf,deltb);

	torque[b][0] -= deltb[0];
	torque[b][1] -= deltb[1];
	torque[b][2] -= deltb[2];

      }

      delta[0] = 0.0;
      delta[1] = 0.0;
      delta[2] = 0.0;
      deltb[0] = 0.0;
      deltb[1] = 0.0;
      deltb[2] = 0.0;

      // pure torques not expressible as r x f 

      // theta1 torque
      if (theta1) {

	tpair = -f1 * df4t1 * f4t4 * factor_lj;
	MathExtra::cross3(ax,bx,t1dir);

	delta[0] += t1dir[0]*tpair;
	delta[1] += t1dir[1]*tpair;
	delta[2] += t1dir[2]*tpair;

	if (newton_pair || b < nlocal) {
	  deltb[0] += t1dir[0]*tpair;
	  deltb[1] += t1dir[1]*tpair;
	  deltb[2] += t1dir[2]*tpair;
	}

      }

      // theta4 torque
      if (theta4) {

	tpair = -f1 * f4t1 * df4t4 * factor_lj;
	MathExtra::cross3(bz,az,t4dir);

	delta[0] += t4dir[0]*tpair;
	delta[1] += t4dir[1]*tpair;
	delta[2] += t4dir[2]*tpair;

	if (newton_pair || b < nlocal) {
	  deltb[0] += t4dir[0]*tpair;
	  deltb[1] += t4dir[1]*tpair;
	  deltb[2] += t4dir[2]*tpair;
	}

      }

      // increment torques

      torque[a][0] += delta[0];
      torque[a][1] += delta[1];
      torque[a][2] += delta[2];

      if (newton_pair || b < nlocal) {

	torque[b][0] -= deltb[0];
	torque[b][1] -= deltb[1];
	torque[b][2] -= deltb[2];

      }

      // increment energy and virial
      if (evflag) ev_tally(a,b,nlocal,newton_pair,evdwl,0.0,fpair,delr_hb[0],delr_hb[1],delr_hb[2]); 

      }
      // end early rejection criteria

    }

  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairOxdnaHbond::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(epsilon_hb,n+1,n+1,"pair:epsilon_hb");
  memory->create(a_hb,n+1,n+1,"pair:a_hb");
  memory->create(cut_hb_0,n+1,n+1,"pair:cut_hb_0");
  memory->create(cut_hb_c,n+1,n+1,"pair:cut_hb_c");
  memory->create(cut_hb_lo,n+1,n+1,"pair:cut_hb_lo");
  memory->create(cut_hb_hi,n+1,n+1,"pair:cut_hb_hi");
  memory->create(cut_hb_lc,n+1,n+1,"pair:cut_hb_lc");
  memory->create(cut_hb_hc,n+1,n+1,"pair:cut_hb_hc");
  memory->create(b_hb_lo,n+1,n+1,"pair:b_hb_lo");
  memory->create(b_hb_hi,n+1,n+1,"pair:b_hb_hi");
  memory->create(shift_hb,n+1,n+1,"pair:shift_hb");
  memory->create(cutsq_hb_hc,n+1,n+1,"pair:cutsq_hb_hc");

  memory->create(a_hb1,n+1,n+1,"pair:a_hb1");
  memory->create(theta_hb1_0,n+1,n+1,"pair:theta_hb1_0");
  memory->create(dtheta_hb1_ast,n+1,n+1,"pair:dtheta_hb1_ast");
  memory->create(b_hb1,n+1,n+1,"pair:b_hb1");
  memory->create(dtheta_hb1_c,n+1,n+1,"pair:dtheta_hb1_c");

  memory->create(a_hb4,n+1,n+1,"pair:a_hb4");
  memory->create(theta_hb4_0,n+1,n+1,"pair:theta_hb4_0");
  memory->create(dtheta_hb4_ast,n+1,n+1,"pair:dtheta_hb4_ast");
  memory->create(b_hb4,n+1,n+1,"pair:b_hb4");
  memory->create(dtheta_hb4_c,n+1,n+1,"pair:dtheta_hb4_c");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairOxdnaHbond::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairOxdnaHbond::coeff(int narg, char **arg)
{
  int count;

  if (narg != 14) error->all(FLERR,"Incorrect args for pair coefficients in oxdna_hbond");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  // h-bonding interaction
  count = 0;

  double epsilon_hb_one, a_hb_one, cut_hb_0_one, cut_hb_c_one, cut_hb_lo_one, cut_hb_hi_one;
  double b_hb_lo_one, b_hb_hi_one, cut_hb_lc_one, cut_hb_hc_one, tmp, shift_hb_one;

  double a_hb1_one, theta_hb1_0_one, dtheta_hb1_ast_one;
  double b_hb1_one, dtheta_hb1_c_one;

  double a_hb4_one, theta_hb4_0_one, dtheta_hb4_ast_one;
  double b_hb4_one, dtheta_hb4_c_one;

  epsilon_hb_one = force->numeric(FLERR,arg[2]);
  a_hb_one = force->numeric(FLERR,arg[3]);
  cut_hb_0_one = force->numeric(FLERR,arg[4]);
  cut_hb_c_one = force->numeric(FLERR,arg[5]);
  cut_hb_lo_one = force->numeric(FLERR,arg[6]);
  cut_hb_hi_one = force->numeric(FLERR,arg[7]);

  a_hb1_one = force->numeric(FLERR,arg[8]);
  theta_hb1_0_one = force->numeric(FLERR,arg[9]);
  dtheta_hb1_ast_one = force->numeric(FLERR,arg[10]);

  a_hb4_one = force->numeric(FLERR,arg[11]);
  theta_hb4_0_one = force->numeric(FLERR,arg[12]);
  dtheta_hb4_ast_one = force->numeric(FLERR,arg[13]);

  b_hb_lo_one = 2*a_hb_one*exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one))*
	2*a_hb_one*exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one))*
	(1-exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one)))*
	(1-exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one)))/
	(4*((1-exp(-a_hb_one*(cut_hb_lo_one -cut_hb_0_one)))*
	(1-exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one)))-
	(1-exp(-a_hb_one*(cut_hb_c_one -cut_hb_0_one)))*
	(1-exp(-a_hb_one*(cut_hb_c_one-cut_hb_0_one)))));

  cut_hb_lc_one = cut_hb_lo_one - a_hb_one*exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one))*
	(1-exp(-a_hb_one*(cut_hb_lo_one-cut_hb_0_one)))/b_hb_lo_one;

  b_hb_hi_one = 2*a_hb_one*exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one))*
	2*a_hb_one*exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one))*
	(1-exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one)))*
	(1-exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one)))/
	(4*((1-exp(-a_hb_one*(cut_hb_hi_one -cut_hb_0_one)))*
	(1-exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one)))-
	(1-exp(-a_hb_one*(cut_hb_c_one -cut_hb_0_one)))*
	(1-exp(-a_hb_one*(cut_hb_c_one-cut_hb_0_one)))));

  cut_hb_hc_one = cut_hb_hi_one - a_hb_one*exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one))*
	(1-exp(-a_hb_one*(cut_hb_hi_one-cut_hb_0_one)))/b_hb_hi_one;

  tmp = 1 - exp(-(cut_hb_c_one-cut_hb_0_one) * a_hb_one);
  shift_hb_one = epsilon_hb_one * tmp * tmp;

  b_hb1_one = a_hb1_one*a_hb1_one*dtheta_hb1_ast_one*dtheta_hb1_ast_one/(1-a_hb1_one*dtheta_hb1_ast_one*dtheta_hb1_ast_one);
  dtheta_hb1_c_one = 1/(a_hb1_one*dtheta_hb1_ast_one);

  b_hb4_one = a_hb4_one*a_hb4_one*dtheta_hb4_ast_one*dtheta_hb4_ast_one/(1-a_hb4_one*dtheta_hb4_ast_one*dtheta_hb4_ast_one);
  dtheta_hb4_c_one = 1/(a_hb4_one*dtheta_hb4_ast_one);

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {

      epsilon_hb[i][j] = epsilon_hb_one;
      a_hb[i][j] = a_hb_one;
      cut_hb_0[i][j] = cut_hb_0_one;
      cut_hb_c[i][j] = cut_hb_c_one;
      cut_hb_lo[i][j] = cut_hb_lo_one;
      cut_hb_hi[i][j] = cut_hb_hi_one;
      cut_hb_lc[i][j] = cut_hb_lc_one;
      cut_hb_hc[i][j] = cut_hb_hc_one;
      b_hb_lo[i][j] = b_hb_lo_one;
      b_hb_hi[i][j] = b_hb_hi_one;
      shift_hb[i][j] = shift_hb_one;

      a_hb1[i][j] = a_hb1_one;
      theta_hb1_0[i][j] = theta_hb1_0_one;
      dtheta_hb1_ast[i][j] = dtheta_hb1_ast_one;
      b_hb1[i][j] = b_hb1_one;
      dtheta_hb1_c[i][j] = dtheta_hb1_c_one;

      a_hb4[i][j] = a_hb4_one;
      theta_hb4_0[i][j] = theta_hb4_0_one;
      dtheta_hb4_ast[i][j] = dtheta_hb4_ast_one;
      b_hb4[i][j] = b_hb4_one;
      dtheta_hb4_c[i][j] = dtheta_hb4_c_one;

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna_hbond");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairOxdnaHbond::init_style()
{
  int irequest;

  // request regular neighbor lists 

  irequest = neighbor->request(this,instance_me);

}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use regular
------------------------------------------------------------------------- */

void PairOxdnaHbond::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  if (id  > 0) error->all(FLERR,"Respa not supported");

}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairOxdnaHbond::init_one(int i, int j)
{

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Coefficient mixing not defined in oxDNA");
  }
  if (offset_flag) {
    error->all(FLERR,"Offset not supported in oxDNA");
  } 

  epsilon_hb[j][i] = epsilon_hb[i][j];
  a_hb[j][i] = a_hb[i][j];
  cut_hb_0[j][i] = cut_hb_0[i][j];
  cut_hb_c[j][i] = cut_hb_c[i][j];
  cut_hb_lo[j][i] = cut_hb_lo[i][j];
  cut_hb_hi[j][i] = cut_hb_hi[i][j];
  b_hb_lo[j][i] = b_hb_lo[i][j];
  b_hb_hi[j][i] = b_hb_hi[i][j];
  cut_hb_lc[j][i] = cut_hb_lc[i][j];
  cut_hb_hc[j][i] = cut_hb_hc[i][j];
  shift_hb[j][i] = shift_hb[i][j];

  a_hb1[j][i] = a_hb1[i][j];
  theta_hb1_0[j][i] = theta_hb1_0[i][j];
  dtheta_hb1_ast[j][i] = dtheta_hb1_ast[i][j];
  b_hb1[j][i] = b_hb1[i][j];
  dtheta_hb1_c[j][i] = dtheta_hb1_c[i][j];

  a_hb4[j][i] = a_hb4[i][j];
  theta_hb4_0[j][i] = theta_hb4_0[i][j];
  dtheta_hb4_ast[j][i] = dtheta_hb4_ast[i][j];
  b_hb4[j][i] = b_hb4[i][j];
  dtheta_hb4_c[j][i] = dtheta_hb4_c[i][j];

  cutsq_hb_hc[i][j] = cut_hb_hc[i][j]*cut_hb_hc[i][j];
  cutsq_hb_hc[j][i] = cutsq_hb_hc[i][j];

  // set the master list distance cutoff
  return cut_hb_hc[i][j];

}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdnaHbond::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {

        fwrite(&epsilon_hb[i][j],sizeof(double),1,fp);
        fwrite(&a_hb[i][j],sizeof(double),1,fp);
        fwrite(&cut_hb_0[i][j],sizeof(double),1,fp);
        fwrite(&cut_hb_c[i][j],sizeof(double),1,fp);
        fwrite(&cut_hb_lo[i][j],sizeof(double),1,fp);
        fwrite(&cut_hb_hi[i][j],sizeof(double),1,fp);
        fwrite(&cut_hb_lc[i][j],sizeof(double),1,fp);
        fwrite(&cut_hb_hc[i][j],sizeof(double),1,fp);
        fwrite(&b_hb_lo[i][j],sizeof(double),1,fp);
        fwrite(&b_hb_hi[i][j],sizeof(double),1,fp);
        fwrite(&shift_hb[i][j],sizeof(double),1,fp);

        fwrite(&a_hb1[i][j],sizeof(double),1,fp);
        fwrite(&theta_hb1_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb1_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_hb1[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb1_c[i][j],sizeof(double),1,fp);

        fwrite(&a_hb4[i][j],sizeof(double),1,fp);
        fwrite(&theta_hb4_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb4_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_hb4[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_hb4_c[i][j],sizeof(double),1,fp);

    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdnaHbond::read_restart(FILE *fp)
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

          fread(&epsilon_hb[i][j],sizeof(double),1,fp);
          fread(&a_hb[i][j],sizeof(double),1,fp);
          fread(&cut_hb_0[i][j],sizeof(double),1,fp);
          fread(&cut_hb_c[i][j],sizeof(double),1,fp);
          fread(&cut_hb_lo[i][j],sizeof(double),1,fp);
          fread(&cut_hb_hi[i][j],sizeof(double),1,fp);
          fread(&cut_hb_lc[i][j],sizeof(double),1,fp);
          fread(&cut_hb_hc[i][j],sizeof(double),1,fp);
          fread(&b_hb_lo[i][j],sizeof(double),1,fp);
          fread(&b_hb_hi[i][j],sizeof(double),1,fp);
          fread(&shift_hb[i][j],sizeof(double),1,fp);

          fread(&a_hb1[i][j],sizeof(double),1,fp);
          fread(&theta_hb1_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb1_ast[i][j],sizeof(double),1,fp);
          fread(&b_hb1[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb1_c[i][j],sizeof(double),1,fp);

          fread(&a_hb4[i][j],sizeof(double),1,fp);
          fread(&theta_hb4_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb4_ast[i][j],sizeof(double),1,fp);
          fread(&b_hb4[i][j],sizeof(double),1,fp);
          fread(&dtheta_hb4_c[i][j],sizeof(double),1,fp);

        }

        MPI_Bcast(&epsilon_hb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a_hb[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_hb_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_hb_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_hb_lo[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_hb_hi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_hb_lc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_hb_hc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_hb_lo[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_hb_hi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&shift_hb[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_hb1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_hb1_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb1_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_hb1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb1_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_hb4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_hb4_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb4_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_hb4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_hb4_c[i][j],1,MPI_DOUBLE,0,world);

      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdnaHbond::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdnaHbond::read_restart_settings(FILE *fp)
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

void PairOxdnaHbond::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d\
	 %g %g %g %g %g %g\
	 %g %g %g %g %g\
	 %g %g %g %g %g\
	 %g %g %g %g %g\
	 \n",i,
	epsilon_hb[i][i],a_hb[i][i],cut_hb_0[i][i],cut_hb_c[i][i],cut_hb_lo[i][i],cut_hb_hi[i][i],
	cut_hb_lc[i][i],cut_hb_hc[i][i],b_hb_lo[i][i],b_hb_hi[i][i],shift_hb[i][i],
        a_hb1[i][i],theta_hb1_0[i][i],dtheta_hb1_ast[i][i],b_hb1[i][i],dtheta_hb1_c[i][i],
        a_hb4[i][i],theta_hb4_0[i][i],dtheta_hb4_ast[i][i],b_hb4[i][i],dtheta_hb4_c[i][i]);

}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairOxdnaHbond::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d\
	 %g %g %g %g %g %g\
	 %g %g %g %g %g\
	 %g %g %g %g %g\
	 %g %g %g %g %g\
	 \n",i,j,
	epsilon_hb[i][j],a_hb[i][j],cut_hb_0[i][j],cut_hb_c[i][j],cut_hb_lo[i][j],cut_hb_hi[i][j],
	cut_hb_lc[i][j],cut_hb_hc[i][j],b_hb_lo[i][j],b_hb_hi[i][j],shift_hb[i][j],
        a_hb1[i][j],theta_hb1_0[i][j],dtheta_hb1_ast[i][j],b_hb1[i][j],dtheta_hb1_c[i][j],
        a_hb4[i][j],theta_hb4_0[i][j],dtheta_hb4_ast[i][j],b_hb4[i][j],dtheta_hb4_c[i][j]);

}

/* ---------------------------------------------------------------------- */

void *PairOxdnaHbond::extract(const char *str, int &dim)
{
  dim = 2;
  
  if (strcmp(str,"epsilon_hb") == 0) return (void *) epsilon_hb;
  if (strcmp(str,"a_hb") == 0) return (void *) a_hb;
  if (strcmp(str,"cut_hb_0") == 0) return (void *) cut_hb_0;
  if (strcmp(str,"cut_hb_c") == 0) return (void *) cut_hb_c;
  if (strcmp(str,"cut_hb_lo") == 0) return (void *) cut_hb_lo;
  if (strcmp(str,"cut_hb_hi") == 0) return (void *) cut_hb_hi;
  if (strcmp(str,"cut_hb_lc") == 0) return (void *) cut_hb_lc;
  if (strcmp(str,"cut_hb_hc") == 0) return (void *) cut_hb_hc;
  if (strcmp(str,"b_hb_lo") == 0) return (void *) b_hb_lo;
  if (strcmp(str,"b_hb_hi") == 0) return (void *) b_hb_hi;
  if (strcmp(str,"shift_hb") == 0) return (void *) shift_hb;

  if (strcmp(str,"a_hb1") == 0) return (void *) a_hb1;
  if (strcmp(str,"theta_hb1_0") == 0) return (void *) theta_hb1_0;
  if (strcmp(str,"dtheta_hb1_ast") == 0) return (void *) dtheta_hb1_ast;
  if (strcmp(str,"b_hb1") == 0) return (void *) b_hb1;
  if (strcmp(str,"dtheta_hb1_c") == 0) return (void *) dtheta_hb1_c;

  if (strcmp(str,"a_hb4") == 0) return (void *) a_hb4;
  if (strcmp(str,"theta_hb4_0") == 0) return (void *) theta_hb4_0;
  if (strcmp(str,"dtheta_hb4_ast") == 0) return (void *) dtheta_hb4_ast;
  if (strcmp(str,"b_hb4") == 0) return (void *) b_hb4;
  if (strcmp(str,"dtheta_hb4_c") == 0) return (void *) dtheta_hb4_c;

  return NULL;
}
