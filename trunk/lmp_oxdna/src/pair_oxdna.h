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

#ifdef PAIR_CLASS

PairStyle(oxdna,PairOxdna)

#else

#ifndef LMP_PAIR_OXDNA_H
#define LMP_PAIR_OXDNA_H

#include "pair.h"
#include "update.h"

namespace LAMMPS_NS {

class PairOxdna : public Pair {
 public:
  PairOxdna(class LAMMPS *);
  virtual ~PairOxdna();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  void *extract(const char *, int &);

 protected:
  // s=sugar-phosphate backbone site, b=base site, st=stacking site

  // excluded volume interaction
  double **epsilon_ss, **sigma_ss, **cut_ss_ast, **cutsq_ss_ast; 
  double **lj1_ss, **lj2_ss, **b_ss, **cut_ss_c, **cutsq_ss_c;
  double **epsilon_sb, **sigma_sb, **cut_sb_ast, **cutsq_sb_ast; 
  double **lj1_sb, **lj2_sb, **b_sb, **cut_sb_c, **cutsq_sb_c;
  double **epsilon_bb, **sigma_bb, **cut_bb_ast, **cutsq_bb_ast; 
  double **lj1_bb, **lj2_bb, **b_bb, **cut_bb_c, **cutsq_bb_c;

  // stacking interaction
  double **epsilon_st, **a_st, **cut_st_0, **cut_st_c;
  double **cut_st_lo, **cut_st_hi;
  double **cut_st_lc, **cut_st_hc, **b_st1_lo, **b_st1_hi, **shift_st;
  double **cutsq_st_hc;
  double **a_st4, **theta_st4_0, **dtheta_st4_ast;
  double **b_st4, **dtheta_st4_c;

  virtual void allocate();

  // modulation factors
  inline double F3(double, double, double, double, double, double, double, double &);
  inline double F1(double, double, double, double, double, double, double, double, double, double, double);
  inline double DF1(double, double, double, double, double, double, double, double, double, double);
  inline double F4(double, double, double, double, double, double);
  inline double DF4(double, double, double, double, double, double);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
