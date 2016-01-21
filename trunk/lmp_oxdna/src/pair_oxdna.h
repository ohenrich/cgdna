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
  double **epsilon_ss, **sigma_ss, **cut_ss_lj, **cutsq_ss_lj; 
  double **b_ss, **cut_ss_sm, **cutsq_ss_sm;
  double **lj1_ss, **lj2_ss, **lj3_ss, **lj4_ss, **offset_ss;
  double **epsilon_sb, **sigma_sb, **cut_sb_lj, **cutsq_sb_lj; 
  double **b_sb, **cut_sb_sm, **cutsq_sb_sm;
  double **lj1_sb, **lj2_sb, **lj3_sb, **lj4_sb, **offset_sb;
  double **epsilon_bb, **sigma_bb, **cut_bb_lj, **cutsq_bb_lj; 
  double **b_bb, **cut_bb_sm, **cutsq_bb_sm;
  double **lj1_bb, **lj2_bb, **lj3_bb, **lj4_bb, **offset_bb;

  virtual void allocate();
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
