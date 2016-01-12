/* -*- c++ -*- ----------------------------------------------------------
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
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  double **epsilon_ss, **sigma_ss, **cut_ss_lj, **b_ss, **cut_ss_sm;
  double **cutsq_ss_lj, **cutsq_ss_sm;
  double **lj1,**lj2,**lj3,**lj4,**offset;

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
