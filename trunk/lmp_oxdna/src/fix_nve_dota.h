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

#ifdef FIX_CLASS

FixStyle(nve/dota,FixNVEDota)

#else

#ifndef LMP_FIX_NVE_DOTA_H
#define LMP_FIX_NVE_DOTA_H

#include "fix_nve.h"
#include "util_oxdna.h"

namespace LAMMPS_NS {

class FixNVEDota : public FixNVE {
 public:
  FixNVEDota(class LAMMPS *, int, char **);
  void init();
  void initial_integrate(int);
  void final_integrate();

 private:
  double dt,dthlf,dthlfm;
  class AtomVecEllipsoid *avec;
};

}
#endif
#endif

/* ERROR/WARNING messages:

E: Compute nve/dota requires atom style ellipsoid

Self-explanatory.

E: Fix nve/dota requires extended particles

This fix can only be used for particles with a shape setting.

*/
