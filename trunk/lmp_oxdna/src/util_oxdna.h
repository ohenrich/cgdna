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

#ifndef UTIL_OXDNA_H
#define UTIL_OXDNA_H

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "error.h"

namespace UtilOxdna {

  inline void vec3_to_vec4(const double *, const double *, double *);
  inline void vec4_to_vec3(const double *, const double *, double *);

}

/* ----------------------------------------------------------------------
  Converts a 3-vector in space frame into a 4-vector in body frame 
------------------------------------------------------------------------- */

inline void UtilOxdna::vec3_to_vec4(const double * q, const double * v3, double * v4) 
{
  v4[0] = -q[1]*v3[0] - q[2]*v3[1] - q[3]*v3[2]; 
  v4[1] =  q[0]*v3[0] + q[3]*v3[1] - q[2]*v3[2]; 
  v4[2] = -q[3]*v3[0] + q[0]*v3[1] + q[1]*v3[2]; 
  v4[3] =  q[2]*v3[0] - q[1]*v3[1] + q[0]*v3[2]; 
}

/* ----------------------------------------------------------------------
  Converts a 4-vector in body frame into a 3-vector in space frame 
------------------------------------------------------------------------- */

inline void UtilOxdna::vec4_to_vec3(const double * q, const double * v4, double * v3) 
{
  v3[0] = -q[1]*v4[0] + q[0]*v4[1] - q[3]*v4[2] + q[2]*v4[3];
  v3[1] = -q[2]*v4[0] + q[3]*v4[1] + q[0]*v4[2] - q[1]*v4[3];
  v3[2] = -q[3]*v4[0] - q[2]*v4[1] + q[1]*v4[2] + q[0]*v4[3];
}

#endif    
