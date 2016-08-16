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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "fix_nve_dota.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

#define INERTIA 0.2          // moment of inertia prefactor for ellipsoid

/* ---------------------------------------------------------------------- */

FixNVEDota::FixNVEDota(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixNVEDota::init()
{
  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec)
    error->all(FLERR,"Compute nve/dota requires atom style ellipsoid");

  // check that all particles are finite-size ellipsoids
  // no point particles allowed, spherical is OK

  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0)
        error->one(FLERR,"Fix nve/dota requires extended particles");

  FixNVE::init();
}

/* ---------------------------------------------------------------------- */

void FixNVEDota::initial_integrate(int vflag)
{
  double *shape,*quat;
  double fquat[4],conjqm[4],inertia[3];

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  dt = update->dt;
  dthlf = 0.5 * dt;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      dthlfm = dthlf / rmass[i];
      quat = bonus[ellipsoid[i]].quat;
      shape = bonus[ellipsoid[i]].shape;

      // update momentum by 1/2 step
      v[i][0] += dthlfm * f[i][0];
      v[i][1] += dthlfm * f[i][1];
      v[i][2] += dthlfm * f[i][2];

      // update position by full step
      x[i][0] += dt * v[i][0];
      x[i][1] += dt * v[i][1];
      x[i][2] += dt * v[i][2];

      // convert angular momentum and torque in space frame into 
      // quaternion 4-momentum and 1/2 of 4-torque in body frame
      vec3_to_vec4(quat,angmom[i],conjqm);
      conjqm[0] *= 2.0;
      conjqm[1] *= 2.0;
      conjqm[2] *= 2.0;
      conjqm[3] *= 2.0;
      vec3_to_vec4(quat,torque[i],fquat);

      // update quaternion 4-momentum by 1/2 step
      conjqm[0] += dt * fquat[0];
      conjqm[1] += dt * fquat[1];
      conjqm[2] += dt * fquat[2];
      conjqm[3] += dt * fquat[3];

      // principal moments of inertia
      inertia[0] = INERTIA*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
      inertia[1] = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
      inertia[2] = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);

      // rotate quaternion and quaternion 4-momentum by full step
      no_squish_rotate(3,conjqm,quat,inertia,dthlf);
      no_squish_rotate(2,conjqm,quat,inertia,dthlf);
      no_squish_rotate(1,conjqm,quat,inertia,dt);
      no_squish_rotate(2,conjqm,quat,inertia,dthlf);
      no_squish_rotate(3,conjqm,quat,inertia,dthlf);

/* Noise

// ran = Gaussian(sigma=1, mean=0)

// v[i][0] = v[i][0] * exp(-gamma * update->dt) + sqrt(kBT/rmass[i] * (1-exp(-2*gamma*update->dt))) * ran0 
// v[i][1] = v[i][1] * exp(-gamma * update->dt) + sqrt(kBT/rmass[i] * (1-exp(-2*gamma*update->dt))) * ran1 
// v[i][2] = v[i][2] * exp(-gamma * update->dt) + sqrt(kBT/rmass[i] * (1-exp(-2*gamma*update->dt))) * ran2 
§
// M = 4 / (1/I_1 + 1/I_2 + 1/I_3) 

// Phi_1 = Pi^T S_1 q
// Phi_2 = Pi^T S_2 q
// Phi_3 = Pi^T S_3 q

// Pi += S_1 q * (exp(-Gamma * M * update->dt / (4 * I_1)) * Phi_1 + sqrt(4 * kBT * I_1 * (1 - exp(-Gamma * M * update->dt / (2 * I_1))))) * ran3
// Pi += S_2 q * (exp(-Gamma * M * update->dt / (4 * I_2)) * Phi_2 + sqrt(4 * kBT * I_2 * (1 - exp(-Gamma * M * update->dt / (2 * I_2))))) * ran4
// Pi += S_3 q * (exp(-Gamma * M * update->dt / (4 * I_3)) * Phi_3 + sqrt(4 * kBT * I_3 * (1 - exp(-Gamma * M * update->dt / (2 * I_3))))) * ran5

end Noise */

      // convert quaternion 4-momentum in body frame back to angular momentum in space frame
      vec4_to_vec3(quat,conjqm,angmom[i]);

      angmom[i][0] *= 0.5;
      angmom[i][1] *= 0.5;
      angmom[i][2] *= 0.5;

    }
}

/* ---------------------------------------------------------------------- */

void FixNVEDota::final_integrate()
{

  double *shape,*quat;
  double fquat[4],conjqm[4];
  double conjqm_dot_quat;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  dt = update->dt;
  dthlf = 0.5 * dt;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      dthlfm = dthlf / rmass[i];
      quat = bonus[ellipsoid[i]].quat;
      shape = bonus[ellipsoid[i]].shape;

      // update momentum
      v[i][0] += dthlfm * f[i][0];
      v[i][1] += dthlfm * f[i][1];
      v[i][2] += dthlfm * f[i][2];

      // convert angular momentum and torque in space frame into 
      // quaternion 4-momentum and 1/2 of 4-torque in body frame
      vec3_to_vec4(quat,angmom[i],conjqm);
      conjqm[0] *= 2.0;
      conjqm[1] *= 2.0;
      conjqm[2] *= 2.0;
      conjqm[3] *= 2.0;
      vec3_to_vec4(quat,torque[i],fquat);

      // update quaternion 4-momentum by 1/2 step
      conjqm[0] += dt * fquat[0];
      conjqm[1] += dt * fquat[1];
      conjqm[2] += dt * fquat[2];
      conjqm[3] += dt * fquat[3];

      // subtract component parallel to quaternion for improved numerical accuracy
      conjqm_dot_quat = conjqm[0]*quat[0] + conjqm[1]*quat[1] + conjqm[2]*quat[2] + conjqm[3]*quat[3];

      conjqm[0] -= conjqm_dot_quat * quat[0]; 
      conjqm[1] -= conjqm_dot_quat * quat[1]; 
      conjqm[2] -= conjqm_dot_quat * quat[2]; 
      conjqm[3] -= conjqm_dot_quat * quat[3]; 

      // convert quaternion 4-momentum in body frame back to angular momentum in space frame
      vec4_to_vec3(quat,conjqm,angmom[i]);

      angmom[i][0] *= 0.5;
      angmom[i][1] *= 0.5;
      angmom[i][2] *= 0.5;

    }
}
