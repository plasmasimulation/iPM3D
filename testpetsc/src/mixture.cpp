/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "ctype.h"
#include "mixture.h"

#include "particle.h"

#include "memory.h"




#define DELTA 8

/* ---------------------------------------------------------------------- */

Mixture::Mixture(SPARTA *sparta, char *userid) : Pointers(sparta)
{
  // mixture ID must be all alphanumeric chars or underscores

  // int n = strlen(userid) + 1;
  // id = new char[n];
  // strcpy(id,userid);

  // for (int i = 0; i < n-1; i++)
  //   if (!isalnum(id[i]) && id[i] != '_')
      // error->all(FLERR,
      //            "Mixture ID must be alphanumeric or underscore characters");

  // special default mixtures

  // all_default = species_default = 0;
  // if (strcmp(id,"all") == 0) all_default = 1;
  // if (strcmp(id,"species") == 0) species_default = 1;

  // initialize mixture values

  nspecies = maxspecies = 0;
  species = NULL;

  nrho_flag = 0;
  vstream_flag = 0;
  temp_thermal_flag = 0;
  temp_rot_flag = 0;
  temp_vib_flag = 0;

  // fraction = NULL;
  // fraction_user = NULL;
  // fraction_flag = NULL;
  // cummulative = NULL;

  // ngroup = maxgroup = 0;
  // groups = NULL;
  // groupsize = NULL;
  // groupspecies = NULL;

  // mix2group = NULL;
  // species2group = NULL;
  // species2species = NULL;

  vscale = NULL;
  active = NULL;

  allocate();
}

/* ---------------------------------------------------------------------- */

Mixture::~Mixture()
{
  // delete [] id;

  memory->destroy(species);
  // memory->destroy(fraction);
  // memory->destroy(fraction_user);
  // memory->destroy(fraction_flag);
  // memory->destroy(cummulative);

  // delete_groups();
  // memory->sfree(groups);
  // delete [] groupsize;
  // memory->destroy(groupspecies);

  // memory->destroy(mix2group);
  // memory->destroy(species2group);
  // memory->destroy(species2species);

  memory->destroy(vscale);
  memory->destroy(active);
}

/* ----------------------------------------------------------------------
   this empty mixture becomes clone of mixture old
------------------------------------------------------------------------- */

void Mixture::copy(Mixture *old)
{
  nrho = old->nrho;
  nrho_flag = old->nrho_flag;
  nrho_user = old->nrho_user;
  vstream[0] = old->vstream[0];
  vstream[1] = old->vstream[1];
  vstream[2] = old->vstream[2];
  vstream_flag = old->vstream_flag;
  vstream_user[0] = old->vstream_user[0];
  vstream_user[1] = old->vstream_user[1];
  vstream_user[2] = old->vstream_user[2];
  temp_thermal = old->temp_thermal;
  temp_thermal_flag = old->temp_thermal_flag;
  temp_thermal_user = old->temp_thermal_user;

  nspecies = maxspecies = old->nspecies;
  allocate();

  // for (int i = 0; i < nspecies; i++) {
  //   species[i] = old->species[i];
  //   fraction[i] = old->fraction[i];
  //   fraction_flag[i] = old->fraction_flag[i];
  //   fraction_user[i] = old->fraction_user[i];
  //   mix2group[i] = old->mix2group[i];
  // }

  // for (int i = 0; i < old->ngroup; i++)
  //   add_group(old->groups[i]);
}

/* ----------------------------------------------------------------------
   process args of a mixture command
   zero of more species IDs followed by optional args
------------------------------------------------------------------------- */


void Mixture::init()
{
  // global attributes

  if (nrho_flag) nrho = nrho_user;
  else nrho = update->nrho;
  if (vstream_flag) {
    vstream[0] = vstream_user[0];
    vstream[1] = vstream_user[1];
    vstream[2] = vstream_user[2];
  } else {
    vstream[0] = update->vstream[0];
    vstream[1] = update->vstream[1];
    vstream[2] = update->vstream[2];
  }
  if (temp_thermal_flag) temp_thermal = temp_thermal_user;
  else temp_thermal = update->temp_thermal;
  if (temp_rot_flag) temp_rot = temp_rot_user;
  else temp_rot = temp_thermal;
  if (temp_vib_flag) temp_vib = temp_vib_user;
  else temp_vib = temp_thermal;

  // mixture temperarate cannot be 0.0 if streaming velocity is 0.0

  if (temp_thermal == 0.0 &&
      vstream[0] == 0.0 && vstream[1] == 0.0 && vstream[2] == 0.0)
    // error->all(FLERR,"Mixture streaming velocity and "
    //            "temperature cannot both be zero");

  // initialize all per-species fraction and cummulative values
  // account for both explicitly and implicitly set fractions

  // int err = init_fraction(fraction_flag,fraction_user,fraction,cummulative);

  // if (err) {
  //   char str[128];
  //   sprintf(str,"Mixture %s fractions exceed 1.0",id);
  //   error->all(FLERR,str);
  // }

  // vscale = factor to scale Gaussian unit variance by
  //          to get thermal distribution of velocities
  // per-species value since includes species mass

  for (int i = 0; i < nspecies; i++) {
    int index = species[i];
    vscale[i] = sqrt(2.0 * update->boltz * temp_thermal /
                     particle->species[index].mass);
  }

  // setup species2group and species2species

  // memory->destroy(species2group);
  // memory->create(species2group,particle->nspecies,"mixture:species2group");
  // for (int i = 0; i < particle->nspecies; i++) species2group[i] = -1;
  // for (int i = 0; i < nspecies; i++) species2group[species[i]] = mix2group[i];

  // memory->destroy(species2species);
  // memory->create(species2species,particle->nspecies,"mixture:species2group");
  // for (int i = 0; i < particle->nspecies; i++) species2species[i] = -1;
  // for (int i = 0; i < nspecies; i++) species2species[species[i]] = i;

  // setup groupsize

  // delete [] groupsize;
  // groupsize = new int[ngroup];
  // for (int i = 0; i < ngroup; i++) groupsize[i] = 0;
  // for (int i = 0; i < nspecies; i++) groupsize[mix2group[i]]++;

  // setup groupspecies

  // memory->destroy(groupspecies);
  // memory->create_ragged(groupspecies,ngroup,groupsize,"mixture:groupspecies");

  // for (int i = 0; i < ngroup; i++) groupsize[i] = 0;
  // for (int i = 0; i < nspecies; i++)
  //   groupspecies[mix2group[i]][groupsize[mix2group[i]]++] = species[i];
}

/* ----------------------------------------------------------------------
   set f = fraction and c = cummulative fraction for each species in mixture
   fflag[I] = 1 if species I fraction is set by fuser[I]
   otherwise species I fraction is an implicit value
   implicit value = 1/nimplicit of unset remainder,
     where nimplicit = # of unset species
   return 0 for success
   return 1 for error if sum of specified fractions > 1.0
   called by init() and also by FixInflowFile::interpolate()
------------------------------------------------------------------------- */



void Mixture::add_species(int narg, char **arg)
{
  int i,j,index;

  // activeflag = 1 if species are listed
  // active[i] = 0 if current mixture species I is not in the list
  // active[i] = 1 if species I is in the list but already existed in mixture
  // active[i] = 2 if species I was just added b/c it was in the list
  // active[i] = 3 if species is being removed from mixture via delete keyword
  //             this flag is set in params()

  if (narg) activeflag = 1;
  else activeflag = 0;
  for (i = 0; i < nspecies; i++) active[i] = 0;

  for (i = 0; i < narg; i++) {
    index = particle->find_species(arg[i]);
    // if (index < 0) error->all(FLERR,"Mixture species is not defined");
    for (j = 0; j < nspecies; j++)
      if (species[j] == index) break;
    if (j < nspecies) active[j] = 1;
    else {
      if (all_default || species_default)
        // error->all(FLERR,"Cannot add new species to mixture all or species");
      if (nspecies == maxspecies) allocate();
      active[nspecies] = 2;
      species[nspecies++] = index;
    }
  }
}


void Mixture::allocate()
{
  int old = maxspecies;
  maxspecies += DELTA;

  memory->grow(species,maxspecies,"mixture:species");
  // memory->grow(fraction,maxspecies,"mixture:species");
  // memory->grow(fraction_flag,maxspecies,"mixture:fraction_flag");
  // memory->grow(fraction_user,maxspecies,"mixture:fraction_user");
  // memory->grow(cummulative,maxspecies,"mixture:cummulative");
  // memory->grow(mix2group,maxspecies,"mixture:cummulative");
  memory->grow(vscale,maxspecies,"mixture:vscale");
  // memory->grow(active,maxspecies,"mixture:active");

  // for (int i = old; i < maxspecies; i++) {
  //   fraction_flag[i] = 0;
  //   fraction_user[i] = 0.0;
  // }
}

