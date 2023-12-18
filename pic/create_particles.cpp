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
#include "stdlib.h"
#include "string.h"
#include "create_particles.h"
 #include "particle.h"
#include "math_const.h"
#include"stdio.h"
#include"time.h"
#include<random>




using namespace MathConst;

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};   // same as Grid

#define MAXATTEMPT 1024      // max attempts to insert a particle into cut/split cell
#define EPSZERO 1.0e-14
extern  "C" {
void InitParticleOne(double *v1,int *ispecies);}
/* ---------------------------------------------------------------------- */

 CreateParticles::CreateParticles(){
  np=0;
 }

/* ---------------------------------------------------------------------- */



/* ----------------------------------------------------------------------
   create particles in parallel
   every proc creates fraction of particles for cells it owns
   cutflag determines whether to insert in all cells or only ones uncut by surfs
   account for cell weighting
   attributes of created particle depend on number of procs
------------------------------------------------------------------------- */

void CreateParticles::create_local(Particle* particle ,double lo[3],double hi[3])
{
  // int dimension = domain->dimension;

  // int me = comm->me;

  // int nglocal = grid->nlocal;
  int nglocal =particle->nlocal;
srand((unsigned)time(NULL));
std::random_device e;
std::uniform_real_distribution<double> randu(0, 1);
  // srand(1);
   double random_num1 =randu(e); 
    double random_num2 = randu(e); 
    double random_num3 = randu(e); 
  // loop over cells I own
  // only add particles to cells eligible for insertion
  // ntarget = floating point # of particles to create in one cell
  // npercell = integer # of particles to create in one cell
  // basing ntarget on accumulated volume and nprev insures Nme total creations
  // particle species = random value based on mixture fractions
  // particle velocity = stream velocity + thermal velocity

  // int *species = particle->mixture[imix]->species;
  // double *cummulative = particle->mixture[imix]->cummulative;
  // double *vstream = particle->mixture[imix]->vstream;
  // double *vscale = particle->mixture[imix]->vscale;
  // int nspecies = particle->mixture[imix]->nspecies;
  // double temp_thermal = particle->mixture[imix]->temp_thermal;
  // double temp_rot = particle->mixture[imix]->temp_rot;
  // double temp_vib = particle->mixture[imix]->temp_vib;

  int npercell,ncreate,isp,ispecies,id,pflag,subcell,icell;
  double x[3],v[3],xcell[3],vstream_variable[3];
  double ntarget,scale,rn,vn,vr,theta1,theta2,erot,evib;
  // double *lo,*hi;

  double tempscale = 1.0;
  double sqrttempscale = 1.0;
  ispecies=0;

  double volsum = 0.0;
  bigint nprev = 0;
  icell=0;
  ncreate=10000;
  // lo=new double[3];
  // hi=new double[3];
  // lo[0]=0;
  // lo[1]=0;
  // lo[2]=0;
  // hi[0]=1;
  // hi[1]=1;
  // hi[2]=1;


    // if surfs in cell, use xcell for all created particle attempts

    for (int m = 0; m < ncreate; m++) {

      // generate random position X for new particle

      // x[0] = lo[0] + random_num1  * (hi[0]-lo[0]);
      // x[1] = lo[1] + random_num2 * (hi[1]-lo[1]);
      // x[2] = lo[2] + random_num3  * (hi[2]-lo[2]);
      // printf("%f,%f,suijishu",random_num1,random_num2);
    //     random_num1 =; 
    //  random_num2 = randu(e); 
    //  random_num3 = randu(e); 

      x[0] = lo[0]+randu(e)*(hi[0]-lo[0]);
      x[1] =lo[1]+randu(e)*(hi[1]-lo[1]);
      x[2] =lo[2]+randu(e)*(hi[2]-lo[2]);
        
      // rn = random->uniform();
      // rn =randu(e); 

      // isp = 0;
      // while (cummulative[isp] < rn) isp++;
      // ispecies = species[isp];

      // if (speciesflag) {
      //   isp = species_variable(x) - 1;
      //   if (isp < 0 || isp >= nspecies) continue;
      //   ispecies = species[isp];
      // }

      // if (tempflag) {
      //   tempscale = temperature_variable(x);
      //   sqrttempscale = sqrt(tempscale);
      // }

      // vn = vscale[isp] * sqrttempscale * sqrt(rand() / RAND_MAX;);
      // vr = vscale[isp] * sqrttempscale * sqrt(rand() / RAND_MAX;);
      // theta1 = MY_2PI *(rand() / RAND_MAX);
      // theta2 = MY_2PI *(rand() / RAND_MAX;);

      //   v[0] = vstream[0] + vn*cos(theta1);
      //   v[1] = vstream[1] + vr*cos(theta2);
      //   v[2] = vstream[2] + vr*sin(theta2);
      // v[0]=rand() / RAND_MAX +10;
      // v[1]=rand() / RAND_MAX +10;
      // v[2]=rand() / RAND_MAX +10;
      //  v[0]=10000000;
      //  v[1]=-100000;
      //  v[2]=1000088;

      // erot = particle->erot(ispecies,temp_rot*tempscale,random);
      // evib = particle->evib(ispecies,temp_vib*tempscale,random);

       id++;
    InitParticleOne(v,&ispecies);//调用Fortran函数初始化粒子速度
    //  if(id<100)
    //  printf("%f,%f,%f",v[0],v[1],v[2]);
      particle->add_particle(id,ispecies,icell,x,v,erot,evib);

      
    }
}




int CreateParticles::outside_region(int dim, double *lo, double *hi)
{
  // int flag = 1;
  // if (hi[0] > region->extent_xlo &&
  //     lo[0] < region->extent_xhi) flag = 0;
  // if (hi[1] > region->extent_ylo &&
  //     lo[1] < region->extent_yhi) flag = 0;
  // if (dim == 3) {
  //   if (hi[2] > region->extent_zlo &&
  //       lo[2] < region->extent_zhi) flag = 0;
  // }
  // return flag;
}

