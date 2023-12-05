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

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "ctype.h"
#include "particle.h"
// #include "mixture.h"

// #include "memory.h"





enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};  // several files


Particle::Particle()
{
  // MPI_Comm_rank(world,&me);

  // exist = sorted = 0;
  nglobal = 0;
  nlocal = maxlocal = 0;
  particles = NULL;

  nspecies = maxspecies = 0;
  species = NULL;
  maxvibmode = 0;
  // lo=new double[3];
  // hi=new double[3];
  
  //maxgrid = 0;
  //cellcount = NULL;
  //first = NULL;
  maxsort = 0;
  next = NULL;

  // create two default mixtures

  nmixture = maxmixture = 0;
  // mixture = NULL;

  char **newarg = new char*[1];
  newarg[0] = (char *) "all";
  // add_mixture(1,newarg);
  newarg[0] = (char *) "species";
  // add_mixture(1,newarg);
  delete [] newarg;

  // custom per-particle vectors/arrays

  ncustom = 0;
  ename = NULL;
  etype = esize = ewhich = NULL;

  ncustom_ivec = ncustom_iarray = 0;
  icustom_ivec = icustom_iarray = NULL;
  eivec = NULL;
  eiarray = NULL;
  eicol = NULL;

  ncustom_dvec = ncustom_darray = 0;
  icustom_dvec = icustom_darray = NULL;
  edvec = NULL;
  edarray = NULL;
  edcol = NULL;

  custom_restart_flag = NULL;

  // RNG for particle weighting

  // wrandom = NULL;

  copy = copymode = 0;


}

/* ---------------------------------------------------------------------- */

Particle::~Particle()
{
  if (copy || copymode) return;

  // memory->sfree(species);
  // for (int i = 0; i < nmixture; i++) delete mixture[i];
  // memory->sfree(mixture);

  // memory->sfree(particles);
  // //memory->destroy(cellcount);
  // //memory->destroy(first);
  // memory->destroy(next);

  // for (int i = 0; i < ncustom; i++) delete [] ename[i];
  // memory->sfree(ename);
  // memory->destroy(etype);
  // memory->destroy(esize);
  // memory->destroy(ewhich);

  // for (int i = 0; i < ncustom_ivec; i++) memory->destroy(eivec[i]);
  // for (int i = 0; i < ncustom_iarray; i++) memory->destroy(eiarray[i]);
  // for (int i = 0; i < ncustom_dvec; i++) memory->destroy(edvec[i]);
  // for (int i = 0; i < ncustom_darray; i++) memory->destroy(edarray[i]);

  // memory->destroy(icustom_ivec);
  // memory->destroy(icustom_iarray);
  // memory->sfree(eivec);
  // memory->sfree(eiarray);
  // memory->destroy(eicol);
  // memory->destroy(icustom_dvec);
  // memory->destroy(icustom_darray);
  // memory->sfree(edvec);
  // memory->sfree(edarray);
  // memory->destroy(edcol);

  // delete wrandom;
}

/* ---------------------------------------------------------------------- */

void Particle::init(double lo[3],double hi[3])
{
  this->lo[0]=lo[0];
  this->lo[1]=lo[1];
  this->lo[2]=lo[2];
  this->hi[0]=hi[0];
  this->hi[1]=hi[1];
  this->hi[2]=hi[2];

}



/* ----------------------------------------------------------------------
   set the initial weight of each particle
   called by Update before particle move
   only called if particle weighting is enabled
   only grid-based weighting is currently implemented
------------------------------------------------------------------------- */

void Particle::pre_weight()
{
  // int icell;
  // Grid::ChildInfo *cinfo = grid->cinfo;

  // for (int i = 0; i < nlocal; i++) {
  //   icell = particles[i].icell;
  //   particles[i].weight = cinfo[icell].weight;
  // }
}

/* ----------------------------------------------------------------------
   clone/delete each particle based on ratio of its initial/final weights
   called by Update after particle move and migration
   only called if particle weighting is enabled
   only grid-based weighting is currently implemented
------------------------------------------------------------------------- */

void Particle::post_weight()
{
  // int m,icell,nclone;
  // double ratio,fraction;

  // int nbytes = sizeof(OnePart);
  // Grid::ChildInfo *cinfo = grid->cinfo;

  // // nlocal_original-1 = index of last original particle

  // int nlocal_original = nlocal;
  // int i = 0;

  // while (i < nlocal_original) {
  //   icell = particles[i].icell;

  //   // next particle will be an original particle
  //   // skip it if no weight change

  //   if (particles[i].weight == cinfo[icell].weight) {
  //     i++;
  //     continue;
    // }

    // // ratio < 1.0 is candidate for deletion
    // // if deleted and particle that takes its place is cloned (Nloc > Norig)
    // //   then skip it via i++, else will examine it on next iteration

    // ratio = particles[i].weight / cinfo[icell].weight;

    // if (ratio < 1.0) {
    //   if (wrandom->uniform() > ratio) {
    //     memcpy(&particles[i],&particles[nlocal-1],nbytes);
    //     if (ncustom) copy_custom(i,nlocal-1);
    //     if (nlocal > nlocal_original) i++;
    //     else nlocal_original--;
    //     nlocal--;
    //   } else i++;
    //   continue;
    // }

    // ratio > 1.0 is candidate for cloning
    // create Nclone new particles each with unique ID

  //   nclone = static_cast<int> (ratio);
  //   fraction = ratio - nclone;
  //   nclone--;
  //   if (wrandom->uniform() < fraction) nclone++;

  //   for (m = 0; m < nclone; m++) {
  //     clone_particle(i);
  //     particles[nlocal-1].id = MAXSMALLINT*wrandom->uniform();
  //   }
  //   i++;
  // }
}

/* ----------------------------------------------------------------------
   insure particle list can hold nextra new particles
   if defined, also grow custom particle arrays and initialize with zeroes
------------------------------------------------------------------------- */

void Particle::grow(int nextra)
{
  bigint target = (bigint) nlocal + nextra;
  if (target <= maxlocal) return;

  int oldmax = maxlocal;
  bigint newmax = maxlocal;
  while (newmax < target) newmax += 1234;
//DELTA
  // if (newmax > MAXSMALLINT)
  //   error->one(FLERR,"Per-processor particle count is too big");

  maxlocal = newmax;
  particles = (OnePart *)
    // memory->srealloc(particles,maxlocal*sizeof(OnePart),
    //                  "particle:particles",SPARTA_GET_ALIGN(OnePart));
  realloc(particles,maxlocal*sizeof(OnePart));
  memset(&particles[oldmax],0,(maxlocal-oldmax)*sizeof(OnePart));

  // if (ncustom == 0) return;

  // for (int i = 0; i < ncustom; i++) {
  //   if (ename[i] == NULL) continue;
  //   grow_custom(i,oldmax,maxlocal);
  // }
}

/* ----------------------------------------------------------------------
   insure species list can hold maxspecies species
   assumes that maxspecies has already been increased
------------------------------------------------------------------------- */

void Particle::grow_species()
{
  species = (Species *)
    // memory->srealloc(species,maxspecies*sizeof(Species),"particle:species");
      realloc(species,maxspecies*sizeof(Species));
}

/* ----------------------------------------------------------------------
   grow next list if more particles now exist than there is room for
   called from Grid::unpack_particles_adapt() when grid adaptation
     takes place and acquire particles from other procs due to coarsening
   unlike sort(), this requires next list be grown, not destroy/create
     b/c sorted particle list is maintained during adaptation
------------------------------------------------------------------------- */

void Particle::grow_next()
{
  // compare maxsort (length of next) to new particle count (nlocal)
  // grow to maxlocal (max length of particles) to avoid frequent re-grow

  // if (maxsort < nlocal) {
  //   maxsort = maxlocal;
  //   memory->grow(next,maxsort,"particle:next");
  // }
}

/* ----------------------------------------------------------------------
   add a particle to particle list
   return 1 if particle array was reallocated, else 0
------------------------------------------------------------------------- */

int Particle::add_particle(int id, int ispecies, int icell,
                           double *x, double *v, double erot, double evib)
{
  int reallocflag = 0;
  if (nlocal == maxlocal) {
    grow(1);
    reallocflag = 1;
  }

  OnePart *p = &particles[nlocal];

  p->id = id;
  p->ispecies = ispecies;
  p->icell = icell;
  p->x[0] = x[0];
  p->x[1] = x[1];
  p->x[2] = x[2];
  p->v[0] = v[0];
  p->v[1] = v[1];
  p->v[2] = v[2];
  p->erot = erot;
  p->evib = evib;
  p->flag = PKEEP;

  //p->dtremain = 0.0;    not needed due to memset in grow() ??
  //p->weight = 1.0;      not needed due to memset in grow() ??

  nlocal++;
  return reallocflag;
}

/* ----------------------------------------------------------------------
   add an empty particle to particle list, caller will fill it
   return 1 if particle array was reallocated, else 0
------------------------------------------------------------------------- */

int Particle::add_particle()
{
  int reallocflag = 0;
  if (nlocal == maxlocal) {
    grow(1);
    reallocflag = 1;
  }

  nlocal++;
  return reallocflag;
}

/* ----------------------------------------------------------------------
   clone particle Index and add to end of particle list
   return 1 if particle array was reallocated, else 0
------------------------------------------------------------------------- */

int Particle::clone_particle(int index)
{
  int reallocflag = 0;
  if (nlocal == maxlocal) {
    grow(1);
    reallocflag = 1;
  }

  memcpy(&particles[nlocal],&particles[index],sizeof(OnePart));
  // if (ncustom) copy_custom(nlocal,index);

  nlocal++;
  return reallocflag;
}

/* ----------------------------------------------------------------------
   return index of ID in list of species IDs
   return -1 if not found
------------------------------------------------------------------------- */

int Particle::find_species(char *id)
{
  for (int i = 0; i < nspecies; i++)
    if (strcmp(id,species[i].id) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   return index of ID in list of species IDs
   return -1 if not found
------------------------------------------------------------------------- */

void Particle::particle_move_comm( ){

double dt;
// local=particle->nlocal;
dt=0.5;
char *psend=NULL,*precv=NULL;
int sendnum,recnum,xyz_np[3],rank,index;
int **plist,l[6],indexmax[6];
plist=new int*[6];
int particle_number_send[2],particle_number_recv[2];
int negb_rank[6];
int boundarytest,temp1,temp2,test[2];
int offset = 0;
int nbytes_particle = sizeof(Particle::OnePart);
for(int i = 0; i < 6; i++) {  
    indexmax[i]=100;
    plist[i] = new int[indexmax[i]]; 
    l[i]=0;
    
}
xyz_np[0]=2;
xyz_np[1]=2;
xyz_np[2]=2;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  

for (int i = 0; i<nlocal ; i++) {
   if(particles[i].flag==1){
    while(particles[nlocal-1].flag==1&&nlocal>1)
    nlocal--;
     memcpy(&particles[i],&particles[nlocal-1],sizeof(OnePart));
     nlocal--;
     continue;
   }
particles[i].x[0]+=dt*particles[i].v[0];
particles[i].x[1]+=dt*particles[i].v[1];
particles[i].x[2]+=dt*particles[i].v[2];

 index=particle_domain_index(&particles[i]);
// index=-1;
if(index>=0)
{
 if(l[index]==indexmax[index])
 { delete[] plist[index];
 indexmax[index]=indexmax[index]*2;
 plist[index]=new int[indexmax[index]];
 }
  plist[index][l[index]]=i;
  l[index]++;
}

}
// if(nlocal==0)
// return ;
for(int dim=0;dim<=2;dim++)
  { test[0]=0;
    test[1]=0;
    particle_number_send[0]=0;
    particle_number_send[1]=0;
    particle_number_recv[0]=0;
    particle_number_recv[1]=0;
    switch(dim){
       case 0:
                    
             negb_rank[0]=rank-1;
             negb_rank[1]=rank+1;
             boundarytest=1;
             break;
       case 1:
             negb_rank[0]=rank-xyz_np[0];
             negb_rank[1]=rank+xyz_np[0];
             boundarytest=xyz_np[0];
             break;
      case 2:
             negb_rank[0]=rank-xyz_np[0]*xyz_np[1];
             negb_rank[1]=rank+xyz_np[0]*xyz_np[1];
             boundarytest=xyz_np[0]*xyz_np[1];
              break;
       default:}
                
   temp1=(rank/boundarytest)%xyz_np[dim];
  for(int i = 0; i<=1;i++)      
     {temp2=(negb_rank[i]/boundarytest)%xyz_np[dim];   
      if(i==0) 
          {if (temp1-temp2==1&&temp2>=0) test[i]=1;}//标记检测结果
      else
          {if (temp2-temp1==1) test[i]=1;}}

                        // !此部分利用取余操作判断了边界上的网格所包含的邻居节点个数，为避免重复判断将结果保存在test数组中
                       
 // 发送消息给其他进程
    int message = 3;

    // int dest = (world_rank + 1) % world_size;
    for(int i=0;i<=1;i++)
    {if(test[i]==1) 
    {                
     particle_number_send[i]=l[i+dim*2];
                       
    MPI_Send(&particle_number_send[i], 1, MPI_INT,negb_rank[i], 0, MPI_COMM_WORLD);
    MPI_Recv(&particle_number_recv[i], 1, MPI_INT,negb_rank[i], 0, MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
  }}
  if(particle_number_recv[0]>1)
  printf("particle_number_recv,%d",particle_number_recv[0]);
    
     for(int i=0;i<=1;i++)
{if(test[i]==1) {
  if(particle_number_send[i]>0)
  
  psend = (char *) malloc(nbytes_particle*particle_number_send[i]);
  if(particle_number_recv[i]>0)
  precv = (char *) malloc(nbytes_particle*particle_number_recv[i]);
}}
  for(int i = 0;i<=1;i++)
   {if(test[i]==1)
    {if(l[i+dim*2]>=1) 
  {offset=0;
  for(int j=0;j<l[i+dim*2];j++)
   { //  { psend[j]=particles[plist[i+dim*2][j]]; 
      memcpy(&psend[offset],&particles[plist[i+dim*2][j]],nbytes_particle);
      offset += nbytes_particle;
     particles[plist[i+dim*2][j]].flag=1;
         printf("testx[0],%.2f  \t", particles[plist[i+dim*2][j]].x[0]); 
         }        
 }}}
 for(int i=0;i<=1;i++)
 { if(test[i]==1&&l[i+dim*2]>=1) {
  offset=0;
    if (particle_number_send[i] > 0) 
     MPI_Send(&psend[offset], particle_number_send[i]*nbytes_particle, MPI_CHAR, negb_rank[i],0, MPI_COMM_WORLD);

     if (particle_number_recv[i] > 0) 
    MPI_Recv(&precv[offset], particle_number_recv[i]*nbytes_particle, MPI_CHAR, negb_rank[i], 0,  MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
 } }                          
// }
 

for(int i=0;i<=1;i++)
{int start=nlocal;
offset=0;
if(test[i]==1)
  {
  for(int j=0;j<particle_number_recv[i];j++)
{
  // particles[nlocal]=precv[j];
  
  memcpy(&particles[nlocal],&precv[offset],nbytes_particle);
      offset += nbytes_particle;
      // nlocal++;
       add_particle();
}
  for(int j=0;j<particle_number_recv[i];j++)
{index=particle_domain_index(&particles[j+start]);
// index=-1;
if(index>=0)
{
 if(l[index]==indexmax[index])
 { delete[] plist[index];
 indexmax[index]=indexmax[index]*2;
 plist[index]=new int[indexmax[index]];
 }
  plist[index][l[index]]=i;
  l[index]++;
}


  }}
  // if (psend != NULL)
  // free(psend);
  // if (precv != NULL)
  // free(precv);
}}

 
// int message = world_rank;
//     int dest = (world_rank + 1) % world_size;
//     MPI_Send(&message, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
}

int Particle::particle_domain_index(OnePart* particle){
if (particle->x[0] < lo[0])  return 0 ; 
else if (particle->x[0] > hi[0])  return 1 ;
else if (particle->x[1] < lo[1])  return 2 ; 
else if (particle->x[1] > hi[1])  return 3 ;
else if (particle->x[2] < lo[2])  return 4 ; 
else if (particle->x[2] > hi[2])  return 5 ;
return -1;
}
