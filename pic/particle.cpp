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
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "ctype.h"
#include "particle.h"
#include "math_const.h"
#include <cstring> 
#include<random>
#include<malloc.h>
// #include "mixture.h"

// #include "memory.h"
using namespace MathConst;
extern "C" {
void getE(double* x, double* y,double* z);
void MCC(double*x1,double*v1,double*x2,double*v2,double*x3,double*v3,int*flag);
void weighting(double*x,double*y,double*z,int *weight);

}

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

void Particle::init(double lo[3],double hi[3],double* dx,double dt)
{
  this->lo[0]=lo[0];
  this->lo[1]=lo[1];
  this->lo[2]=lo[2];
  this->hi[0]=hi[0];
  this->hi[1]=hi[1];
  this->hi[2]=hi[2];
  this->dx=dx[0];
  this->dy=dx[1];
  this->dz=dx[2];
  this->dt=dt;

}

  void  Particle:: add_species(int num, char **name,double *coll_ratio)
  {
    
    species=(Species *)
    // memory->srealloc(particles,maxlocal*sizeof(OnePart),
    //                  "particle:particles",SPARTA_GET_ALIGN(OnePart));
  realloc(species,2*sizeof(Species));
    for(int i=0;i<num;i++)
    {
      
       strcpy(species[i].id,name[i]);
       species[i].coll_ratio=coll_ratio[i];
      //  printf("%s",species[i].id);
      // printf("%f,collisionratio\n", species[i].coll_ratio);
  }
  species[0].charge=-1.6022e-19*6103.514;
  species[0].mass=9.1095e-31*6103.514;
  species[1].charge=1.6022e-19*6103.514;
  species[1].mass=39.948*AtomicMass*6103.514;


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
  // maxlocal=100000;
  if (target <= maxlocal) return;

  int oldmax = maxlocal;
  bigint newmax = maxlocal;
  while (newmax < target) newmax += 42344;
//DELTA
  // if (newmax > MAXSMALLINT)
  //   error->one(FLERR,"Per-processor particle count is too big");

  maxlocal = newmax;
  // particles = (OnePart *)  
  // realloc(particles,maxlocal*sizeof(OnePart));
  // memory->srealloc(particles,maxlocal*sizeof(OnePart),
    //                  "particle:particles",SPARTA_GET_ALIGN(OnePart));
    // bigint nbytes=maxlocal*sizeof(OnePart);
    // int align=64;

   particles = (OnePart *)
    srealloc(particles,maxlocal*sizeof(OnePart),
                     "particle:particles",64);
//  ptr = realloc(ptr,nbytes);

  // if (ptr == NULL) {
  //   char str[128];
  //   sprintf(str,"Failed to reallocate " BIGINT_FORMAT " bytes for array %s",
  //           nbytes,name);
    // error->one(FLERR,str);
  // }



  memset(&particles[oldmax],0,(maxlocal-oldmax)*sizeof(OnePart));

  // if (ncustom == 0) return;

  // for (int i = 0; i < ncustom; i++) {
  //   if (ename[i] == NULL) continue;
  //   grow_custom(i,oldmax,maxlocal);
  // }
}

void *Particle::smalloc(bigint nbytes, const char *name, int align)
{
  if (nbytes == 0) return NULL;

  void *ptr;

  if (align) {
    int retval = posix_memalign(&ptr, align, nbytes);
    if (retval) ptr = NULL;
  } else {
    ptr = malloc(nbytes);
  }

  if (ptr == NULL) {
    char str[128];
    // printf(str,"Failed to allocate " BIGINT_FORMAT " bytes for array %s", nbytes,name);
    // error->one(FLERR,str);
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   safe realloc
------------------------------------------------------------------------- */

void *Particle::srealloc(void *ptr, bigint nbytes, const char *name, int align)
{
  if (nbytes == 0) {
    // destroy(ptr);
    return NULL;
  }

  if (align) {
    ptr = realloc(ptr, nbytes);
    uintptr_t offset = ((uintptr_t)(const void *)(ptr)) % align;
    if (offset) {
      void *optr = ptr;
      ptr = smalloc(nbytes, name, align);
#if defined(__APPLE__)
      memcpy(ptr, optr, MIN(nbytes,malloc_size(optr)));
#else
      memcpy(ptr, optr, MIN(nbytes,malloc_usable_size(optr)));
#endif
      free(optr);
    }
  } else
    ptr = realloc(ptr,nbytes);

  if (ptr == NULL) {
    char str[128];
    // sprintf(str,"Failed to reallocate " BIGINT_FORMAT " bytes for array %s",
    //         nbytes,name);
    // error->one(FLERR,str);
  }
  return ptr;
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
  p->v[0] = v[0]*dt/dx;
  p->v[1] = v[1]*dt/dy;
  p->v[2] = v[2]*dt/dz;
  p->a[0] =0;
  p->a[1] =0;
  p->a[2] =0;
  // p->erot = erot;
  // p->evib = evib;
  p->flag = PKEEP;
  // weighting(&x[0],&x[1],&x[2],&ispecies);

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

// double dt;
// local=particle->nlocal;
// dt=0.5;
char *psend[2]={NULL,NULL},*precv[2]={NULL,NULL};
int sendnum,recnum,xyz_np[3],rank,index;
int **plist,l[6],indexmax[6];
plist=new int*[6];
int particle_number_send[2],particle_number_recv[2];
int negb_rank[6],speciesid;
int boundarytest,temp1,temp2,test[2];
int offset = 0;
int nbytes_particle = sizeof(OnePart);
double Ex,Ey,Ez;
int int_x,int_y,int_z;
double double_x,double_y,double_z ; 
double EFactor[2];
std::random_device rd;
std::mt19937 e(rd());
std::uniform_real_distribution<double> randu(0, 1);
for(int i=0;i<2;i++)
{EFactor[i]=species[i].charge/species[i].mass*dt/(dx/dt)/dx;}
for(int i = 0; i < 6; i++) {  
    indexmax[i]=10000;
    plist[i]=(int*)malloc(indexmax[i]*sizeof(int));
  memset(plist[i],0,indexmax[i]*sizeof(int));
    // plist[i] = new int[indexmax[i]]; 
    l[i]=0;
    
}
xyz_np[0]=2;
xyz_np[1]=2;
xyz_np[2]=2;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  

for (int i = 0; i<nlocal ; i++) {
   




//move
   Ex=particles[i].x[0]; //假定x,y,z是真实电子坐标，设定粒子相对网格位置；
   Ey=particles[i].x[1];
   Ez=particles[i].x[2];
   //传给E需要的粒子相对网格位置,返回对应的电场大小；
  //  Ex=lo[0];
  //  Ey=lo[1];
  //  Ez=lo[2];
  //  printf("Ez,%f",Ez);
  // index=particle_domain_index(&particles[i]);
  // if(index<0)
  //   getE(&Ex,&Ey,&Ez);
  //   else 
  //   printf("not 0");
   getE(&Ex,&Ey,&Ez);
speciesid=particles[i].ispecies;
  //  double EFactor= species[particles[i].ispecies].charge* ECharge/EMass*dt/(dx/dt);
   Ex*=EFactor[speciesid];
   Ey*=EFactor[speciesid];
   Ez*=EFactor[speciesid];
  
   particles[i].v[0]+=Ex;
   particles[i].v[1]+=Ey;
   particles[i].v[2]+=Ez;

   particles[i].x[0]+=Ex;
   particles[i].x[1]+=Ey;
   particles[i].x[2]+=Ez;
   
   particles[i].a[0]=0.5*(particles[i].a[0]+Ex);
   particles[i].a[1]=0.5*(particles[i].a[1]+Ey);
   particles[i].a[2]=0.5*(particles[i].a[2]+Ez);
  
   particles[i].v[0]+=particles[i].a[0];
   particles[i].v[1]+=particles[i].a[1];
   particles[i].v[2]+=particles[i].a[2];
particles[i].x[0]+=particles[i].v[0];
particles[i].x[1]+=particles[i].v[1];
particles[i].x[2]+=particles[i].v[2];
//  if(i<6)
//    printf("%f  %f %f\n",particles[i].x[0],particles[i].v[0],Ez*1e6);
index=particle_domain_index(&particles[i]);
//particle comm
//  index=particle_domain_index(&particles[i]);
// index=-1;
if(index>=0)
{
 if(l[index]==indexmax[index])
{//  { delete[] plist[index];

//  printf("(sizeof(int))*indexmax[index]%d",(sizeof(int))*indexmax[index]);
//  plist[index]=new int[indexmax[index]];
//  plist[index]=(int*)malloc((sizeof(int))*indexmax[index]);
  plist[index]=(int*)realloc(plist[index],(sizeof(int))*indexmax[index]*2);
  memset(&plist[index][indexmax[index]],0,indexmax[index]*sizeof(int));
  indexmax[index]=indexmax[index]*2;

 }
  // realloc(particles,maxlocal*sizeof(OnePart));
  // memset(&particles[oldmax],0,(maxlocal-oldmax)*sizeof(OnePart));
  plist[index][l[index]]=i;
  l[index]++;
}

}
// if(nlocal==0)
// return ;
for(int dim=0;dim<=2;dim++)
  { 
   
    test[0]=0;
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
    // int message = 3;

    // int dest = (world_rank + 1) % world_size;
    for(int i=0;i<=1;i++)
    {if(test[i]==1) 
    {                
     particle_number_send[i]=l[i+dim*2];
                       
    MPI_Send(&particle_number_send[i], 1, MPI_INT,negb_rank[i], 0, MPI_COMM_WORLD);
    MPI_Recv(&particle_number_recv[i], 1, MPI_INT,negb_rank[i], 0, MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
  //  if(particle_number_recv[i]>=1)
  // printf("particle_number_recv,%d",particle_number_recv[i]);
  //  printf("particle_number_recv,%d",particle_number_recv[i]);
  }}
 
    
     for(int i=0;i<=1;i++)
{if(test[i]==1) {
  if(particle_number_send[i]>0)
  
  psend[i] = (char *) malloc(nbytes_particle*particle_number_send[i]);
  if(particle_number_recv[i]>0)
  precv[i] = (char *) malloc(nbytes_particle*particle_number_recv[i]);
}}
  for(int i = 0;i<=1;i++)
   {if(test[i]==1)
    {if(l[i+dim*2]>=1) 
  {offset=0;
  for(int j=0;j<l[i+dim*2];j++)
   { //  { psend[j]=particles[plist[i+dim*2][j]]; 
      memcpy(&(psend[i][offset]),&particles[plist[i+dim*2][j]],nbytes_particle);
      offset += nbytes_particle;
     particles[plist[i+dim*2][j]].flag=1;
        //  printf("plist[i+dim*2][j],%.2f  \t", particles[plist[i+dim*2][j]].v[0]); 
          }        
 }}
 else
 {//dealing with boundary  specular reflection
//  {if(l[i+dim*2]>=1) 
//   for(int j=0;j<l[i+dim*2];j++)
//   particles[plist[i+dim*2][j]].flag=1;
//   }
//    if(particles[plist[i+dim*2][j]].x[0]<lo[0]||particles[plist[i+dim*2][j]].x[0]>hi[0])
//    particles[plist[i+dim*2][j]].x[0]=lo[0]+randu(e)*(hi[0]-lo[0]);
//    if(particles[plist[i+dim*2][j]].x[1]<lo[1]||particles[plist[i+dim*2][j]].x[1]>hi[1])
// particles[plist[i+dim*2][j]].x[1]=lo[1]+randu(e)*(hi[1]-lo[1]);
//  if(particles[plist[i+dim*2][j]].x[1]<lo[1]||particles[plist[i+dim*2][j]].x[1]>hi[1])
// particles[plist[i+dim*2][j]].x[2]=lo[2]+randu(e)*(hi[2]-lo[2]);

//  }
//  l[i+dim*2]=0;

//  }
  }
 
 }
 for(int i=0;i<=1;i++)
 { if(test[i]==1) {
  //  printf("recv[0]");
  offset=0;
  double x=9;
     if (particle_number_send[i] > 0) 
     MPI_Send(&(psend[i][offset]), (particle_number_send[i]*nbytes_particle), MPI_CHAR, negb_rank[i],0, MPI_COMM_WORLD);

 
      if (particle_number_recv[i] > 0) 
  { MPI_Recv(&(precv[i][offset]), (particle_number_recv[i]*nbytes_particle), MPI_CHAR, negb_rank[i], 0,  MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
//  memcpy(&x,&(precv[i][offset+4*sizeof(int)]),sizeof(double));
//      printf("recv[0],%d ,%f \t",nbytes_particle, x);
     }
 } }                          
// }
 

for(int i=0;i<=1;i++)
{int start=nlocal;

if(test[i]==1)
  {if(particle_number_recv[i]>=1)
  {
    int offset=0;
    for(int j=0;j<particle_number_recv[i];j++)
{
  // particles[nlocal]=precv[j];
    //  printf("recv[0],%.2f  \t", particles[nlocal].x[0]); 
 
     double x=9;
     int reallocflag = 0;
  if (nlocal == maxlocal) {
    grow(1);
    reallocflag = 1;
  }
     
      //  memcpy(&precv[i][offset+4*sizeof(int)],&x,sizeof(double));
     memcpy(&particles[nlocal],&(precv[i][offset]),nbytes_particle);
    //  printf("recv[0],%d ,%f \t",nbytes_particle,particles[nlocal].x[0]);


      offset += nbytes_particle;
       nlocal++;
      // char*ptr=&precv[i][offset];
      //  ptr +=5*sizeof(int);
      // memcpy(&particles[nlocal].v[0],ptr,sizeof(double));
       
     
    }
     
          
      //  add_particle();
}
  for(int j=0;j<particle_number_recv[i];j++)
{index=particle_domain_index(&particles[j+start]);
// index=-1;
if(index>=0)
{
 if(l[index]==indexmax[index])
 {
  //  delete[] plist[index];
//  indexmax[index]=indexmax[index]*2;
// //  plist[index]=new int[indexmax[index]];
//  realloc(plist[index],(sizeof(int))*indexmax[index]);
//  memset(plist[index],0,indexmax[index]*sizeof(int));
//  printf("(sizeof(int))*indexmax[index]%d",(sizeof(int))*indexmax[index]);
//  int a=sizeof(int)*indexmax[index];
//  plist[index]=new int[indexmax[index]];
//  plist[index]=(int*)malloc((sizeof(int))*indexmax[index]);
  plist[index]=(int*)realloc(plist[index],(sizeof(int))*indexmax[index]*2);
//  plist[index]=(int*) realloc(plist[index],(sizeof(int))*indexmax[index]);
  memset(&plist[index][indexmax[index]],0,indexmax[index]*sizeof(int));
  indexmax[index]=indexmax[index]*2;
 }
  plist[index][l[index]]=j+start;
  l[index]++;
  // particles[j+start].flag=1;
}}


  }
  
  }
for(int i=0;i<=1;i++)
  if(test[i]!=1)
  {
    if(l[i+dim*2]>=1) 
  for(int j=0;j<l[i+dim*2];j++)
  particles[plist[i+dim*2][j]].flag=1;
  }
   for(int h=0;h<2;h++)
    {
      if(psend[h]!=NULL)
      free(psend[h]);
      psend[h]=NULL;
      if(precv[h]!=NULL)
      free(precv[h]);
      precv[h]=NULL;

    }

  }

  
int weight=1;
double nx,ny,nz; 
double x0[3],v0[3],x1[3],v1[3],x2[3],v2[3],x3[3],v3[3];
int ispecies[3];
for (int i = 0; i<nlocal ; i++) {
if(particles[i].flag==1){
    while(particles[nlocal-1].flag==1&&nlocal>1)
    nlocal--;
     memcpy(&particles[i],&particles[nlocal-1],sizeof(OnePart));
     nlocal--;
      continue;
   }
  //  continue;
//         if(particles[i].x[0]<lo[0]||particles[i].x[0]>hi[0])
//    {particles[i].x[0]=lo[0]+randu(e)*(hi[0]-lo[0]);
//    printf("%dflag\t",particles[i].flag);}
//    if(particles[i].x[1]<lo[1]||particles[i].x[1]>hi[1])
// particles[i].x[1]=lo[1]+randu(e)*(hi[1]-lo[1]);
//  if(particles[i].x[2]<lo[2]||particles[i].x[2]>hi[2])
// particles[i].x[2]=lo[2]+randu(e)*(hi[2]-lo[2]);
   //weighting
   weight=particles[i].ispecies;
  //  weight=0;
    nx=particles[i].x[0];
 ny=particles[i].x[1];
 nz=particles[i].x[2];
  // printf("%f ,%f ,%f /n",nx,ny,nz);
//  printf("%f,%f,%f/n",nx,ny,nz);
 weighting(&nx,&ny,&nz,&weight);

  //  printf("weighting over");
   

ispecies[0]=particles[i].ispecies;
ispecies[1]=-1;
ispecies[2]=-1;

if(randu(e)<species[ispecies[0]].coll_ratio)
{
  //printf("%f,x3",particles[i].x[2]);
// continue;
for(int k=0;k<3;k++)
{x0[k]=particles[i].x[k];
v0[k]=particles[i].v[k];}
  // MCC(particles[i].x,particles[i].v,x2,v2,x3,v3,ispecies);
   MCC(x0,v0,x2,v2,x3,v3,ispecies);
  if(ispecies[1]>=0)
  {
    if (nlocal == maxlocal) {
    grow(1);
    // reallocflag = 1;
  }
    particles[nlocal].x[0]=particles[i].x[0];
    particles[nlocal].x[1]=particles[i].x[1];
    particles[nlocal].x[2]=particles[i].x[2];
    particles[nlocal].v[0]=v2[0]*dt/dx;
    particles[nlocal].v[1]=v2[1]*dt/dy;
    particles[nlocal].v[2]=v2[2]*dt/dz;
    particles[nlocal].a[0]=0;
    particles[nlocal].a[1]=0;
    particles[nlocal].a[2]=0;
    particles[nlocal].ispecies=ispecies[1];
    particles[nlocal].flag=0;
    nlocal++;

  }
  if(ispecies[2]>=0)
  {
    if (nlocal == maxlocal) {
    grow(1);
    // reallocflag = 1;
  }
    particles[nlocal].x[0]=particles[i].x[0];
    particles[nlocal].x[1]=particles[i].x[1];
    particles[nlocal].x[2]=particles[i].x[2];
    particles[nlocal].v[0]=v3[0]*dt/dx;
    particles[nlocal].v[1]=v3[1]*dt/dy;
    particles[nlocal].v[2]=v3[2]*dt/dz;
    particles[nlocal].a[0]=0;
    particles[nlocal].a[1]=0;
    particles[nlocal].a[2]=0;
    particles[nlocal].ispecies=ispecies[2];
    particles[nlocal].flag=0;
    nlocal++;

  }
  }
// int_x=static_cast<int>(std::floor(nx));
// int_y=static_cast<int>(std::floor(ny));
// int_z=static_cast<int>(std::floor(nz));
// double_x=nx-int_x;
// double_y=ny-int_y;
// double_z=nz-int_z;

// barray[int_z][int_y][int_x]+=1*double_x*double_y*double_z;
// barray[int_z][int_y][int_x+1]+=1*(1-double_x)*double_y*double_z;
// barray[int_z][int_y+1][int_x]+=1*double_x*(1-double_y)*double_z;
// barray[int_z+1][int_y][int_x]+=1*double_x*double_y*(1-double_z);
// barray[int_z][int_y+1][int_x+1]+=1*(1-double_x)*(1-double_y)*double_z;
// barray[int_z+1][int_y][int_x+1]+=1*(1-double_x)*double_y*(1-double_z);
// barray[int_z+1][int_y+1][int_x]+=1*double_x*(1-double_y)*(1-double_z);
// barray[int_z+1][int_y+1][int_x+1]+=1*(1-double_x)*(1-double_y)*(1-double_z);












}

for(int i=0;i<6;i++)
{
  if(plist[i]!=NULL)
  free(plist[i]);

}

}
  
 
// int message = world_rank;
//     int dest = (world_rank + 1) % world_size;
//     MPI_Send(&message, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);


int Particle::particle_domain_index(OnePart* particle){
if (particle->x[0] < lo[0])  return 0 ; 
else if (particle->x[0] > hi[0])  return 1 ;
else if (particle->x[1] < lo[1])  return 2 ; 
else if (particle->x[1] > hi[1])  return 3 ;
else if (particle->x[2] < lo[2])  return 4 ; 
else if (particle->x[2] > hi[2])  return 5 ;
return -1;
}
