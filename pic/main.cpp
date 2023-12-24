#include <iostream>
 #include <mpi.h>

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscdmda.h>
#include<fstream>
#include"material.h"
#include"fieldsolver.h"
#include "hdf5.h" 
#include"particle.h"
#include"create_particles.h"
#include"domain.h"
#include <cstring> 
#include <ctime>
//  #include"create_particles.h"
// #include"sparta.h"
//  #include"pointers.h"
 
/*
MCC 相关初始化信息在MCCinterface.F90中
*/
using namespace std;
 extern "C" {
  void MCCBundleInit(double *coll_ratio,double *dx,double *dt);
}
  

  
int main(int argc, char** argv) {
 
  PetscInt Mx, My, Mz;
  int rank;
  
  int*** data2;
  int *xyz_np=new int[3]{2,2,2};
  double dx[3],dt;
  double *coll_ratio=new double[2];
  // xyz_np={2,2,2};

   Particle* particle=new Particle();
   CreateParticles* createparticles=new CreateParticles(); 
   FieldSolver* fieldsolver=new FieldSolver();
     MPI_Init(nullptr, nullptr);
     PetscInitialize(NULL, NULL, (char *)0, NULL);
     MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   
  
    // CreateParticles createparticles(sparta);
    // 初始化MPI环境
    // parameter set
    Mx = 65;
    My = 65;
    Mz = 65; //包含边界
    dx[0]=3.125e-4; //1mm
    dx[1]=3.125e-4;
    dx[2]=3.125e-4;
    dt=1e-10;
    int num =2;
    char**name;
    name=new char*[2];
    name[0]=new char[16];
    name[1]=new char[16];
    strcpy(name[0],"electron");
    strcpy(name[1],"Ar+");
    int ncreate=1638400; //每个进程初始电子数

 
  Domain* domain=new Domain(xyz_np,Mx,My,Mz,dx,rank);
 MCCBundleInit(coll_ratio,dx,&dt);

   fieldsolver->initpetsc(Mx, My, Mz, xyz_np,dx,domain->lo,domain->hi);
    particle->init(domain->lo,domain->hi,domain->dx,dt);
    particle->add_species(num,name,coll_ratio);
    // particle->grow(50);
    // particle->add_particle(id,ispecies,icell,x,v,erot,evib);
     createparticles->create_local(particle,domain->lo,domain->hi,ncreate);
int step=1;//循环5个周期。
int peroid=1/dt/13.56e6;

     for(int i=0;i<step;i++)
     { for(int j=0;j<40;j++)
      {clock_t start = clock();
      fieldsolver->fieldsolve(j*dt);
     
        clock_t end1   = clock();
        // cout<<"fieldsolve time cost "<<(double)(end1-start) / CLOCKS_PER_SEC<<endl;
         particle->particle_move_comm();
       clock_t end2  = clock();
       if(rank==6){
 cout<<"fieldsolve time cost "<<(double)(end1-start) / CLOCKS_PER_SEC<<endl;
 cout<<"particlemove time cost"<<(double)(end2-end1) / CLOCKS_PER_SEC<<endl;
 cout<<j<<"step"<<endl;
 cout<<"particles"<<particle->nlocal<<endl;
       }}
       
      

       
      }
    //   cout<<"first"<<endl;
    //   particle->particle_move_comm();
    //      cout<<"second"<<endl;
    // particle->particle_move_comm();
    //  fieldsolver->fieldsolve();
   
  // cout<<"ok";
  // create_material_file();//创建介质文件，有了这个介质文件才可以修改网格数大小


 
      //  fieldsolver->fieldsolve();
      PetscFinalize(); 

     MPI_Finalize();
     delete particle;
     delete createparticles;
     delete fieldsolver;
   return 0;
}
