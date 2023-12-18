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
//  #include"create_particles.h"
// #include"sparta.h"
//  #include"pointers.h"
 
/*本测试测试petsc三维场求解，整体网格为5x5x5(包含边界),电荷密度设置为0，边界条件
设置z方向上边界电势为1，下边界电势为0，然后中间区域插值获得电势，petsc求解电势后
传电势回去求解出电场并输出，此部分结果包含电场数据以及电势数据，
保存在run/solve文件夹内*/
using namespace std;
 extern "C" {
  void MCCBundleInit(double *coll_ratio,double *dx,double *dt);
}
  

  
int main(int argc, char** argv) {
 
  PetscInt Mx, My, Mz;
  int rank;
  int id=1; int ispecies=1;int icell=1;
  double *x; double *v; double erot=1; double evib=1;
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
   
  
  //  x=new double[3];
  //  v=new double[3];
  //  x[2]=2;
  //  x[0]=1;
  //   x[1]=2;
  //   v[1]=1;
  //   v[2]=2;
  //   v[0]=1;
    // CreateParticles createparticles(sparta);
    // 初始化MPI环境
    // parameter set
    Mx = 95;
    My = 95;
    Mz = 95; //包含边界
    dx[0]=1e-3; //1mm
    dx[1]=1e-3;
    dx[2]=1e-3;
    dt=1e-10;
    int num =2;
    char**name;
    name=new char*[2];
    name[0]=new char[16];
    name[1]=new char[16];
    strcpy(name[0],"electron");
    strcpy(name[1],"Ar+");

 
  Domain* domain=new Domain(xyz_np,Mx,My,Mz,dx,rank);
 MCCBundleInit(coll_ratio,dx,&dt);

   fieldsolver->initpetsc(Mx, My, Mz, xyz_np,dx,domain->lo,domain->hi);
    particle->init(domain->lo,domain->hi,domain->dx,dt);
    particle->add_species(num,name,coll_ratio);
    particle->grow(50);
    // particle->add_particle(id,ispecies,icell,x,v,erot,evib);
     createparticles->create_local(particle,domain->lo,domain->hi);
int step=1 ;//循环100次。
     for(int i=0;i<step;i++)
     {fieldsolver->fieldsolve();
       particle->particle_move_comm();
      //  cout<<i<<endl;
      }
    //   cout<<"first"<<endl;
    //   particle->particle_move_comm();
    //      cout<<"second"<<endl;
    // particle->particle_move_comm();
    //  fieldsolver->fieldsolve();
   
  // cout<<"ok";
//  create_material_file();//创建介质文件，有了这个介质文件才可以修改网格数大小


 
      //  fieldsolver->fieldsolve();
      PetscFinalize(); 

     MPI_Finalize();
     delete particle;
     delete createparticles;
     delete fieldsolver;
   return 0;
}
