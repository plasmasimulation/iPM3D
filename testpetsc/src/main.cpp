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
//  #include"create_particles.h"
// #include"sparta.h"
//  #include"pointers.h"
 
/*本测试测试petsc三维场求解，整体网格为5x5x5(包含边界),电荷密度设置为0，边界条件
设置z方向上边界电势为1，下边界电势为0，然后中间区域插值获得电势，petsc求解电势后
传电势回去求解出电场并输出，此部分结果包含电场数据以及电势数据，
保存在run/solve文件夹内*/
using namespace std;
 extern "C" {
  void MCCBundleInit();
}
  

  
int main(int argc, char** argv) {
 
   PetscInt Mx, My, Mz;
 int rank;
 int id=1; int ispecies=1;int icell=1;
   double *x; double *v; double erot=1; double evib=1;
  int data2[5][5][5];
  int *xyz_np=new int[3]{2,2,2};
  int dx,dy,dz,dt;
  // xyz_np={2,2,2};
  int data[64][64][64];
   Particle* particle=new Particle();
   CreateParticles* createparticles=new CreateParticles(); 
   FieldSolver* fieldsolver=new FieldSolver();
     MPI_Init(nullptr, nullptr);
     PetscInitialize(NULL, NULL, (char *)0, NULL);
     MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   
  
   x=new double[3];
   v=new double[3];
   x[2]=2;
   x[0]=1;
    x[1]=2;
    v[1]=1;
    v[2]=2;
    v[0]=1;
    // CreateParticles createparticles(sparta);
    // 初始化MPI环境
    // parameter set
    Mx = 5;
    My = 5;
    Mz = 5; //包含边界
    dx=1;
    dy=1;
    dz=1;
    dt=1;
 
  Domain* domain=new Domain(xyz_np,Mx,My,Mz,dx,dy,dz,rank);
 MCCBundleInit();
 load_material(data2);
   fieldsolver->initpetsc(Mx, My, Mz, data2,xyz_np);
    particle->init(domain->lo,domain->hi,domain->dx,domain->dy,domain->dz,dt);
    particle->grow(50);
    particle->add_particle(id,ispecies,icell,x,v,erot,evib);
     createparticles->create_local(particle,domain->lo,domain->hi);
     fieldsolver->fieldsolve();
      particle->particle_move_comm(fieldsolver->barray);
    //   cout<<"first"<<endl;
    //   particle->particle_move_comm();
    //      cout<<"second"<<endl;
    // particle->particle_move_comm();
   
  // cout<<"ok";
 


  // 创建HDF5文件  
if(rank==0)
   {// 创建HDF5文件  
    hid_t file_id = H5Fcreate("material.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);  
  
    // 创建数据集  
    hsize_t dims[3] = {64,64,64}; // 定义数据集的大小为10  
    hid_t dataspace_id = H5Screate_simple(3, dims, NULL);  
    hid_t dataset_id = H5Dcreate2(file_id, "my_dataset", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
  
    for(int i=0;i<64;i++)
    for(int j=0;j<64;j++)
    for(int k=0;k<64;k++)
    {if(i==0||i==63)
      data[i][j][k]=1;
      else
      data[i][j][k]=0;
    }
    // 写入数据  
    // int data[5][5][5] =     {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    //                       }, {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2
    //                       }, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    //                       }, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    //                       }, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};  
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_INT, dataspace_id, dataspace_id, H5P_DEFAULT, data);  
  
    // 关闭数据集、数据空间和文件  
    status = H5Dclose(dataset_id);  
    status = H5Sclose(dataspace_id);  
    status = H5Fclose(file_id);  
  
    if (status < 0) {  
        cerr << "Failed to write data or close file/dataset/dataspace" << endl;  
        return -1;  
    }  
  
    cout << "Data written to " <<"example.h5" << " successfully!" << endl;  
  
   }
    
    
    
       fieldsolver->fieldsolve();
      PetscFinalize(); 

     MPI_Finalize();
     delete particle;
     delete createparticles;
     delete fieldsolver;
   return 0;
}
