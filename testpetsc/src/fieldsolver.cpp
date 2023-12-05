/*-----------------------------------------------------------------------------
3维泊松方程并行求解

方程:
                        - Laplacian u = 0      0 < x < 5, 0 < y < 5，0 < z < 5

边界：设置z方向上边界电势为1，下边界电势为0，然后中间区域插值获得电势

离散：

    - u(i-1, j, k) + 2u(i, j, k) - u(i+1, j, k )
    - u(i, j-1, k)  + 2u(i, j, k) - u(i, j+1, k) 
    - u(i, j, k-1)  + 2u(i, j, k) - u(i, j+1, k+1)== 0
-----------------------------------------------------------------------------*/

#include <iostream>

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscdmda.h>
#include<fstream>
#include"material.h"
#include"fieldsolver.h"
using namespace std;

extern "C" {
  void  GetRho(int coord_x,int coord_y, int coord_z, int width_x, int width_y,int width_z, float *a);
  void  SendPhi(int coord_x,int coord_y, int coord_z, int width_x, int width_y,int width_z, float*a);
  void  Finalize();
}


int FieldSolver::initpetsc( PetscInt Mx,  PetscInt My, PetscInt Mz, int data [5][5][5]){
    
    
    // PetscErrorCode ierr = 0;
    // Mat A;
    // Vec b, x,d;
    // KSP ksp;
    // DM dm;
    // PetscMPIInt size, rank;
    // PetscScalar v[7];
    // MatStencil row, col[7];
    //  double ***barray = nullptr,***array=nullptr,***geometry=nullptr;
    //  PetscScalar
    
    double*charge;
    float *rho,*phi,*rho1;

    Material *material[3];
    Material vacuum;
    Metal metal;
    Dielectric dielectric;

    //赋值
    material[0]=&vacuum;
    material[1]=&metal;
    material[2]=&dielectric;
    for(int i=0;i<3;i++)
    {
        // cout<<"epsilon"<<material[i]->epsilon<<endl;
        //  cout<<"epsilon"<<vacuum.epsilon<<endl;
    }

   
   
    // // init petsc
    // PetscInitialize(NULL, NULL, (char *)0, NULL);
  

    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // shape

    DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_STAR, Mx, My, Mz, 2, 2, 2,
                 1, 1, NULL, NULL, NULL, &dm);
    DMSetFromOptions(dm);
    DMSetUp(dm);

    // create mat
    DMCreateMatrix(dm, &A);
    MatSetFromOptions(A);

    //此部分设置数据存储格式，可以提升速度
     MatMPIAIJSetPreallocation(A, 7, NULL, 7, NULL);
     MatSeqAIJSetPreallocation(A, 7, NULL);
     MatSeqSBAIJSetPreallocation(A, 1, 7, NULL);
     MatMPISBAIJSetPreallocation(A, 1, 7, NULL, 7, NULL);
     MatMPISELLSetPreallocation(A, 7, NULL, 7, NULL);
     MatSeqSELLSetPreallocation(A, 7, NULL);

    // create vec
    DMCreateGlobalVector(dm, &b);
    DMCreateGlobalVector(dm, &x);

    // 此部分获取该进程x,y,z起始坐标以及各方向的宽度
   
    DMDAGetCorners(dm, &coord_x, &coord_y, &coord_z,
                   &width_x, &width_y, &width_z);
    rho=new float[width_x*width_y*width_z];
    // for(i=0,i<)
    //  GetRho(coord_x, coord_y, coord_z, width_x,width_y,width_z,rho);
    // lo=new double[3];
    // hi=new double[3];
    // lo[0]=coord_x;
    // hi[0]=coord_x+width_x;
    
   std:: cout << rank << ": " << coord_x << " " << coord_y << " "
         << coord_z << " " << width_x << " " << width_y << " " << width_z << std::endl;

    for (auto i =  coord_x; i < coord_x + width_x; i++) {
        for (auto j = coord_y; j < coord_y +width_y; j++) {
            for (auto k = coord_z; k < coord_z + width_z; k++) {
                row.i = i; row.j = j; row.k = k;
                if (0 == i || Mx-1 == i || 0 == j || My-1 == j||0==k||Mz-1==k) {
                    v[0] = 1.0;
                     //设置边界条件为0，此部分v[0]=1是保证只有1个未知数，即v[0]*x(i,j,k)=0得x(i,j,k)=0
                    MatSetValuesStencil(A, 1, &row, 1, &row, v, INSERT_VALUES);
                }
                else {
                    col[0].i = i-1; col[0].j = j;   col[0].k = k;
                    col[1].i = i+1; col[1].j = j;   col[1].k = k;
                    col[2].i = i;   col[2].j = j-1; col[2].k = k;
                    col[3].i = i;   col[3].j = j+1; col[3].k = k;
                    col[4].i = i;   col[4].j = j;   col[4].k = k-1;
                    col[5].i = i;   col[5].j = j;   col[5].k = k+1;
                    col[6].i = i;   col[6].j = j;   col[6].k = k;

                    v[0] = -1*material[data[i-1][j][k]]->epsilon;
                    v[1] = -1*material[data[i+1][j][k]]->epsilon;
                    v[2] = -1*material[data[i][j-1][k]]->epsilon;
                    v[3] = -1*material[data[i][j+1][k]]->epsilon;
                    v[4] = -1*material[data[i][j][k-1]]->epsilon;
                    v[5] = -1*material[data[i][j][k+1]]->epsilon;
                    v[6] =  -(v[0]+v[1]+v[2]+v[3]+v[4]+v[5]);
                    // v[0] = -1;
                    // v[1] = -1;
                    // v[2] = -1;
                    // v[3] = -1;
                    // v[4] = -1;
                    // v[5] = -1;
                    // v[6] = 6;
                     cout<<" "<<data[i][j][k]<<' '<<v[6];
                    MatSetValuesStencil(A, 1, &row, 7, col, v, INSERT_VALUES);
                    //将泊松方程置入，其中col的i,j,k为变量坐标，row设置第几个方程
                }
            }
        }
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    GetRho(coord_x, coord_y, coord_z, width_x,width_y,width_z,rho);
    DMDAVecGetArray(dm, b, &barray);
     for (auto i =  coord_x; i < coord_x + width_x; i++) {
        for (auto j = coord_y; j < coord_y +width_y; j++) {
            for (auto k = coord_z; k < coord_z + width_z; k++) {
                if (0 == k)
                {
                    barray[k][j][i] = 0;
                }
                else if (Mz - 1 == k)
                {
                    barray[k][j][i] = 1;
                }
                else if (0 == i || Mx - 1 == i || 0 == j || My - 1 == j)
                {
                    barray[k][j][i] = k / (Mz - 1.0);
                }
                else {
                      barray[k][j][i] = rho[(i-coord_x)*width_y*width_z+(j-coord_y)*width_z+k-coord_z];
                    //   barray[k][j][i]=0;
                }
            }
        }
    }

    DMDAVecRestoreArray(dm, b, &barray);

    // create ksp
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);

    // solve once
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
  
 }
 int FieldSolver::fieldsolve(){
   double*charge;
    float *rho,*phi,*rho1;
  ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);
 
  
    DMDAVecGetArray(dm,x,&array);

   phi=new float[width_x*width_y*width_z];
     std::ofstream log("./solve/petscphi"+std::to_string(rank)+".txt", std::ios_base::app | std::ios_base::out);
     log<< width_x << " " << width_y << " " << width_z << " "<<std::endl;
    //  log<<"电势 存储顺序：先遍历z,再遍历x和y"<<std::endl;
     for (auto i =  coord_x; i < coord_x + width_x; i++) {
        for (auto j = coord_y; j < coord_y +width_y; j++) {
            for (auto k = coord_z; k < coord_z + width_z; k++) {
                phi[(i-coord_x)*width_y*width_z+(j-coord_y)*width_z+k-coord_z]=array[k][j][i];
            log<<array[k][j][i]<<" ";}
        log<<std::endl;}}
  
    DMDAVecRestoreArray(dm,x,&array);
     
    // SendPhi(coord_x, coord_y, coord_z,width_x,width_y,width_z,phi);
    // Finalize();
//     PetscFinalize(); 
 }


