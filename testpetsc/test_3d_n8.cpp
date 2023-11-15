/*-----------------------------------------------------------------------------
3维泊松方程并行求解

方程:
                        - Laplacian u = 1      0 < x < 32, 0 < y < 32，0 < z < 32

边界：
                            边界为 0

离散：

    - u(i-1, j, k) + 2u(i, j, k) - u(i+1, j, k )
    - u(i, j-1, k)  + 2u(i, j, k) - u(i, j+1, k) 
    - u(i, j, k-1)  + 2u(i, j, k) - u(i, j+1, k+1)== 10
-----------------------------------------------------------------------------*/

#include <iostream>

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscdmda.h>
#include <time.h>
#include<fstream>

extern "C" {
  void print_1d(int n, int*a);
}

int testfieldcomm()
{
    // test for 2d deps
    PetscErrorCode ierr = 0;
    Mat A;
    Vec b, x,d;
    KSP ksp;
    DM dm;
    PetscMPIInt size, rank;
    PetscScalar v[7];
    MatStencil row, col[7];
    PetscScalar ***barray = nullptr,***array=nullptr,***x_array=nullptr;

    PetscInt Mx, My, Mz;
    long start,end;//定义clock_t变量
    start = clock();
 

   int*a;
 a=new int[3];
  print_1d(3, a);
  for(int i=0;i<3;i++)
   std::cout<<"C++read"<<a[i]<<std::endl;

    // init petsc
    PetscInitialize(NULL, NULL, (char *)0, NULL);

    // parameter set
    Mx = 4;
    My = 4;
    Mz = 4;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // shape

    DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_STAR, Mx, My, Mz, 2, 2, 2,
                 1, 2, NULL, NULL, NULL, &dm);
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
    PetscInt coord_x, coord_y, coord_z, width_x, width_y, width_z;
    DMDAGetCorners(dm, &coord_x, &coord_y, &coord_z,
                   &width_x, &width_y, &width_z);

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
                    col[1].i = i;   col[1].j = j;   col[1].k = k;
                    col[2].i = i+1; col[2].j = j;   col[2].k = k;
                    col[3].i = i;   col[3].j = j-1; col[3].k = k;
                    col[4].i = i;   col[4].j = j+1; col[4].k = k;
                    col[5].i = i;   col[5].j = j;   col[5].k = k-1;
                    col[6].i = i;   col[6].j = j;   col[6].k = k+1;

                    v[0] = -1;
                    v[1] =  6;
                    v[2] = -1;
                    v[3] = -1;
                    v[4] = -1;
                    v[5] = -1;
                    v[6] = -1;
                    MatSetValuesStencil(A, 1, &row, 7, col, v, INSERT_VALUES);
                    //将泊松方程置入，其中col的i,j,k为变量坐标，row设置第几个方程
                }
            }
        }
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    // set b
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
                    barray[k][j][i] = -10;
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
    ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);
    // KSPGetSolution(ksp,&d);
    
// dk=array.begin();
//    std::cout<<dk<<std::endl;
DMDAVecGetArray(dm,x,&array);

 end = clock();   //结束时间
std::cout<<"time = "<<double(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;  //输出时间（单位：ｓ）

    // std::cout<<dk<<std::endl;
     std::ofstream log("logfile"+std::to_string(rank)+".txt", std::ios_base::app | std::ios_base::out);
     log<< width_x << " " << width_y << " " << width_z << std::endl;
     for (auto i =  coord_x; i < coord_x + width_x; i++) {
        for (auto j = coord_y; j < coord_y +width_y; j++) {
            for (auto k = coord_z; k < coord_z + width_z; k++) {
    //  log<<array[i][orrd_y+width_y-1][orrd_x+width_x-1]<<std::endl;
 log<<array[k][j][i]<<" ";
}
log<<std::endl;}}

 DMDAVecRestoreArray(dm,x,&array);

   DMDAGetGhostCorners(dm, &coord_x, &coord_y, &coord_z,
                     &width_x, &width_y, &width_z);
 DMGetLocalVector(dm, &d);
    DMGlobalToLocalBegin( dm, x, INSERT_VALUES, d);
   
      DMGlobalToLocalEnd(dm, x, INSERT_VALUES, d);
      DMDAVecGetArray(dm,d, &x_array);

   std:: cout <<"afterlocal"<< rank << ": " << coord_x << " " << coord_y << " "
          << coord_z << " " << width_x << " " << width_y << " " << width_z << std::endl;
//  DMDAVecGetArrayDOF(dm,d,&x);
    // output
     // std::cout<<dk<<std::endl;
     std::ofstream log2("logyyfile"+std::to_string(rank)+".txt", std::ios_base::app | std::ios_base::out);
     log<< width_x << " " << width_y << " " << width_z << std::endl;
     for (auto i =  coord_x; i < coord_x + width_x; i++) {
        for (auto j = coord_y; j < coord_y +width_y; j++) {
            for (auto k = coord_z; k < coord_z + width_z; k++) {
    //  log<<array[i][orrd_y+width_y-1][orrd_x+width_x-1]<<std::endl;
 log2<<x_array[k][j][i]<<" ";
}
log2<<std::endl;}}
     PetscViewer myViewer;
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, "phi.dat", &myViewer);
     VecView(x, myViewer);
    // VecView(d, myViewer);

    // PetscViewerPushFormat(myViewer,PETSC_VIEWER_STDOUT_COMMON)
    
    // if (!(rank))
    // {
    //     ostringstream os;
    //     os << "../../" << "plotVec.py" << " " << "phi.dat" << " " << Mx << " " << My;
    //     string cmd = os.str();
    //     system(cmd.c_str());
    // }
   

    PetscFinalize();
   
}

int testpetsc(){
    
    
    PetscErrorCode ierr = 0;
    Mat A;
    Vec b, x,d;
    KSP ksp;
    DM dm;
    PetscMPIInt size, rank;
    PetscScalar v[7];
    MatStencil row, col[7];
    PetscScalar ***barray = nullptr,***array=nullptr,***x_array=nullptr;
    PetscInt Mx, My, Mz;
    double*charge;

    long start,end;//定义clock_t变量
    start = clock();
    // init petsc
    PetscInitialize(NULL, NULL, (char *)0, NULL);
    // parameter set
    Mx = 80;
    My = 80;
    Mz = 80;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // shape

    DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_STAR, Mx, My, Mz, 2, 2, 2,
                 1, 2, NULL, NULL, NULL, &dm);
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
    PetscInt coord_x, coord_y, coord_z, width_x, width_y, width_z;
    DMDAGetCorners(dm, &coord_x, &coord_y, &coord_z,
                   &width_x, &width_y, &width_z);

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
                    col[1].i = i;   col[1].j = j;   col[1].k = k;
                    col[2].i = i+1; col[2].j = j;   col[2].k = k;
                    col[3].i = i;   col[3].j = j-1; col[3].k = k;
                    col[4].i = i;   col[4].j = j+1; col[4].k = k;
                    col[5].i = i;   col[5].j = j;   col[5].k = k-1;
                    col[6].i = i;   col[6].j = j;   col[6].k = k+1;

                    v[0] = -1;
                    v[1] =  6;
                    v[2] = -1;
                    v[3] = -1;
                    v[4] = -1;
                    v[5] = -1;
                    v[6] = -1;
                    MatSetValuesStencil(A, 1, &row, 7, col, v, INSERT_VALUES);
                    //将泊松方程置入，其中col的i,j,k为变量坐标，row设置第几个方程
                }
            }
        }
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);


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
                    barray[k][j][i] = -10;
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
    ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);
  

    end = clock();   //结束时间
    std::cout<<"time = "<<double(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;  //输出时间（单位：ｓ）
    DMDAGetGhostCorners(dm, &coord_x, &coord_y, &coord_z, &width_x, &width_y, &width_z);
    DMGetLocalVector(dm, &d);
    DMGlobalToLocalBegin( dm, x, INSERT_VALUES, d);
    DMGlobalToLocalEnd(dm, x, INSERT_VALUES, d);
    DMDAVecGetArray(dm,d, &x_array);

    std:: cout <<"afterlocal"<< rank << ": " << coord_x << " " << coord_y << " "
          << coord_z << " " << width_x << " " << width_y << " " << width_z << std::endl;

     std::ofstream logg("loggfile"+std::to_string(rank)+".txt", std::ios_base::app | std::ios_base::out);
     logg<< width_x << " " << width_y << " " << width_z << std::endl;
     for (auto i =  coord_x; i < coord_x + width_x; i++) {
        for (auto j = coord_y; j < coord_y +width_y; j++) {
            for (auto k = coord_z; k < coord_z + width_z; k++) {

                    logg<<x_array[k][j][i]<<" ";}
                    logg<<std::endl;}}
    PetscFinalize(); 
}
