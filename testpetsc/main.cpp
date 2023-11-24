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
/*本测试测试petsc三维场求解，整体网格为5x5x5(包含边界),电荷密度设置为0，边界条件
设置z方向上边界电势为1，下边界电势为0，然后中间区域插值获得电势，petsc求解电势后
传电势回去求解出电场并输出，此部分结果包含电场数据以及电势数据，
保存在run/solve文件夹内*/
using namespace std;
 
  

  
int main(int argc, char** argv) {
 
  PetscInt Mx, My, Mz;
 int rank;
    // 初始化MPI环境
    // parameter set
    Mx = 5;
    My = 5;
    Mz = 5; //包含边界
  MPI_Init(nullptr, nullptr);
  PetscInitialize(NULL, NULL, (char *)0, NULL);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  testpetsc(Mx, My, Mz);

  PetscFinalize(); 
  MPI_Finalize();
  // 创建HDF5文件  
if(rank==0)
   {// 创建HDF5文件  
    hid_t file_id = H5Fcreate("example.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);  
  
    // 创建数据集  
    hsize_t dims[1] = {10}; // 定义数据集的大小为10  
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);  
    hid_t dataset_id = H5Dcreate2(file_id, "my_dataset", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
  
    // 写入数据  
    float data[10] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f};  
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, dataspace_id, dataspace_id, H5P_DEFAULT, data);  
  
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
  return 0;
}
