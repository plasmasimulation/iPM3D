#include"material.h"
#include "hdf5.h" 
#include<iostream>
// #include "H5Cpp.h" 

#define FILE_NAME "../input/material.h5"  
#define DATASET_NAME "my_dataset"  
#define NUM_THREADS 8  
using namespace std;

Material::Material(){
    epsilon=1;
}

Metal::Metal(){
    epsilon=1;
    charge=0;
    charge_onestep=0;
    capacitance=0;
    voltage=0;
}


Dielectric::Dielectric(){
    epsilon=2; //假定介质相对介电常数为2
    permittivity=0;
    permeability=0;
    conductivity=0;
}

int create_material_file(int rank){
     int data[64][64][64];
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
    
    
    
}
 void load_material( int data[5][5][5]){

      int rank, size;  
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
    MPI_Comm_size(MPI_COMM_WORLD, &size);  
  
    hid_t file_id, dataset_id, memspace_id;  
    herr_t status;  
    hsize_t start[3], count[3];  
    //  float data[5][5][5];  
    hsize_t dims[3] = {5,5,5};
     
  
    /* Open HDF5 file */  
    file_id = H5Fopen(FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);  
  
    /* Open dataset */  
    dataset_id = H5Dopen2(file_id, DATASET_NAME, H5P_DEFAULT);  
  
    /* Create memory dataspace */  
    // memspace_id = H5Screate_simple(3,dims, NULL);  
   memspace_id = H5Dget_space(dataset_id);  
    /* Set the start and count arrays for each process */  
    start[0] = 0;  
    start[1] = 0;  
    start[2] = 0;
    count[0] = 5;  
    count[1] = 5; 
    count[2] = 5; 
  
    /* Select hyperslab in the memory dataspace */  
    status = H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, start, NULL, count, NULL);  
  
    /* Read data in parallel */  
     status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);  
   cout<<" "<<data[1][1][1]<<" hdf5 ";
    /* Process the data as needed */  
    // ...  
  
    /* Close dataset and file */  
    status = H5Dclose(dataset_id);  
    status = H5Fclose(file_id);  
  
    
 }

