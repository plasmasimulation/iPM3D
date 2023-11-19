#include <iostream>
 #include <mpi.h>
#include"fieldsolver.h"
/*本测试测试petsc三维场求解，整体网格为5x5x5(包含边界),电荷密度设置为0，边界条件
设置z方向上边界电势为1，下边界电势为0，然后中间区域插值获得电势，petsc求解电势后
传电势回去求解出电场并输出，此部分结果包含电场数据以及电势数据，
保存在run/solve文件夹内*/



int main(int argc, char** argv) {

    // 初始化MPI环境
  MPI_Init(nullptr, nullptr);
  testpetsc();
  MPI_Finalize();
  
  return 0;
}
