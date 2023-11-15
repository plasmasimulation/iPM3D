#include <iostream>
 #include <mpi.h>
#include"test_2d_n4.h"
/*本测试分为两部分，第一部分是测试petsc获取邻居场信息，
此部分为简化特意选取了4x4x4的小区域，此部分结果保存在run/com文件夹内，
有两组文件对比未获取和获取邻居节点信息的结果
第二部分是petsc场求解，整体网格为80x80x80,电荷密度设置为0，边界条件
设置z方向上边界电势为1，下边界电势为0，然后中间区域插值获得电势，petsc求解电势后
传电势回去求解出电场并输出，此部分结果包含电场数据以及电势切面图，
保存在run/solve文件夹内*/


int main(int argc, char** argv) {

    // 初始化MPI环境
  MPI_Init(nullptr, nullptr);
  testfieldcomm();
  testpetsc();
  MPI_Finalize();
  //   // 获取进程总数
  //   int world_size;
  //   MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  //   // 获取当前进程的排名
  //   int world_rank;
  //   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  //   // 获取当前进程的名称
  //   char processor_name[MPI_MAX_PROCESSOR_NAME];
  //   int name_len;
  //   MPI_Get_processor_name(processor_name, &name_len);

  //   // 打印当前进程的信息
  //   std::cout << "Hello world from processor " << processor_name << ", rank " << world_rank << " out of " << world_size << " processors\n";

  //   // 发送消息给其他进程
  //   int message = world_rank;
  //   int dest = (world_rank + 1) % world_size;
  //   MPI_Send(&message, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);

  //   // 接收消息来自其他进程
  //   int message2;
  //   if(world_rank>0)
  //  { MPI_Recv(&message2, 1, MPI_INT, world_rank-1 , 0, MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
  //   std::cout << "Processor " << processor_name << " received message " << message2 << "\n";}

  //   // 清理MPI环境
   
   
  //   std::cout << "Processor " << processor_name << " received message " << message2 << "\n";



  return 0;
}
