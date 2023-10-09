cmake_minimum_required(VERSION 2.8)
project(iPM3D)#指定工程名
....
add_subdirectory(test)#递归进入test目录
enable_testing()#见附录1 可以开启项目的测试功能。一定要在根目录下的CMakeLists.txt中开启
#add_test(NAME mytest COMMAND Test)#见附录1