cmake_minimum_required(VERSION 3.10.0)

set(CMAKE_C_COMPILER icc)
set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_Fortran_COMPILER ifort)
project(main0 LANGUAGES CXX Fortran)
set(CMAKE_CXX_STANDARD 17)

# 设置头文件路径
# include_directories(${PROJECT_SOURCE_DIR}/include/)
include_directories($ENV{PETSC_ENV_PATH}/include/)
include_directories($ENV{PATH_HDF5}/include/)

# 设置源文件和头文件路径
# file(GLOB_RECURSE SRCS ${PROJECT_SOURCE_DIR}/source/*.cpp)


# build lib
# add_library(solver STATIC ${SRCS})

# 链接库
find_library(PETSc_LIB petsc $ENV{PETSC_ENV_PATH}/lib)
find_library(HDF5_LIB hdf5 $ENV{PATH_HDF5}/lib)
# target_link_libraries(solver PUBLIC ${PETSc_LIB})
# target_link_libraries(solver PUBLIC MPI::MPI_CXX)

# build lib
find_package(MPI REQUIRED)
add_library(cxx-call-f-f ${CMAKE_CURRENT_SOURCE_DIR}/FieldInterface.f90 FieldEM.f90 FieldOne.f90 PICFieldCom3d.f90)
set_property(TARGET cxx-call-f-f PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(cxx-call-f-f MPI::MPI_Fortran)

# add_executable(cxx-call-f-c ${CMAKE_CURRENT_SOURCE_DIR}/cxx-call-f.c)
# target_link_libraries(cxx-call-f-c cxx-call-f-f)
# set_property(TARGET cxx-call-f-c PROPERTY LINKER_LANGUAGE C)
include_directories ("${MPI_CXX_INCLUDE_DIRS}")

add_executable(cxx-call-f-cxx ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp fieldsolver.cpp material.cpp )
target_link_libraries(cxx-call-f-cxx cxx-call-f-f)
target_link_libraries(cxx-call-f-cxx ${MPI_LIBRARIES})
target_link_libraries(cxx-call-f-cxx ${PETSc_LIB})
target_link_libraries(cxx-call-f-cxx ${HDF5_LIB})
# target_link_libraries(cxx-call-f-cxx ${PETSC_LIBRARIES})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/run)