# CTest
cmake_minimum_required(VERSION 2.8)
project(executablefile LANGUAGES CXX Fortran)

    enable_testing()
    set(PATH_PROJECT ${PROJECT_SOURCE_DIR} CACHE STRING "project path")
    # 06 test for domain
    set(testname "06_test_domain")
    find_package(MPI REQUIRED)
    add_executable(${PROJECT_NAME} testdomain3d.f90 Domain3d.f90 FileName.F90 ParallelDump.F90 M_strings_oop.f90 M_strings.f90)
#add_definitions(${MPI_Fortran_COMPILE_FLAGS})
#include_directories(${MPI_Fortran_INCLUDE_PATH})
#link_directories(${MPI_Fortran_LIBRARIES})
set(PATH_HDF5 $ENV{HDF5_ENV_PATH} CACHE STRING "hdf5 lib path")
# include
target_include_directories(${PROJECT_NAME} PUBLIC ${PATH_HDF5}/include)
# link libraries
find_library(HDF5_LIB_HDF5 hdf5 ${PATH_HDF5}/lib)
target_link_libraries(${PROJECT_NAME} ${HDF5_LIB_HDF5})

find_library(HDF5_LIB_HDF5_FORTRAN hdf5_fortran ${PATH_HDF5}/lib)
target_link_libraries(${PROJECT_NAME} ${HDF5_LIB_HDF5_FORTRAN})

find_library(HDF5_LIB_HDF5_HL hdf5_hl ${PATH_HDF5}/lib)
target_link_libraries(${PROJECT_NAME} ${HDF5_LIB_HDF5_HL})

find_library(HDF5_LIB_HDF5HL_FORTRAN hdf5hl_fortran ${PATH_HDF5}/lib)
target_link_libraries(${PROJECT_NAME} ${HDF5_LIB_HDF5HL_FORTRAN})

#target_link_libraries(${NAME_FLIB} ${NAME_SOLVER})

target_link_libraries(${NAME_FLIB}MPI::MPI_Fortran)

#target_link_libraries(${PROJECT_NAME} otherlibs ${MPI_Fortran_LIBRARIES})
   # link flib to master
#target_link_libraries(${PROJECT_NAME}  ${NAME_FLIB})
   

    add_test(NAME ${PROJECT_NAME} COMMAND "mpirun" "-np" "8" $<TARGET_FILE:${PROJECT_NAME}>)

    

#endif()
