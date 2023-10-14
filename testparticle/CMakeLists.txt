# CTest
cmake_minimum_required(VERSION 2.8)
project(particlecomm LANGUAGES CXX Fortran)

    enable_testing()
    set(PATH_PROJECT ${PROJECT_SOURCE_DIR} CACHE STRING "project path")
    # 06 test for domain
    set(testname "06_test_domain")
    find_package(MPI REQUIRED)
    add_executable(${PROJECT_NAME} test_particle_com3d.F90 ParticleBundle.f90 ParticleOne.f90 Constants.F90 PICParticleCom3d.f90)
#add_definitions(${MPI_Fortran_COMPILE_FLAGS})
#include_directories(${MPI_Fortran_INCLUDE_PATH})
#link_directories(${MPI_Fortran_LIBRARIES})

#target_link_libraries(${NAME_FLIB} ${NAME_SOLVER})

target_link_libraries(${NAME_FLIB}MPI::MPI_Fortran)

#target_link_libraries(${PROJECT_NAME} otherlibs ${MPI_Fortran_LIBRARIES})
   # link flib to master
#target_link_libraries(${PROJECT_NAME}  ${NAME_FLIB})
   

    add_test(NAME ${PROJECT_NAME} COMMAND "mpirun" "-np" "8" $<TARGET_FILE:${PROJECT_NAME}>)

    

#endif()
