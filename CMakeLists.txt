cmake_minimum_required(VERSION 3.10.0)

set(CMAKE_C_COMPILER icc)
set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_Fortran_COMPILER ifort)
project(main0 LANGUAGES CXX Fortran)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS "-O3") 
set(CMAKE_CXX_FLAGS "-O3") 
set(CMAKE_Fortran_FLAGS "-O3")

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
file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/mcc/*)  
# add_library(cxx-call-f-f ${CMAKE_CURRENT_SOURCE_DIR}/mcc/FieldInterface.f90 ${CMAKE_CURRENT_SOURCE_DIR}/mcc/FieldEM.f90 ${CMAKE_CURRENT_SOURCE_DIR}/mcc/FieldOne.f90 ${CMAKE_CURRENT_SOURCE_DIR}/mcc/PICFieldCom3d.f90
# ${CMAKE_CURRENT_SOURCE_DIR}/mcc/MCCinterface.f90 ${CMAKE_CURRENT_SOURCE_DIR}/mcc/SpecyOne.F90 ${CMAKE_CURRENT_SOURCE_DIR}/mcc/Constants.F90 ${CMAKE_CURRENT_SOURCE_DIR}/mcc/ReactionOnePegasus.f90 ${CMAKE_CURRENT_SOURCE_DIR}/mcc/Numrical.F90  ${CMAKE_CURRENT_SOURCE_DIR}/mcc/MCCSigma.f90 
# ${CMAKE_CURRENT_SOURCE_DIR}/mcc/TypeMCC.f90 ${CMAKE_CURRENT_SOURCE_DIR}/mcc/ParticleOne.f90
# ${CMAKE_CURRENT_SOURCE_DIR}/mcc/MCCInit.F90 ${CMAKE_CURRENT_SOURCE_DIR}/mcc/ParticleOneIndex.f90 ${CMAKE_CURRENT_SOURCE_DIR}/mcc/MCCChi.f90 ${CMAKE_CURRENT_SOURCE_DIR}/mcc/MCCElectron.F90 ${CMAKE_CURRENT_SOURCE_DIR}/mcc/MCCIon.F90) 
add_library(cxx-call-f-f ${SOURCES} )
set_property(TARGET cxx-call-f-f PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(cxx-call-f-f MPI::MPI_Fortran)

# add_executable(cxx-call-f-c ${CMAKE_CURRENT_SOURCE_DIR}/cxx-call-f.c)
# target_link_libraries(cxx-call-f-c cxx-call-f-f)main.cpp fieldsolver.cpp material.cpp
# set_property(TARGET cxx-call-f-c PROPERTY LINKER_LANGUAGE C) ${SOURCES}
include_directories ("${MPI_CXX_INCLUDE_DIRS}")
file(GLOB_RECURSE SOURCESCXX ${CMAKE_CURRENT_SOURCE_DIR}/pic/*)  
# file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/ main.cpp fieldsolver.cpp material.cpp )  
# add_executable(cxx-call-f-cxx ${CMAKE_CURRENT_SOURCE_DIR}/pic/main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/pic/particle.cpp 
# ${CMAKE_CURRENT_SOURCE_DIR}/pic/fieldsolver.cpp ${CMAKE_CURRENT_SOURCE_DIR}/pic/material.cpp 
#   ${CMAKE_CURRENT_SOURCE_DIR}/pic/create_particles.cpp ${CMAKE_CURRENT_SOURCE_DIR}/pic/domain.cpp  )
# ${CMAKE_CURRENT_SOURCE_DIR}/src/particle.cpp 
# ${CMAKE_CURRENT_SOURCE_DIR}/src/mixture.cpp  )
add_executable(cxx-call-f-cxx ${SOURCESCXX})
target_link_libraries(cxx-call-f-cxx cxx-call-f-f)
target_link_libraries(cxx-call-f-cxx ${MPI_LIBRARIES})
target_link_libraries(cxx-call-f-cxx ${PETSc_LIB})
target_link_libraries(cxx-call-f-cxx ${HDF5_LIB})
# target_link_libraries(cxx-call-f-cxx ${PETSC_LIBRARIES})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/run)
file(COPY
${CMAKE_CURRENT_SOURCE_DIR}/input/
    DESTINATION
    ${PROJECT_SOURCE_DIR}/run/input)