rm -rf build


export FC=mpiifort
export F9X=mpiifort
export CC=mpiicc
export CXX=mpiicpc

# export FC=ifort
# export F9X=ifort
# export CC=icc
# export CXX=icx

cmake -S . -B build
cmake --build build -j 8
