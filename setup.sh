# g++ -c -fpic   test.f90 
# g++ -shared test.o  -o test.so
#g++ -shared thermal.o  -o thermal.so
rm -rf build
rm -rf run
export FC=mpiifort
export F9X=mpiifort
export CC=mpiicc 
export CXX=mpiicpc 
cmake -S . -B build
cmake --build build -j 8

mkdir run/solve
mkdir run/field

