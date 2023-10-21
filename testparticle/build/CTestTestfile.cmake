# CMake generated Testfile for 
# Source directory: /home/sapphire/sparta/iPM3D/testparticle
# Build directory: /home/sapphire/sparta/iPM3D/testparticle/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(particlecomm "mpirun" "-np" "8" "/home/sapphire/sparta/iPM3D/testparticle/build/particlecomm")
set_tests_properties(particlecomm PROPERTIES  _BACKTRACE_TRIPLES "/home/sapphire/sparta/iPM3D/testparticle/CMakeLists.txt;24;add_test;/home/sapphire/sparta/iPM3D/testparticle/CMakeLists.txt;0;")
