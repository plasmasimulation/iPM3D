

// #include "spatype.h"
#include "mpi.h"
#include "stdio.h"
#ifndef sparta.h
#define sparta.h
//  #include "sparta.h"
//  #include "memory.h"
//  #include "particle.h"

// universal defines inside namespace

#define FLERR __FILE__,__LINE__

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

// roundup a char ptr to 8-byte boundary
// roundup an int to multiple of 8

#define ROUNDUP(A) (char *) (((uint64_t) (A) + 7) & ~7)
#define IROUNDUP(A) ((((int) (A) + 7) / 8) * 8)
#define BIROUNDUP(A) ((((bigint) (A) + 7) / 8) * 8)

class Pointers {
 public:
  Pointers(SPARTA *ptr) :
    sparta(ptr),
    memory(ptr->memory){}
   
    // particle(ptr->particle){}
   
  virtual ~Pointers() {}

 protected:
  SPARTA *sparta;
  Memory *&memory;
  // Particle *&particle;


};

#endif