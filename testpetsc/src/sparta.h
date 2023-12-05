

// #include "stdio.h"

class SPARTA {
 public:

  // fundamental SPARTA classes

   class Memory *memory;          // memory allocation functions

  // class Particle *particle;      // particles


 
  // other top-level SPARTA classes and variables

  // SPARTA(int, char **, MPI_Comm);
  SPARTA();
  ~SPARTA();
 
  void init();
  void destroy();

};



