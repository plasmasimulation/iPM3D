

#include "mpi.h"
#include "string.h"
#include "ctype.h"
#include "sparta.h"
// #include "particle.h"
// #include "memory.h"


/* ----------------------------------------------------------------------
   start up SPARTA
   allocate fundamental classes (memory, error, universe, input)
   parse input switches
   initialize communicators, screen & logfile output
   input is allocated at end after MPI info is setup
------------------------------------------------------------------------- */

SPARTA::SPARTA()
{
//   memory = new Memory(this);
 

  
}

/* ----------------------------------------------------------------------
   shutdown SPARTA
   delete top-level classes
   delete fundamental classes
   close screen and log files in world and universe
   output files were already closed in destroy()
   delete fundamental classes
------------------------------------------------------------------------- */

SPARTA::~SPARTA()
{
//    delete particle;
//   delete memory;
}

/* ----------------------------------------------------------------------
   initialize top-level classes
------------------------------------------------------------------------- */

void SPARTA::init()
{
  
}

/* ----------------------------------------------------------------------
   delete single instance of top-level classes
   fundamental classes are deleted in destructor
------------------------------------------------------------------------- */

void SPARTA::destroy()
{
  
//   delete particle;
 
}


