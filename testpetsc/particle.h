#define SPARTA_ALIGN(n)
#include "stdio.h"
#include "pointers.h"
#include"spatype.h"



class Particle  {
 public:
  int exist;                // 1 if particles exist
  int sorted;               // 1 if particles are sorted by grid cell

  enum{MAXVIBMODE=4};       // increase value if species need more vib modes

  struct Species {          // info on each particle species, read from file
    char id[16];            // species ID
    double molwt;           // molecular weight
    double mass;            // molecular mass
    double specwt;          // species weight
    double charge;          // multiple of electron charge
    double rotrel;          // inverse rotational relaxation number
    double rottemp[3];      // rotational temperature(s)
    double vibtemp[MAXVIBMODE];   // vibrational temperature(s)
    double vibrel[MAXVIBMODE];    // inverse vibrational relaxation number(s)
    int vibdegen[MAXVIBMODE];     // vibrational mode degeneracies
    int rotdof,vibdof;      // rotational/vibrational DOF
    int nrottemp,nvibmode;  // # of rotational/vibrational temps/modes defined
    int internaldof;        // 1 if either rotdof or vibdof != 0
    int vibdiscrete_read;   // 1 if species.vib file read for this species
    double magmoment;       // magnetic moment, set by species_modify command
  };

  struct RotFile {          // extra rotation info read from rotfile
    char id[16];
    double rottemp[4];
    int ntemp;
  };

  struct VibFile {          // extra vibration info read from vibfile
    char id[16];
    double vibrel[MAXVIBMODE];
    double vibtemp[MAXVIBMODE];
    int vibdegen[MAXVIBMODE];
    int nmode;
  };

  Species *species;         // list of particle species info
  int nspecies;             // # of defined species
  int maxvibmode;           // max vibmode of any species (mode = dof/2)

  class Mixture **mixture;
  int nmixture;
  int maxmixture;

  struct SPARTA_ALIGN(64) OnePart {
    int id;                 // particle ID
    int ispecies;           // particle species index
    int icell;              // which local Grid::cells the particle is in
    double x[3];            // particle position
    double v[3];            // particle velocity
    double erot;            // rotational energy
    double evib;            // vibrational energy
    int flag;               // used for migration status
    double dtremain;        // portion of move timestep remaining
    double weight;          // particle or cell weight, if weighting enabled
  };

  struct OnePartRestart {
    int id;                 // particle ID
    int ispecies;           // particle species index
    cellint icell;          // cell ID the particle is in
    int nsplit;             // 1 for unsplit cell
                            // else neg of sub cell index (0 to Nsplit-1)
    double x[3];            // particle position
    double v[3];            // particle velocity
    double erot;            // rotational energy
    double evib;            // vibrational energy
  };

  bigint nglobal;           // global # of particles
  int nlocal;               // # of particles I own
  int maxlocal;             // max # particles list can hold
  OnePart *particles;       // list of particles I own

  // currently stored in grid.h for every cell, whether I own it or not
  // not sure why storing it here is slower

  //int *cellcount;           // count of particles in each grid cell I own
  //int *first;               // index of first particle in each grid cell

  int *next;                // index of next particle in each grid cell

  // extra custom vectors/arrays for per-particle data
  // ncustom > 0 if there are any extra arrays
  // custom attributes are created by various commands
  // these variables are public, others below are private

  int ncustom;              // # of custom attributes, some may be deleted
  int *etype;               // type = INT/DOUBLE of each attribute
  int *esize;               // size = 0 for vector, N for array columns
  int *ewhich;              // index into eivec,eiarray,edvec,edarray for data

  int **eivec;              // pointer to each integer vector
  int ***eiarray;           // pointer to each integer array
  double **edvec;           // pointer to each double vector
  double ***edarray;        // pointer to each double array

  // restart buffers, filled by read_restart

  int nlocal_restart;
  char *particle_restart;

  // Kokkos settings

  int copy,copymode;        // 1 if copy of class (prevents deallocation of
                            //  base class when child copy is destroyed)

  // methods

  Particle(class SPARTA *);
  virtual ~Particle();
  void init();
  virtual void compress_migrate(int, int *);
  void compress_rebalance();
  void compress_rebalance_sorted();
  void compress_reactions(int, int *);
  void sort();
  void sort_allocate();
  void remove_all_from_cell(int);
  virtual void grow(int);
  virtual void grow_species();
  void grow_next();
  virtual void pre_weight();
  virtual void post_weight();

  virtual int add_particle(int, int, int, double *, double *, double, double);
  virtual int add_particle();
  int clone_particle(int);
  void add_species(int, char **);
  int find_species(char *);
  void species_modify(int, char **);
  void add_mixture(int, char **);
  int find_mixture(char *);
  double erot(int, double, class RanKnuth *);
  double evib(int, double, class RanKnuth *);

  void write_restart_species(FILE *fp);
  void read_restart_species(FILE *fp);
  void write_restart_mixture(FILE *fp);
  void read_restart_mixture(FILE *fp);

  int size_restart();
  bigint size_restart_big();
  int pack_restart(char *);
  void pack_restart(char *, int, int);
  int unpack_restart(char *);
  void unpack_restart(char *, int &, int, int);

  int find_custom(char *);
  void error_custom();
  virtual int add_custom(char *, int, int);
  virtual void grow_custom(int, int, int);
  virtual void remove_custom(int);
  virtual void copy_custom(int, int);
  int sizeof_custom();
  void write_restart_custom(FILE *fp);
  void read_restart_custom(FILE *fp);
  virtual void pack_custom(int, char *);
  virtual void unpack_custom(char *, int);

  bigint memory_usage();

 protected:
  int me;
  int maxgrid;              // max # of indices first can hold
  int maxsort;              // max # of particles next can hold
  int maxspecies;           // max size of species list

  FILE *fp;                 // file pointer for species, rotation, vibration
  int nfile;                // # of species read from file
  int maxfile;              // max size of file list

  Species *filespecies;     // list of species read from file
  RotFile *filerot;         // list of species rotation info read from file
  VibFile *filevib;         // list of species vibration info read from file

  class RanKnuth *wrandom;   // RNG for particle weighting

  // extra custom vectors/arrays for per-particle data
  // ncustom > 0 if there are any extra arrays
  // these variables are private, others above are public

  char **ename;             // name of each attribute

  int ncustom_ivec;         // # of integer vector attributes
  int ncustom_iarray;       // # of integer array attributes
  int *icustom_ivec;        // index into ncustom for each integer vector
  int *icustom_iarray;      // index into ncustom for each integer array
  int *eicol;               // # of columns in each integer array (esize)

  int ncustom_dvec;         // # of double vector attributes
  int ncustom_darray;       // # of double array attributes
  int *icustom_dvec;        // index into ncustom for each double vector
  int *icustom_darray;      // index into ncustom for each double array
  int *edcol;               // # of columns in each double array (esize)

  int *custom_restart_flag; // flag on each custom vec/array read from restart
                            // used to delete them if not redefined in
                            // restart script

  // private methods

  void read_species_file();
  void read_rotation_file();
  void read_vibration_file();
  int wordcount(char *, char **);
};


  subroutine VelocityMaxwellianInitializationParticleOne(PO,Mass,Temperature)
            Class(ParticleOne), intent(inout) :: PO
            real(8), intent(in) :: Mass,Temperature
            real(8) :: V,Beta,FuncA,FuncB
            real(8) :: Theta,CosTheta,SinTheta,Phi

            Beta=1.d0/(kB*Temperature)
            FuncA=1.d0
            FuncB=0.d0
            do while(FuncA>FuncB)
                call RANDOM_NUMBER(R)
                FuncA=R*R
                call RANDOM_NUMBER(R)
                FuncB=-exp*R*Dlog(R)
            end do
            V=DSqrt(-3.d0*Dlog(R)/Beta/Mass)
            call VelocityRandomInitializationParticleOne(PO,V)

