#include"particle.h"

class CreateParticles  {

 public:
//   CreateParticles(class SPARTA *);
CreateParticles();
//   void command(int, char **);
//   int evib(int);
//   double erot(int);

 
  int imix,single,cutflag,mspecies,twopass;
  int np;
  double xp,yp,zp,vx,vy,vz;
//   class Region *region;

  int speciesflag,densflag,velflag,tempflag,normflag;
  char *sstr,*sxstr,*systr,*szstr;
  char *dstr,*dxstr,*dystr,*dzstr;
  char *tstr,*txstr,*tystr,*tzstr;
  char *vxstr,*vystr,*vzstr,*vstrx,*vstry,*vstrz;
  int svar,sxvar,syvar,szvar;
  int dvar,dxvar,dyvar,dzvar;
  int tvar,txvar,tyvar,tzvar;
  int vxvar,vyvar,vzvar,vvarx,vvary,vvarz;
  char *sxstr_copy,*systr_copy,*szstr_copy;
  char *dxstr_copy,*dystr_copy,*dzstr_copy;
  char *txstr_copy,*tystr_copy,*tzstr_copy;
  char *vstrx_copy,*vstry_copy,*vstrz_copy;

   void create_single();
   void create_local(Particle* particle);
  //  void particle_move(Particle* particle);
  void create_local_twopass();
  int species_variable(double *);
  double density_variable(double *, double *);
  double temperature_variable(double *);
  void velocity_variable(double *, double *, double *);
  int outside_region(int, double *, double *);
};

