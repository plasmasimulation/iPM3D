// int testpetsc(  PetscInt Mx,  PetscInt  My,   PetscInt Mz ,int data[5][5][5]);
class FieldSolver{
    public:
    PetscErrorCode ierr = 0;
    Mat A;
    Vec b, x,m;
    KSP ksp;
    DM dm;
    PetscMPIInt size, rank;
    PetscScalar v[7];
    MatStencil row, col[7];
     PetscInt coord_x, coord_y, coord_z, width_x, width_y, width_z,Mx,My,Mz;
     double *lo;
     double *hi;
    double ***array=nullptr,***geometry=nullptr,***barray = nullptr;
    // int ***data=nullptr;
    double dx[3];
    
    int initpetsc(PetscInt Mx,  PetscInt My, PetscInt Mz,int* xyz_np,double dx[3],double *lo,double*hi);
     int fieldsolve(double wt); 
    // int petscfinalsize();


};