// int testpetsc(  PetscInt Mx,  PetscInt  My,   PetscInt Mz ,int data[5][5][5]);
class FieldSolver{
    public:
    PetscErrorCode ierr = 0;
    Mat A;
    Vec b, x,d;
    KSP ksp;
    DM dm;
    PetscMPIInt size, rank;
    PetscScalar v[7];
    MatStencil row, col[7];
     PetscInt coord_x, coord_y, coord_z, width_x, width_y, width_z,Mx,My,Mz;
     double *lo;
     double *hi;
    double ***barray = nullptr,***array=nullptr,***geometry=nullptr;
    int initpetsc(PetscInt Mx,  PetscInt My, PetscInt Mz,int data [5][5][5],int* xyz_np);
     int fieldsolve(); 
    // int petscfinalsize();


};