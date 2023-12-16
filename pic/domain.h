class Domain{
public:
int xyz_np[3];
double lo[3];
double hi[3];
int rank;
int col,row,layer;
int Mx,My,Mz;
double dx[3];
Domain(int* xyz_np ,int Mx,int My, int Mz,double *dx ,int rank);

};