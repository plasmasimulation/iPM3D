#include"domain.h"
#include"stdio.h"
Domain::Domain(int* xyz_np ,int Mx,int My, int Mz,double *dx  ,int rank){
   this->xyz_np[0]=xyz_np[0];
   this->xyz_np[1]=xyz_np[1];
   this->xyz_np[2]=xyz_np[2];
   this->Mx=Mx;
   this->My=My;
   this->Mz=Mz;
   this->rank=rank;
   this->col=(rank%(xyz_np[0]*xyz_np[1]))%xyz_np[0];
   this->row=(rank%(xyz_np[0]*xyz_np[1]))/xyz_np[0];
   this->layer=rank/(xyz_np[0]*xyz_np[1]);
   lo[0]=(Mx/xyz_np[0])*col;
   lo[1]=(My/xyz_np[1])*row;
   lo[2]=(Mz/xyz_np[2])*layer;
   hi[0]=lo[0]+Mx/xyz_np[0];
   hi[1]=lo[1]+My/xyz_np[1];
   hi[2]=lo[2]+Mz/xyz_np[2];
   this->dx[0]=dx[0];
   this->dx[1]=dx[1];
   this->dx[2]=dx[2];
    printf("zuobiao,%d,%f,%f,%f,%f,\n",rank,lo[0],lo[1],lo[2],hi[0]);

   
   
}