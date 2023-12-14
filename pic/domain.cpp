#include"domain.h"
#include"stdio.h"
Domain::Domain(int* xyz_np ,int Mx,int My, int Mz,double dx,double dy, double dz ,int rank){
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
   lo[0]=col;
   lo[1]=row;
   lo[2]=layer;
   hi[0]=col+1;
   hi[1]=row+1;
   hi[2]=layer+1;
   this->dx=dx;
   this->dy=dy;
   this->dz=dz;
   // printf("zuobiao,%f,%f,%f,%f",lo[0],lo[1],lo[2],hi[0]);

   
   
}