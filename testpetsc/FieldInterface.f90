module FieldInterface
use mpi
use ModuleFieldOne
use ModuleFieldEM
use ModulePICFieldCom
use iso_c_binding
implicit none
    ! real(8),allocatable:: Field(:,:,:,:)
    integer(8):: dimension=3,i,j,k,p,q,m
    real(8):: dx=1,dy=1,dz=1
    type(FieldOne) FO
    type(FieldEM) FG
    
    type(PICCom3D) :: mycom
    real(8), allocatable :: array(:, :,:),rho(:,:,:)
    integer(4) :: xstart, xend, ystart, yend,zstart,zend
   

    
     
contains

subroutine GetRho(xstart, ystart,zstart,  width_x,  width_y, width_z,  rho1)bind(C, name="GetRho")
    
    integer(C_INT), value :: xstart, ystart,zstart,  width_x,  width_y, width_z
    real(C_Double) rho1(width_z,width_y,width_x)

    integer(8) xend,yend,zend
    xend=xstart-1+width_x
    yend=ystart-1+width_y
    zend=zstart-1+width_z
    allocate(FO%RhoOne(xstart:xend,ystart:yend,zstart:zend))
    ! allocate( rho1(xstart:xend))
    FO%RhoOne=0
 !设置电荷密度,后期耦合粒子后添加，暂时设置为0
    !  rho1=0
    write(*,"(3f9.3,/)")rho1
    ! do  i = xstart, xend
    !     do j=ystart, yend
    !         do k=zstart,zend
    !             rho(width_x*(i-xstart)+width_y*(j-ystart)+k-zstart+1)=FO%RhoOne(i,j,k)
    !         end do
    !     end do
    ! end do    
end subroutine GetRho

subroutine SendPhi( coor_x, coor_y,coor_z,  width_x,  width_y, width_z, phi)bind(C, name="SendPhi")

    integer(C_INT), value ::  coor_x, coor_y,coor_z,width_x,  width_y, width_z
    real(C_FLOAT) phi(width_x*width_y*width_z)
    integer(4) :: lx = 2, ly = 2,lz=2 !边的个数,两条边三个点 
    integer(4) :: xyz_np(3)=[2,2,2],com_type_nonblock
    ! integer(4) :: xstart, xend, ystart, yend,zstart,zend
    integer(4) :: size, rank, ierr
    character(len=99) :: file_name
    call mycom%init(lx, ly, lz,xyz_np, com_type_nonblock)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    
    xstart = mycom%col * lx
    xend   = (mycom%col+1) * lx
    ystart = mycom%row * ly
    yend   = (mycom%row+1) * ly
    zstart=mycom%layer*lz
    zend=(mycom%layer+1)*lz 
   
    allocate(array(xstart-1: xend+1, ystart-1:yend+1,zstart-1:zend+1))
    array=0
    !获取petsc数据，保存在array数组中
    do  i = coor_x,coor_x+width_x-1
        do j=coor_y, coor_y+width_y-1
            do k=coor_z,coor_z+width_z-1
                p=i-coor_x
                q=j-coor_y
                m=k-coor_z
                array(i,j,k)=phi(p*width_y*width_z+q*width_z+m+1)
   
            end do
                
        end do 
    end do

    !保存场通讯前后的数据，petsc场通讯无法获取角点数据，因此还是采用了自己编写的场通讯
    
    write(file_name, '(i1)') rank
    open(13, file="./solve/initphi"//trim(file_name)//".txt")
    write(13,*)"通讯前petsc求解出的电势,场通讯之前数据并不完整&
    每个只有一小部分,此处因为petsc包含了边界,但边界不是每个网格都有,因此又多&
    了外面一圈"
    do i=zstart-1,zend+1
        write(13,*)" "
        write(13,"(5f9.3,/)")array(:,:,i)
    end do
  

    close(13)
    call mycom%comfext(array,xstart, xend, ystart, yend,zstart,zend)

    open(13, file="./solve/finalphi"//trim(file_name)//".txt")
    write(13,*)"通讯后的电势,此处因为petsc包含了边界,但边界不是每个网格都有，因此又多&
    了外面一圈,这一圈在有边界的地方是无效电势,值为初值0,后面会舍去,电势存储顺序是&
    先遍历x和y,再遍历z,为便于观察每一层都用空格隔开"
    do i=zstart-1,zend+1
    write(13,*)" "
    write(13,"(5f9.3,/)")array(:,:,i)
    end do
   
    close(13)

!根据电势求解电场
    allocate(FG%Ex(xstart:xend+1,ystart:yend+1,zstart:zend+1))
    allocate(FG%Ey(xstart:xend+1,ystart:yend+1,zstart:zend+1))
    allocate(FG%Ez(xstart:xend+1,ystart:yend+1,zstart:zend+1))
    FG%Ex(xstart:xend+1,ystart:yend+1,zstart:zend+1)=&
    array(xstart:xend+1,ystart:yend+1,zstart:zend+1)-&
    array(xstart-1:xend,ystart:yend+1,zstart:zend+1)
    FG%Ey(xstart:xend+1,ystart:yend+1,zstart:zend+1)=&
    array(xstart:xend+1,ystart:yend+1,zstart:zend+1)-&
    array(xstart:xend+1,ystart-1:yend,zstart:zend+1)
    FG%Ez(xstart:xend+1,ystart:yend+1,zstart:zend+1)=&
    array(xstart:xend+1,ystart:yend+1,zstart:zend+1)-&
    array(xstart:xend+1,ystart:yend+1,zstart-1:zend)

    !下面是处理冗余边界，petsc已将边界保存在网格中，在场通讯时在
    !边界外又多留了一格，有冗余
     if(mycom%col==0) then
    xstart=xstart+1
    end if
    if(mycom%row==0) then
    ystart=ystart+1
    end if
    if(mycom%layer==0)  then
    zstart=zstart+1
    end if
    if(mycom%col==1) then
    xend=xend-1
    end if
    if(mycom%row==1) then
    yend=yend-1
    end if
    if(mycom%layer==1) then
    zend=zend-1
    end if

end subroutine SendPhi
subroutine Finalize()bind(C, name="Finalize")
    call FO%destroy()
    call FG%dump(xstart, xend, ystart, yend,zstart,zend)
    call FG%destroy()
end subroutine
end module FieldInterface