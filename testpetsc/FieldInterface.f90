module FieldInterface
use mpi
use ModuleFieldOne
use ModuleParticleBundle
use iso_c_binding
implicit none
    real(8),allocatable:: Field(:,:,:,:)
    integer(8) dimension=3,i,j,k
    real(8) dx,dy,dz
    type(FieldOne) FO
    type(ParticleBundle) PB

     
contains
subroutine FieldOneInitial()bind(C, name="ieldOneInitial")
PB%Init()
FO%Init()
end subroutine
subroutine GetRho(xstart, ystart,zstart,  xend,  yend, zend,  rho)bind(C, name="GetRho")
    
    integer(C_INT), value :: xstart, ystart,zstart,  xend,  yend, zend
    integer(C_INT) rho((xend-xstart+1)*(yend-ystart+1)*(zend-zstart+1))
    rho=0 !设置电荷密度
end subroutine GetRho

subroutine SendPhi( xstart, ystart,zstart,  xend,  yend, zend, phi)bind(C, name="SendPhi")

    integer(C_INT), value ::  xstart, ystart,zstart,  xend,  yend, zend
    integer(C_INT) phi((xend-xstart+1)*(yend-ystart+1)*(zend-zstart+1))
    !根据电势求出电场强度
    
    allocate(Field(dimension,xstart:xend,ystart:yend,zstart:zend))
    do  i = zstart, zend
        do j=ystart, yend
            do z=zstart,zend
                
            end do
        end do
    end do    
    ! do i = zstart, zend
    !     Field(1,xstart:xend, ystart:yend,i) &
    !     = (Phi(1,xstart:xend, ystart:yend,i) - Phi(1,xstart:xend, ystart:yend,i-1))/dz
    ! end do
    ! do j = ystart, yend
    !     Field(1,xstart:xend,i, zstart:zend) &
    !     = (Phi(1,xstart:xend, i,zstart:zend) - Phi(1,xstart:xend,i-1, zstart:zend))/dy
    ! end do
    ! do k = xstart, xend
    !     Field(1,xstart:xend,i, zstart:zend) &
    !     = (Phi(1,xstart:xend, i,zstart:zend) - Phi(1,xstart:xend,i-1, zstart:zend))/dz
    ! end do
    deallocate(Field)
end subroutine SendPhi
end module FieldInterface