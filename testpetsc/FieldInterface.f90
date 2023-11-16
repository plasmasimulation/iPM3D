
module mpi_field_test
    use mpi
   
    use iso_c_binding
    implicit none
    real(8),allocatable:: Field(:,:,:,:)
    integer(8) dimension=3
    
    
contains


subroutine GetRho( coord_x, coord_y,coord_z,  width_x,  width_y, width_z, rho)bind(C, name="GetRho")

    integer(C_INT), value :: coord_x, coord_y, coord_z, width_x, width_y, width_z
    integer(C_INT) rho(width_x*width_y*width_z)
    rho=0 !设置电荷密度
end subroutine GetRho
subroutine SendPhi( coord_x, coord_y,coord_z,  width_x,  width_y, width_z, phi)bind(C, name="SendPhi")

    integer(C_INT), value :: coord_x, coord_y, coord_z, width_x, width_y, width_z
    integer(C_INT) phi(width_x*width_y*width_z)
    !根据电势求出电场强度
    
    allocate(Field(dimension,coord_x:(coord_x+width_x),coord_y:(coord_y+width_y),coord_z:(coord_z+width_z)))
    do i = zstart-1, zend
        Ez(i, :) = Phi(i, rstart:rend) - Phi(i+1, rstart:rend)
    end do
    Ez = Ez / dz

    do i = rstart-1, rend
        Er(:, i) = Phi(zstart:zend, i) - Phi(zstart:zend, i+1)
    end do
    Er = Er / dr
end subroutine SendPhi
end module mpi_field_test