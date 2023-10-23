program fortran_mpi
    use mpi
    use ModulePICFieldCom
    implicit none

    integer(4) :: size, rank, ierr, i,j, k
    type(PICCom3D) :: mycom
    real(8), allocatable :: array(:, :,:), array_ext(:, :)
    real(8), allocatable :: array_out(:, :,:)
    character(len=99) :: file_name

    integer(4) :: lx = 5, ly = 5,lz=5
    integer(4) :: xyz_np(3)=[2,2,2],com_type_nonblock, com_field_opt_sum
    integer(4) :: xstart, xend, ystart, yend,zstart,zend

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

     call mycom%init(lx, ly, lz,xyz_np, com_type_nonblock)

    xstart = mycom%col * lx + 1
    xend   = (mycom%col+1) * lx
    ystart = mycom%row * ly + 1
    yend   = (mycom%row+1) * ly
    zstart=mycom%layer*lz+1
    zend=(mycom%layer+1)*lz

    allocate(array(xstart:xend, ystart:yend,zstart:zend))
    ! allocate(array_ext(xstart-1:xend+1, ystart-1:yend+1))
    ! allocate(array_out(lx*px, ly*py))

    do i = xstart, xend
        do j = ystart, yend
            do k=zstart,zend
            array(i, j,k) = 1
        end do
    end do
    end do

    ! array_ext = -1.d0
    ! array_ext(xstart:xend, ystart:yend) = array(xstart:xend, ystart:yend)

    ! mycom%left_ext_type = com_field_ext_type_symmetry
    write(file_name, '(i1)') rank
    open(10, file="raw_field_"//trim(file_name)//".txt")
        write(10,*)array
    close(10)

     call mycom%comf(array, com_field_opt_sum, xstart, xend, ystart, yend,zstart,zend)
    ! call mycom%comf(array_ext, com_field_opt_ext, xstart, xend, ystart, yend)

    !  call mycom%gather(array, array_out, xstart, xend, ystart, yend, lx*px, ly*py)

    
     ! dump
    write(file_name, '(i1)') rank
    open(10, file="final_field_"//trim(file_name)//".txt")
        write(10,*)array
    close(10)

    ! write(file_name, '(i2)') rank
    ! open(10, file="./result_field_com_"//trim(file_name)//".txt")
    !     do i = ystart-1, yend+1
    !         write(10, '(*(f10.4, 1x))') array_ext(:, i)
    !     end do
    ! close(10)

    deallocate(array)
    ! deallocate(array_ext)
    ! deallocate(array_out)

    ! call mycom%destroy()
     call MPI_FINALIZE(ierr)

end program fortran_mpi