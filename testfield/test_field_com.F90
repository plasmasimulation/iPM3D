program fortran_mpi
    use mpi
    use ModulePICCommunication
    implicit none

    integer(4) :: size, rank, ierr, i, j
    type(PICCom2D) :: mycom
    real(8), allocatable :: array(:, :), array_ext(:, :)
    real(8), allocatable :: array_out(:, :)
    character(len=99) :: file_name

    integer(4) :: lx = 4, ly = 5
    integer(4) :: px = 4, py = 2
    integer(4) :: xstart, xend, ystart, yend

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    call mycom%init(lx, ly, lz,xyz_np, com_type_nonblock)

    xstart = dble(mycom%col)
    xend = dble(mycom%col+1)
    ystart = dble(mycom%row)
    yend = dble(mycom%row+1)
    zstart=dble(mycom%layer)
    zend=dble(mycom%layer+1)

    allocate(array(xstart:xend, ystart:yend))
    allocate(array_ext(xstart-1:xend+1, ystart-1:yend+1))
    allocate(array_out(lx*px, ly*py))

    do i = xstart, xend
        do j = ystart, yend
            array(i, j) = i + j
        end do
    end do

    array_ext = -1.d0
    array_ext(xstart:xend, ystart:yend) = array(xstart:xend, ystart:yend)

    mycom%left_ext_type = com_field_ext_type_symmetry

    call mycom%comf(array, com_field_opt_sum, xstart, xend, ystart, yend)
    ! call mycom%comf(array_ext, com_field_opt_ext, xstart, xend, ystart, yend)

    ! call mycom%gather(array, array_out, xstart, xend, ystart, yend, lx*px, ly*py)

    if (0 == rank) then
        open(10, file="./result_field_com.txt")
            do i = 1, ly*py
                write(10, '(*(f10.4, 1x))') array_out(:, i)
            end do
        close(10)
    end if

    write(file_name, '(i2)') rank
    open(10, file="./result_field_com_"//trim(file_name)//".txt")
        do i = ystart-1, yend+1
            write(10, '(*(f10.4, 1x))') array_ext(:, i)
        end do
    close(10)

    deallocate(array)
    deallocate(array_ext)
    deallocate(array_out)

    call mycom%destroy()
    call MPI_FINALIZE(ierr)

end program fortran_mpi