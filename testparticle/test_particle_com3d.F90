program fortran_mpi
    use mpi
    use ModulePICParticleCom3d
    implicit none

    integer(4) :: size, rank, ierr, i, j
    type(PICCom3D) :: mycom
    type(ParticleBundle) :: pb
    type(ParticleOne) :: one
    character(len=99) :: file_name
    real(8) :: RR
    real(8) :: xstart, xend, ystart, yend, zstart, zend
    real(8) :: x_lb, x_ub, y_lb, y_ub,z_lb,z_ub

    integer(4) :: lx = 1, ly = 1,lz=1
    integer(4) :: xyz_np(3)=[2,2,2]

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    call mycom%init(lx, ly, lz,xyz_np, com_type_nonblock)
    call pb%init(1000, 10, 8000, 0.25d0, 0.5d0, 1.5d0)

    xstart = dble(mycom%col)
    xend = dble(mycom%col+1)
    ystart = dble(mycom%row)
    yend = dble(mycom%row+1)
    zstart=dble(mycom%layer)
    zend=dble(mycom%layer+1)
    !write(*,*)mycom%rank,mycom%col,mycom%row,mycom%layer
    x_lb = 0.d0
    x_ub = dble(xyz_np(1))
    y_lb = 0.d0
    y_ub = dble(xyz_np(2))
    z_lb = 0.d0
    z_ub = dble(xyz_np(3))

    ! 生成粒子
    do i = 1, pb%size
        call random_number(RR)
        one%X = xstart - 1.d0 + (xend - xstart + 2.d0) * RR

        call random_number(RR)
        one%Y = ystart - 1.d0 + (yend - ystart + 2.d0) * RR

        call random_number(RR)
        one%Z = zstart - 1.d0 + (zend - zstart + 2.d0) * RR

        call pb%addone(one)
    end do

    ! 边界处理
    do i = pb%npar, 1, -1
        if (pb%PO(i)%X <= x_lb .or. &
            pb%PO(i)%X >= x_ub .or. &
            pb%PO(i)%Y <= y_lb .or. &
            pb%PO(i)%Y >= y_ub .or. &
            pb%PO(i)%Z <= z_lb .or. &
            pb%PO(i)%Z >= z_ub) then

            call pb%delone(i)
        end if
    end do

    ! dump
    write(file_name, '(i1)') rank
    open(10, file="raw_par_"//trim(file_name)//".txt")
        do i = 1, pb%npar
            write(10, '(*(f10.4, 1x))') pb%PO(i)%X, pb%PO(i)%Y,pb%PO(i)%Z
        end do
    close(10)
    ! write(*,*)"before",pb%npar
    ! com
    ! write(*,*)" xstart, xend, ystart, yend,zstart,zend", xstart, xend, ystart, yend,zstart,zend
    call mycom%comp(pb, xstart, xend, ystart, yend,zstart,zend)
    ! write(*,*)"after",pb%NPar
    ! dump
    write(file_name, '(i1)') rank
    open(10, file="final_par_"//trim(file_name)//".txt")
        do i = 1, pb%npar
            write(10, '(*(f10.4, 1x))') pb%PO(i)%X, pb%PO(i)%Y,pb%PO(i)%Z
        end do
    close(10)

    !call mycom%destroy()
    call pb%destroy()
    call MPI_FINALIZE(ierr)

end program fortran_mpi