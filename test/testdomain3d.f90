program fortran_mpi
    use mpi
    use ModuleDomain3d

    implicit none

    integer(4) :: size, rank, ierr, i
    type(Domain) :: dm

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

   ! call InitFileName()

    call dm%Init(17, 33, 33, 0.d0, 1.d0, 0.d0, 2.d0,0.d0,2.d0, 2, 2, 2)
    call dm%Specify([0.2d0, 0.8d0], [0.5d0, 0.5d0],[0.5d0, 0.5d0])
    call dm%ShowAll()
    call dm%Show()
   ! call dm%Save()
    !call dm%Dump()

    call MPI_FINALIZE(ierr)

end program fortran_mpi