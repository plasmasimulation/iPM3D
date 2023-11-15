program fortran_mpi
    use ModuleControlFlow
    use ModuleFieldSolver
    use ModulePICCommunication

    implicit none

    type(ControlFlow) :: CF
    type(FieldSolver) :: FS
    integer(4) :: size, rank, ierr, i, j
    type(PICCom2D) :: mycom

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    call InitFileName()
    
    call CF%Init()
    call CF%DM%Init(17, 33, 0.d0, 1.d0, 0.d0, 2.d0, 2, 4)
    call FS%Init(CF, 1)
    call mycom%init(CF%DM%LocalShape(1), CF%DM%LocalShape(2), 2, 4, com_type=com_type_nonblock)

    associate(Nz => CF%DM%GlobalShape(1), Nr => CF%DM%GlobalShape(2), &
                zstart => CF%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                zend => CF%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                rstart => CF%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                rend => CF%DM%CornerIndex(BOUNDARY_RIGHT, 2), &
                dz => CF%DM%SpaceStep(1), dr => CF%DM%SpaceStep(2), &
                coeA => FS%CoeA, coeB => FS%CoeB, coeC => FS%CoeC, &
                coeD => FS%CoeD, coeE => FS%CoeE, source => FS%Source)

        do j = rstart, rend                                 ! set coe and source for RZ
            do i = zstart, zend
                if (Nr == j .or. 1 == i .or. Nz == i) then  ! 边界
                    source(i, j) = 0.d0
                
                    coeA(i, j) = 0.d0                       ! the value is not use
                    coeB(i, j) = 0.d0
                    coeC(i, j) = 0.d0
                    coeD(i, j) = 0.d0
                    coeE(i, j) = 0.d0

                elseif (1 == j) then                       ! 中轴
                    source(i, j) = -10.d0

                    coeA(i, j) =  1.d0
                    coeC(i, j) =  1.d0
                    coeD(i, j) =  0.d0
                    coeE(i, j) =  4.d0
                    coeB(i, j) = -(coeA(i, j) + coeC(i, j) + coeD(i, j) + coeE(i, j))

                else
                    source(i, j) = -10.d0

                    coeA(i, j) =  1.d0
                    coeC(i, j) =  1.d0
                    coeD(i, j) =  (j - 1.5d0) / (j - 1.d0)
                    coeE(i, j) =  (j - 0.5d0) / (j - 1.d0)
                    coeB(i, j) = -(coeA(i, j) + coeC(i, j) + coeD(i, j) + coeE(i, j))

                end if
            end do
        end do

        ! solve
        call FS%SolveOne()
        call mycom%comf(FS%Solve, com_field_opt_sum, zstart, zend, rstart, rend)
        call FS%Diag('result')

    end associate

    ! do j = 0, size-1
    !     if (rank == j) then
    !         write(*, *) ""
    !         write(*, '(*(1x, i))') CF%DM%DomainIndex
            
    !         do i = CF%DM%CornerIndex(BOUNDARY_LEFT, 2), CF%DM%CornerIndex(BOUNDARY_RIGHT, 2)
    !             write(*, '(*(f10.6, 1x))') FS%Solve(:, i)
    !         end do
            
    !         call execute_command_line(" ")
    !     end if
    !     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    ! end do

    call CF%Destroy()
    call FS%Destroy()
    call mycom%destroy()
    call MPI_FINALIZE(ierr)

end program fortran_mpi