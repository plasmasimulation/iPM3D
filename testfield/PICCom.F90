module ModulePICCommunication
    use mpi
    use ModuleParticleBundle
    implicit none

    integer(4), parameter :: com_negb_type_boundary = 0
    integer(4), parameter :: com_negb_type_domain   = 1

    integer(4), parameter :: com_negb_left_bottom   = 1
    integer(4), parameter :: com_negb_bottom        = 2
    integer(4), parameter :: com_negb_right_bottom  = 3
    integer(4), parameter :: com_negb_right         = 4
    integer(4), parameter :: com_negb_right_top     = 5
    integer(4), parameter :: com_negb_top           = 6
    integer(4), parameter :: com_negb_left_top      = 7
    integer(4), parameter :: com_negb_left          = 8

    integer(4), parameter :: com_negb_num           = 8
    integer(4), parameter :: com_corner_num         = 4

    integer(4), parameter :: com_type_nonblock      = 0
    integer(4), parameter :: com_type_onesided      = 1

    integer(4), parameter :: com_field_opt_sum      = 0
    integer(4), parameter :: com_field_opt_ext      = 1

    integer(4), parameter :: com_field_ext_type_intep    = 0
    integer(4), parameter :: com_field_ext_type_symmetry = 1


    type :: PICCom2D
        integer(4) :: lx, ly
        integer(4) :: dpx, dpy
        integer(4) :: size, rank
        integer(4) :: com_type

        logical :: is_init = .False.
        integer(4) :: col, row

        integer(4) :: left_ext_type   = com_field_ext_type_intep
        integer(4) :: right_ext_type  = com_field_ext_type_intep
        integer(4) :: bottom_ext_type = com_field_ext_type_intep
        integer(4) :: top_ext_type    = com_field_ext_type_intep

        ! for non-block field com
        integer(4) :: reqs_left(2), reqs_right(2), reqs_top(2), reqs_bottom(2)
        integer(4) :: reqs_corner(8)

        integer(4) :: status_left(MPI_STATUS_SIZE, 2), status_right(MPI_STATUS_SIZE, 2)
        integer(4) :: status_top(MPI_STATUS_SIZE, 2), status_bottom(MPI_STATUS_SIZE, 2)
        integer(4) :: status_corner(MPI_STATUS_SIZE, 8)
        integer(4) :: recv_count, recv_count_corner(4)

        ! for one-sided field com
        integer(4) :: win_left, win_right, win_top, win_bottom
        integer(kind=MPI_ADDRESS_KIND) :: win_size, disp_aint

        ! negb info
        integer(4) :: negb_type(com_negb_num)
        integer(4) :: negb_rank(com_negb_num)

        ! field buffer
        real(8), allocatable :: send_buff_edge_left(:)
        real(8), allocatable :: send_buff_edge_right(:)
        real(8), allocatable :: send_buff_edge_bottom(:)
        real(8), allocatable :: send_buff_edge_top(:)
        real(8) :: send_buff_corner(com_corner_num)

        real(8), allocatable :: recv_buff_edge_left(:)
        real(8), allocatable :: recv_buff_edge_right(:)
        real(8), allocatable :: recv_buff_edge_bottom(:)
        real(8), allocatable :: recv_buff_edge_top(:)
        real(8) :: recv_buff_corner(com_corner_num)

        ! for non-block particle com
        integer(4) :: reqs_send(com_negb_num), reqs_recv(com_negb_num)
        integer(4) :: status_send(MPI_STATUS_SIZE, com_negb_num), status_recv(MPI_STATUS_SIZE, com_negb_num)

        ! particleone datatype
        integer(4) :: mpi_type_particle_one

    contains

        procedure :: init       => initPICCom2D
        procedure :: destroy    => destroyPICCom2D
        procedure :: comf       => communicationPICFieldCom2D
        procedure :: comp       => communicationPICParticleCom2D
        procedure :: gather     => gatherPICFieldCom2D

    end type PICCom2D


    contains

        subroutine initPICCom2D(this, x_local_size, y_local_size, deps_x, deps_y, com_type)
            class(PICCom2D), intent(inout) :: this
            integer(4), intent(in) :: x_local_size, y_local_size, deps_x, deps_y
            integer(4), optional, intent(in) :: com_type
            integer(4) :: ierr, i

            call MPI_COMM_SIZE(MPI_COMM_WORLD, this%size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, this%rank, ierr)

            if (x_local_size > 0 .and. y_local_size > 0 .and. deps_x > 0 .and. deps_y > 0) then
                this%lx = x_local_size
                this%ly = y_local_size
                this%dpx = deps_x
                this%dpy = deps_y
                this%com_type = com_type_nonblock

                if (present(com_type)) then
                    if (com_type == com_type_nonblock) this%com_type = com_type_nonblock
                    if (com_type == com_type_onesided) this%com_type = com_type_onesided
                end if

                ! 二维数组沿行和列分解，必须恰好整数
                if (this%dpx * this%dpy /= this%size) then
                    write(*, *) "The domian decomposition is error."
                    stop
                end if

                this%row = this%rank / this%dpx
                this%col = mod(this%rank, this%dpx)

                ! do i = 0, this%size-1
                !     if (this%rank == i) then
                !         write(*, '(a, i2, a, i2, a, i2, i2)') "Start: ", this%rank, "/", this%size, ", ", this%row, this%col
                !     end if

                !     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                ! end do
                
                ! 邻近domain的rank和type
                this%negb_rank = -1
                this%negb_type = com_negb_type_boundary

                if (this%row - 1 >= 0 .and. this%col - 1 >= 0) then                 ! pos 1
                    this%negb_rank(com_negb_left_bottom) = (this%row-1) * this%dpx + this%col - 1
                end if

                if (this%row - 1 >= 0) then                                         ! pos 2
                    this%negb_rank(com_negb_bottom) = (this%row-1) * this%dpx + this%col
                end if

                if (this%row - 1 >= 0 .and. this%col + 1 < this%dpx) then           ! pos 3
                    this%negb_rank(com_negb_right_bottom) = (this%row-1) * this%dpx + this%col + 1
                end if

                if (this%col + 1 < this%dpx) then                                   ! pos 4
                    this%negb_rank(com_negb_right) = this%row * this%dpx + this%col + 1
                end if

                if (this%row + 1 < this%dpy .and. this%col + 1 < this%dpx) then     ! pos 5
                    this%negb_rank(com_negb_right_top) = (this%row+1) * this%dpx + this%col + 1
                end if

                if (this%row + 1 < this%dpy) then                                   ! pos 6
                    this%negb_rank(com_negb_top) = (this%row+1) * this%dpx + this%col
                end if

                if (this%row + 1 < this%dpy .and. this%col - 1 >= 0) then           ! pos 7
                    this%negb_rank(com_negb_left_top) = (this%row+1) * this%dpx + this%col - 1
                end if

                if (this%col - 1 >= 0) then                                         ! pos 8
                    this%negb_rank(com_negb_left) = this%row * this%dpx + this%col - 1
                end if

                do i = 1, com_negb_num
                    if (this%negb_rank(i) >= 0) then
                        this%negb_type(i) = com_negb_type_domain
                    end if
                end do

                ! do i = 0, this%size-1
                !     if (this%rank == i) then
                !         write(*, '(i2, a, 8(i3), a, 8(i3))') this%rank, ": ", this%negb_rank, " | ", this%negb_type
                !     end if

                !     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                ! end do

                call this%destroy()
                allocate(this%send_buff_edge_left(this%ly))
                allocate(this%send_buff_edge_right(this%ly))
                allocate(this%send_buff_edge_bottom(this%lx))
                allocate(this%send_buff_edge_top(this%lx))

                allocate(this%recv_buff_edge_left(this%ly))
                allocate(this%recv_buff_edge_right(this%ly))
                allocate(this%recv_buff_edge_bottom(this%lx))
                allocate(this%recv_buff_edge_top(this%lx))

                ! 创建新的数据类型
                call MPI_TYPE_CONTIGUOUS(11, MPI_DOUBLE, this%mpi_type_particle_one, ierr)
                call MPI_TYPE_COMMIT(this%mpi_type_particle_one, ierr)

                this%is_init = .True.

                ! 创建远程访问窗口
                if (com_type_onesided == this%com_type) then
                    this%win_size = 8 * this%ly
                    call MPI_WIN_CREATE(this%send_buff_edge_left, this%win_size, 8, MPI_INFO_NULL, MPI_COMM_WORLD, this%win_left, ierr)

                    this%win_size = 8 * this%ly
                    call MPI_WIN_CREATE(this%send_buff_edge_right, this%win_size, 8, MPI_INFO_NULL, MPI_COMM_WORLD, this%win_right, ierr)

                    this%win_size = 8 * this%lx
                    call MPI_WIN_CREATE(this%send_buff_edge_bottom, this%win_size, 8, MPI_INFO_NULL, MPI_COMM_WORLD, this%win_bottom, ierr)

                    this%win_size = 8 * this%lx
                    call MPI_WIN_CREATE(this%send_buff_edge_top, this%win_size, 8, MPI_INFO_NULL, MPI_COMM_WORLD, this%win_top, ierr)
                end if

            else
                if (0 == this%rank) then
                    write(*, *) "The input parameters of the FieldCom2D are invalid."
                end if
                stop
            end if

        end subroutine initPICCom2D


        subroutine destroyPICCom2D(this)
            class(PICCom2D), intent(inout) :: this
            integer(4) :: ierr, i

            if (allocated(this%send_buff_edge_left))    deallocate(this%send_buff_edge_left)
            if (allocated(this%send_buff_edge_right))   deallocate(this%send_buff_edge_right)
            if (allocated(this%send_buff_edge_bottom))  deallocate(this%send_buff_edge_bottom)
            if (allocated(this%send_buff_edge_top))     deallocate(this%send_buff_edge_top)

            if (allocated(this%recv_buff_edge_left))    deallocate(this%recv_buff_edge_left)
            if (allocated(this%recv_buff_edge_right))   deallocate(this%recv_buff_edge_right)
            if (allocated(this%recv_buff_edge_bottom))  deallocate(this%recv_buff_edge_bottom)
            if (allocated(this%recv_buff_edge_top))     deallocate(this%recv_buff_edge_top)

            if (com_type_onesided == this%com_type .and. this%is_init) then
                call MPI_WIN_FREE(this%win_left, ierr)
                call MPI_WIN_FREE(this%win_right, ierr)
                call MPI_WIN_FREE(this%win_top, ierr)
                call MPI_WIN_FREE(this%win_bottom, ierr)
            end if

            if (this%is_init) call MPI_TYPE_FREE(this%mpi_type_particle_one, ierr)

            this%is_init = .False.

        end subroutine destroyPICCom2D


        subroutine communicationPICFieldCom2D(this, array2d, opt_type, xstart, xend, ystart, yend)
            class(PICCom2D), intent(inout) :: this
            real(8), allocatable, intent(inout) :: array2d(:, :)
            integer(4), intent(in) :: opt_type
            integer(4), intent(in) :: xstart, xend, ystart, yend
            integer(4) :: ierr

            if (this%is_init) then
                if (xend-xstart+1 /= this%lx .or. yend-ystart+1 /= this%ly) then
                    if (0 == this%rank) then
                        write(*, *) "The input parameters (xstart, xend, ystart, yend) are invalid."
                    end if
                    stop
                end if

                if (com_field_opt_sum == opt_type) then
                    if (com_type_nonblock == this%com_type) then
                        this%send_buff_edge_left    = array2d(xstart, ystart:yend)
                        this%send_buff_edge_right   = array2d(xend, ystart:yend)
                        this%send_buff_edge_bottom  = array2d(xstart:xend, ystart)
                        this%send_buff_edge_top     = array2d(xstart:xend, yend)

                        this%send_buff_corner(1)    = array2d(xstart, ystart)
                        this%send_buff_corner(2)    = array2d(xend, ystart)
                        this%send_buff_corner(3)    = array2d(xend, yend)
                        this%send_buff_corner(4)    = array2d(xstart, yend)

                        this%recv_buff_edge_left    = 0
                        this%recv_buff_edge_right   = 0
                        this%recv_buff_edge_bottom  = 0
                        this%recv_buff_edge_top     = 0
                        this%recv_buff_corner       = 0

                        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

                        ! send and recv edge
                        if (this%negb_type(com_negb_left) == com_negb_type_domain) then         ! left
                            call MPI_ISEND(this%send_buff_edge_left, this%ly, MPI_DOUBLE, this%negb_rank(com_negb_left), 1, MPI_COMM_WORLD, this%reqs_left(1), ierr)
                            call MPI_IRECV(this%recv_buff_edge_left, this%ly, MPI_DOUBLE, this%negb_rank(com_negb_left), 1, MPI_COMM_WORLD, this%reqs_left(2), ierr)
                        end if

                        if (this%negb_type(com_negb_right) == com_negb_type_domain) then        ! right
                            call MPI_ISEND(this%send_buff_edge_right, this%ly, MPI_DOUBLE, this%negb_rank(com_negb_right), 1, MPI_COMM_WORLD, this%reqs_right(1), ierr)
                            call MPI_IRECV(this%recv_buff_edge_right, this%ly, MPI_DOUBLE, this%negb_rank(com_negb_right), 1, MPI_COMM_WORLD, this%reqs_right(2), ierr)
                        end if

                        if (this%negb_type(com_negb_bottom) == com_negb_type_domain) then       ! bottom
                            call MPI_ISEND(this%send_buff_edge_bottom, this%lx, MPI_DOUBLE, this%negb_rank(com_negb_bottom), 1, MPI_COMM_WORLD, this%reqs_bottom(1), ierr)
                            call MPI_IRECV(this%recv_buff_edge_bottom, this%lx, MPI_DOUBLE, this%negb_rank(com_negb_bottom), 1, MPI_COMM_WORLD, this%reqs_bottom(2), ierr)
                        end if

                        if (this%negb_type(com_negb_top) == com_negb_type_domain) then          ! top
                            call MPI_ISEND(this%send_buff_edge_top, this%lx, MPI_DOUBLE, this%negb_rank(com_negb_top), 1, MPI_COMM_WORLD, this%reqs_top(1), ierr)
                            call MPI_IRECV(this%recv_buff_edge_top, this%lx, MPI_DOUBLE, this%negb_rank(com_negb_top), 1, MPI_COMM_WORLD, this%reqs_top(2), ierr)
                        end if

                        ! send and recv corner
                        if (this%negb_type(com_negb_left_bottom) == com_negb_type_domain) then  ! corner left bottom
                            call MPI_ISEND(this%send_buff_corner(1), 1, MPI_DOUBLE, this%negb_rank(com_negb_left_bottom), 1, MPI_COMM_WORLD, this%reqs_corner(1), ierr)
                            call MPI_IRECV(this%recv_buff_corner(1), 1, MPI_DOUBLE, this%negb_rank(com_negb_left_bottom), 1, MPI_COMM_WORLD, this%reqs_corner(2), ierr)
                        end if

                        if (this%negb_type(com_negb_right_bottom) == com_negb_type_domain) then ! corner right bottom
                            call MPI_ISEND(this%send_buff_corner(2), 1, MPI_DOUBLE, this%negb_rank(com_negb_right_bottom), 1, MPI_COMM_WORLD, this%reqs_corner(3), ierr)
                            call MPI_IRECV(this%recv_buff_corner(2), 1, MPI_DOUBLE, this%negb_rank(com_negb_right_bottom), 1, MPI_COMM_WORLD, this%reqs_corner(4), ierr)
                        end if

                        if (this%negb_type(com_negb_right_top) == com_negb_type_domain) then    ! corner right top
                            call MPI_ISEND(this%send_buff_corner(3), 1, MPI_DOUBLE, this%negb_rank(com_negb_right_top), 1, MPI_COMM_WORLD, this%reqs_corner(5), ierr)
                            call MPI_IRECV(this%recv_buff_corner(3), 1, MPI_DOUBLE, this%negb_rank(com_negb_right_top), 1, MPI_COMM_WORLD, this%reqs_corner(6), ierr)
                        end if

                        if (this%negb_type(com_negb_left_top) == com_negb_type_domain) then     ! corner left top
                            call MPI_ISEND(this%send_buff_corner(4), 1, MPI_DOUBLE, this%negb_rank(com_negb_left_top), 1, MPI_COMM_WORLD, this%reqs_corner(7), ierr)
                            call MPI_IRECV(this%recv_buff_corner(4), 1, MPI_DOUBLE, this%negb_rank(com_negb_left_top), 1, MPI_COMM_WORLD, this%reqs_corner(8), ierr)
                        end if

                        ! wait and accumulation
                        this%recv_count_corner = 0
                        if (this%negb_type(com_negb_left) == com_negb_type_domain) then         ! left
                            call MPI_WAITALL(2, this%reqs_left, this%status_left, ierr)
                            array2d(xstart, ystart:yend) = array2d(xstart, ystart:yend) + this%recv_buff_edge_left
                            this%recv_count_corner(1) = this%recv_count_corner(1) + 1
                            this%recv_count_corner(4) = this%recv_count_corner(4) + 1
                        end if

                        if (this%negb_type(com_negb_right) == com_negb_type_domain) then        ! right
                            call MPI_WAITALL(2, this%reqs_right, this%status_right, ierr)
                            array2d(xend, ystart:yend) = array2d(xend, ystart:yend) + this%recv_buff_edge_right
                            this%recv_count_corner(2) = this%recv_count_corner(2) + 1
                            this%recv_count_corner(3) = this%recv_count_corner(3) + 1
                        end if

                        if (this%negb_type(com_negb_bottom) == com_negb_type_domain) then       ! bottom
                            call MPI_WAITALL(2, this%reqs_bottom, this%status_bottom, ierr)
                            array2d(xstart:xend, ystart) = array2d(xstart:xend, ystart) + this%recv_buff_edge_bottom
                            this%recv_count_corner(1) = this%recv_count_corner(1) + 1
                            this%recv_count_corner(2) = this%recv_count_corner(2) + 1
                        end if

                        if (this%negb_type(com_negb_top) == com_negb_type_domain) then          ! top
                            call MPI_WAITALL(2, this%reqs_top, this%status_top, ierr)
                            array2d(xstart:xend, yend) = array2d(xstart:xend, yend) + this%recv_buff_edge_top
                            this%recv_count_corner(4) = this%recv_count_corner(4) + 1
                            this%recv_count_corner(3) = this%recv_count_corner(3) + 1
                        end if

                        if (this%negb_type(com_negb_left_bottom) == com_negb_type_domain) then  ! left bottom
                            call MPI_WAITALL(2, this%reqs_corner(1:2), this%status_corner(:, 1:2), ierr)
                            array2d(xstart, ystart) = array2d(xstart, ystart) + this%recv_buff_corner(1)
                            this%recv_count_corner(1) = this%recv_count_corner(1) + 1
                        end if
                        
                        if (this%negb_type(com_negb_right_bottom) == com_negb_type_domain) then ! right bottom
                            call MPI_WAITALL(2, this%reqs_corner(3:4), this%status_corner(:, 3:4), ierr)
                            array2d(xend, ystart) = array2d(xend, ystart) + this%recv_buff_corner(2)
                            this%recv_count_corner(2) = this%recv_count_corner(2) + 1
                        end if

                        if (this%negb_type(com_negb_right_top) == com_negb_type_domain) then    ! right top
                            call MPI_WAITALL(2, this%reqs_corner(5:6), this%status_corner(:, 5:6), ierr)
                            array2d(xend, yend) = array2d(xend, yend) + this%recv_buff_corner(3)
                            this%recv_count_corner(3) = this%recv_count_corner(3) + 1
                        end if

                        if (this%negb_type(com_negb_left_top) == com_negb_type_domain) then     ! left top
                            call MPI_WAITALL(2, this%reqs_corner(7:8), this%status_corner(:, 7:8), ierr)
                            array2d(xstart, yend) = array2d(xstart, yend) + this%recv_buff_corner(4)
                            this%recv_count_corner(4) = this%recv_count_corner(4) + 1
                        end if

                    else if (com_type_onesided == this%com_type) then
                        this%send_buff_edge_left    = array2d(xstart, ystart:yend)
                        this%send_buff_edge_right   = array2d(xend, ystart:yend)
                        this%send_buff_edge_bottom  = array2d(xstart:xend, ystart)
                        this%send_buff_edge_top     = array2d(xstart:xend, yend)

                        this%recv_buff_edge_left    = 0
                        this%recv_buff_edge_right   = 0
                        this%recv_buff_edge_bottom  = 0
                        this%recv_buff_edge_top     = 0
                        this%recv_buff_corner       = 0

                        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

                        ! send and recv edge
                        if (this%negb_type(com_negb_left) == com_negb_type_domain) then         ! left
                            this%disp_aint = 0
                            call MPI_WIN_LOCK(MPI_LOCK_SHARED, this%negb_rank(com_negb_left), 0, this%win_right, ierr)
                            call MPI_GET(this%recv_buff_edge_left, this%ly, MPI_DOUBLE, this%negb_rank(com_negb_left), this%disp_aint, this%ly, MPI_DOUBLE, this%win_right, ierr)
                            call MPI_WIN_UNLOCK(this%negb_rank(com_negb_left), this%win_right, ierr)

                            array2d(xstart, ystart:yend) = array2d(xstart, ystart:yend) + this%recv_buff_edge_left
                        end if

                        if (this%negb_type(com_negb_right) == com_negb_type_domain) then        ! right
                            this%disp_aint = 0
                            call MPI_WIN_LOCK(MPI_LOCK_SHARED, this%negb_rank(com_negb_right), 0, this%win_left, ierr)
                            call MPI_GET(this%recv_buff_edge_right, this%ly, MPI_DOUBLE, this%negb_rank(com_negb_right), this%disp_aint, this%ly, MPI_DOUBLE, this%win_left, ierr)
                            call MPI_WIN_UNLOCK(this%negb_rank(com_negb_right), this%win_left, ierr)

                            array2d(xend, ystart:yend) = array2d(xend, ystart:yend) + this%recv_buff_edge_right
                        end if

                        if (this%negb_type(com_negb_bottom) == com_negb_type_domain) then       ! bottom
                            this%disp_aint = 0
                            call MPI_WIN_LOCK(MPI_LOCK_SHARED, this%negb_rank(com_negb_bottom), 0, this%win_top, ierr)
                            call MPI_GET(this%recv_buff_edge_bottom, this%lx, MPI_DOUBLE, this%negb_rank(com_negb_bottom), this%disp_aint, this%lx, MPI_DOUBLE, this%win_top, ierr)
                            call MPI_WIN_UNLOCK(this%negb_rank(com_negb_bottom), this%win_top, ierr)

                            array2d(xstart:xend, ystart) = array2d(xstart:xend, ystart) + this%recv_buff_edge_bottom
                        end if

                        if (this%negb_type(com_negb_top) == com_negb_type_domain) then          ! top
                            this%disp_aint = 0
                            call MPI_WIN_LOCK(MPI_LOCK_SHARED, this%negb_rank(com_negb_top), 0, this%win_bottom, ierr)
                            call MPI_GET(this%recv_buff_edge_top, this%lx, MPI_DOUBLE, this%negb_rank(com_negb_top), this%disp_aint, this%lx, MPI_DOUBLE, this%win_bottom, ierr)
                            call MPI_WIN_UNLOCK(this%negb_rank(com_negb_top), this%win_bottom, ierr)

                            array2d(xstart:xend, yend) = array2d(xstart:xend, yend) + this%recv_buff_edge_top
                        end if

                        ! send and recv corner
                        if (this%negb_type(com_negb_left_bottom) == com_negb_type_domain) then  ! corner left bottom
                            this%disp_aint = this%lx-1
                            call MPI_WIN_LOCK(MPI_LOCK_SHARED, this%negb_rank(com_negb_left_bottom), 0, this%win_top, ierr)
                            call MPI_GET(this%recv_buff_corner(1), 1, MPI_DOUBLE, this%negb_rank(com_negb_left_bottom), this%disp_aint, 1, MPI_DOUBLE, this%win_top, ierr)
                            call MPI_WIN_UNLOCK(this%negb_rank(com_negb_left_bottom), this%win_top, ierr)

                            array2d(xstart, ystart) = array2d(xstart, ystart) + this%recv_buff_corner(1)
                        end if

                        if (this%negb_type(com_negb_right_bottom) == com_negb_type_domain) then ! corner right bottom
                            this%disp_aint = 0
                            call MPI_WIN_LOCK(MPI_LOCK_SHARED, this%negb_rank(com_negb_right_bottom), 0, this%win_top, ierr)
                            call MPI_GET(this%recv_buff_corner(2), 1, MPI_DOUBLE, this%negb_rank(com_negb_right_bottom), this%disp_aint, 1, MPI_DOUBLE, this%win_top, ierr)
                            call MPI_WIN_UNLOCK(this%negb_rank(com_negb_right_bottom), this%win_top, ierr)

                            array2d(xend, ystart) = array2d(xend, ystart) + this%recv_buff_corner(2)
                        end if

                        if (this%negb_type(com_negb_right_top) == com_negb_type_domain) then    ! corner right top
                            this%disp_aint = 0
                            call MPI_WIN_LOCK(MPI_LOCK_SHARED, this%negb_rank(com_negb_right_top), 0, this%win_bottom, ierr)
                            call MPI_GET(this%recv_buff_corner(3), 1, MPI_DOUBLE, this%negb_rank(com_negb_right_top), this%disp_aint, 1, MPI_DOUBLE, this%win_bottom, ierr)
                            call MPI_WIN_UNLOCK(this%negb_rank(com_negb_right_top), this%win_bottom, ierr)

                            array2d(xend, yend) = array2d(xend, yend) + this%recv_buff_corner(3)
                        end if

                        if (this%negb_type(com_negb_left_top) == com_negb_type_domain) then     ! corner left top
                            this%disp_aint = this%lx-1
                            call MPI_WIN_LOCK(MPI_LOCK_SHARED, this%negb_rank(com_negb_left_top), 0, this%win_bottom, ierr)
                            call MPI_GET(this%recv_buff_corner(4), 1, MPI_DOUBLE, this%negb_rank(com_negb_left_top), this%disp_aint, 1, MPI_DOUBLE, this%win_bottom, ierr)
                            call MPI_WIN_UNLOCK(this%negb_rank(com_negb_left_top), this%win_bottom, ierr)

                            array2d(xstart, yend) = array2d(xstart, yend) + this%recv_buff_corner(4)
                        end if
                        
                    end if

                else if (com_field_opt_ext == opt_type) then
                    if (com_type_nonblock == this%com_type) then
                        this%send_buff_edge_left    = array2d(xstart+1, ystart:yend)
                        this%send_buff_edge_right   = array2d(xend-1, ystart:yend)
                        this%send_buff_edge_bottom  = array2d(xstart:xend, ystart+1)
                        this%send_buff_edge_top     = array2d(xstart:xend, yend-1)

                        this%recv_buff_edge_left    = 0
                        this%recv_buff_edge_right   = 0
                        this%recv_buff_edge_bottom  = 0
                        this%recv_buff_edge_top     = 0

                        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

                        ! send and recv edge
                        if (this%negb_type(com_negb_left) == com_negb_type_domain) then         ! left
                            call MPI_ISEND(this%send_buff_edge_left, this%ly, MPI_DOUBLE, this%negb_rank(com_negb_left), 1, MPI_COMM_WORLD, this%reqs_left(1), ierr)
                            call MPI_IRECV(this%recv_buff_edge_left, this%ly, MPI_DOUBLE, this%negb_rank(com_negb_left), 1, MPI_COMM_WORLD, this%reqs_left(2), ierr)
                        end if

                        if (this%negb_type(com_negb_right) == com_negb_type_domain) then        ! right
                            call MPI_ISEND(this%send_buff_edge_right, this%ly, MPI_DOUBLE, this%negb_rank(com_negb_right), 1, MPI_COMM_WORLD, this%reqs_right(1), ierr)
                            call MPI_IRECV(this%recv_buff_edge_right, this%ly, MPI_DOUBLE, this%negb_rank(com_negb_right), 1, MPI_COMM_WORLD, this%reqs_right(2), ierr)
                        end if

                        if (this%negb_type(com_negb_bottom) == com_negb_type_domain) then       ! bottom
                            call MPI_ISEND(this%send_buff_edge_bottom, this%lx, MPI_DOUBLE, this%negb_rank(com_negb_bottom), 1, MPI_COMM_WORLD, this%reqs_bottom(1), ierr)
                            call MPI_IRECV(this%recv_buff_edge_bottom, this%lx, MPI_DOUBLE, this%negb_rank(com_negb_bottom), 1, MPI_COMM_WORLD, this%reqs_bottom(2), ierr)
                        end if

                        if (this%negb_type(com_negb_top) == com_negb_type_domain) then          ! top
                            call MPI_ISEND(this%send_buff_edge_top, this%lx, MPI_DOUBLE, this%negb_rank(com_negb_top), 1, MPI_COMM_WORLD, this%reqs_top(1), ierr)
                            call MPI_IRECV(this%recv_buff_edge_top, this%lx, MPI_DOUBLE, this%negb_rank(com_negb_top), 1, MPI_COMM_WORLD, this%reqs_top(2), ierr)
                        end if

                        ! wait and accumulation
                        if (this%negb_type(com_negb_left) == com_negb_type_domain) then         ! left
                            call MPI_WAITALL(2, this%reqs_left, this%status_left, ierr)
                            array2d(xstart-1, ystart:yend) = this%recv_buff_edge_left
                            array2d(xstart-1, ystart-1)    = 2.d0 * array2d(xstart-1, ystart) - array2d(xstart-1, ystart+1)
                            array2d(xstart-1, yend+1)      = 2.d0 * array2d(xstart-1, yend) - array2d(xstart-1, yend-1)

                        else if (this%negb_type(com_negb_left) == com_negb_type_boundary) then
                            if (this%left_ext_type == com_field_ext_type_intep) then
                                array2d(xstart-1, ystart:yend) = 2.d0 * array2d(xstart, ystart:yend) - array2d(xstart+1, ystart:yend)
                                array2d(xstart-1, ystart-1)    = 2.d0 * array2d(xstart-1, ystart) - array2d(xstart-1, ystart+1)
                                array2d(xstart-1, yend+1)      = 2.d0 * array2d(xstart-1, yend) - array2d(xstart-1, yend-1)
                            
                            else if (this%left_ext_type == com_field_ext_type_symmetry) then
                                array2d(xstart-1, ystart:yend) = array2d(xstart+1, ystart:yend)
                                array2d(xstart-1, ystart-1)    = 2.d0 * array2d(xstart-1, ystart) - array2d(xstart-1, ystart+1)
                                array2d(xstart-1, yend+1)      = 2.d0 * array2d(xstart-1, yend) - array2d(xstart-1, yend-1)

                            end if
                        end if

                        if (this%negb_type(com_negb_right) == com_negb_type_domain) then        ! right
                            call MPI_WAITALL(2, this%reqs_right, this%status_right, ierr)
                            array2d(xend+1, ystart:yend)   = this%recv_buff_edge_right
                            array2d(xend+1, ystart-1)      = 2.d0 * array2d(xend+1, ystart) - array2d(xend+1, ystart+1)
                            array2d(xend+1, yend+1)        = 2.d0 * array2d(xend+1, yend) - array2d(xend+1, yend-1)

                        else if (this%negb_type(com_negb_right) == com_negb_type_boundary) then
                            if (this%right_ext_type == com_field_ext_type_intep) then
                                array2d(xend+1, ystart:yend)   = 2.d0 * array2d(xend, ystart:yend) - array2d(xend-1, ystart:yend)
                                array2d(xend+1, ystart-1)      = 2.d0 * array2d(xend+1, ystart) - array2d(xend+1, ystart+1)
                                array2d(xend+1, yend+1)        = 2.d0 * array2d(xend+1, yend) - array2d(xend+1, yend-1)

                            else if (this%right_ext_type == com_field_ext_type_symmetry) then
                                array2d(xend+1, ystart:yend)   = array2d(xend-1, ystart:yend)
                                array2d(xend+1, ystart-1)      = 2.d0 * array2d(xend+1, ystart) - array2d(xend+1, ystart+1)
                                array2d(xend+1, yend+1)        = 2.d0 * array2d(xend+1, yend) - array2d(xend+1, yend-1)

                            end if
                        end if

                        if (this%negb_type(com_negb_bottom) == com_negb_type_domain) then       ! bottom
                            call MPI_WAITALL(2, this%reqs_bottom, this%status_bottom, ierr)
                            array2d(xstart:xend, ystart-1) = this%recv_buff_edge_bottom
                            array2d(xstart-1, ystart-1)    = 2.d0 * array2d(xstart, ystart-1) - array2d(xstart+1, ystart-1)
                            array2d(xend+1, ystart-1)      = 2.d0 * array2d(xend, ystart-1) - array2d(xend-1, ystart-1)
                        
                        else if (this%negb_type(com_negb_bottom) == com_negb_type_boundary) then
                            if (this%bottom_ext_type == com_field_ext_type_intep) then
                                array2d(xstart:xend, ystart-1) = 2.d0 * array2d(xstart:xend, ystart) - array2d(xstart:xend, ystart+1)
                                array2d(xstart-1, ystart-1)    = 2.d0 * array2d(xstart, ystart-1) - array2d(xstart+1, ystart-1)
                                array2d(xend+1, ystart-1)      = 2.d0 * array2d(xend, ystart-1) - array2d(xend-1, ystart-1)
                            else if (this%bottom_ext_type == com_field_ext_type_symmetry) then
                                array2d(xstart:xend, ystart-1) = array2d(xstart:xend, ystart+1)
                                array2d(xstart-1, ystart-1)    = 2.d0 * array2d(xstart, ystart-1) - array2d(xstart+1, ystart-1)
                                array2d(xend+1, ystart-1)      = 2.d0 * array2d(xend, ystart-1) - array2d(xend-1, ystart-1)

                            end if
                        end if

                        if (this%negb_type(com_negb_top) == com_negb_type_domain) then          ! top
                            call MPI_WAITALL(2, this%reqs_top, this%status_top, ierr)
                            array2d(xstart:xend, yend+1)   = this%recv_buff_edge_top
                            array2d(xstart-1, yend+1)      = 2.d0 * array2d(xstart, yend+1) - array2d(xstart+1, yend+1)
                            array2d(xend+1, yend+1)        = 2.d0 * array2d(xend, yend+1) - array2d(xend-1, yend+1)
                        
                        else if (this%negb_type(com_negb_top) == com_negb_type_boundary) then
                            if (this%top_ext_type == com_field_ext_type_intep) then
                                array2d(xstart:xend, yend+1)   = 2.d0 * array2d(xstart:xend, yend) - array2d(xstart:xend, yend-1)
                                array2d(xstart-1, yend+1)      = 2.d0 * array2d(xstart, yend+1) - array2d(xstart+1, yend+1)
                                array2d(xend+1, yend+1)        = 2.d0 * array2d(xend, yend+1) - array2d(xend-1, yend+1)

                            else if (this%top_ext_type == com_field_ext_type_symmetry) then
                                array2d(xstart:xend, yend+1)   = array2d(xstart:xend, yend-1)
                                array2d(xstart-1, yend+1)      = 2.d0 * array2d(xstart, yend+1) - array2d(xstart+1, yend+1)
                                array2d(xend+1, yend+1)        = 2.d0 * array2d(xend, yend+1) - array2d(xend-1, yend+1)

                            end if
                        end if

                    else if (com_type_onesided == this%com_type) then
                        this%send_buff_edge_left    = array2d(xstart+1, ystart:yend)
                        this%send_buff_edge_right   = array2d(xend-1, ystart:yend)
                        this%send_buff_edge_bottom  = array2d(xstart:xend, ystart+1)
                        this%send_buff_edge_top     = array2d(xstart:xend, yend-1)

                        this%recv_buff_edge_left    = 0
                        this%recv_buff_edge_right   = 0
                        this%recv_buff_edge_bottom  = 0
                        this%recv_buff_edge_top     = 0

                        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

                        ! send and recv edge
                        if (this%negb_type(com_negb_left) == com_negb_type_domain) then         ! left
                            this%disp_aint = 0
                            call MPI_WIN_LOCK(MPI_LOCK_SHARED, this%negb_rank(com_negb_left), 0, this%win_right, ierr)
                            call MPI_GET(this%recv_buff_edge_left, this%ly, MPI_DOUBLE, this%negb_rank(com_negb_left), this%disp_aint, this%ly, MPI_DOUBLE, this%win_right, ierr)
                            call MPI_WIN_UNLOCK(this%negb_rank(com_negb_left), this%win_right, ierr)

                            array2d(xstart-1, ystart:yend) = this%recv_buff_edge_left
                            array2d(xstart-1, ystart-1)    = 2.d0 * array2d(xstart-1, ystart) - array2d(xstart-1, ystart+1)
                            array2d(xstart-1, yend+1)      = 2.d0 * array2d(xstart-1, yend) - array2d(xstart-1, yend-1)

                        else if (this%negb_type(com_negb_left) == com_negb_type_boundary) then
                            if (this%left_ext_type == com_field_ext_type_intep) then
                                array2d(xstart-1, ystart:yend) = 2.d0 * array2d(xstart, ystart:yend) - array2d(xstart+1, ystart:yend)
                                array2d(xstart-1, ystart-1)    = 2.d0 * array2d(xstart-1, ystart) - array2d(xstart-1, ystart+1)
                                array2d(xstart-1, yend+1)      = 2.d0 * array2d(xstart-1, yend) - array2d(xstart-1, yend-1)
                            
                            else if (this%left_ext_type == com_field_ext_type_symmetry) then
                                array2d(xstart-1, ystart:yend) = array2d(xstart+1, ystart:yend)
                                array2d(xstart-1, ystart-1)    = 2.d0 * array2d(xstart-1, ystart) - array2d(xstart-1, ystart+1)
                                array2d(xstart-1, yend+1)      = 2.d0 * array2d(xstart-1, yend) - array2d(xstart-1, yend-1)

                            end if
                        end if

                        if (this%negb_type(com_negb_right) == com_negb_type_domain) then        ! right
                            this%disp_aint = 0
                            call MPI_WIN_LOCK(MPI_LOCK_SHARED, this%negb_rank(com_negb_right), 0, this%win_left, ierr)
                            call MPI_GET(this%recv_buff_edge_right, this%ly, MPI_DOUBLE, this%negb_rank(com_negb_right), this%disp_aint, this%ly, MPI_DOUBLE, this%win_left, ierr)
                            call MPI_WIN_UNLOCK(this%negb_rank(com_negb_right), this%win_left, ierr)

                            array2d(xend+1, ystart:yend)   = this%recv_buff_edge_right
                            array2d(xend+1, ystart-1)      = 2.d0 * array2d(xend+1, ystart) - array2d(xend+1, ystart+1)
                            array2d(xend+1, yend+1)        = 2.d0 * array2d(xend+1, yend) - array2d(xend+1, yend-1)

                        else if (this%negb_type(com_negb_right) == com_negb_type_boundary) then
                            if (this%right_ext_type == com_field_ext_type_intep) then
                                array2d(xend+1, ystart:yend)   = 2.d0 * array2d(xend, ystart:yend) - array2d(xend-1, ystart:yend)
                                array2d(xend+1, ystart-1)      = 2.d0 * array2d(xend+1, ystart) - array2d(xend+1, ystart+1)
                                array2d(xend+1, yend+1)        = 2.d0 * array2d(xend+1, yend) - array2d(xend+1, yend-1)

                            else if (this%right_ext_type == com_field_ext_type_symmetry) then
                                array2d(xend+1, ystart:yend)   = array2d(xend-1, ystart:yend)
                                array2d(xend+1, ystart-1)      = 2.d0 * array2d(xend+1, ystart) - array2d(xend+1, ystart+1)
                                array2d(xend+1, yend+1)        = 2.d0 * array2d(xend+1, yend) - array2d(xend+1, yend-1)

                            end if
                        end if

                        if (this%negb_type(com_negb_bottom) == com_negb_type_domain) then       ! bottom
                            this%disp_aint = 0
                            call MPI_WIN_LOCK(MPI_LOCK_SHARED, this%negb_rank(com_negb_bottom), 0, this%win_top, ierr)
                            call MPI_GET(this%recv_buff_edge_bottom, this%lx, MPI_DOUBLE, this%negb_rank(com_negb_bottom), this%disp_aint, this%lx, MPI_DOUBLE, this%win_top, ierr)
                            call MPI_WIN_UNLOCK(this%negb_rank(com_negb_bottom), this%win_top, ierr)

                            array2d(xstart:xend, ystart-1) = this%recv_buff_edge_bottom
                            array2d(xstart-1, ystart-1)    = 2.d0 * array2d(xstart, ystart-1) - array2d(xstart+1, ystart-1)
                            array2d(xend+1, ystart-1)      = 2.d0 * array2d(xend, ystart-1) - array2d(xend-1, ystart-1)
                        
                        else if (this%negb_type(com_negb_bottom) == com_negb_type_boundary) then
                            if (this%bottom_ext_type == com_field_ext_type_intep) then
                                array2d(xstart:xend, ystart-1) = 2.d0 * array2d(xstart:xend, ystart) - array2d(xstart:xend, ystart+1)
                                array2d(xstart-1, ystart-1)    = 2.d0 * array2d(xstart, ystart-1) - array2d(xstart+1, ystart-1)
                                array2d(xend+1, ystart-1)      = 2.d0 * array2d(xend, ystart-1) - array2d(xend-1, ystart-1)
                            else if (this%bottom_ext_type == com_field_ext_type_symmetry) then
                                array2d(xstart:xend, ystart-1) = array2d(xstart:xend, ystart+1)
                                array2d(xstart-1, ystart-1)    = 2.d0 * array2d(xstart, ystart-1) - array2d(xstart+1, ystart-1)
                                array2d(xend+1, ystart-1)      = 2.d0 * array2d(xend, ystart-1) - array2d(xend-1, ystart-1)

                            end if
                        end if

                        if (this%negb_type(com_negb_top) == com_negb_type_domain) then          ! top
                            this%disp_aint = 0
                            call MPI_WIN_LOCK(MPI_LOCK_SHARED, this%negb_rank(com_negb_top), 0, this%win_bottom, ierr)
                            call MPI_GET(this%recv_buff_edge_top, this%lx, MPI_DOUBLE, this%negb_rank(com_negb_top), this%disp_aint, this%lx, MPI_DOUBLE, this%win_bottom, ierr)
                            call MPI_WIN_UNLOCK(this%negb_rank(com_negb_top), this%win_bottom, ierr)

                            array2d(xstart:xend, yend+1)   = this%recv_buff_edge_top
                            array2d(xstart-1, yend+1)      = 2.d0 * array2d(xstart, yend+1) - array2d(xstart+1, yend+1)
                            array2d(xend+1, yend+1)        = 2.d0 * array2d(xend, yend+1) - array2d(xend-1, yend+1)
                        
                        else if (this%negb_type(com_negb_top) == com_negb_type_boundary) then
                            if (this%top_ext_type == com_field_ext_type_intep) then
                                array2d(xstart:xend, yend+1)   = 2.d0 * array2d(xstart:xend, yend) - array2d(xstart:xend, yend-1)
                                array2d(xstart-1, yend+1)      = 2.d0 * array2d(xstart, yend+1) - array2d(xstart+1, yend+1)
                                array2d(xend+1, yend+1)        = 2.d0 * array2d(xend, yend+1) - array2d(xend-1, yend+1)

                            else if (this%top_ext_type == com_field_ext_type_symmetry) then
                                array2d(xstart:xend, yend+1)   = array2d(xstart:xend, yend-1)
                                array2d(xstart-1, yend+1)      = 2.d0 * array2d(xstart, yend+1) - array2d(xstart+1, yend+1)
                                array2d(xend+1, yend+1)        = 2.d0 * array2d(xend, yend+1) - array2d(xend-1, yend+1)

                            end if
                        end if

                    end if

                else
                    if (0 == this%rank) then
                        write(*, *) "The input parameter of the opt_type is invalid."
                    end if
                    stop
                end if

            else
                write(*, *) "The FieldCom2D is not initialized."
                stop
            end if

        end subroutine communicationPICFieldCom2D


        subroutine communicationPICParticleCom2D(this, pb, x_lb, x_ub, y_lb, y_ub)
            class(PICCom2D), intent(inout) :: this
            type(ParticleBundle), intent(inout) :: pb
            real(8), intent(in) :: x_lb, x_ub, y_lb, y_ub
            integer(4) :: ierr, index, i
            integer(4) :: particle_number_send(com_negb_num)
            integer(4) :: particle_number_recv(com_negb_num)
            type(ParticleBundle) :: particle_send_buff(com_negb_num)
            type(ParticleBundle) :: particle_recv_buff(com_negb_num)

            ! 处理超过domain的粒子
            particle_number_send = 0
            particle_number_recv = 0

            if (pb%NPar > 0) then
                do i = 1, pb%NPar
                    index = getParticleDomainIndex(pb%PO(i), x_lb, x_ub, y_lb, y_ub)
                    if (index > 0) particle_number_send(index) = particle_number_send(index) + 1
                end do
            end if

            do i = 1, com_negb_num
                if (this%negb_type(i) == com_negb_type_domain) then
                    call MPI_ISEND(particle_number_send(i), 1, &
                                   MPI_INTEGER, this%negb_rank(i), &
                                   1, MPI_COMM_WORLD, this%reqs_send(i), ierr)
                    call MPI_IRECV(particle_number_recv(i), 1, &
                                   MPI_INTEGER, this%negb_rank(i), 1, &
                                   MPI_COMM_WORLD, this%reqs_recv(i), ierr)
                end if
            end do

            do i = 1, com_negb_num
                if (this%negb_type(i) == com_negb_type_domain) then
                    call MPI_WAIT(this%reqs_send(i), this%status_send(:, i), ierr)
                    call MPI_WAIT(this%reqs_recv(i), this%status_recv(:, i), ierr)
                end if
            end do

            ! send and recv buffer
            do i = 1, com_negb_num
                call particle_send_buff(i)%init(particle_number_send(i))
                call particle_recv_buff(i)%init(particle_number_recv(i))
            end do

            do i = pb%NPar, 1, -1
                index = getParticleDomainIndex(pb%PO(i), x_lb, x_ub, y_lb, y_ub)

                if (index > 0) then
                    call particle_send_buff(index)%addone(pb%PO(i))
                    call pb%delone(i)
                end if
            end do

            ! send and recv
            do i = 1, com_negb_num
                if (this%negb_type(i) == com_negb_type_domain) then
                    if (particle_number_send(i) > 0) &
                    call MPI_ISEND(particle_send_buff(i)%PO, particle_number_send(i), &
                                   this%mpi_type_particle_one, this%negb_rank(i), &
                                   1, MPI_COMM_WORLD, this%reqs_send(i), ierr)

                    if (particle_number_recv(i) > 0) &
                    call MPI_IRECV(particle_recv_buff(i)%PO, particle_number_recv(i), &
                                   this%mpi_type_particle_one, this%negb_rank(i), 1, &
                                   MPI_COMM_WORLD, this%reqs_recv(i), ierr)

                end if
            end do

            ! wait
            do i = 1, com_negb_num
                if (this%negb_type(i) == com_negb_type_domain .and. particle_number_send(i) > 0) then
                    call MPI_WAIT(this%reqs_send(i), this%status_send(:, i), ierr)
                end if

                if (this%negb_type(i) == com_negb_type_domain .and. particle_number_recv(i) > 0) then
                    call MPI_WAIT(this%reqs_recv(i), this%status_recv(:, i), ierr)
                end if
            end do

            do i = 1, com_negb_num
                if (this%negb_type(i) == com_negb_type_domain .and. particle_number_recv(i) > 0) then
                    particle_recv_buff(i)%NPar = particle_number_recv(i)
                    call pb%addbun(particle_recv_buff(i))
                end if
            end do

            do i = 1, com_negb_num
                call particle_send_buff(i)%destroy()
                call particle_recv_buff(i)%destroy()
            end do

            do i = pb%NPar, 1, -1
                index = getParticleDomainIndex(pb%PO(i), x_lb, x_ub, y_lb, y_ub)

                if (index > 0) then
                    call pb%delone(i)
                end if
            end do

        end subroutine communicationPICParticleCom2D


        subroutine gatherPICFieldCom2D(this, array2d_local, array2d_global, xstart, xend, ystart, yend, x_global_size, y_global_size)
            class(PICCom2D), intent(inout) :: this
            real(8), allocatable, intent(inout) :: array2d_local(:, :)
            real(8), allocatable, intent(inout) :: array2d_global(:, :)
            integer(4), intent(in) :: xstart, xend, ystart, yend, x_global_size, y_global_size
            real(8) :: array2d_tmp(x_global_size, y_global_size)
            integer(4) :: ierr

            if (this%is_init) then
                array2d_tmp = 0
                array2d_tmp(xstart:xend, ystart:yend) = array2d_local(xstart:xend, ystart:yend)
                
                call MPI_ALLREDUCE(array2d_tmp, array2d_global, this%lx*this%ly*this%size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)

            else
                write(*, *) "The FieldCom2D is not initialized."
                stop
            end if

        end subroutine gatherPICFieldCom2D


        function getParticleDomainIndex(one, x_lb, x_ub, y_lb, y_ub)
            type(ParticleOne), intent(in) :: one
            real(8), intent(in) :: x_lb, x_ub, y_lb, y_ub
            integer(4) :: getParticleDomainIndex

            getParticleDomainIndex = 0
            if (one%Z <= x_lb) then
                if (one%R <= y_lb) then
                    getParticleDomainIndex = com_negb_left_bottom
                else if (one%R >= y_ub) then
                    getParticleDomainIndex = com_negb_left_top
                else
                    getParticleDomainIndex = com_negb_left
                end if

            else if (one%Z >= x_ub) then
                if (one%R <= y_lb) then
                    getParticleDomainIndex = com_negb_right_bottom
                else if (one%R >= y_ub) then
                    getParticleDomainIndex = com_negb_right_top
                else
                    getParticleDomainIndex = com_negb_right
                end if

            else
                if (one%R <= y_lb) then
                    getParticleDomainIndex = com_negb_bottom
                else if (one%R >= y_ub) then
                    getParticleDomainIndex = com_negb_top
                else
                    getParticleDomainIndex = 0
                end if
            end if
        end

end module ModulePICCommunication
