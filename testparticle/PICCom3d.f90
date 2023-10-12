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

    type :: PICCom2D
        integer(4) :: lx, ly, lz
        integer(4) :: xyz_np(3)
        integer(4) :: size, rank
        integer(4) :: com_type

        logical :: is_init = .False.
        integer(4) :: col, row

        

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

        subroutine initPICCom2D(this, x_local_size, y_local_size, z_local_size, xyz_np(3), com_type)
            class(PICCom2D), intent(inout) :: this
            integer(4), intent(in) :: x_local_size, y_local_size,z_local_size, xyz_np(3)
            integer(4), optional, intent(in) :: com_type
            integer(4) :: ierr, i

            call MPI_COMM_SIZE(MPI_COMM_WORLD, this%size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, this%rank, ierr)

            if (x_local_size > 0 .and. y_local_size > 0.and. z_local_size > 0.and. deps_x > 0 .and. deps_y > 0.and. deps_z > 0) then
                this%lx = x_local_size
                this%ly = y_local_size
                this%lz = z_local_size
                this%xyz_np = xyz_np
        
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

                call this%destroy()
                

                ! 创建新的数据类型
                call MPI_TYPE_CONTIGUOUS(11, MPI_DOUBLE, this%mpi_type_particle_one, ierr)
                call MPI_TYPE_COMMIT(this%mpi_type_particle_one, ierr)

                this%is_init = .True.
            end if
        end subroutine initPICCom2D


        subroutine destroyPICCom2D(this)
            class(PICCom2D), intent(inout) :: this
            integer(4) :: ierr, i

            if (com_type_onesided == this%com_type .and. this%is_init) then
                call MPI_WIN_FREE(this%win_left, ierr)
                call MPI_WIN_FREE(this%win_right, ierr)
                call MPI_WIN_FREE(this%win_top, ierr)
                call MPI_WIN_FREE(this%win_bottom, ierr)
            end if

            if (this%is_init) call MPI_TYPE_FREE(this%mpi_type_particle_one, ierr)

            this%is_init = .False.

        end subroutine destroyPICCom2D

    

        subroutine communicationPICParticleCom3D(this, pb, x_lb, x_ub, y_lb, y_ub,z_lb,z_ub)
            class(PICCom2D), intent(inout) :: this
            type(ParticleBundle), intent(inout) :: pb
            real(8), intent(in) :: x_lb, x_ub, y_lb, y_ub,z_lb,z_ub
            integer(4) :: ierr, index, i,Indx(6,pb%NPar),l(6),j
            integer(4) :: particle_number_send(com_negb_num)
            integer(4) :: particle_number_recv(com_negb_num),negb_rank(2),boundarytest,temp1,temp2
            type(ParticleBundle) :: particle_send_buff(com_negb_num)
            type(ParticleBundle) :: particle_recv_buff(com_negb_num)

            ! 处理超过domain的粒子
            particle_number_send = 0
            particle_number_recv = 0

            if (pb%NPar > 0) then
                k=0,l=0,m=0
                do i = 1, pb%NPar
                    index = getParticleDomainRank(this,pb%PO(i), x_lb, x_ub, y_lb, y_ub,z_lb,z_ub)    
                    Indx(index,l(index))=i
                    l(index)=l(index)+1
                end do
            end if
            do dim=1,3,-1
                select case (dim)
                case (1)
                    orrd=one%X
                    negb_rank(1)=this%rank-1
                    negb_rank(2)=this%rank+1
                    boundarytest=1
                case (2)
                    negb_rank(1)=this%rank-this%xyz_np(1)
                    negb_rank(2)=this%rank+this%xyz_np(1)
                    boundarytest=this%xyz_np(1)
                case (3)
                    negb_rank(1)=this%rank-this%dpx*this%xyz_np(1)*xyz_np(2)
                    negb_rank(2)=this%rank+this%dpx*this%xyz_np(1)*xyz_np(2)
                    boundarytest=this%xyz_np(1)*xyz_np(2)
                case default
                end select
                    temp1=mod(this%rank/boundarytest,xyz_np(dim)) 
                do i = 1, 2       
                    temp2=mod(neighbor(dim)/boundarytest,xyz_np(dim))
                    if (abs(temp2-temp1)==1) then!检测邻居节点是否存在
                        call MPI_ISEND(particle_number_send(i), 1, &
                                   MPI_INTEGER, negb_rank(i), &
                                   1, MPI_COMM_WORLD, this%reqs_send(i), ierr)
                        call MPI_IRECV(particle_number_recv(i), 1, &
                                   MPI_INTEGER, negb_rank(i), 1, &
                                   MPI_COMM_WORLD, this%reqs_recv(i), ierr)
                    end if
                end do

                do i = dim, dim+1
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

                do i = k, 1, -1
                    !index = getParticleDomainIndex(pb%PO(i), x_lb, x_ub, y_lb, y_ub)

                ! if (index > 0) then
                    call particle_send_buff(Indx(k))%addone(pb%PO(NPar(k)))
                   ! call pb%delone(NPar(k))
                ! end if
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
                        if(particle_recv_buff(i)%NPar>0)
                            !particle_recv_buff(i)%NPar = particle_number_recv(i)
                            do j=pb%NPar,pb%NPar+particle_recv_buff(i)%NPar
                                index = getParticleDomainRank(this,particle_recv_buff(i)%PO(j), x_lb, x_ub, y_lb, y_ub,z_lb,z_ub)    
                                Indx(index,l(index))=i
                                l(index)=l(index)+1
                            end do  particle_recv_buff(i)%PO(j)
                            call pb%addbun(particle_recv_buff(i))
                        end if

                    end if
                end do

                do i = 1, com_negb_num
                    call particle_send_buff(i)%destroy()
                    call particle_recv_buff(i)%destroy()
                end do
            end do
            do i=1,6
                Indx(index,l(index))=i
           ! do i = pb%NPar, 1, -1
           !     index = getParticleDomainIndex(pb%PO(i), x_lb, x_ub, y_lb, y_ub)

           !     if (index > 0) then
           !         call pb%delone(i)
           !     end if
           ! end do

        end subroutine communicationPICParticleCom2D



        function getParticleDomainRank(this,one, x_lb, x_ub, y_lb, y_ub,z_lb,z_ub)
            type(ParticleOne), intent(in) :: one
            class(PICCom2D), intent(inout) :: this
            real(8), intent(in) :: x_lb, x_ub, y_lb, y_ub,z_lb,z_ub
            integer(4) :: getParticleDomainIndex

            getParticleDomainIndex = 0
            if (one%X <= x_lb) then
                    getParticleDomainIndex = this%rank-1
                if (one%X >= x_lb) then
                    getParticleDomainIndex = this%rank+1
                    if (one%Y <= y_lb) then
                        getParticleDomainIndex = this%rank-this%xyz_np(1)
                        if (one%Y >= y_lb) then
                            getParticleDomainIndex = this%rank+this%xyz_np(1)
                            if (one%Z <= z_lb) then
                                getParticleDomainIndex = this%rank-this%dpx*this%xyz_np(1)*xyz_np(2)
                                if (one%Z >= z_lb) then
                                    getParticleDomainIndex =this%rank+this%dpx*this%xyz_np(1)*xyz_np(2)
                                else
                                    getParticleDomainIndex =this%rank
                                end if
                            end if
                        end if
                    end if            
                end if
            end if       
        end
        function get_ParticleDomainIndex(one, l_b, u_b,dim)
            type(ParticleOne), intent(in) :: one
            real(8), intent(in) :: l_b,u_b
            integer(4) :: get_ParticleDomainIndex,dim
            real(8)::orrd
            select case (dim)
            case (1)
                orrd=one%X
            case (2)
                orrd=one%Y
            case (3)
                orrd=one%Z
            case default
            end select
            
            get_ParticleDomainIndex = 0
            select case (orrd)
            case (u_b:)
                get_ParticleDomainIndex = com_negb_right
            case (:l_b)
                get_ParticleDomainIndex = com_negb_left
            case default
                get_ParticleDomainIndex = 0  
            end select
            
        end

end module ModulePICCommunication
