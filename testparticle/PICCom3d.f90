module ModulePICCommunication
    use mpi
    use ModuleParticleBundle
    implicit none

    integer com_negb_num =2 !分别是上下 左右 前后 两个方向，分三次通信

    type :: PICCom3D
        integer(4) :: lx, ly, lz
        integer(4) :: xyz_np(3)
        integer(4) :: size, rank
        logical :: is_init = .False.
       

        ! for non-block particle com
        integer(4) :: reqs_send(com_negb_num), reqs_recv(com_negb_num)
        integer(4) :: status_send(MPI_STATUS_SIZE, com_negb_num), status_recv(MPI_STATUS_SIZE, com_negb_num)

        ! particleone datatype
        integer(4) :: mpi_type_particle_one

    contains

        procedure :: init       => initPICCom3D
        procedure :: destroy    => destroyPICCom3D
        procedure :: comp       => communicationPICParticleCom3D

    end type PICCom3D


    contains

        subroutine initPICCom3D(this, x_local_size, y_local_size, z_local_size, xyz_np(3), com_type)
            class(PICCom2D), intent(inout) :: this
            integer(4), intent(in) :: x_local_size, y_local_size,z_local_size, xyz_np(3)
            integer(4), optional, intent(in) :: com_type
            integer(4) :: ierr, i

            call MPI_COMM_SIZE(MPI_COMM_WORLD, this%size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, this%rank, ierr)

            if (x_local_size > 0 .and. y_local_size > 0.and. z_local_size > 0.and. xyz_np(1) > 0 .and. xyz_np(2) > 0.and. xyz_np(3) > 0) then
                this%lx = x_local_size
                this%ly = y_local_size
                this%lz = z_local_size
                this%xyz_np = xyz_np
        

                ! 进程数分解到三个维度，必须恰好等于总进程数
                if (this%xyz_np(1) * this%xyz_np(2)*this%xyz_np(3) /= this%size) then
                    write(*, *) "The domian decomposition is error."
                    stop
                end if
                ! 创建新的数据类型
                call MPI_TYPE_CONTIGUOUS(11, MPI_DOUBLE, this%mpi_type_particle_one, ierr)
                call MPI_TYPE_COMMIT(this%mpi_type_particle_one, ierr)

                this%is_init = .True.
            end if
        end subroutine initPICCom3D


        subroutine destroyPICCom3D(this)
            class(PICCom2D), intent(inout) :: this
            integer(4) :: ierr, i

            

        end subroutine destroyPICCom3D

    

        subroutine communicationPICParticleCom3D(this, pb, x_lb, x_ub, y_lb, y_ub,z_lb,z_ub)
            class(PICCom2D), intent(inout) :: this
            type(ParticleBundle), intent(inout) :: pb
            real(8), intent(in) :: x_lb, x_ub, y_lb, y_ub,z_lb,z_ub
            integer(4) :: ierr, index, i,Indx(6,pb%NPar+5000),l(6),j,test(com_negb_num)
            integer(4) :: particle_number_send(com_negb_num)
            integer(4) :: particle_number_recv(com_negb_num),negb_rank(com_negb_num),boundarytest,temp1,temp2
            type(ParticleBundle) :: particle_send_buff(com_negb_num)
            type(ParticleBundle) :: particle_recv_buff(com_negb_num)


            ! 处理超过domain的粒子
            particle_number_send = 0
            particle_number_recv = 0
            test=0

            if (pb%NPar > 0) then
                do i = 1, pb%NPar
                    index = getParticleDomainIndex(this,pb%PO(i), x_lb, x_ub, y_lb, y_ub,z_lb,z_ub)    
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
                        test(i)=1!标记检测结果
                        call MPI_ISEND(particle_number_send(i), 1, &
                                   MPI_INTEGER, negb_rank(i), &
                                   1, MPI_COMM_WORLD, this%reqs_send(i), ierr)
                        call MPI_IRECV(particle_number_recv(i), 1, &
                                   MPI_INTEGER, negb_rank(i), 1, &
                                   MPI_COMM_WORLD, this%reqs_recv(i), ierr)
                    end if
                end do

                do i = 1, 2
                    if (test(i)==1) then!检测邻居节点是否存在
                        call MPI_WAIT(this%reqs_send(i), this%status_send(:, i), ierr)
                        call MPI_WAIT(this%reqs_recv(i), this%status_recv(:, i), ierr)
                    end if
                end do

                ! send and recv buffer
                do i = 1, 2
                    call particle_send_buff(i)%init(particle_number_send(i))
                    call particle_recv_buff(i)%init(particle_number_recv(i))
                end do

                do i = 1, 2
                    
                    do j=1,l(i+dim*2-2)
                
                    call particle_send_buff(i)%addone(pb%PO(Indx(i+dim*2-2,j)))
                    end do
                ! end if
                end do

                ! send and recv
                do i = 1, 2    
                    if (test(i)==1) then!检测邻居节点是否存在
                        if (particle_number_send(i) > 0) &
                        call MPI_ISEND(particle_send_buff(i)%PO, particle_number_send(i), &
                                   this%mpi_type_particle_one, negb_rank(i), &
                                   1, MPI_COMM_WORLD, this%reqs_send(i), ierr)

                    if (particle_number_recv(i) > 0) &
                    call MPI_IRECV(particle_recv_buff(i)%PO, particle_number_recv(i), &
                                   this%mpi_type_particle_one, negb_rank(i), 1, &
                                   MPI_COMM_WORLD, this%reqs_recv(i), ierr)
                     end if
                end do

                ! wait
                do i = 1, 2
                    if (test(i)==1 .and. particle_number_send(i) > 0) then
                        call MPI_WAIT(this%reqs_send(i), this%status_send(:, i), ierr)
                    end if

                    if (test(i)==1 .and. particle_number_recv(i) > 0) then
                        call MPI_WAIT(this%reqs_recv(i), this%status_recv(:, i), ierr)
                    end if
                end do

                do i = 1, com_negb_num
                    if (test(i)==1 .and. particle_number_recv(i) > 0) then
                        if(particle_recv_buff(i)%NPar>0)
                            do j=pb%NPar,pb%NPar+particle_recv_buff(i)%NPar
                                index = getParticleDomainRank(this,particle_recv_buff(i)%PO(j-pb%NPar+1), x_lb, x_ub, y_lb, y_ub,z_lb,z_ub)    
                                if(index>0)  then
                                Indx(index,l(index))=j 
                                l(index)=l(index)+1
                                end if
                            end do  
                            call pb%addbun(particle_recv_buff(i))
                        end if

                    end if
                end do

                do i = 1, 2
                    call particle_send_buff(i)%destroy()
                    call particle_recv_buff(i)%destroy()
                end do
            end do
            do i=1,6
                do j=l(i),1
                pb%DelOne(Indx(i,j))
                end do
            end do    

        end subroutine communicationPICParticleCom2D



        function getParticleDomainIndex(this,one, x_lb, x_ub, y_lb, y_ub,z_lb,z_ub)
            type(ParticleOne), intent(in) :: one
            class(PICCom2D), intent(inout) :: this
            real(8), intent(in) :: x_lb, x_ub, y_lb, y_ub,z_lb,z_ub
            integer(4) :: getParticleDomainIndex

            getParticleDomainIndex = 0
            if (one%X <= x_lb) then
                    getParticleDomainIndex = 1
                if (one%X >= x_lb) then
                    getParticleDomainIndex = 2
                    if (one%Y <= y_lb) then
                        getParticleDomainIndex = 3
                        if (one%Y >= y_lb) then
                            getParticleDomainIndex = 4
                            if (one%Z <= z_lb) then
                                getParticleDomainIndex = 5
                                if (one%Z >= z_lb) then
                                    getParticleDomainIndex =6
                                else
                                    getParticleDomainIndex =0
                                end if
                            end if
                        end if
                    end if            
                end if
            end if       
        end

end module ModulePICCommunication
