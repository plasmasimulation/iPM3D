module ModulePICParticleCom3d
    use mpi
    use ModuleParticleBundle
    implicit none

    integer(4) , parameter :: com_negb_num =2 !分别是上下 左右 前后 两个方向，分三次通信
    integer(4), parameter :: com_type_nonblock      = 0

    type :: PICCom3D
        integer(4) :: lx, ly, lz
        integer(4) :: xyz_np(3)
        integer(4) :: size, rank
        integer(4) :: col, row,layer
        logical :: is_init = .False.
       

        ! for non-block particle com
        integer(4) :: reqs_send(com_negb_num), reqs_recv(com_negb_num)
        integer(4) :: status_send(MPI_STATUS_SIZE, com_negb_num), status_recv(MPI_STATUS_SIZE, com_negb_num)

        ! particleone datatype
        integer(4) :: mpi_type_particle_one

    contains

        procedure :: init       => initPICCom3D
       ! procedure :: destroy    => destroyPICCom3D
        procedure :: comp       => communicationPICParticleCom3D

    end type PICCom3D


    contains

        subroutine initPICCom3D(this, x_local_size, y_local_size, z_local_size, xyz_np, com_type)
            class(PICCom3D), intent(inout) :: this
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
        
                this%col = mod(mod(this%rank, this%xyz_np(1)*this%xyz_np(2)), this%xyz_np(1))
                this%row = mod(this%rank, this%xyz_np(1)*this%xyz_np(2))/ this%xyz_np(1) 
                this%layer=this%rank/(this%xyz_np(1)*this%xyz_np(2))

                ! 进程数分解到三个维度，必须恰好等于总进程数
                if (this%xyz_np(1) * this%xyz_np(2)*this%xyz_np(3) /= this%size) then
                    write(*, *) "The domian decomposition is error."
                    stop
                end if
                ! 创建新的数据类型
                call MPI_TYPE_CONTIGUOUS(12, MPI_DOUBLE, this%mpi_type_particle_one, ierr)
                call MPI_TYPE_COMMIT(this%mpi_type_particle_one, ierr)

                this%is_init = .True.
            end if
        end subroutine initPICCom3D


        subroutine destroyPICCom3D(this)
            class(PICCom3D), intent(inout) :: this
            integer(4) :: ierr, i

            

        end subroutine destroyPICCom3D

    

        subroutine communicationPICParticleCom3D(this, pb,  xstart, xend, ystart, yend,zstart,zend)
            class(PICCom3D), intent(inout) :: this
            type(ParticleBundle), intent(inout) :: pb
            real(8), intent(in) :: xstart, xend, ystart, yend,zstart,zend
            integer(4) :: ierr, index, i,Indx(6,pb%NPar+5000),l(6),j,test(com_negb_num),dim
            integer(4) :: particle_number_send(com_negb_num)
            integer(4) :: particle_number_recv(com_negb_num),negb_rank(com_negb_num),boundarytest,temp1,temp2
            type(ParticleBundle) :: particle_send_buff(com_negb_num)
            type(ParticleBundle) :: particle_recv_buff(com_negb_num)


            ! 处理超过domain的粒子
            particle_number_send = 0
            particle_number_recv = 0
            
            l=1
           
            if (pb%NPar > 0) then
                do i = 1, pb%NPar
                    index = getParticleDomainIndex(this,pb%PO(i), xstart, xend, ystart, yend,zstart,zend)    
                    if(index>0) then
                    Indx(index,l(index))=i
                    l(index)=l(index)+1
                    end if
                end do
            end if
            !此部分标记了所有不在domain的粒子，粒子编号保存在index数组中，l(index)记录了应传送到每个邻居domain的粒子个数
           
            do dim=1,3
               
                test=0
                select case (dim)
                case (1)
                    
                    negb_rank(1)=this%rank-1
                    negb_rank(2)=this%rank+1
                    boundarytest=1
                case (2)
                    negb_rank(1)=this%rank-this%xyz_np(1)
                    negb_rank(2)=this%rank+this%xyz_np(1)
                    boundarytest=this%xyz_np(1)
                case (3)
                    negb_rank(1)=this%rank-this%xyz_np(1)*this%xyz_np(2)
                    negb_rank(2)=this%rank+this%xyz_np(1)*this%xyz_np(2)
                    boundarytest=this%xyz_np(1)*this%xyz_np(2)
                case default
                end select
                    temp1=mod((this%rank)/boundarytest,this%xyz_np(dim)) 
                do i = 1, 2       
                    temp2=mod((negb_rank(i))/boundarytest,this%xyz_np(dim))
                    
                    if(i==1) then

                        if (temp1-temp2==1.and.temp2>=0) test(i)=1!标记检测结果
                    else
                        if (temp2-temp1==1) test(i)=1
                        !此部分利用取余操作判断了边界上的网格所包含的邻居节点个数，为避免重复判断将结果保存在test数组中
                    end if    
                    if(test(i)==1) then
                      
                        particle_number_send(i)=l(i+dim*2-2)-1
                       
                       
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
                    ! if (test(i)==1) then!检测邻居节点是否存在
                    call particle_send_buff(i)%init(particle_number_send(i))
                    call particle_recv_buff(i)%init(particle_number_recv(i))
                ! end if

                    
                end do
               
                do i = 1, 2
                    if(l(i+dim*2-2)>1) then
                    do j=1,(l(i+dim*2-2)-1)
                
                    call particle_send_buff(i)%addone(pb%PO(Indx(i+dim*2-2,j)))
                    end do
                    
                    !  write(*,*)"particle_send_buff(i)%PO%X" ,particle_send_buff(i)%PO(particle_number_send(i)-1)%X,particle_number_send(i)
                 end if
                end do

                ! send and recv
                do i = 1, 2    
                    if (test(i)==1) then!检测邻居节点是否存在
                       
                        if (particle_number_send(i) > 0) &
                        call MPI_SEND(particle_send_buff(i)%PO, particle_number_send(i), &
                                   this%mpi_type_particle_one, negb_rank(i), &
                                   1, MPI_COMM_WORLD, this%reqs_send(i), ierr)

                    if (particle_number_recv(i) > 0) &
                    call MPI_RECV(particle_recv_buff(i)%PO, particle_number_recv(i), &
                                   this%mpi_type_particle_one, negb_rank(i), 1, &
                                   MPI_COMM_WORLD, this%reqs_recv(i), ierr)
                     end if
                end do

                do i = 1, com_negb_num
                    if (test(i)==1 .and. particle_number_recv(i) > 0) then
                         particle_recv_buff(i)%NPar=particle_number_recv(i)
                        
                        
                             do j=1,particle_recv_buff(i)%NPar
                                 index = getParticleDomainIndex(this,particle_recv_buff(i)%PO(j), xstart, xend, ystart, yend,zstart,zend)    
                               
                                 if(index>0) then
                                 Indx(index,l(index))=j+pb%NPar 
                                 l(index)=l(index)+1
                                 end if
                                
                                 !此部分将传送来的粒子进行标记，index为0表明传送粒子在该domain中，否则存入index数组进行下一步传送
                                 
                             end do  
                            call pb%addbun(particle_recv_buff(i))
                      

                    end if
                end do

                do i = 1, 2
                    call particle_send_buff(i)%destroy()
                    call particle_recv_buff(i)%destroy()
                end do
            end do
            do i=1,6
                
                if(l(i)>1) then
                    do j=1,l(i)-1
                    call pb%DelOne2(Indx(i,j))
                    ! 为了减少重复检测粒子在哪个domain，重写了删除粒子函数，只对需删除粒子进行标记操作，
                    !然后在DelOneready() 函数中将标记粒子删除
                    
                    
                    end do
                end if
            end do  
            call pb%DelOneready()  

        end subroutine communicationPICParticleCom3D



        function getParticleDomainIndex(this,one,  xstart, xend, ystart, yend,zstart,zend)
            type(ParticleOne), intent(in) :: one
            class(PICCom3D), intent(inout) :: this
            real(8), intent(in) ::  xstart, xend, ystart, yend,zstart,zend
            integer(4) :: getParticleDomainIndex
            if (one%X < xstart) then
                    getParticleDomainIndex = 1  
            else if (one%X > xend) then
                    getParticleDomainIndex = 2
            else if (one%Y < ystart) then
                    getParticleDomainIndex = 3
            else if (one%Y > yend) then
                    getParticleDomainIndex = 4
            else if (one%Z < zstart) then
                     getParticleDomainIndex = 5
            else if (one%Z > zend) then
                     getParticleDomainIndex =6 
            else
                    getParticleDomainIndex =0
          end if
                            
        end

end module ModulePICParticleCom3d
