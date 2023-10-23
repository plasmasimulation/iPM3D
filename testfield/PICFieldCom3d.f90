module ModulePICFieldCom
    use mpi
integer(4), parameter :: com_field_negh=2

    type :: PICCom3D
        integer(4) :: negb_type(com_field_negh)
        integer(4) :: negb_rank(com_field_negh)
        integer(4) :: lx, ly,lz
        integer(4) ::  xyz_np(3)
        integer(4) :: size, rank
        integer(4) :: com_type
        logical :: is_init = .False.
        integer(4) :: col, row,layer
        integer(4) :: reqs_send(2),reqs_recv(2)
        ! field buffer
        real(8), allocatable :: send_buff(:,:,:)
        real(8), allocatable :: recv_buff(:,:,:)
         integer(4), allocatable :: orrd_x(:,:),orrd_y(:,:),orrd_z(:,:)
         integer(4) :: status_send(MPI_STATUS_SIZE, com_field_negh), status_recv(MPI_STATUS_SIZE,com_field_negh)

    contains
        procedure :: init       => initPICCom3D
        procedure :: comf => communicationPICFieldCom3D
        procedure :: destroy=>destroyFieldCom3D
        
    end type PICCom3D
contains
subroutine destroyFieldCom3D(this)
    class(PICCom3D), intent(inout) :: this
  
        if (allocated(this%send_buff))deallocate(this%send_buff)
        if (allocated(this%recv_buff))deallocate(this%recv_buff)
        if (allocated(this%orrd_x))deallocate(this%orrd_x)
        if (allocated(this%orrd_y))deallocate(this%orrd_y)
        if (allocated(this%orrd_z))deallocate(this%orrd_z)
       
    end subroutine destroyFieldCom3D
subroutine initPICCom3D(this, x_local_size, y_local_size, z_local_size, xyz_np, com_type)
    class(PICCom3D), intent(inout) :: this
    integer(4), intent(in) :: x_local_size, y_local_size,z_local_size, xyz_np(3)
    integer(4), optional, intent(in) :: com_type
    integer(4) :: ierr

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
       
        this%is_init = .True.
    end if
end subroutine initPICCom3D

    subroutine communicationPICFieldCom3D(this, array3d, opt_type, xstart, xend, ystart, yend,zstart,zend)
        class(PICCom3D), intent(inout) :: this
        real(8), allocatable, intent(inout) :: array3d(:, :,:)
        integer(4), intent(in) :: opt_type
        integer(4), intent(in) :: xstart, xend, ystart, yend,zstart,zend
        integer(4) :: size,negb_rank(com_field_negh),boundarytest,temp1,temp2,boundarysize,test(2)
        integer(4) :: ierr,dim,arrysize(3)
        integer(4):: send,recv
        integer(4) :: reqs_send(6),reqs_recv(6)
        send=1
        recv=0
        
       
            if (this%is_init) then
                if (xend-xstart+1 /= this%lx .or. yend-ystart+1 /= this%ly) then
                    if (0 == this%rank) then
                        write(*, *) "The input parameters (xstart, xend, ystart, yend) are invalid."
                    end if
                    stop
                end if
            end if    
            do dim=1,3
              if (allocated(this%send_buff))deallocate(this%send_buff)
              if (allocated(this%recv_buff))deallocate(this%recv_buff)
              if (allocated(this%orrd_x))deallocate(this%orrd_x)
              if (allocated(this%orrd_y))deallocate(this%orrd_y)
              if (allocated(this%orrd_z))deallocate(this%orrd_z)
       
                select case(dim)
                case(1)
                   
                    negb_rank(1)=this%rank-1
                    negb_rank(2)=this%rank+1
                    boundarytest=1
                    boundarysize=this%ly*this%lz
                    allocate(this%send_buff(2,this%ly,this%lz))
                    allocate(this%recv_buff(2,this%ly,this%lz))
                     allocate(this%orrd_x(2,1))
                     allocate(this%orrd_y(2,this%ly))
                     allocate(this%orrd_z(2,this%lz))
                    this%orrd_x(:,1)=[xstart,xend]
                    this%orrd_y(1,:)=(/(i,i=ystart,yend)/)
                    this%orrd_y(2,:)=(/(i,i=ystart,yend)/)
                    this%orrd_z(1,:)=(/(i,i=zstart,zend)/)
                    this%orrd_z(2,:)=(/(i,i=zstart,zend)/)
                    arrysize=(/1,this%ly,this%lz/)

                    
                case(2)
                    
                    negb_rank(1)=this%rank-this%xyz_np(1)
                    negb_rank(2)=this%rank+this%xyz_np(1)
                    boundarytest=this%xyz_np(1)
                    boundarysize=this%lx*this%lz
                    allocate(this%send_buff(2,this%lx,this%lz))
                    allocate(this%recv_buff(2,this%lx,this%lz))
                    allocate(this%orrd_x(2,this%lx))
                    allocate(this%orrd_y(2,1))
                    allocate(this%orrd_z(2,this%lz))
                    this%orrd_y(:,1)=(/ystart,yend/)
                    this%orrd_x(1,:)=(/(i,i=xstart,xend)/)
                    this%orrd_x(2,:)=(/(i,i=xstart,xend)/)
                    this%orrd_z(1,:)=(/(i,i=zstart,zend)/)
                    this%orrd_z(2,:)=(/(i,i=zstart,zend)/)
                    arrysize=(/this%lx,1,this%lz/)

                case(3)
                    negb_rank(1)=this%rank-this%xyz_np(1)*this%xyz_np(2)
                    negb_rank(2)=this%rank+this%xyz_np(1)*this%xyz_np(2)
                    boundarytest=this%xyz_np(1)*this%xyz_np(2)
                    boundarysize=this%lx*this%ly
                    allocate(this%send_buff(2,this%lx,this%ly))
                    allocate(this%recv_buff(2,this%lx,this%ly))
                    allocate(this%orrd_x(2,this%lx))
                    allocate(this%orrd_y(2,this%ly))
                    allocate(this%orrd_z(2,1))
                    this%orrd_z(:,1)=(/zstart,zend/)
                    this%orrd_y(1,:)=(/(i,i=ystart,yend)/)
                    this%orrd_y(2,:)=(/(i,i=ystart,yend)/)
                    this%orrd_x(1,:)=(/(i,i=xstart,xend)/)
                    this%orrd_x(2,:)=(/(i,i=xstart,xend)/)
                    arrysize=(/this%ly,this%lx,1/)
                case default    
                end select    
                test=0
                this%recv_buff = 0
                temp1=mod((this%rank)/boundarytest,this%xyz_np(dim)) 
                do i = 1, 2       
                    temp2=mod((negb_rank(i))/boundarytest,this%xyz_np(dim))    
                    if(i==1) then
                        if (temp1-temp2==1.and.temp2>=0) test(i)=1!标记检测结果
                    else 
                         if (temp2-temp1==1) test(i)=1    !此部分利用取余操作判断了边界上的网格所包含的邻居节点个数，为避免重复判断将结果保存在test数组中
                    end if    
                    if(test(i)==1) then
                           this%send_buff(i,:,:) = reshape(array3d(this%orrd_x(i,:), this%orrd_y(i,:),this%orrd_z(i,:)),(/5,5/))
                    !        if (0 == this%rank) then
                    !         write(*,*)this%send_buff(i,:),this%rank,"this%send_buff(k(i,:),:,:)"
                    !        write(*,*)array3d(this%orrd_x(i,:), this%orrd_y(i,:),this%orrd_z(i,:)),this%rank,"array3d(this%orrd_x(i,:), this%orrd_y(i,:),this%orrd_z(i,:))"
                    !       write(*,*)boundarysize
                    !   end if
                           
                              call MPI_SEND(this%send_buff(i,:,:), boundarysize, MPI_DOUBLE, negb_rank(i), 1, MPI_COMM_WORLD, reqs_send(i+dim*2-2), ierr)
                              call MPI_RECV(this%recv_buff(i,:,:), boundarysize, MPI_DOUBLE, negb_rank(i), 1, MPI_COMM_WORLD, reqs_recv(i+dim*2-2), ierr)
                            !   if(i==1) then
                            !     array3d(this%orrd_x(i,:), this%orrd_y(i,:),this%orrd_z(i,:)) = array3d(this%orrd_x(i,:), this%orrd_y(i,:),this%orrd_z(i,:)) + reshape(this%recv_buff(i,:,:),(/arrysize(1),arrysize(2),arrysize(3)/))
                            !   else
                            !     array3d(this%orrd_x(i,:), this%orrd_y(i,:),this%orrd_z(i,:)) =  reshape(this%recv_buff(i,:,:),(/arrysize(1),arrysize(2),arrysize(3)/))
                            !   end if
                   
                    end if 
                end do
                do i=1,2
                    if(test(i)==1) then
                        array3d(this%orrd_x(i,:), this%orrd_y(i,:),this%orrd_z(i,:)) = array3d(this%orrd_x(i,:), this%orrd_y(i,:),this%orrd_z(i,:)) + reshape(this%recv_buff(i,:,:),(/arrysize(1),arrysize(2),arrysize(3)/))
                     end if   
                end do
                !  do i = 1, 2
                !      if (test(i)==1) then!检测邻居节点是否存在
                !          call MPI_WAIT(this%send_buff(i,:,:), this%status_send(:, i), ierr)
                !          call MPI_WAIT(this%recv_buff(i,:,:), this%status_recv(:, i), ierr)
                !         !  array3d(this%orrd_x(i,:), this%orrd_y(i,:),this%orrd_z(i,:)) = array3d(this%orrd_x(i,:), this%orrd_y(i,:),this%orrd_z(i,:)) + reshape(this%recv_buff(i,:,:),(/arrysize(1),arrysize(2),arrysize(3)/))
                !      end if
                !  end do
            !    call this%destroy()
            end do

    end subroutine communicationPICFieldCom3D
end module ModulePICFieldCom