module ModulePICFieldCom
integer(4), parameter :: com_field_negh=2

    type :: PICCom3D
        integer(4) :: negb_type(com_field_negh)
        integer(4) :: negb_rank(com_field_negh)
        integer(4) :: lx, ly,lz
        integer(4) ::  xyz_np(3)
        integer(4) :: size, rank
        integer(4) :: com_type
        logical :: is_init = .False.
        integer(4) :: col, row
        integer(4) :: reqs(2)
        ! field buffer
        real(8), allocatable :: send_buff(:,:,:)
        real(8), allocatable :: recv_buff(:,:,:)
        integer(4), allocatable :: orrd_x(:),orrd_y(:),orrd_z(:)

    contains
        procedure :: init       => initPICCom3D
        procedure :: comf => communicationPICFieldCom3D
        
    end type PICCom3D
contains
    subroutine initPICCom3D(this, x_local_size, y_local_size, deps_x, deps_y, com_type)
        class(PICCom3D), intent(inout) :: this
            integer(4), intent(in) :: x_local_size, y_local_size, deps_x, deps_y
            integer(4), optional, intent(in) :: com_type
            integer(4) :: ierr, i
    end subroutine initPICCom3D

    subroutine communicationPICFieldCom3D(this, array3d, opt_type, xstart, xend, ystart, yend,zstart,zend)
        class(PICCom3D), intent(inout) :: this
        real(8), allocatable, intent(inout) :: array3d(:, :,:)
        integer(4), intent(in) :: opt_type
        integer(4), intent(in) :: xstart, xend, ystart, yend,zstart,zend
        integer(4) :: particle_number_recv(com_field_negh),negb_rank(com_field_negh),boundarytest,temp1,temp2,boundarysize,test(2)
        integer(4) :: ierr,dim
       
            if (this%is_init) then
                if (xend-xstart+1 /= this%lx .or. yend-ystart+1 /= this%ly) then
                    if (0 == this%rank) then
                        write(*, *) "The input parameters (xstart, xend, ystart, yend) are invalid."
                    end if
                    stop
                end if
            end if    
            do dim=1,3
                select case(dim)
                case(1)
                    negb_rank(1)=this%rank-1
                    negb_rank(2)=this%rank+1
                    boundarytest=1
                    boundarysize=this%ly*this%lz
                    allocate(this%send_buff(2,this%ly,this%lz))
                    allocate(this%recv_buff(2,this%ly,this%lz))
                    allocate(this%orrd_x(2))
                    allocate(this%orrd_y(ystart:yend))
                    allocate(this%orrd_z(zstart:zend))
                    this%orrd_x(:)=(/xstart,xend/)
                    
                case(2)
                    negb_rank(1)=this%rank-this%xyz_np(1)
                    negb_rank(2)=this%rank+this%xyz_np(1)
                    boundarytest=this%xyz_np(1)
                    boundarysize=this%lx*this%lz
                    allocate(this%send_buff(2,this%lx,this%lz))
                    allocate(this%recv_buff(2,this%lx,this%lz))
                    allocate(this%orrd_x(xstart:xend))
                    allocate(this%orrd_y(2))
                    allocate(this%orrd_z(zstart:zend))
                    this%orrd_y(:)=(/ystart,yend/)

                case(3)
                    negb_rank(1)=this%rank-this%xyz_np(1)*this%xyz_np(2)
                    negb_rank(2)=this%rank+this%xyz_np(1)*this%xyz_np(2)
                    boundarytest=this%xyz_np(1)*this%xyz_np(2)
                    boundarysize=this%lx*this%ly
                    allocate(this%send_buff(2,this%lx,this%ly))
                    allocate(this%recv_buff(2,this%lx,this%ly))
                    allocate(this%orrd_x(xstart:xend))
                    allocate(this%orrd_y(ystart:yend))
                    allocate(this%orrd_z(2))
                    this%orrd_z(:)=(/zstart,zend/)
                case default    
                end select    
            
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
                        this%send_buff(i,:,:) = array3d(this%orrd_x(i), this%orrd_y,this%orrd_z)
                        call MPI_SEND(this%send_buff(i,:,:), boundarysize, MPI_DOUBLE, negb_rank(i), 1, MPI_COMM_WORLD, this%reqs(1), ierr)
                        call MPI_RECV(this%recv_buff(i,:,:), boundarysize, MPI_DOUBLE, negb_rank(i), 1, MPI_COMM_WORLD, this%reqs(2), ierr)
                        array3d(this%orrd_x(i), this%orrd_y,this%orrd_z) = array3d(this%orrd_x(i), this%orrd_y,this%orrd_z) + this%recv_buff
                    end if 
                end do
                if (allocated(this%send_buff))    deallocate(this%send_buff)
                if (allocated(this%recv_buff))    deallocate(this%recv_buff)
                if (allocated(this%orrd_x))    deallocate(this%orrd_x)
                if (allocated(this%orrd_y))    deallocate(this%orrd_y)
                if (allocated(this%orrd_z))    deallocate(this%orrd_z)
            end do

    end subroutine communicationPICFieldCom3D
end module ModulePICFieldCom