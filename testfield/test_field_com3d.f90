! 测试分两部分，第一部分是求和测试，此部分将所有值赋值为1，场通讯求和后会表现为
!有一个邻居的值为2，两个邻居的值为3，对于角点，有七个邻居，值为8，此部分结果保存在
!sum文件夹下

!另一个测试为电势测试，测试区域分为8个方块，
!上下左右前后分别两个方块，建立三维坐标系，设置电势沿x，y,z方向均匀增加，坐标范围(0,0,0)
!到(1,1,1),其中点(0,0,0)处电势为0
!(1,1,1)处电势为3，以此类推，边界电势为0
!初始时刻只对xstart到xend,ystart到yend,zstart到zend赋值，而通过场通讯获得 
!（xstart-1,:,:）,(xend+1,:,:)等处的值,在场通讯完后计算出电场分布(Ex,Ey,Ez)并绘图，
!此部分结果保存在ext文件夹下,在结果正确的情况下，可以预计电场为匀强电场。

program fortran_mpi
    use mpi
    use ModulePICFieldCom
    implicit none
    
    integer(4) :: size, rank, ierr, i,j, k
    type(PICCom3D) :: mycom
    real(8), allocatable :: array(:, :,:), array_ext(:, :,:)
    real(8), allocatable :: array_out(:, :,:),electric_field_x(:, :,:),electric_field_y(:, :,:),electric_field_z(:, :,:)
    character(len=99) :: file_name

    integer(4) :: lx = 3, ly = 3,lz=3
    integer(4) :: xyz_np(3)=[2,2,2],com_type_nonblock, com_field_opt_sum
    integer(4) :: xstart, xend, ystart, yend,zstart,zend

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

     call mycom%init(lx, ly, lz,xyz_np, com_type_nonblock)

    xstart = mycom%col * lx + 1
    xend   = (mycom%col+1) * lx
    ystart = mycom%row * ly + 1
    yend   = (mycom%row+1) * ly
    zstart=mycom%layer*lz+1
    zend=(mycom%layer+1)*lz
!生成的是2x2的网格，两条边有三个节点
    allocate(array(xstart:xend, ystart:yend,zstart:zend))
    allocate(array_ext(xstart-1:xend+1, ystart-1:yend+1,zstart-1:zend+1))
    allocate(electric_field_x(xstart:xend+1, ystart:yend+1,zstart:zend+1))
    allocate(electric_field_y(xstart:xend+1, ystart:yend+1,zstart:zend+1))
    allocate(electric_field_z(xstart:xend+1, ystart:yend+1,zstart:zend+1))
    !电场插值到半格点上，实际位置为xstart-0.5：xend+0.5,以此类推，和数组索引有0.5的差值
    ! allocate(array_out(lx*px, ly*py))
    array_ext=0
    do i = xstart, xend
        do j = ystart, yend
            do k=zstart,zend
            array(i, j,k) = 1 
             array_ext(i,j,k)=(i-xstart*1.0)/(lx-1)+xstart/lx+(j-ystart*1.0)/(ly-1)+ystart/ly+(k-zstart*1.0)/(lz-1)+zstart/lz !电势赋值
            ! array_ext(i,j,k)=1
        end do
    end do
    end do

    ! array_ext = -1.d0
    ! array_ext(xstart:xend, ystart:yend) = array(xstart:xend, ystart:yend)

    ! mycom%left_ext_type = com_field_ext_type_symmetry
    write(file_name, '(i1)') rank
    open(10, file="./sum/raw_field_"//trim(file_name)//".txt")
        write(10,*)array !保存求和场 write数组是先遍历列，即先x,再y,最后z
        !所保存的先是底层，再是中层，然后是顶层
    close(10)
    write(file_name, '(i1)') rank
    open(11, file="./ext/potential/raw_potential"//trim(file_name)//".txt")
        write(11,"(5f9.3,/)")array_ext
    close(11)
    
     call mycom%comfsum(array,xstart, xend, ystart, yend,zstart,zend)
     
      call mycom%comfext(array_ext,xstart, xend, ystart, yend,zstart,zend)
    
!求解电场
      do i = xstart, xend+1
        do j = ystart, yend+1
            do k=zstart,zend+1
                electric_field_x(i,j,k)=(array_ext(i,j,k)-array_ext(i-1,j,k))/1
                electric_field_y(i,j,k)=(array_ext(i,j,k)-array_ext(i,j-1,k))/1
                electric_field_z(i,j,k)=(array_ext(i,j,k)-array_ext(i,j,k-1))/1
            end do
        end do     
      end do
                !/1为除以网格间距dx,这里为方便设为dx=1
    !  call mycom%gather(array, array_out, xstart, xend, ystart, yend, lx*px, ly*py)

    
     ! dump
    write(file_name, '(i1)') rank
    open(12, file="./sum/final_field_"//trim(file_name)//".txt")
        write(12,*)array
    close(12)
    write(file_name, '(i1)') rank
    open(13, file="./ext/potential/final_potential"//trim(file_name)//".txt")
        write(13,"(5f9.3,/)")array_ext
    close(13)

    !将所有电场写入3个文件，分别为x方向电场，y方向电场，z方向电场
    write(file_name, '(i1)') rank
    open(14, file="./ext/field/fieldx"//trim(file_name)//".txt")
        write(14,"(4f9.3,/)")electric_field_x
    close(14)
    write(file_name, '(i1)') rank
    open(15, file="./ext/field/fieldy"//trim(file_name)//".txt")
        write(15,"(4f9.3,/)")electric_field_y
    close(15)
    write(file_name, '(i1)') rank
    open(16, file="./ext/field/fieldz"//trim(file_name)//".txt")
        write(16,"(4f9.3,/)")electric_field_z
    close(16)
       


    
    ! write(file_name, '(i2)') rank
    ! open(10, file="./result_field_com_"//trim(file_name)//".txt")
    !     do i = ystart-1, yend+1
    !         write(10, '(*(f10.4, 1x))') array_ext(:, i)
    !     end do
    ! close(10)

    deallocate(array)
    deallocate(array_ext)
    deallocate(electric_field_x)
    deallocate(electric_field_y)
    deallocate(electric_field_z)
    ! deallocate(array_ext)
    ! deallocate(array_out)

    ! call mycom%destroy()
     call MPI_FINALIZE(ierr)

end program fortran_mpi