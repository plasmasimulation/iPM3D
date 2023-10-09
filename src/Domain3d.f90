Module ModuleDomain3d
    use mpi
    use ModuleFileName
    use ModuleParallelDump

    implicit none

    integer(4), parameter :: DIM_DOMAIN = 3         !三维
    integer(4), parameter :: BOUNDARY_PER_DIM = 3   !三维
    integer(4), parameter :: BOUNDARY_LEFT = 1      !uncertain
    integer(4), parameter :: BOUNDARY_RIGHT = 2     !uncertain
    integer(4), parameter :: BUFFER_SIZE = BOUNDARY_PER_DIM * DIM_DOMAIN
    integer(4), parameter :: NEIGHBOR_TYPE_BOUNDARY = 11  !uncertain
    integer(4), parameter :: NEIGHBOR_TYPE_DOMAIN = 12
    integer(4), parameter :: NEIGHBOR_TYPE_MID = 13

    integer(4), parameter :: deps_type_uniform = 0 !deposit particle
    integer(4), parameter :: deps_type_specify = 1
    integer(4), parameter :: deps_type_balance = 2

    integer(4), save :: deps_type = deps_type_uniform

    type Domain
        type(FileName)  :: IOName
        integer(4)  :: GlobalShape(DIM_DOMAIN) !全局形态，xyz总分布
        integer(4)  :: LocalShape(DIM_DOMAIN) !单个进程形态，部分x,y,z
        real(8)     :: SpaceStep(DIM_DOMAIN) !网格对应空间大小

        integer(4)  :: MyId = 0, master = 1
        integer(4)  :: DomainIndex(DIM_DOMAIN) !区间编号，三维数组，作用未知
        integer(4)  :: NProcess = 0         !processor number
        integer(4)  :: ImageShape(DIM_DOMAIN)  !镜像空间形状
        
        integer(4)  :: NeighborDomainIndex(BOUNDARY_PER_DIM, DIM_DOMAIN, DIM_DOMAIN) !邻居空间三维邻居6个
        integer(4)  :: NeighborImageId(BOUNDARY_PER_DIM, DIM_DOMAIN) !需要区分id和domainindex的区别
        integer(4)  :: NeighborType(BOUNDARY_PER_DIM, DIM_DOMAIN) !uncetrain

        integer(4)  :: CornerIndex(BOUNDARY_PER_DIM, DIM_DOMAIN) !neighbor和corner分别代表什么
        real(8)     :: CoordinateOrigin(DIM_DOMAIN)         !坐标
        real(8)     :: CornerCoordinateValue(BOUNDARY_PER_DIM, DIM_DOMAIN)

        integer(4), allocatable :: LocalShapeList(:, :)

    contains

        procedure :: Init => InitDomain3D
        procedure :: Destroy => DestroyDomain

        procedure :: Uniform => DepsUniform
        procedure :: Specify => DepsSpecify
        procedure :: Balance => DepsBalance

        procedure :: Dump => DumpDomain
        procedure :: Load => LoadDomain

        procedure :: ShowAll => ShowDomainAllInfo
        procedure :: Show => ShowDomain
        procedure :: Save => SaveDomainDepsInfo

        procedure, private :: SetLocalShape

        procedure, private :: LoadDomainHDF5
        procedure, private :: LoadDomainDAT

        procedure, private :: DumpDomainHDF5
        procedure, private :: DumpDomainDAT

    end type Domain

    type(HDF5_PDump), save, private :: hdf5DomainDump

    contains

        subroutine InitDomain3D(this, Nx, Ny,Nz, x_start, x_end, y_start, y_end,z_start,z_end,px, py, pz, dmname)
            class(Domain), intent(inout) :: this
            integer(4), intent(in) :: Nx, Ny,Nz !事先规定的网格数
            real(8), intent(in) :: x_start, x_end, y_start, y_end,z_start,z_end !所有进程模拟的区域
            integer(4), intent(in) :: px, py,pz !x,x,z方向上的process number
            character(*), optional, intent(in) :: dmname
            integer(4) :: size, rank, ierr

            call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

            if (Nx > 1 .and. Ny > 1 .and. Nz>1) then
                this%GlobalShape = [Nx, Ny,Nz]
                this%CoordinateOrigin = [x_start, y_start,z_start]
                this%SpaceStep = [x_end - x_start, y_end - y_start,z_end-z_start] / (this%GlobalShape - 1)!每个网格的宽度，长度，高度

                this%MyId = rank
                this%NProcess = size

                if (px*py*pz == this%NProcess) then
                    this%ImageShape = [px, py,pz]
                else
                    this%ImageShape = [1, 1,this%NProcess]!单个进程处理的区域，待定
                end if

                this%DomainIndex = getDomainIndex(this%MyId, this%ImageShape, DIM_DOMAIN)

                if (present(dmname)) then
                    call this%IOName%Init(dmname, RESTART_FILE_NAME)
                else
                    call this%IOName%Init("Domain", RESTART_FILE_NAME)
                end if
                
                call this%Uniform()
            else
                if (0 == this%MyId) write(*, '(a)') "The Domain init failed, please check input parameter."
                stop
            end if

        end subroutine InitDomain3D


        subroutine DepsUniform(this) !确定每个进程处理的区域
            class(Domain), intent(inout) :: this
            integer(4) :: i, j, k

            ! set list of image size
            call this%Destroy()
            allocate(this%LocalShapeList(DIM_DOMAIN, this%NProcess))

            this%LocalShapeList = 0
            do i = 1, DIM_DOMAIN
                do j = 1, this%ImageShape(i)
                    this%LocalShapeList(i, j) = int((this%GlobalShape(i)-1) / this%ImageShape(i)) + 1 !网格数除进程数得到每个进程处理的网格数
                    if (j <= (this%GlobalShape(i)-1 - int((this%GlobalShape(i)-1) / this%ImageShape(i)) * this%ImageShape(i))) &
                        this%LocalShapeList(i, j) = this%LocalShapeList(i, j) + 1 !操作未知
                end do
            end do
            
            call this%SetLocalShape()

            ! mid-axis
            if (this%DomainIndex(2) == 0) then
                this%NeighborImageId(BOUNDARY_LEFT, 2) = -1
                this%NeighborType(BOUNDARY_LEFT, 2) = NEIGHBOR_TYPE_MID
            end if

            call this%Save()

        end subroutine DepsUniform


        subroutine DepsSpecify(this, deps_x, deps_y,deps_z)
            class(Domain), intent(inout) :: this
            real(8), intent(in) :: deps_x(this%ImageShape(1)), deps_y(this%ImageShape(2)),deps_z(this%ImageShape(3))
            integer(4) :: dp_x(this%ImageShape(1)), dp_y(this%ImageShape(2)),dp_z(this%ImageShape(3))
            integer(4) :: i, j
            real(8) :: tmp_sum

            ! norm
            if (this%ImageShape(1) == 1) then
                dp_x(1) = this%GlobalShape(1)

            else
                tmp_sum = sum(deps_x)
                do j = 1, this%ImageShape(1)
                    dp_x(j) = int(deps_x(j) / tmp_sum * dble(this%GlobalShape(1)-1) + 0.5d0) + 1
                end do

                if (sum(dp_x-1) /= this%GlobalShape(1)-1) then
                    dp_x(this%ImageShape(1)) = this%GlobalShape(1) - sum(dp_x(1:this%ImageShape(1)-1)-1)
                end if
            end if

            if (this%ImageShape(2) == 1) then
                dp_y(1) = this%GlobalShape(2)

            else
                tmp_sum = sum(deps_y)
                do j = 1, this%ImageShape(2)
                    dp_y(j) = int(deps_y(j) / tmp_sum * dble(this%GlobalShape(2)-1) + 0.5d0) + 1
                end do

                if (sum(dp_y-1) /= this%GlobalShape(2)-1) then
                    dp_y(this%ImageShape(2)) = this%GlobalShape(2) - sum(dp_y(1:this%ImageShape(2)-1)-1)
                end if
            end if

            if (this%ImageShape(3) == 1) then
                dp_z(1) = this%GlobalShape(3)

            else
                tmp_sum = sum(deps_z)
                do j = 1, this%ImageShape(2)
                    dp_z(j) = int(deps_z(j) / tmp_sum * dble(this%GlobalShape(2)-1) + 0.5d0) + 1
                end do

                if (sum(dp_z-1) /= this%GlobalShape(2)-1) then
                    dp_z(this%ImageShape(2)) = this%GlobalShape(2) - sum(dp_z(1:this%ImageShape(2)-1)-1)
                end if
            end if

            ! set list of image size
            call this%Destroy()
            allocate(this%LocalShapeList(DIM_DOMAIN, this%NProcess))

            this%LocalShapeList = 0
            this%LocalShapeList(1, 1:this%ImageShape(1)) = dp_x
            this%LocalShapeList(2, 1:this%ImageShape(2)) = dp_y
            this%LocalShapeList(3, 1:this%ImageShape(1)) = dp_z

            call this%SetLocalShape()

            ! mid-axis
            if (this%DomainIndex(2) == 0) then
                this%NeighborImageId(BOUNDARY_LEFT, 2) = -1
                this%NeighborType(BOUNDARY_LEFT, 2) = NEIGHBOR_TYPE_MID
            end if

            call this%Save()

        end subroutine DepsSpecify


        subroutine DepsBalance(this)
            class(Domain), intent(inout) :: this

        end subroutine DepsBalance


        subroutine SetLocalShape(this)
            class(Domain), intent(inout) :: this
            integer(4) :: i, j
            integer(4) :: tempDomainIndex(DIM_DOMAIN)

            ! local size
            do i = 1, DIM_DOMAIN
                this%LocalShape(i) = this%LocalShapeList(i, this%DomainIndex(i)+1)
            end do

            ! set index for local image
            do i = 1, DIM_DOMAIN
                this%CornerIndex(BOUNDARY_LEFT, i) = 1 + sum(this%LocalShapeList(i, 1:this%DomainIndex(i)) - 1)
                this%CornerIndex(BOUNDARY_RIGHT, i) = 1 + sum(this%LocalShapeList(i, 1:this%DomainIndex(i)+1) - 1)
            end do

            ! set coordinate values for local image
            do i = 1, DIM_DOMAIN
                this%CornerCoordinateValue(BOUNDARY_LEFT, i) = this%CoordinateOrigin(i) + &
                                                             dble((this%CornerIndex(BOUNDARY_LEFT, i) - 1)) * this%SpaceStep(i)

                this%CornerCoordinateValue(BOUNDARY_RIGHT, i) = this%CoordinateOrigin(i) + &
                                                             dble((this%CornerIndex(BOUNDARY_RIGHT, i) - 1)) * this%SpaceStep(i)

            end do

            ! set neighbor image index
            this%NeighborType = NEIGHBOR_TYPE_DOMAIN
            do i = 1, DIM_DOMAIN
                tempDomainIndex = this%DomainIndex
                tempDomainIndex(i) = tempDomainIndex(i) - 1
                this%NeighborDomainIndex(BOUNDARY_LEFT, i, 1:DIM_DOMAIN) = tempDomainIndex
                this%NeighborImageId(BOUNDARY_LEFT, i) = getImageIndex(tempDomainIndex, this%ImageShape, DIM_DOMAIN)

                if (tempDomainIndex(i) < 0 .or. &
                    tempDomainIndex(i) > this%ImageShape(i)-1) then

                    this%NeighborImageId(BOUNDARY_LEFT, i) = -1
                    this%NeighborType(BOUNDARY_LEFT, i) = NEIGHBOR_TYPE_BOUNDARY
                end if

                tempDomainIndex = this%DomainIndex
                tempDomainIndex(i) = tempDomainIndex(i) + 1
                this%NeighborDomainIndex(BOUNDARY_RIGHT, i, 1:DIM_DOMAIN) = tempDomainIndex
                this%NeighborImageId(BOUNDARY_RIGHT, i) = getImageIndex(tempDomainIndex, this%ImageShape, DIM_DOMAIN)

                if (tempDomainIndex(i) < 0 .or. &
                    tempDomainIndex(i) > this%ImageShape(i)-1) then

                    this%NeighborImageId(BOUNDARY_RIGHT, i) = -1
                    this%NeighborType(BOUNDARY_RIGHT, i) = NEIGHBOR_TYPE_BOUNDARY
                end if
            end do

        end subroutine SetLocalShape


        subroutine DestroyDomain(this)
            class(Domain), intent(inout) :: this

            if (allocated(this%LocalShapeList)) deallocate(this%LocalShapeList)

        end subroutine DestroyDomain


        subroutine ShowDomainAllInfo(this)
            class(Domain), intent(inout) :: this
            integer(4) :: i, j, k
            integer(4) :: size, rank, ierr

            call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

            do k = 0, size-1
                if (rank == k) then
                    write(*, *) ""
                    write(*, '(a30, i)') "Dimension: ", DIM_DOMAIN
                    write(*, '(a30, *(1x, i))') "Global shape: ", this%GlobalShape
                    write(*, '(a30, *(1x, i))') "Local shape: ", this%LocalShape
                    write(*, '(a30, *(1x, es))') "Space step: ", this%SpaceStep
                    write(*, '(a30, i, i)') "Domain id: ", this%MyId, this%NProcess
                    write(*, '(a30, *(1x, i))') "Domain index: ", this%DomainIndex
                    write(*, '(a30, *(1x, i))') "Image shape: ", this%ImageShape
                    write(*, '(a30, *(1x, es))') "Coordinate Origin: ", this%CoordinateOrigin
                    write(*, '(a30, *(1x, i))') "Corner Start: ", this%CornerIndex(BOUNDARY_LEFT, :)
                    write(*, '(a30, *(1x, i))') "Corner End: ", this%CornerIndex(BOUNDARY_RIGHT, :)
                    write(*, '(a30, *(1x, es))') "Corner Start Value: ", this%CornerCoordinateValue(BOUNDARY_LEFT, :)
                    write(*, '(a30, *(1x, es))') "Corner End Value: ", this%CornerCoordinateValue(BOUNDARY_RIGHT, :)

                    do i = 1, DIM_DOMAIN
                        write(*, '(a30, i)') "Info in dimension-", i
                        do j = 1, BOUNDARY_PER_DIM
                            write(*, '(a30, *(1x, i))') "Neighbor Domain index: ", this%NeighborDomainIndex(j, i, :)
                            write(*, '(a30, *(1x, i))') "Neighbor Domain id: ", this%NeighborImageId(j, i)
                            write(*, '(a30, *(1x, i))') "Neighbor Domain type: ", this%NeighborType(j, i)
                        end do

                        write(*, '(a30, *(1x, i))') "The list of local shape: ", this%LocalShapeList(i, 1:this%ImageShape(i))
                    end do
                    write(*, *) ""

                end if

                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            end do

        end subroutine ShowDomainAllInfo


        subroutine ShowDomain(this)
            class(Domain), intent(inout) :: this
            integer(4) :: i, j

            if(0 == this%MyId) then
                write(*, '(a30, i)') "Dimension:", DIM_DOMAIN
                write(*, '(a30, *(1x, i))') "Global shape:", this%GlobalShape
                write(*, '(a30, *(1x, es12.4))') "Space step:", this%SpaceStep
                write(*, '(a30, *(1x, i))') "Image shape:", this%ImageShape

                do i = 1, DIM_DOMAIN
                    write(*, '(a30, *(1x, i))') "The list of local shape:", this%LocalShapeList(i, 1:this%ImageShape(i))
                end do
                write(*, *) ""

            end if

        end subroutine ShowDomain


        subroutine SaveDomainDepsInfo(this)
            class(Domain), intent(inout) :: this
            type(FileName) :: name
            integer :: k

            if (0 == this%MyId) then
                call name%Init('domain', CHECK_FILE_NAME)

                open(10, file=name%FullName%str, position='append')

                    write(10, '(*(i16, 1x))') DIM_DOMAIN, this%GlobalShape, this%ImageShape, &
                    (this%LocalShapeList(k, 1:this%ImageShape(k)), k=1, DIM_DOMAIN)

                close(10)
            end if

        end subroutine SaveDomainDepsInfo


        subroutine DumpDomain(this)
            class(Domain), intent(inout) :: this

            if (FILE_EXTENSION_MODE_H5 == this%IOName%ExtensionMode) then
                call this%DumpDomainHDF5()

            else if (FILE_EXTENSION_MODE_DAT == this%IOName%ExtensionMode) then
                call this%DumpDomainDAT()
            
            else
                call this%DumpDomainHDF5()

            end if

        endsubroutine DumpDomain


        subroutine LoadDomain(this)
            class(Domain), intent(inout) :: this
            logical :: alive

            Inquire(file=this%IOName%FullName%str, exist=alive)
            if (alive) then
                if (FILE_EXTENSION_MODE_H5 == this%IOName%ExtensionMode) then
                    call this%LoadDomainHDF5()

                else if (FILE_EXTENSION_MODE_DAT == this%IOName%ExtensionMode) then
                    call this%LoadDomainDAT()
                
                end if
            end if

        end subroutine LoadDomain


        subroutine LoadDomainHDF5(this)
            class(Domain), intent(inout) :: this
            integer(4), allocatable :: TempInt1D(:)
            integer(4), allocatable :: TempInt2D(:,:)
            real(8), allocatable :: TempReal1D(:)
            real(8), allocatable :: TempReal2D(:,:)
            real(8), allocatable :: TempReal3D(:,:,:)
            integer(4) :: size, rank, ierr

            call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

            if (0 == rank) write(*, '(a)') "The domain is load from hdf5 file."

            call hdf5DomainDump%init(filename=this%IOName%FullName%str, mode='read')
            call hdf5DomainDump%open()

            call hdf5DomainDump%readattr('/', 'MyId', this%MyId)
            call hdf5DomainDump%readattr('/', 'master', this%master)
            call hdf5DomainDump%readattr('/', 'NProcess', this%NProcess)

            if (allocated(TempInt1D)) deallocate(TempInt1D)
            allocate(TempInt1D(1:DIM_DOMAIN))
            call hdf5DomainDump%read('GlobalShape', TempInt1D)
            this%GlobalShape = TempInt1D

            call hdf5DomainDump%read('LocalShape', TempInt1D)
            this%LocalShape = TempInt1D

            call hdf5DomainDump%read('DomainIndex', TempInt1D)
            this%DomainIndex = TempInt1D

            call hdf5DomainDump%read('ImageShape', TempInt1D)
            this%ImageShape = TempInt1D

            if (allocated(TempReal1D)) deallocate(TempReal1D)
            allocate(TempReal1D(1:DIM_DOMAIN))
            call hdf5DomainDump%read('SpaceStep', TempReal1D)
            this%SpaceStep = TempReal1D

            call hdf5DomainDump%read('CoordinateOrigin', TempReal1D)
            this%CoordinateOrigin = TempReal1D

            if (allocated(TempInt2D)) deallocate(TempInt2D)
            allocate(TempInt2D(BOUNDARY_PER_DIM, DIM_DOMAIN))
            call hdf5DomainDump%read('NeighborImageId', TempInt2D)
            this%NeighborImageId = TempInt2D

            call hdf5DomainDump%read('NeighborType', TempInt2D)
            this%NeighborType = TempInt2D

            call hdf5DomainDump%read('CornerIndex', TempInt2D)
            this%CornerIndex = TempInt2D

            if (allocated(TempReal2D)) deallocate(TempReal2D)
            allocate(TempReal2D(BOUNDARY_PER_DIM, DIM_DOMAIN))
            call hdf5DomainDump%read('CornerCoordinateValue', TempReal2D)
            this%CornerCoordinateValue = TempReal2D

            if (allocated(TempReal3D)) deallocate(TempReal3D)
            allocate(TempReal3D(BOUNDARY_PER_DIM, DIM_DOMAIN, DIM_DOMAIN))
            call hdf5DomainDump%read('NeighborDomainIndex', TempReal3D)
            this%NeighborDomainIndex = TempReal3D

            call this%Destroy()
            allocate(this%LocalShapeList(DIM_DOMAIN, this%NProcess))

            call hdf5DomainDump%read('LocalShapeList', this%LocalShapeList)

            call hdf5DomainDump%close()

            if (allocated(TempInt1D)) deallocate(TempInt1D)
            if (allocated(TempInt2D)) deallocate(TempInt2D)
            if (allocated(TempReal1D)) deallocate(TempReal1D)
            if (allocated(TempReal2D)) deallocate(TempReal2D)
            if (allocated(TempReal3D)) deallocate(TempReal3D)

        end subroutine LoadDomainHDF5


        subroutine LoadDomainDAT(this)
            class(Domain), intent(inout) :: this
            integer(4) :: size, rank, ierr

            call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

            if (0 == rank) write(*, '(a)') "The domain is load from dat file."
            open(10, file=this%IOName%FullName%str)

                read(10, *) this%GlobalShape
                read(10, *) this%LocalShape
                read(10, *) this%SpaceStep

                read(10, *) this%MyId
                read(10, *) this%master
                read(10, *) this%DomainIndex
                read(10, *) this%NProcess
                read(10, *) this%ImageShape
                read(10, *) this%NeighborDomainIndex
                read(10, *) this%NeighborImageId
                read(10, *) this%NeighborType

                read(10, *) this%CornerIndex
                read(10, *) this%CornerCoordinateValue
                read(10, *) this%CoordinateOrigin

                call this%Destroy()
                allocate(this%LocalShapeList(DIM_DOMAIN, this%NProcess))

                read(10, *) this%LocalShapeList

            close(10)

        end subroutine LoadDomainDAT


        subroutine DumpDomainHDF5(this)
            class(Domain), intent(inout) :: this
            integer(4) :: i, j

            if (0 == this%MyId) write(*, '(a)') "The domain is save as hdf5 file."

            call hdf5DomainDump%init(filename=this%IOName%FullName%str, mode='write', serial=.True.)
            call hdf5DomainDump%open()

            call hdf5DomainDump%write('GlobalShape', this%GlobalShape)
            call hdf5DomainDump%write('LocalShape', this%LocalShape)
            call hdf5DomainDump%write('SpaceStep', this%SpaceStep)

            call hdf5DomainDump%writeattr('/', 'MyId', this%MyId)
            call hdf5DomainDump%writeattr('/', 'master', this%master)
            call hdf5DomainDump%write('DomainIndex', this%DomainIndex)
            call hdf5DomainDump%writeattr('/', 'NProcess', this%NProcess)
            call hdf5DomainDump%write('ImageShape', this%ImageShape)
            call hdf5DomainDump%write('NeighborDomainIndex', this%NeighborDomainIndex)
            call hdf5DomainDump%write('NeighborImageId', this%NeighborImageId)
            call hdf5DomainDump%write('NeighborType', this%NeighborType)

            call hdf5DomainDump%write('CornerIndex', this%CornerIndex)
            call hdf5DomainDump%write('CornerCoordinateValue', this%CornerCoordinateValue)
            call hdf5DomainDump%write('CoordinateOrigin', this%CoordinateOrigin)

            call hdf5DomainDump%write('LocalShapeList', this%LocalShapeList)

            call hdf5DomainDump%close()

        end subroutine DumpDomainHDF5


        subroutine DumpDomainDAT(this)
            class(Domain), intent(inout) :: this

            if (0 == this%MyId) write(*, '(a)') "The domain is save as dat file."
            open(10, file=this%IOName%FullName%str)

                write(10, *) this%GlobalShape
                write(10, *) this%LocalShape
                write(10, *) this%SpaceStep

                write(10, *) this%MyId
                write(10, *) this%master
                write(10, *) this%DomainIndex
                write(10, *) this%NProcess
                write(10, *) this%ImageShape
                write(10, *) this%NeighborDomainIndex
                write(10, *) this%NeighborImageId
                write(10, *) this%NeighborType

                write(10, *) this%CornerIndex
                write(10, *) this%CornerCoordinateValue
                write(10, *) this%CoordinateOrigin

                write(10, *) this%LocalShapeList

            close(10)

        end subroutine DumpDomainDAT


        function getDomainIndex(rank, image_shape, domain_dim)
            integer(4) :: domain_dim ,shapesquare,square
            integer(4) :: getDomainIndex(domain_dim)
            integer(4), intent(in) :: rank 
            integer(4), intent(in) :: image_shape(domain_dim)
            integer(4) :: temp, temp_rank, temp2
            
            if (rank >= 0 .and. rank < product(image_shape)) then
                if (1 == domain_dim) then
                    getDomainIndex = [rank]
                
                else if (2 == domain_dim) then
                    getDomainIndex = [mod(rank, image_shape(1)), rank / image_shape(1)]

                else if (3==domain_dim)then
                    shapesqure=image_shape(1)*image_shape(2)
                    square=mod(rank,shapesquare)
                    getDomainIndex = [mod(square, image_shape(1)), square / image_shape(1),rank/shapesquare]
                else

                    write(*, *) "The domain_dim must be in [1, 2]."

                end if
            end if

            return
        end

        function getImageIndex(domain_index, image_shape, domain_dim)
            integer(4) :: domain_dim
            integer(4) :: getImageIndex
            integer(4), intent(in) :: domain_index(domain_dim)
            integer(4), intent(in) :: image_shape(domain_dim)
            integer(4) :: i

            if (1 == domain_dim) then
                getImageIndex = domain_index(1)

            else if (2 == domain_dim) then
                getImageIndex = domain_index(1) + domain_index(2) * image_shape(1)

            else
                write(*, *) "The domain_dim must be in [1, 2]."

            end if

            return
        end

end Module ModuleDomain3d
