Module ModuleDomain
    use mpi  
    implicit none

    integer(4), parameter :: DIM_DOMAIN = 2
    integer(4), parameter :: BOUNDARY_PER_DIM = 2
    integer(4), parameter :: BOUNDARY_LEFT = 1
    integer(4), parameter :: BOUNDARY_RIGHT = 2
    integer(4), parameter :: BUFFER_SIZE = BOUNDARY_PER_DIM * DIM_DOMAIN
    integer(4), parameter :: NEIGHBOR_TYPE_BOUNDARY = 11
    integer(4), parameter :: NEIGHBOR_TYPE_DOMAIN = 12
    integer(4), parameter :: NEIGHBOR_TYPE_MID = 13

    integer(4), parameter :: deps_type_uniform = 0
    integer(4), parameter :: deps_type_specify = 1
    integer(4), parameter :: deps_type_balance = 2

    integer(4), save :: deps_type = deps_type_uniform

    type Domain
        !type(FileName)  :: IOName
        integer(4)  :: GlobalShape(DIM_DOMAIN)
        integer(4)  :: LocalShape(DIM_DOMAIN)
        real(8)     :: SpaceStep(DIM_DOMAIN)

        integer(4)  :: MyId = 0, master = 1
        integer(4)  :: DomainIndex(DIM_DOMAIN)
        integer(4)  :: NProcess = 0
        integer(4)  :: ImageShape(DIM_DOMAIN)
        integer(4)  :: NeighborDomainIndex(BOUNDARY_PER_DIM, DIM_DOMAIN, DIM_DOMAIN)
        integer(4)  :: NeighborImageId(BOUNDARY_PER_DIM, DIM_DOMAIN)
        integer(4)  :: NeighborType(BOUNDARY_PER_DIM, DIM_DOMAIN)

        integer(4)  :: CornerIndex(BOUNDARY_PER_DIM, DIM_DOMAIN)
        real(8)     :: CoordinateOrigin(DIM_DOMAIN)
        real(8)     :: CornerCoordinateValue(BOUNDARY_PER_DIM, DIM_DOMAIN)

        integer(4), allocatable :: LocalShapeList(:, :)

    contains

        procedure :: Init => InitDomain2DZR
        procedure :: Destroy => DestroyDomain

        procedure :: Uniform => DepsUniform
        procedure :: Specify => DepsSpecify
        procedure :: Balance => DepsBalance

       ! procedure :: Dump => DumpDomain
        !procedure :: Load => LoadDomain

        procedure :: ShowAll => ShowDomainAllInfo
        procedure :: Show => ShowDomain
       ! procedure :: Save => SaveDomainDepsInfo

        procedure, private :: SetLocalShape

        !procedure, private :: LoadDomainHDF5
        !procedure, private :: LoadDomainDAT

        !procedure, private :: DumpDomainHDF5
        !procedure, private :: DumpDomainDAT

    end type Domain

    !type(HDF5_PDump), save, private :: hdf5DomainDump

    contains

        subroutine InitDomain2DZR(this, Nz, Nr, z_start, z_end, r_start, r_end, pz, pr, dmname)
            class(Domain), intent(inout) :: this
            integer(4), intent(in) :: Nz, Nr
            real(8), intent(in) :: z_start, z_end, r_start, r_end
            integer(4), intent(in) :: pz, pr
            character(*), optional, intent(in) :: dmname
            integer(4) :: size, rank, ierr

            call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

            if (Nz > 1 .and. Nr > 1) then
                this%GlobalShape = [Nz, Nr]
                this%CoordinateOrigin = [z_start, r_start]
                this%SpaceStep = [z_end - z_start, r_end - r_start] / (this%GlobalShape - 1)

                this%MyId = rank
                this%NProcess = size

                if (pz * pr == this%NProcess) then
                    this%ImageShape = [pz, pr]
                else
                    this%ImageShape = [1, this%NProcess]
                end if

                this%DomainIndex = getDomainIndex(this%MyId, this%ImageShape, DIM_DOMAIN)

              !  if (present(dmname)) then
               !     call this%IOName%Init(dmname, RESTART_FILE_NAME)
               ! else
               !     call this%IOName%Init("Domain", RESTART_FILE_NAME)
               ! end if
                
                call this%Uniform()
            else
                if (0 == this%MyId) write(*, '(a)') "The Domain init failed, please check input parameter."
                stop
            end if

        end subroutine InitDomain2DZR


        subroutine DepsUniform(this)
            class(Domain), intent(inout) :: this
            integer(4) :: i, j

            ! set list of image size
            call this%Destroy()
            allocate(this%LocalShapeList(DIM_DOMAIN, this%NProcess))

            this%LocalShapeList = 0
            do i = 1, DIM_DOMAIN
                do j = 1, this%ImageShape(i)
                    this%LocalShapeList(i, j) = int((this%GlobalShape(i)-1) / this%ImageShape(i)) + 1
                    if (j <= (this%GlobalShape(i)-1 - int((this%GlobalShape(i)-1) / this%ImageShape(i)) * this%ImageShape(i))) &
                        this%LocalShapeList(i, j) = this%LocalShapeList(i, j) + 1
                end do
            end do
            
            call this%SetLocalShape()

            ! mid-axis
            if (this%DomainIndex(2) == 0) then
                this%NeighborImageId(BOUNDARY_LEFT, 2) = -1
                this%NeighborType(BOUNDARY_LEFT, 2) = NEIGHBOR_TYPE_MID
            end if

           ! call this%Save()

        end subroutine DepsUniform


        subroutine DepsSpecify(this, deps_z, deps_r)
            class(Domain), intent(inout) :: this
            real(8), intent(in) :: deps_z(this%ImageShape(1)), deps_r(this%ImageShape(2))
            integer(4) :: dp_z(this%ImageShape(1)), dp_r(this%ImageShape(2))
            integer(4) :: i, j
            real(8) :: tmp_sum

            ! norm
            if (this%ImageShape(1) == 1) then
                dp_z(1) = this%GlobalShape(1)

            else
                tmp_sum = sum(deps_z)
                do j = 1, this%ImageShape(1)
                    dp_z(j) = int(deps_z(j) / tmp_sum * dble(this%GlobalShape(1)-1) + 0.5d0) + 1
                end do

                if (sum(dp_z-1) /= this%GlobalShape(1)-1) then
                    dp_z(this%ImageShape(1)) = this%GlobalShape(1) - sum(dp_z(1:this%ImageShape(1)-1)-1)
                end if
            end if

            if (this%ImageShape(2) == 1) then
                dp_r(1) = this%GlobalShape(2)

            else
                tmp_sum = sum(deps_r)
                do j = 1, this%ImageShape(2)
                    dp_r(j) = int(deps_r(j) / tmp_sum * dble(this%GlobalShape(2)-1) + 0.5d0) + 1
                end do

                if (sum(dp_r-1) /= this%GlobalShape(2)-1) then
                    dp_r(this%ImageShape(2)) = this%GlobalShape(2) - sum(dp_r(1:this%ImageShape(2)-1)-1)
                end if
            end if

            ! set list of image size
            call this%Destroy()
            allocate(this%LocalShapeList(DIM_DOMAIN, this%NProcess))

            this%LocalShapeList = 0
            this%LocalShapeList(1, 1:this%ImageShape(1)) = dp_z
            this%LocalShapeList(2, 1:this%ImageShape(2)) = dp_r

            call this%SetLocalShape()

            ! mid-axis
            if (this%DomainIndex(2) == 0) then
                this%NeighborImageId(BOUNDARY_LEFT, 2) = -1
                this%NeighborType(BOUNDARY_LEFT, 2) = NEIGHBOR_TYPE_MID
            end if

            !call this%Save()

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


       


      

        

        

       

        function getDomainIndex(rank, image_shape, domain_dim)
            integer(4) :: domain_dim
            integer(4) :: getDomainIndex(domain_dim)
            integer(4), intent(in) :: rank
            integer(4), intent(in) :: image_shape(domain_dim)
            integer(4) :: temp, temp_rank, temp2
            
            if (rank >= 0 .and. rank < product(image_shape)) then
                if (1 == domain_dim) then
                    getDomainIndex = [rank]
                
                else if (2 == domain_dim) then
                    getDomainIndex = [mod(rank, image_shape(1)), rank / image_shape(1)]

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

end Module ModuleDomain
