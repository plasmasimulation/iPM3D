module ModuleMaterials
    use mpi
    use ModuleFileName
    use ModuleParallelDump
    use ModuleGeometry

    implicit none

    ! 三种材料类型
    integer(4), parameter :: material_type_vacuum = 0
    integer(4), parameter :: material_type_metal = 1
    integer(4), parameter :: material_type_dielectric = 2

    ! 介质
    type DielectricOne
        real(8) :: permittivity = 0.d0
        real(8) :: permeability = 0.d0
        real(8) :: conductivity = 0.d0
    end type DielectricOne

    ! 金属
    type MetalOne
        real(8) :: charge       = 0.d0
        real(8) :: charge_one_step = 0.d0
        real(8) :: capacitance  = 0.d0
        real(8) :: voltage      = 0.d0
    end type MetalOne

    ! 材料
    type Materials
        type(FileName) :: IOName
        integer(4) :: dielectric_count                  ! 介质索引为负
        integer(4) :: metal_count                       ! 金属索引为正

        integer(4) :: Nz
        integer(4) :: Nr

        real(8), allocatable :: sigma(:, :)             ! 总的介质表面电荷密度
        real(8), allocatable :: sigma_one_step(:, :)    ! 每个步长累积的介质表面电荷密度
        type(DielectricOne), allocatable :: diels(:)    ! 介质数组
        type(MetalOne), allocatable :: metls(:)         ! 金属数组
        integer(4), allocatable :: material_type(:)     ! 材料类型
        character(99), allocatable :: material_name(:)  ! 材料名称

        integer(4) :: index_length
        integer(4), allocatable :: sigma_index(:, :)

    contains

        procedure :: init => initMaterials              ! 只负责创建数组, 具体参数需要额外设定
        procedure :: zero => zeroMaterials
        procedure :: dump => dumpMaterials
        procedure :: load => loadMaterials
        procedure :: destroy => destroyMaterials

        procedure :: updataIndex => updataMaterialsSigmaIndex

        procedure :: show => showMaterials

        procedure, private :: dumpMaterialsHDF5
        procedure, private :: dumpMaterialsDAT

        procedure, private :: loadMaterialsHDF5
        procedure, private :: loadMaterialsDAT

    end type Materials

    type(HDF5_PDump), private :: hdf5Materials

contains

    subroutine initMaterials(this, Nz, Nr, dielectric_num, metals_num, name)
        class(Materials), intent(inout) :: this
        integer(4), intent(in) :: Nz, Nr
        integer(4), intent(in) :: dielectric_num, metals_num
        character(*), optional, intent(in) :: name
        integer(4) :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        if (Nz > 0 .and. Nr > 0 .and. metals_num >= 0 .and. dielectric_num >= 0) then
            call this%destroy()

            this%Nz = Nz
            this%Nr = Nr
            this%dielectric_count = dielectric_num
            this%metal_count = metals_num

            if (present(name)) then
                call this%IOName%Init(name, RESTART_FILE_NAME)
            else
                call this%IOName%Init("Materials", RESTART_FILE_NAME)
            end if

            allocate(this%sigma(this%Nz, this%Nr))
            allocate(this%sigma_one_step(this%Nz, this%Nr))
            allocate(this%material_type(-this%dielectric_count:this%metal_count))
            allocate(this%material_name(-this%dielectric_count:this%metal_count))

            if (this%dielectric_count > 0) allocate(this%diels(-this%dielectric_count:-1))
            if (this%metal_count > 0) allocate(this%metls(1:this%metal_count))

            call this%zero()
        else
            if (0 == rank) write(*, *) "Invalid input parameters."

        end if

    end subroutine initMaterials


    subroutine zeroMaterials(this)
        class(Materials), intent(inout) :: this
        integer(4) :: i

        if (allocated(this%sigma)) then
            this%sigma = 0.d0
            this%sigma_one_step = 0.d0
        end if

        if (allocated(this%metls) .and. this%metal_count > 0) then
            do i = 1, this%metal_count
                this%metls(i)%charge = 0.d0
                this%metls(i)%charge_one_step = 0.d0
                this%metls(i)%capacitance = 0.d0
                this%metls(i)%voltage = 0.d0

            end do
        end if

    end subroutine zeroMaterials


    subroutine showMaterials(this)
        class(Materials), intent(inout) :: this
        integer(4) :: i, size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        if (0 == rank .and. (allocated(this%material_type))) then
            write(*, *) ''
            write(*, '(a, i8, a, i8, a)') "Materials information: ", this%metal_count, " metals, ", this%dielectric_count, " media"

            do i = -this%dielectric_count, this%metal_count

                if (this%material_type(i) == material_type_dielectric) then
                    write(*, '(i, 2x, a, 2x, a, es12.4, es12.4, es12.4)') i, 'dielectric', &
                                trim(this%material_name(i)), &
                                this%diels(i)%permittivity, &
                                this%diels(i)%permeability, &
                                this%diels(i)%conductivity

                else if (this%material_type(i) == material_type_metal) then
                    write(*, '(i, 2x, a, 2x, a, es12.4, es12.4, es12.4)') i, 'metal', &
                                trim(this%material_name(i)), &
                                this%metls(i)%charge, &
                                this%metls(i)%capacitance, &
                                this%metls(i)%voltage

                else
                    write(*, '(i, 2x, a, 2x, a)') i, 'vacuum', &
                                trim(this%material_name(i))

                end if
            end do
        end if

    end subroutine showMaterials


    subroutine dumpMaterials(this)
        class(Materials), intent(inout) :: this

        if (FILE_EXTENSION_MODE_H5 == this%IOName%ExtensionMode) then
            call this%dumpMaterialsHDF5()
        
        else if (FILE_EXTENSION_MODE_DAT == this%IOName%ExtensionMode) then
            call this%dumpMaterialsDAT()

        end if

    end subroutine dumpMaterials


    subroutine dumpMaterialsHDF5(this)
        class(Materials), intent(inout) :: this
        integer(4) :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        if (0 == rank) write(*, '(a)') "The Materials is save as hdf5 file."

        call hdf5Materials%init(filename=this%IOName%FullName%str, mode='write', serial=.True.)
        call hdf5Materials%open()

        call hdf5Materials%write('sigma', this%sigma)

        if (this%metal_count > 0) then
            call hdf5Materials%write('charge', this%metls%charge)
            call hdf5Materials%write('capacitance', this%metls%capacitance)
            call hdf5Materials%write('voltage', this%metls%voltage)
        end if

        call hdf5Materials%close()

    end subroutine dumpMaterialsHDF5


    subroutine dumpMaterialsDAT(this)
        class(Materials), intent(inout) :: this
        integer(4) :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        if (allocated(this%sigma)) then
            if (0 == rank) write(*, '(a)') "The Materials is save as dat file."

            open(10, file=this%IOName%FullName%str)
                write(10, *) this%sigma

                if (this%metal_count > 0) then
                    write(10, *) this%metal_count
                    write(10, *) this%metls(1:this%metal_count)%charge
                    write(10, *) this%metls(1:this%metal_count)%capacitance
                    write(10, *) this%metls(1:this%metal_count)%voltage
                end if
            close(10)
        end if

    end subroutine dumpMaterialsDAT


    subroutine loadMaterials(this)
        class(Materials), intent(inout) :: this
        logical :: alive

        inquire(file=this%IOName%FullName%str, exist=alive)
        if (alive) then
            if (FILE_EXTENSION_MODE_H5 == this%IOName%ExtensionMode) then
                call this%loadMaterialsHDF5()

            else if (FILE_EXTENSION_MODE_DAT == this%IOName%ExtensionMode) then
                call this%loadMaterialsDAT()

            end if

        end if

    end subroutine loadMaterials


    subroutine loadMaterialsHDF5(this)
        class(Materials), intent(inout) :: this
        real(8), allocatable :: Temp2D(:, :)
        real(8), allocatable :: Temp1D_1(:), Temp1D_2(:), Temp1D_3(:)
        integer(4) :: i, count
        integer(4) :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        if (0 == rank) write(*, '(a)') "The Materials is load from hdf5 file."

        call hdf5Materials%init(filename=this%IOName%FullName%str, mode='read')
        call hdf5Materials%open()

        allocate(Temp2D, mold=this%sigma)
        call hdf5Materials%read('sigma',Temp2D)
        this%sigma=Temp2D

        if (this%metal_count > 0) then
            allocate(Temp1D_1(this%metal_count))
            allocate(Temp1D_2(this%metal_count))
            allocate(Temp1D_3(this%metal_count))
            call hdf5Materials%read('charge',Temp1D_1)
            call hdf5Materials%read('capacitance',Temp1D_2)
            call hdf5Materials%read('voltage',Temp1D_3)

            do i = 1, this%metal_count
                this%metls(i)%charge = Temp1D_1(i)
                this%metls(i)%capacitance = Temp1D_2(i)
                this%metls(i)%voltage = Temp1D_3(i)
            end do
        end if

        call hdf5Materials%close()

        if (allocated(Temp2D)) deallocate(Temp2D)
        if (allocated(Temp1D_1)) deallocate(Temp1D_1)
        if (allocated(Temp1D_2)) deallocate(Temp1D_2)
        if (allocated(Temp1D_3)) deallocate(Temp1D_3)

    end subroutine loadMaterialsHDF5


    subroutine loadMaterialsDAT(this)
        class(Materials), intent(inout) :: this
        integer(4) :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        if (0 == rank) write(*, '(a)') "The Materials is load from dat file."

        open(10, file=this%IOName%FullName%str)
            read(10, *) this%sigma
            read(10, *) this%metal_count
            if (this%metal_count > 0) then
                read(10, *) this%metls(1:this%metal_count)%charge
                read(10, *) this%metls(1:this%metal_count)%capacitance
                read(10, *) this%metls(1:this%metal_count)%voltage
            end if
        close(10)

    end subroutine loadMaterialsDAT


    subroutine destroyMaterials(this)
        class(Materials), intent(inout) :: this

        this%dielectric_count = 0
        this%metal_count = 0
        if (allocated(this%sigma)) deallocate(this%sigma)
        if (allocated(this%sigma_one_step)) deallocate(this%sigma_one_step)
        if (allocated(this%diels)) deallocate(this%diels)
        if (allocated(this%metls)) deallocate(this%metls)
        if (allocated(this%material_type)) deallocate(this%material_type)
        if (allocated(this%material_name)) deallocate(this%material_name)

        this%index_length = 0
        if (allocated(this%sigma_index)) deallocate(this%sigma_index)

    end subroutine destroyMaterials


    subroutine updataMaterialsSigmaIndex(this, Geom)
        class(Materials), intent(inout) :: this
        type(Geometry), intent(in) :: Geom
        integer(4) :: i, j, tmp

        if (this%dielectric_count > 0) then
            this%index_length = 0

            do i = 2, this%Nz-1
                do j = 2, this%Nr-1
                    if (Geom%hedge_type(i, j-1) == edge_type_dielectric_to_vacuum .or. &
                        Geom%hedge_type(i, j-1) == edge_type_vacuum_to_dielectric .or. &
                        Geom%hedge_type(i, j) == edge_type_dielectric_to_vacuum .or. &
                        Geom%hedge_type(i, j) == edge_type_vacuum_to_dielectric .or. &
                        Geom%vedge_type(i-1, j) == edge_type_dielectric_to_vacuum .or. &
                        Geom%vedge_type(i-1, j) == edge_type_vacuum_to_dielectric .or. &
                        Geom%vedge_type(i, j) == edge_type_dielectric_to_vacuum .or. &
                        Geom%vedge_type(i, j) == edge_type_vacuum_to_dielectric) then

                        this%index_length = this%index_length + 1

                    end if
                end do
            end do

            if (this%index_length > 0) then
                if (allocated(this%sigma_index)) deallocate(this%sigma_index)
                allocate(this%sigma_index(2, this%index_length))

                tmp = 0
                do i = 2, this%Nz-1
                    do j = 2, this%Nr-1
                        if (Geom%hedge_type(i, j-1) == edge_type_dielectric_to_vacuum .or. &
                            Geom%hedge_type(i, j-1) == edge_type_vacuum_to_dielectric .or. &
                            Geom%hedge_type(i, j) == edge_type_dielectric_to_vacuum .or. &
                            Geom%hedge_type(i, j) == edge_type_vacuum_to_dielectric .or. &
                            Geom%vedge_type(i-1, j) == edge_type_dielectric_to_vacuum .or. &
                            Geom%vedge_type(i-1, j) == edge_type_vacuum_to_dielectric .or. &
                            Geom%vedge_type(i, j) == edge_type_dielectric_to_vacuum .or. &
                            Geom%vedge_type(i, j) == edge_type_vacuum_to_dielectric) then

                            tmp = tmp + 1
                            this%sigma_index(1, tmp) = i
                            this%sigma_index(2, tmp) = j

                        end if
                    end do
                end do

            end if
        else
            this%index_length = 0
            if (allocated(this%sigma_index)) deallocate(this%sigma_index)

        end if

    end subroutine updataMaterialsSigmaIndex
    
end module ModuleMaterials
