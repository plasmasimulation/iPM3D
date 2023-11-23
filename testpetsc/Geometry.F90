module ModuleGeometry
    use mpi
    use ModuleFileName
    use ModuleParallelDump

    implicit none

    ! 格点类型
    integer(4), parameter :: point_type_dielectric_vacuum_boundary = -2
    integer(4), parameter :: point_type_dielectric_inner = -1
    integer(4), parameter :: point_type_vacuum = 0
    integer(4), parameter :: point_type_metal = 1                           ! 大于0的都是金属

    ! 边类型，从坐标值小的方向指向值大的方向
    integer(4), parameter :: edge_type_vacuum_to_metal = 0
    integer(4), parameter :: edge_type_vacuum_to_dielectric = 1
    integer(4), parameter :: edge_type_metal_to_vacuum = 2
    integer(4), parameter :: edge_type_dielectric_to_vacuum = 3
    integer(4), parameter :: edge_type_other = 4

    type Geometry
        integer(4) :: Nz
        integer(4) :: Nr

        integer(4), allocatable :: cell_material(:, :)      ! 每个网格的材料索引
        real(8), allocatable :: cell_epsilon(:, :)       ! 每个网格的相对介电常数

        integer(4), allocatable :: point_type(:, :)         ! 每个格点的类型
        real(8), allocatable :: point_volume(:, :)          ! 每个格点的有效体积 假设dz=dr
        real(8), allocatable :: point_inv_volume(:, :)      ! 每个格点的有效体积的倒数 假设dz=dr
        real(8), allocatable :: point_area(:, :)            ! 每个格点的表面积 假设dz=dr

        integer(4), allocatable :: hedge_type(:, :)         ! 水平方向的边类型
        integer(4), allocatable :: vedge_type(:, :)         ! 垂直方向的边类型

    contains

        procedure :: init => initGeometry                   ! 只负责创建数组, 具体参数需要额外设定
        procedure :: check => checkGeometry
        procedure :: destroy => destroyGeometry

    end type Geometry

    contains

        subroutine initGeometry(this, Nz, Nr)
            class(Geometry), intent(inout) :: this
            integer, intent(in) :: Nz, Nr
            integer(4) :: size, rank, ierr

            call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

            if (Nz > 1 .and. Nr > 1) then
                call this%destroy()
                this%Nz = Nz
                this%Nr = Nr
                allocate(this%cell_material(this%Nz-1, this%Nr-1))
                allocate(this%cell_epsilon(this%Nz-1, this%Nr-1))

                allocate(this%point_type(this%Nz, this%Nr))
                allocate(this%point_volume(this%Nz, this%Nr))
                allocate(this%point_inv_volume(this%Nz, this%Nr))
                allocate(this%point_area(this%Nz, this%Nr))

                allocate(this%hedge_type(this%Nz, this%Nr-1))
                allocate(this%vedge_type(this%Nz-1, this%Nr))

            else
                if (0 == rank) write(*, *) "Invalid input parameters."

            end if

        end subroutine initGeometry


        subroutine checkGeometry(this, name)
            class(Geometry), intent(inout) :: this
            character(*), intent(in) :: name
            type(FileName) :: IOName
            type(HDF5_PDump) :: h5dumpcheck
            integer(4) :: size, rank, ierr

            call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

            if (0 == rank) then
                call IOName%Init(name, &
                                PathMode=FILE_PATH_MODE_CHECK, &
                                ExtensionMode=FILE_EXTENSION_MODE_H5, &
                                ParallelMode=FILE_PARALLEL_MODE_CLOSE, &
                                DynamicIndex=FILE_DYNAMIC_MODE_CLOSE)

                call h5dumpcheck%init(filename=IOName%FullName%str, mode='write', serial=.True.)
                call h5dumpcheck%open()

                call h5dumpcheck%writeattr('/', 'Nz', this%Nz)
                call h5dumpcheck%writeattr('/', 'Nr', this%Nr)

                call h5dumpcheck%write('cell_material', this%cell_material)
                call h5dumpcheck%write('cell_epsilon', this%cell_epsilon)

                call h5dumpcheck%write('point_type', this%point_type)
                call h5dumpcheck%write('point_volume', this%point_volume)
                call h5dumpcheck%write('point_inv_volume', this%point_inv_volume)
                call h5dumpcheck%write('point_area', this%point_area)

                call h5dumpcheck%write('hedge_type', this%hedge_type)
                call h5dumpcheck%write('vedge_type', this%vedge_type)

                call h5dumpcheck%close()

            end if

        end subroutine checkGeometry


        subroutine destroyGeometry(this)
            class(Geometry), intent(inout) :: this
        
            if (allocated(this%cell_material)) deallocate(this%cell_material)
            if (allocated(this%cell_epsilon)) deallocate(this%cell_epsilon)

            if (allocated(this%point_type)) deallocate(this%point_type)
            if (allocated(this%point_volume)) deallocate(this%point_volume)
            if (allocated(this%point_inv_volume)) deallocate(this%point_inv_volume)
            if (allocated(this%point_area)) deallocate(this%point_area)

            if (allocated(this%hedge_type)) deallocate(this%hedge_type)
            if (allocated(this%vedge_type)) deallocate(this%vedge_type)

        end subroutine destroyGeometry

end module ModuleGeometry
