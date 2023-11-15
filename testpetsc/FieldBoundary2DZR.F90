Module ModuleFieldBoundary2DZR
    use mpi
    use ModuleGeometry
    use ModuleMaterials

    Implicit none

    Type FieldBoundary2DZR
        Integer(4) :: Nz, Nr

        real(8), allocatable :: bound(:, :)
        integer(4), allocatable :: bound_mask(:, :)

    contains

        procedure :: Init    => InitilalizationFieldBounday
        procedure :: Destroy => DestroyFieldBounday
        procedure :: Zero    => ZeroFieldBounday

        procedure :: One     => UpdateOneFieldBounday
        procedure :: Update  => UpdateFieldBounday

        procedure :: intep   => IntepFieldBounday

        procedure :: diag    => diagFieldBounday

    end Type FieldBoundary2DZR

    contains

        subroutine InitilalizationFieldBounday(FB, Nz, Nr)
            Class(FieldBoundary2DZR), intent(inout) :: FB
            Integer(4), intent(in) :: Nz, Nr
            integer(4) :: size, rank, ierr

            call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

            call FB%Destroy()

            if (Nz > 1 .and. Nr > 1) then

                FB%Nz = Nz
                FB%Nr = Nr

                Allocate(FB%bound(Nz, Nr))
                Allocate(FB%bound_mask(Nz, Nr))

                call FB%Zero()

            else
                if (0 == rank) then
                    write(*, '(a)') "Nz and Nz must be greater than 0, please check input parameter."
                    stop
                end if
                
            end if

        end subroutine InitilalizationFieldBounday


        subroutine DestroyFieldBounday(FB)
            Class(FieldBoundary2DZR), intent(inout) :: FB

            if (Allocated(FB%bound)) Deallocate(FB%bound)
            if (Allocated(FB%bound_mask)) Deallocate(FB%bound_mask)

        end subroutine DestroyFieldBounday


        subroutine ZeroFieldBounday(FB)
            Class(FieldBoundary2DZR), intent(inout) :: FB

            FB%bound = 0.d0
            FB%bound_mask = 0

        end subroutine ZeroFieldBounday


        subroutine UpdateOneFieldBounday(FB, MT, Geom, metal_index)
            Class(FieldBoundary2DZR), intent(inout) :: FB
            Type(Materials), intent(in) :: MT
            Type(Geometry), intent(in) :: Geom
            Integer(4), intent(in) :: metal_index

            if (metal_index > 0) then
                FB%bound = merge(MT%metls(metal_index)%voltage, FB%bound, Geom%point_type == metal_index)
                FB%bound_mask = merge(1, FB%bound_mask, Geom%point_type == metal_index)
            end if

        end subroutine UpdateOneFieldBounday


        subroutine UpdateFieldBounday(FB, MT, Geom)
            Class(FieldBoundary2DZR), intent(inout) :: FB
            Type(Materials), intent(inout) :: MT
            Type(Geometry), intent(in) :: Geom
            Integer(4) :: i

            if (MT%metal_count > 0) then
                do i = 1, MT%metal_count
                    FB%bound = merge(MT%metls(i)%voltage, FB%bound, Geom%point_type == i)
                    FB%bound_mask = merge(1, FB%bound_mask, Geom%point_type == i)
                end do
            end if

        end subroutine UpdateFieldBounday


        subroutine IntepFieldBounday(FB, MT, Geom)
            Class(FieldBoundary2DZR), intent(inout) :: FB
            Type(Materials), intent(in) :: MT
            Type(Geometry), intent(in) :: Geom
            integer(4) :: i
            integer(4), allocatable :: mask(:)

            ! R
            if (allocated(mask)) deallocate(mask)
            allocate(mask(FB%Nr))
            do i = 1, FB%Nr
                if (Geom%point_type(1, i) /= point_type_metal) then
                    mask(i) = 1
                end if
            end do
            call intep_log(FB%Nr, FB%bound(1, 1:FB%Nr), mask)
            FB%bound_mask(1, 1:FB%Nr) = 1

            mask = 0
            do i = 1, FB%Nr
                if (Geom%point_type(FB%Nz, i) /= point_type_metal) then
                    mask(i) = 1
                end if
            end do
            call intep_log(FB%Nr, FB%bound(FB%Nz, 1:FB%Nr), mask)
            FB%bound_mask(FB%Nz, 1:FB%Nr) = 1

            ! Z
            if (allocated(mask)) deallocate(mask)
            allocate(mask(FB%Nz))
            do i = 1, FB%Nz
                if (Geom%point_type(i, FB%Nr) /= point_type_metal) then
                    mask(i) = 1
                end if
            end do
            call intep_line(FB%Nz, FB%bound(1:FB%Nz, FB%Nr), mask)
            FB%bound_mask(1:FB%Nz, FB%Nr) = 1

            if (allocated(mask)) deallocate(mask)

        end subroutine IntepFieldBounday


        subroutine intep_log(data_length, data_1d, mask)
            integer(4), intent(in) :: data_length
            real(8), intent(inout) :: data_1d(data_length)        ! 待处理的一维数组
            integer(4), intent(in) :: mask(data_length)           ! 掩码数组, 尺寸与 data_1d 相同, 值0表示不需要处理, 值1表示需要插值
            real(8) :: vm, vn
            integer(4) :: intep_start, intep_stop, m, n, i, j

            if (sum(mask) /= 0) then
                intep_start = 0;
                intep_stop = 0;

                do i = 2, data_length
                    if (intep_start == 0 .and. mask(i) == 1) then
                        intep_start = i
                    end if
                    
                    if (intep_start /= 0 .and. mask(i) == 0) then
                        intep_stop = i-1
                        
                        m = intep_start-1
                        n = intep_stop+1
                        vm = data_1d(m)
                        vn = data_1d(n)
                        do j = intep_start, intep_stop
                            data_1d(j) = vm + log(dble(j)/dble(m))/log(dble(n)/dble(m))*(vn-vm)
                        end do
                        
                        intep_start = 0
                        intep_stop = 0
                        
                    end if
                end do
            end if

        end subroutine intep_log


        subroutine intep_line(data_length, data_1d, mask)
            integer(4), intent(in) :: data_length
            real(8), intent(inout) :: data_1d(data_length)        ! 待处理的一维数组
            integer(4), intent(in) :: mask(data_length)           ! 掩码数组, 尺寸与 data_1d 相同, 值0表示不需要处理, 值1表示需要插值
            real(8) :: vm, vn
            integer(4) :: intep_start, intep_stop, m, n, i, j

            if (sum(mask) /= 0) then
                intep_start = 0;
                intep_stop = 0;

                do i = 2, data_length
                    if (intep_start == 0 .and. mask(i) == 1) then
                        intep_start = i
                    end if
                    
                    if (intep_start /= 0 .and. mask(i) == 0) then
                        intep_stop = i-1
                        
                        m = intep_start-1
                        n = intep_stop+1
                        vm = data_1d(m)
                        vn = data_1d(n)
                        do j = intep_start, intep_stop
                            data_1d(j) = vm + (dble(j)/dble(m))/(dble(n)/dble(m))*(vn-vm)
                        end do
                        
                        intep_start = 0
                        intep_stop = 0
                        
                    end if
                end do
            end if

        end subroutine intep_line


        subroutine diagFieldBounday(FB, name)
            Class(FieldBoundary2DZR), intent(inout) :: FB
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

                call h5dumpcheck%writeattr('/', 'Nz', FB%Nz)
                call h5dumpcheck%writeattr('/', 'Nr', FB%Nr)

                call h5dumpcheck%write('bound', FB%bound)
                call h5dumpcheck%write('mask', FB%bound_mask)

                call h5dumpcheck%close()

            end if

        end subroutine diagFieldBounday

end Module ModuleFieldBoundary2DZR
