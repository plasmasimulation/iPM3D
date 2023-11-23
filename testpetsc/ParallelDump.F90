module ModuleParallelDump
    use mpi
    use HDF5
    implicit none
    type DomainDump
        integer(HSIZE_T), allocatable :: count(:)
        integer(HSSIZE_T), allocatable :: offset(:)
        integer, allocatable :: location(:)
    contains
        procedure,private :: init => init_domain
        procedure,private :: calc_loc => calculate_domain_location
        procedure,private :: calc_off => calculate_domain_offset
        final destroy
    end type

    type :: HDF5_PDump
        character(len=99) :: filename
        type(DomainDump) :: thisdomain
        integer, allocatable :: chunk(:)
        INTEGER(HID_T) :: file_id = 0
        integer :: chunkrank
		logical :: serialmode = .False.
		character(len=99) :: openmode = 'write'
    contains
        procedure init => init_HDF5_PDump
        procedure open => open_HDF5_PDump

        procedure,private :: setchunk => set_parallel_chunk
        procedure mkdir => mkdir_HDF5_PDump
        procedure close => close_HDF5_PDump

        procedure,private :: write_attr_integer, write_attr_real, write_attr_double, write_attr_character, write_attr_logical

		procedure,private :: read_attr_integer, read_attr_real, read_attr_double, read_attr_character, read_attr_logical

        procedure,private :: write_integer_1D, write_integer_2D, write_integer_3D, write_integer_4D, write_integer_5D, write_integer_6D, write_integer_7D, &
                write_real_1D, write_real_2D, write_real_3D, write_real_4D, write_real_5D, write_real_6D, write_real_7D, &
                write_double_1D, write_double_2D, write_double_3D, write_double_4D, write_double_5D, write_double_6D, write_double_7D 

		procedure,private :: read_integer_1D, read_integer_2D, read_integer_3D, read_integer_4D, read_integer_5D, read_integer_6D, read_integer_7D, &
				read_real_1D, read_real_2D, read_real_3D, read_real_4D, read_real_5D, read_real_6D, read_real_7D, &
				read_double_1D, read_double_2D, read_double_3D, read_double_4D, read_double_5D, read_double_6D, read_double_7D

        generic :: write => write_integer_1D, write_integer_2D, write_integer_3D, write_integer_4D, write_integer_5D, write_integer_6D, write_integer_7D, &
                write_real_1D, write_real_2D, write_real_3D, write_real_4D, write_real_5D, write_real_6D, write_real_7D, &
                write_double_1D, write_double_2D, write_double_3D, write_double_4D, write_double_5D, write_double_6D, write_double_7D 

        generic :: writeattr => write_attr_integer, write_attr_real, write_attr_double, write_attr_character, write_attr_logical

		generic :: read => read_integer_1D, read_integer_2D, read_integer_3D, read_integer_4D, read_integer_5D, read_integer_6D, read_integer_7D, &
				read_real_1D, read_real_2D, read_real_3D, read_real_4D, read_real_5D, read_real_6D, read_real_7D, &
				read_double_1D, read_double_2D, read_double_3D, read_double_4D, read_double_5D, read_double_6D, read_double_7D

		generic :: readattr => read_attr_integer, read_attr_real, read_attr_double, read_attr_character, read_attr_logical

    end type

    integer, private :: error
	logical, private :: WARN_MESSAGE = .FALSE.
contains

    subroutine init_HDF5_PDump(thisHDF5file, filename, mode, serial)
        implicit none
        class(HDF5_PDump), intent(inout) :: thisHDF5file
        character(len=*), intent(in) :: filename  ! File name
		character(len=*), intent(in) :: mode
		logical, intent(in), optional :: serial
		if (present(serial) .AND. serial == .True.) then 
			thisHDF5file%serialmode = .True.
			if(WARN_MESSAGE) write(*,*) "=== Warning: Serial mode is ON so that chunk will be ignored ==="
		end if
		thisHDF5file%openmode=TRIM(mode)
        thisHDF5file%filename = TRIM(filename)
        CALL h5open_f(error)

        return
    end subroutine init_HDF5_PDump

    subroutine open_HDF5_PDump(thisHDF5file)
        implicit none
        class(HDF5_PDump), intent(inout) :: thisHDF5file
        INTEGER(HID_T) :: plist_id
		logical :: exists=.False.
        integer :: size_image, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size_image, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
		select case(thisHDF5file%openmode)
		case('write')
			if(thisHDF5file%serialmode == .False.) then
				CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
			end if
			CALL h5pset_fclose_degree_f(plist_id, H5F_CLOSE_SEMI_F, error)
        	CALL h5fcreate_f(thisHDF5file%filename, H5F_ACC_TRUNC_F, thisHDF5file%file_id, error, access_prp=plist_id)
		case('read')
			! CALL h5pset_fclose_degree_f(plist_id, H5F_CLOSE_SEMI_F, error)
			INQUIRE(file=thisHDF5file%filename,exist=exists)
			if(exists) then
				CALL h5fopen_f(thisHDF5file%filename, H5F_ACC_RDONLY_F, thisHDF5file%file_id, error, access_prp=plist_id)
			else
				write(*,*) "=== Error: Image <",rank,"> report that file ",thisHDF5file%filename," does not exist ==="
				stop
			end if
		case default
			write(*,*) "=== Error: Image <",rank,"> report that must have a valid open mode ==="
			stop
		end select
        CALL h5pclose_f(plist_id, error)

        return
    end subroutine open_HDF5_PDump

    subroutine set_parallel_chunk(thisHDF5file, chunkdim)
        implicit none
        class(HDF5_PDump), intent(inout) :: thisHDF5file
        integer, intent(in)  :: chunkdim(:)
        integer :: i, j, k
        integer :: size_image, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size_image, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        if (product(chunkdim) /= size_image) then
            if (0 == rank) write (*, *) "=== Error: Number of Chunks Don't Match Images ==="
            stop
        end if

        thisHDF5file%chunk = chunkdim
        thisHDF5file%chunkrank = size(chunkdim)

        call thisHDF5file%thisdomain%init(thisHDF5file%chunkrank)
        call thisHDF5file%thisdomain%calc_loc(thisHDF5file%chunk)

        return
    end subroutine set_parallel_chunk

    subroutine mkdir_HDF5_PDump(thisHDF5file, groupname)
        implicit none
        class(HDF5_PDump), intent(inout) :: thisHDF5file
        CHARACTER(LEN=*), intent(in) :: groupname
        INTEGER(HID_T) :: group_id
        integer :: slash=0,cursor=1
        logical :: exists
        do
            slash=index(groupname(cursor:),'/')
            if (slash==0) exit
            cursor=cursor+slash
            call h5lexists_f(thisHDF5file%file_id, groupname(:cursor-1), exists, error)
            if (.not.exists) then
                CALL h5gcreate_f(thisHDF5file%file_id, groupname(:cursor-1), group_id, error)
                CALL h5gclose_f(group_id, error)
            end if
        end do
        slash=0;cursor=1
        return
    end subroutine mkdir_HDF5_PDump
    
	subroutine read_integer_1D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		integer(4), dimension(:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(1) :: dims
		integer :: datarank=1
		dtype = H5T_NATIVE_INTEGER
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1)))
			end if
		else
			allocate(data(dims(1)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_integer_1D

	subroutine read_integer_2D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		integer(4), dimension(:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(2) :: dims
		integer :: datarank=2
		dtype = H5T_NATIVE_INTEGER
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2)))
			end if
		else
			allocate(data(dims(1),dims(2)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_integer_2D

	subroutine read_integer_3D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		integer(4), dimension(:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(3) :: dims
		integer :: datarank=3
		dtype = H5T_NATIVE_INTEGER
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_integer_3D

	subroutine read_integer_4D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		integer(4), dimension(:,:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(4) :: dims
		integer :: datarank=4
		dtype = H5T_NATIVE_INTEGER
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3),dims(4)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3),dims(4)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_integer_4D

	subroutine read_integer_5D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		integer(4), dimension(:,:,:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(5) :: dims
		integer :: datarank=5
		dtype = H5T_NATIVE_INTEGER
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_integer_5D

	subroutine read_integer_6D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		integer(4), dimension(:,:,:,:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(6) :: dims
		integer :: datarank=6
		dtype = H5T_NATIVE_INTEGER
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_integer_6D

	subroutine read_integer_7D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		integer(4), dimension(:,:,:,:,:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(7) :: dims
		integer :: datarank=7
		dtype = H5T_NATIVE_INTEGER
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_integer_7D

	subroutine read_real_1D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(4), dimension(:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(1) :: dims
		integer :: datarank=1
		dtype = H5T_NATIVE_REAL
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1)))
			end if
		else
			allocate(data(dims(1)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_real_1D

	subroutine read_real_2D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(4), dimension(:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(2) :: dims
		integer :: datarank=2
		dtype = H5T_NATIVE_REAL
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2)))
			end if
		else
			allocate(data(dims(1),dims(2)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_real_2D

	subroutine read_real_3D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(4), dimension(:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(3) :: dims
		integer :: datarank=3
		dtype = H5T_NATIVE_REAL
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_real_3D

	subroutine read_real_4D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(4), dimension(:,:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(4) :: dims
		integer :: datarank=4
		dtype = H5T_NATIVE_REAL
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3),dims(4)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3),dims(4)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_real_4D

	subroutine read_real_5D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(4), dimension(:,:,:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(5) :: dims
		integer :: datarank=5
		dtype = H5T_NATIVE_REAL
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_real_5D

	subroutine read_real_6D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(4), dimension(:,:,:,:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(6) :: dims
		integer :: datarank=6
		dtype = H5T_NATIVE_REAL
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_real_6D

	subroutine read_real_7D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(4), dimension(:,:,:,:,:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(7) :: dims
		integer :: datarank=7
		dtype = H5T_NATIVE_REAL
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_real_7D

	subroutine read_double_1D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(8), dimension(:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(1) :: dims
		integer :: datarank=1
		dtype = H5T_NATIVE_DOUBLE
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1)))
			end if
		else
			allocate(data(dims(1)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_double_1D

	subroutine read_double_2D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(8), dimension(:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(2) :: dims
		integer :: datarank=2
		dtype = H5T_NATIVE_DOUBLE
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2)))
			end if
		else
			allocate(data(dims(1),dims(2)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_double_2D

	subroutine read_double_3D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(8), dimension(:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(3) :: dims
		integer :: datarank=3
		dtype = H5T_NATIVE_DOUBLE
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_double_3D

	subroutine read_double_4D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(8), dimension(:,:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(4) :: dims
		integer :: datarank=4
		dtype = H5T_NATIVE_DOUBLE
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3),dims(4)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3),dims(4)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_double_4D

	subroutine read_double_5D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(8), dimension(:,:,:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(5) :: dims
		integer :: datarank=5
		dtype = H5T_NATIVE_DOUBLE
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_double_5D

	subroutine read_double_6D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(8), dimension(:,:,:,:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(6) :: dims
		integer :: datarank=6
		dtype = H5T_NATIVE_DOUBLE
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_double_6D

	subroutine read_double_7D(thisHDF5file, dsetname, data)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(8), dimension(:,:,:,:,:,:,:), allocatable, intent(inout) :: data
		INTEGER(HID_T) :: dset_id, dspace_id, plist_id, dtype
		INTEGER(HSIZE_T), DIMENSION(7) :: dims
		integer :: datarank=7
		dtype = H5T_NATIVE_DOUBLE
		call preread(thisHDF5file, dsetname, dset_id, dspace_id, dims, datarank)
		if(allocated(data)) then
			if(any(shape(data) /= dims)) then
				deallocate(data)
				allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))
			end if
		else
			allocate(data(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))
		end if
		CALL h5dread_f(dset_id,dtype,data,dims,error)
		CALL postread(thisHDF5file,dspace_id,dset_id)
	end subroutine read_double_7D

    subroutine print_chunkdim_present_error_info()
		integer :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        write(*,"(a,i0,a)") "=== Error: Image <",rank,"> report that the serial mode is FALSE but chundim is not given ==="

    end subroutine print_chunkdim_present_error_info

	subroutine write_integer_1D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		integer(4), dimension(:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(1) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(1) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(1) :: counttemp, fulldataShape
		integer(4), dimension(:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_INTEGER

        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=1
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_integer_1D

	subroutine write_integer_2D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		integer(4), dimension(:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(2) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(2) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(2) :: counttemp, fulldataShape
		integer(4), dimension(:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_INTEGER
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=2
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_integer_2D

	subroutine write_integer_3D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		integer(4), dimension(:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(3) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(3) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(3) :: counttemp, fulldataShape
		integer(4), dimension(:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_INTEGER
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=3
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_integer_3D

	subroutine write_integer_4D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		integer(4), dimension(:,:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(4) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(4) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(4) :: counttemp, fulldataShape
		integer(4), dimension(:,:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_INTEGER
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=4
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_integer_4D

	subroutine write_integer_5D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		integer(4), dimension(:,:,:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(5) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(5) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(5) :: counttemp, fulldataShape
		integer(4), dimension(:,:,:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_INTEGER
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=5
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_integer_5D

	subroutine write_integer_6D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		integer(4), dimension(:,:,:,:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(6) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(6) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(6) :: counttemp, fulldataShape
		integer(4), dimension(:,:,:,:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_INTEGER
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=6
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_integer_6D

	subroutine write_integer_7D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		integer(4), dimension(:,:,:,:,:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(7) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(7) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(7) :: counttemp, fulldataShape
		integer(4), dimension(:,:,:,:,:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_INTEGER
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=7
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_integer_7D

	subroutine write_real_1D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(4), dimension(:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(1) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(1) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(1) :: counttemp, fulldataShape
		real(4), dimension(:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_REAL
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=1
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_real_1D

	subroutine write_real_2D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(4), dimension(:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(2) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(2) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(2) :: counttemp, fulldataShape
		real(4), dimension(:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_REAL
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=2
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_real_2D

	subroutine write_real_3D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(4), dimension(:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(3) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(3) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(3) :: counttemp, fulldataShape
		real(4), dimension(:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_REAL
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=3
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_real_3D

	subroutine write_real_4D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(4), dimension(:,:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(4) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(4) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(4) :: counttemp, fulldataShape
		real(4), dimension(:,:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_REAL
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=4
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_real_4D

	subroutine write_real_5D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(4), dimension(:,:,:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(5) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(5) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(5) :: counttemp, fulldataShape
		real(4), dimension(:,:,:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_REAL
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=5
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_real_5D

	subroutine write_real_6D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(4), dimension(:,:,:,:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(6) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(6) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(6) :: counttemp, fulldataShape
		real(4), dimension(:,:,:,:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_REAL
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=6
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_real_6D

	subroutine write_real_7D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(4), dimension(:,:,:,:,:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(7) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(7) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(7) :: counttemp, fulldataShape
		real(4), dimension(:,:,:,:,:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_REAL
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=7
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_real_7D

	subroutine write_double_1D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(8), dimension(:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(1) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(1) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(1) :: counttemp, fulldataShape
		real(8), dimension(:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_DOUBLE
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=1
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_double_1D

	subroutine write_double_2D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(8), dimension(:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(2) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(2) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(2) :: counttemp, fulldataShape
		real(8), dimension(:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_DOUBLE
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=2
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_double_2D

	subroutine write_double_3D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(8), dimension(:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(3) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(3) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(3) :: counttemp, fulldataShape
		real(8), dimension(:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_DOUBLE
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=3
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_double_3D

	subroutine write_double_4D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(8), dimension(:,:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(4) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(4) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(4) :: counttemp, fulldataShape
		real(8), dimension(:,:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_DOUBLE
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=4
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_double_4D

	subroutine write_double_5D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(8), dimension(:,:,:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(5) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(5) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(5) :: counttemp, fulldataShape
		real(8), dimension(:,:,:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_DOUBLE
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=5
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_double_5D

	subroutine write_double_6D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(8), dimension(:,:,:,:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(6) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(6) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(6) :: counttemp, fulldataShape
		real(8), dimension(:,:,:,:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_DOUBLE
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=6
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_double_6D

	subroutine write_double_7D(thisHDF5file, dsetname, data, chunkdim)
		implicit none
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: dsetname
		real(8), dimension(:,:,:,:,:,:,:), intent(in) :: data
		integer, intent(in), optional  :: chunkdim(:)
		INTEGER(HSIZE_T), DIMENSION(7) :: dimsf, count
		INTEGER(HSSIZE_T), DIMENSION(7) :: offset
		INTEGER(HID_T) :: dset_id, filespace, memspace, plist_id, dtype
		integer, dimension(7) :: counttemp, fulldataShape
		real(8), dimension(:,:,:,:,:,:,:),allocatable ::datatemp
		integer :: i
		dtype = H5T_NATIVE_DOUBLE
        if(thisHDF5file%serialmode == .False.) then
			if (.NOT. present(chunkdim)) then
				call print_chunkdim_present_error_info()
				stop
			end if
			call thisHDF5file%setchunk(chunkdim)
			allocate (datatemp, mold=data)
			datatemp = data
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%thisdomain%count = shape(datatemp)
			call thisHDF5file%thisdomain%calc_off(thisHDF5file%chunkrank, dimsf)
			count = thisHDF5file%thisdomain%count
			offset = thisHDF5file%thisdomain%offset
		else
			call thisHDF5file%mkdir(dsetname)
			thisHDF5file%chunkrank=7
			dimsf = shape(data)
			count = shape(data)
			offset = 0
		end if
		call prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
		if(thisHDF5file%serialmode == .False.) then
			CALL h5dwrite_f(dset_id, dtype, datatemp, dimsf, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
		else
			CALL h5dwrite_f(dset_id, dtype, data, dimsf, error, file_space_id=filespace, xfer_prp=plist_id)
		end if
		call postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
		if (allocated(datatemp)) deallocate (datatemp)

		return
	end subroutine write_double_7D

	subroutine read_attr_character(thisHDF5File, objname, attr_name, attr_data)
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: objname,attr_name
		character(len=*),intent(INOUT) :: attr_data
		INTEGER(HID_T) :: dset_id, obj_id, aspace_id, attr_id, dtype
		INTEGER(HSIZE_T), DIMENSION(1) :: adims=[1]
		logical :: dset_exists=.False.,attr_exists=.False.
		integer :: slash=0,cursor=1,arank=1
		INTEGER(SIZE_T):: max_size
		integer :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        do
			slash=index(objname(cursor:),'/')
			cursor=cursor+slash
			if (slash==0) then
				call h5lexists_f(thisHDF5file%file_id, objname, dset_exists, error)
				if (.not.dset_exists) then
					if(0 == rank)write(*,*) "object '", objname, "' don't exist!"
					stop
				end if
				exit
			end if
			call h5lexists_f(thisHDF5file%file_id, objname(:cursor-1), dset_exists, error)
			if (.not.dset_exists) then
				if(0 == rank)write(*,*) "object '", objname(:cursor-1), "' don't exist!"
				stop
			end if
		end do
		slash=0;cursor=1
		CALL h5Oopen_f(thisHDF5file%file_id, objname, obj_id, error)
		CALL h5aexists_f(obj_id,attr_name,attr_exists,error)
		if(.NOT. attr_exists) then
			write(*,"(a,i0,a,a,a)") "=== Error: Image <",rank,"> report that attribute '",attr_name,"' does not exist ==="
			stop
		end if
		CALL h5aopen_f(obj_id,attr_name,attr_id,error)
		CALL h5aget_storage_size_f(attr_id,max_size,error)
		CALL h5aget_type_f(attr_id,dtype,error)
		CALL h5aread_f(attr_id,dtype,attr_data,adims,error)
		CALL h5aclose_f(attr_id, error)
		CALL h5Oclose_f(obj_id,error)
	end subroutine read_attr_character

	subroutine read_attr_integer(thisHDF5File, objname, attr_name, attr_data)
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: objname,attr_name
		integer(4), intent(inout) :: attr_data
		INTEGER(HID_T) :: dset_id, obj_id, aspace_id, attr_id, dtype
		INTEGER(HSIZE_T), DIMENSION(1) :: adims=[1]
		logical :: dset_exists=.False.,attr_exists=.False.
		integer :: slash=0,cursor=1,arank=1
		integer :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        do
			slash=index(objname(cursor:),'/')
			cursor=cursor+slash
			if (slash==0) then
				call h5lexists_f(thisHDF5file%file_id, objname, dset_exists, error)
				if (.not.dset_exists) then
					if(0 == rank)write(*,*) "object '", objname, "' don't exist!"
					stop
				end if
				exit
			end if
			call h5lexists_f(thisHDF5file%file_id, objname(:cursor-1), dset_exists, error)
			if (.not.dset_exists) then
				if(0 == rank)write(*,*) "object '", objname(:cursor-1), "' don't exist!"
				stop
			end if
		end do
		slash=0;cursor=1
		dtype = H5T_NATIVE_INTEGER
		CALL h5Oopen_f(thisHDF5file%file_id, objname, obj_id, error)
		CALL h5aexists_f(obj_id,attr_name,attr_exists,error)
		if(.NOT. attr_exists) then
			write(*,"(a,i0,a,a,a)") "=== Error: Image <",rank,"> report that attribute '",attr_name,"' does not exist ==="
			stop
		end if
		CALL h5aopen_f(obj_id,attr_name,attr_id,error)
		CALL h5aread_f(attr_id,dtype,attr_data,adims,error)
		CALL h5aclose_f(attr_id, error)
		CALL h5Oclose_f(obj_id,error)
	end subroutine read_attr_integer

	subroutine read_attr_real(thisHDF5File, objname, attr_name, attr_data)
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: objname,attr_name
		real(4), intent(inout) :: attr_data
		INTEGER(HID_T) :: dset_id, obj_id, aspace_id, attr_id, dtype
		INTEGER(HSIZE_T), DIMENSION(1) :: adims=[1]
		logical :: dset_exists=.False.,attr_exists=.False.
		integer :: slash=0,cursor=1,arank=1
		integer :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        do
			slash=index(objname(cursor:),'/')
			cursor=cursor+slash
			if (slash==0) then
				call h5lexists_f(thisHDF5file%file_id, objname, dset_exists, error)
				if (.not.dset_exists) then
					if(0 == rank)write(*,*) "object '", objname, "' don't exist!"
					stop
				end if
				exit
			end if
			call h5lexists_f(thisHDF5file%file_id, objname(:cursor-1), dset_exists, error)
			if (.not.dset_exists) then
				if(0 == rank)write(*,*) "object '", objname(:cursor-1), "' don't exist!"
				stop
			end if
		end do
		slash=0;cursor=1
		dtype = H5T_NATIVE_REAL
		CALL h5Oopen_f(thisHDF5file%file_id, objname, obj_id, error)
		CALL h5aexists_f(obj_id,attr_name,attr_exists,error)
		if(.NOT. attr_exists) then
			write(*,"(a,i0,a,a,a)") "=== Error: Image <",rank,"> report that attribute '",attr_name,"' does not exist ==="
			stop
		end if
		CALL h5aopen_f(obj_id,attr_name,attr_id,error)
		CALL h5aread_f(attr_id,dtype,attr_data,adims,error)
		CALL h5aclose_f(attr_id, error)
		CALL h5Oclose_f(obj_id,error)
	end subroutine read_attr_real

	subroutine read_attr_double(thisHDF5File, objname, attr_name, attr_data)
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: objname,attr_name
		real(8), intent(inout) :: attr_data
		INTEGER(HID_T) :: dset_id, obj_id, aspace_id, attr_id, dtype
		INTEGER(HSIZE_T), DIMENSION(1) :: adims=[1]
		logical :: dset_exists=.False.,attr_exists=.False.
		integer :: slash=0,cursor=1,arank=1
		integer :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        do
			slash=index(objname(cursor:),'/')
			cursor=cursor+slash
			if (slash==0) then
				call h5lexists_f(thisHDF5file%file_id, objname, dset_exists, error)
				if (.not.dset_exists) then
					if(0 == rank)write(*,*) "object '", objname, "' don't exist!"
					stop
				end if
				exit
			end if
			call h5lexists_f(thisHDF5file%file_id, objname(:cursor-1), dset_exists, error)
			if (.not.dset_exists) then
				if(0 == rank)write(*,*) "object '", objname(:cursor-1), "' don't exist!"
				stop
			end if
		end do
		slash=0;cursor=1
		dtype = H5T_NATIVE_DOUBLE
		CALL h5Oopen_f(thisHDF5file%file_id, objname, obj_id, error)
		CALL h5aexists_f(obj_id,attr_name,attr_exists,error)
		if(.NOT. attr_exists) then
			write(*,"(a,i0,a,a,a)") "=== Error: Image <",rank,"> report that attribute '",attr_name,"' does not exist ==="
			stop
		end if
		CALL h5aopen_f(obj_id,attr_name,attr_id,error)
		CALL h5aread_f(attr_id,dtype,attr_data,adims,error)
		CALL h5aclose_f(attr_id, error)
		CALL h5Oclose_f(obj_id,error)
	end subroutine read_attr_double

	subroutine read_attr_logical(thisHDF5File, objname, attr_name, attr_data)
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: objname,attr_name
		logical,intent(INOUT) :: attr_data
		character :: logicread
		INTEGER(HID_T) :: dset_id, obj_id, aspace_id, attr_id, dtype
		INTEGER(HSIZE_T), DIMENSION(1) :: adims=[1]
		logical :: dset_exists=.False.,attr_exists=.False.
		integer :: slash=0,cursor=1,arank=1
		INTEGER(SIZE_T):: max_size
		integer :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        do
			slash=index(objname(cursor:),'/')
			cursor=cursor+slash
			if (slash==0) then
				call h5lexists_f(thisHDF5file%file_id, objname, dset_exists, error)
				if (.not.dset_exists) then
					if(0 == rank)write(*,*) "object '", objname, "' don't exist!"
					stop
				end if
				exit
			end if
			call h5lexists_f(thisHDF5file%file_id, objname(:cursor-1), dset_exists, error)
			if (.not.dset_exists) then
				if(0 == rank)write(*,*) "object '", objname(:cursor-1), "' don't exist!"
				stop
			end if
		end do
		slash=0;cursor=1
		CALL h5Oopen_f(thisHDF5file%file_id, objname, obj_id, error)
		CALL h5aexists_f(obj_id,attr_name,attr_exists,error)
		if(.NOT. attr_exists) then
			write(*,"(a,i0,a,a,a)") "=== Error: Image <",rank,"> report that attribute '",attr_name,"' does not exist ==="
			stop
		end if
		CALL h5aopen_f(obj_id,attr_name,attr_id,error)
		CALL h5aget_storage_size_f(attr_id,max_size,error)
		CALL h5aread_f(attr_id,H5T_NATIVE_CHARACTER,logicread,adims,error)
		CALL h5aclose_f(attr_id, error)
		CALL h5Oclose_f(obj_id,error)
		if(logicread=='T') then
			attr_data=.True.
		elseif(logicread=='F') then
			attr_data=.False.
		else
			write(*,"(a,i0,a,a,a)") "=== Error: Image <",rank,"> report that '",logicread,"' is not a logical value"
			stop
		end if
	end subroutine read_attr_logical

    subroutine write_attr_character(thisHDF5File, objname, attr_name, attr_data)
        class(HDF5_PDump), intent(inout) :: thisHDF5file
        character(len=*), intent(in) :: objname,attr_name
        character(len=*) :: attr_data
        INTEGER(HID_T) :: dset_id, obj_id, aspace_id, attr_id, dtype
        INTEGER(HSIZE_T), DIMENSION(1) :: adims=[1]
        logical :: dset_exists=.False.
        integer :: slash=0,cursor=1,arank=1
        integer(size_t) :: attrlen
		integer :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        do
            slash=index(objname(cursor:),'/')
            cursor=cursor+slash
            if (slash==0) then
                call h5lexists_f(thisHDF5file%file_id, objname, dset_exists, error)
                if (.not.dset_exists) then
                    if(0 == rank)write(*,*) "object '", objname, "' don't exist!"
                    stop
                end if
                exit
            end if
            call h5lexists_f(thisHDF5file%file_id, objname(:cursor-1), dset_exists, error)
            if (.not.dset_exists) then
                if(0 == rank)write(*,*) "object '", objname(:cursor-1), "' don't exist!"
                stop
            end if
        end do
        slash=0;cursor=1

        ! dtype = H5T_NATIVE_INTEGER
        ! attr_data=trim(attr_data)
        ! write(*,*) attr_data
        attrlen=len(attr_data)

        CALL h5Oopen_f(thisHDF5file%file_id, objname, obj_id, error)
        CALL h5screate_simple_f(arank, adims, aspace_id, error)
        CALL h5tcopy_f(H5T_NATIVE_CHARACTER, dtype, error)
        CALL h5tset_size_f(dtype, attrlen, error)
        CALL h5acreate_f(obj_id, attr_name, dtype, aspace_id,attr_id,error)
        CALL h5awrite_f(attr_id, dtype, attr_data, adims, error)
        CALL h5aclose_f(attr_id, error)
        CALL h5tclose_f(dtype, error)
        CALL h5sclose_f(aspace_id, error)
        CALL h5Oclose_f(obj_id,error)
    end subroutine write_attr_character

    subroutine write_attr_integer(thisHDF5File, objname, attr_name, attr_data)
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: objname,attr_name
        integer(4), intent(in) :: attr_data
		INTEGER(HID_T) :: dset_id, obj_id, aspace_id, attr_id, dtype
		INTEGER(HSIZE_T), DIMENSION(1) :: adims=[1]
		logical :: dset_exists=.False.
		integer :: slash=0,cursor=1,arank=1
		integer :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        do
			slash=index(objname(cursor:),'/')
			cursor=cursor+slash
			if (slash==0) then
				call h5lexists_f(thisHDF5file%file_id, objname, dset_exists, error)
				if (.not.dset_exists) then
					if(0 == rank)write(*,*) "object '", objname, "' don't exist!"
					stop
				end if
				exit
			end if
			call h5lexists_f(thisHDF5file%file_id, objname(:cursor-1), dset_exists, error)
			if (.not.dset_exists) then
				if(0 == rank)write(*,*) "object '", objname(:cursor-1), "' don't exist!"
				stop
			end if
		end do
		slash=0;cursor=1
        dtype = H5T_NATIVE_INTEGER
		CALL h5Oopen_f(thisHDF5file%file_id, objname, obj_id, error)
		CALL h5screate_simple_f(arank, adims, aspace_id, error)
		CALL h5acreate_f(obj_id, attr_name, dtype, aspace_id,attr_id,error)
		CALL h5awrite_f(attr_id, dtype, attr_data, adims, error)
		CALL h5aclose_f(attr_id, error)
		CALL h5sclose_f(aspace_id, error)
		CALL h5Oclose_f(obj_id,error)
    end subroutine write_attr_integer

    subroutine write_attr_real(thisHDF5File, objname, attr_name, attr_data)
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: objname,attr_name
        real(4), intent(in) :: attr_data
		INTEGER(HID_T) :: dset_id, obj_id, aspace_id, attr_id, dtype
		INTEGER(HSIZE_T), DIMENSION(1) :: adims=[1]
		logical :: dset_exists=.False.
		integer :: slash=0,cursor=1,arank=1
		integer :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        do
			slash=index(objname(cursor:),'/')
			cursor=cursor+slash
			if (slash==0) then
				call h5lexists_f(thisHDF5file%file_id, objname, dset_exists, error)
				if (.not.dset_exists) then
					if(0 == rank)write(*,*) "object '", objname, "' don't exist!"
					stop
				end if
				exit
			end if
			call h5lexists_f(thisHDF5file%file_id, objname(:cursor-1), dset_exists, error)
			if (.not.dset_exists) then
				if(0 == rank)write(*,*) "object '", objname(:cursor-1), "' don't exist!"
				stop
			end if
		end do
		slash=0;cursor=1
        dtype = H5T_NATIVE_REAL
		CALL h5Oopen_f(thisHDF5file%file_id, objname, obj_id, error)
		CALL h5screate_simple_f(arank, adims, aspace_id, error)
		CALL h5acreate_f(obj_id, attr_name, dtype, aspace_id,attr_id,error)
		CALL h5awrite_f(attr_id, dtype, attr_data, adims, error)
		CALL h5aclose_f(attr_id, error)
		CALL h5sclose_f(aspace_id, error)
		CALL h5Oclose_f(obj_id,error)
    end subroutine write_attr_real

    subroutine write_attr_double(thisHDF5File, objname, attr_name, attr_data)
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: objname,attr_name
        real(8), intent(in) :: attr_data
		INTEGER(HID_T) :: dset_id, obj_id, aspace_id, attr_id, dtype
		INTEGER(HSIZE_T), DIMENSION(1) :: adims=[1]
		logical :: dset_exists=.False.
		integer :: slash=0,cursor=1,arank=1
		integer :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        do
			slash=index(objname(cursor:),'/')
			cursor=cursor+slash
			if (slash==0) then
				call h5lexists_f(thisHDF5file%file_id, objname, dset_exists, error)
				if (.not.dset_exists) then
					if(0 == rank)write(*,*) "object '", objname, "' don't exist!"
					stop
				end if
				exit
			end if
			call h5lexists_f(thisHDF5file%file_id, objname(:cursor-1), dset_exists, error)
			if (.not.dset_exists) then
				if(0 == rank)write(*,*) "object '", objname(:cursor-1), "' don't exist!"
				stop
			end if
		end do
		slash=0;cursor=1
        dtype = H5T_NATIVE_DOUBLE
		CALL h5Oopen_f(thisHDF5file%file_id, objname, obj_id, error)
		CALL h5screate_simple_f(arank, adims, aspace_id, error)
		CALL h5acreate_f(obj_id, attr_name, dtype, aspace_id,attr_id,error)
		CALL h5awrite_f(attr_id, dtype, attr_data, adims, error)
		CALL h5aclose_f(attr_id, error)
		CALL h5sclose_f(aspace_id, error)
		CALL h5Oclose_f(obj_id,error)
    end subroutine write_attr_double

	subroutine write_attr_logical(thisHDF5File, objname, attr_name, attr_data)
		class(HDF5_PDump), intent(inout) :: thisHDF5file
		character(len=*), intent(in) :: objname,attr_name
        logical, intent(in) :: attr_data
		character(len=1) :: attr_save
		INTEGER(HID_T) :: dset_id, obj_id, aspace_id, attr_id, dtype
		INTEGER(HSIZE_T), DIMENSION(1) :: adims=[1]
		logical :: dset_exists=.False.
		integer :: slash=0,cursor=1,arank=1
		integer :: size, rank, ierr

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        do
			slash=index(objname(cursor:),'/')
			cursor=cursor+slash
			if (slash==0) then
				call h5lexists_f(thisHDF5file%file_id, objname, dset_exists, error)
				if (.not.dset_exists) then
					if(0 == rank)write(*,*) "object '", objname, "' don't exist!"
					stop
				end if
				exit
			end if
			call h5lexists_f(thisHDF5file%file_id, objname(:cursor-1), dset_exists, error)
			if (.not.dset_exists) then
				if(0 == rank)write(*,*) "object '", objname(:cursor-1), "' don't exist!"
				stop
			end if
		end do
		slash=0;cursor=1
		attr_save=merge('T','F',attr_data)
        dtype = H5T_NATIVE_CHARACTER
		CALL h5Oopen_f(thisHDF5file%file_id, objname, obj_id, error)
		CALL h5screate_simple_f(arank, adims, aspace_id, error)
		CALL h5acreate_f(obj_id, attr_name, dtype, aspace_id,attr_id,error)
		CALL h5awrite_f(attr_id, dtype, attr_save, adims, error)
		CALL h5aclose_f(attr_id, error)
		CALL h5sclose_f(aspace_id, error)
		CALL h5Oclose_f(obj_id,error)
    end subroutine write_attr_logical

    subroutine prewrite(thisHDF5file, dsetname, dimsf, count, offset, dset_id, filespace, memspace, plist_id, dtype)
        implicit none
        class(HDF5_PDump), intent(inout) :: thisHDF5file
        CHARACTER(LEN=*), intent(in) :: dsetname
        INTEGER(HSIZE_T), intent(inout) :: dimsf(thisHDF5file%chunkrank)
        INTEGER(HSIZE_T), DIMENSION(thisHDF5file%chunkrank), intent(inout) :: count
        INTEGER(HSSIZE_T), DIMENSION(thisHDF5file%chunkrank), intent(inout) :: offset
        INTEGER(HID_T), intent(inout) :: dset_id, filespace, memspace, plist_id, dtype

        !
        ! Create the data space for the  dataset.
        !
        CALL h5screate_simple_f(size(dimsf), dimsf, filespace, error)

        !
        ! Create the dataset with default properties.
        !
        CALL h5dcreate_f(thisHDF5file%file_id, dsetname, dtype, filespace, dset_id, error)
		if(thisHDF5file%serialmode ==.FALSE.) then
			CALL h5sclose_f(filespace, error)

			!
			! Each process defines dataset in memory and writes it to the hyperslab
			! in the file.
			!
			CALL h5screate_simple_f(size(dimsf), count, memspace, error)

			!
			! Select hyperslab in the file.
			!
			CALL h5dget_space_f(dset_id, filespace, error)
			CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)
		end if

        !
        ! Create property list for collective dataset write
        !
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
		if(thisHDF5file%serialmode ==.FALSE.) then
        	CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
		end if
        return
    end subroutine prewrite

    subroutine postwrite(thisHDF5file, dset_id, filespace, memspace, plist_id)
        implicit none
        class(HDF5_PDump), intent(inout) :: thisHDF5file
        INTEGER(HID_T), intent(inout) :: dset_id, filespace, memspace, plist_id
        !
        ! Close dataspaces.
        !
        CALL h5sclose_f(filespace, error)
		if(thisHDF5file%serialmode ==.FALSE.) then
        	CALL h5sclose_f(memspace, error)
		end if
        !
        ! Close the dataset and property list.
        !
        CALL h5dclose_f(dset_id, error)
        CALL h5pclose_f(plist_id, error)

        return
    end subroutine postwrite

	subroutine preread(thisHDF5file,dsetname,dset_id,dspace_id,dims,rank)
		implicit none
        class(HDF5_PDump), intent(inout) :: thisHDF5file
		INTEGER(HID_T), intent(inout) :: dset_id, dspace_id
		character(len=*),intent(in) :: dsetname
		integer(HSIZE_T),dimension(:),intent(inout) :: dims
		integer,intent(in) :: rank
		integer(HSIZE_T),dimension(rank) :: maxdims
		INTEGER :: ndims
		CALL h5dopen_f(thisHDF5File%file_id, dsetname, dset_id, error)
		CALL h5dget_space_f(dset_id, dspace_id, error)
		CALL h5sget_simple_extent_ndims_f(dspace_id,ndims,error)
		if(ndims /= rank) then
			write(*,"(a,i0,a)") "=== Error: report that the container array does not match data shape ==="
			stop
		end if
		CALL h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,error)
	end subroutine preread

	subroutine postread(thisHDF5File, dspace_id, dset_id)
		implicit none
        class(HDF5_PDump), intent(inout) :: thisHDF5file
		INTEGER(HID_T), intent(inout) :: dset_id, dspace_id
		CALL h5sclose_f(dspace_id,error)
		CALL h5dclose_f(dset_id,error)
	end subroutine postread


    subroutine close_HDF5_PDump(thisHDF5file)
        implicit none
        class(HDF5_PDump), intent(inout) :: thisHDF5file

        !
        ! Close the file.
        !
        CALL h5fclose_f(thisHDF5file%file_id, error)
        !
        ! Close FORTRAN predefined datatypes.
        !
        CALL h5close_f(error)

        return
    end subroutine close_HDF5_PDump

    subroutine init_domain(selfdomain, domainrank)
        implicit none
        class(DomainDump), intent(inout) :: selfdomain
        integer, intent(in) :: domainrank

        if (allocated(selfdomain%count)) deallocate (selfdomain%count)
        if (allocated(selfdomain%offset)) deallocate (selfdomain%offset)
        if (allocated(selfdomain%location)) deallocate (selfdomain%location)

        allocate (selfdomain%count(domainrank))
        allocate (selfdomain%offset(domainrank))
        allocate (selfdomain%location(domainrank))

        return
    end subroutine init_domain

    subroutine destroy(selfdomain)
        implicit none
        type(DomainDump), intent(inout) :: selfdomain

        if (allocated(selfdomain%count)) deallocate (selfdomain%count)
        if (allocated(selfdomain%offset)) deallocate (selfdomain%offset)
        if (allocated(selfdomain%location)) deallocate (selfdomain%location)

        return
    end subroutine destroy

    subroutine calculate_domain_location(selfdomain, chunkdim)
        class(DomainDump), intent(inout) :: selfdomain
        integer, intent(in)  :: chunkdim(:)
        integer :: i, rank, size_image, ierr
        integer, allocatable :: images(:)

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size_image, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        allocate(images(size_image))
        images = [(i, i=0, size_image-1)]

        select case (size(chunkdim))
        case (1)
            selfdomain%location = findloc(reshape(images, chunkdim(1:1)), rank)
        case (2)
            selfdomain%location = findloc(reshape(images, chunkdim(1:2)), rank)
        case (3)
            selfdomain%location = findloc(reshape(images, chunkdim(1:3)), rank)
        case (4)
            selfdomain%location = findloc(reshape(images, chunkdim(1:4)), rank)
        case (5)
            selfdomain%location = findloc(reshape(images, chunkdim(1:5)), rank)
        case (6)
            selfdomain%location = findloc(reshape(images, chunkdim(1:6)), rank)
        case (7)
            selfdomain%location = findloc(reshape(images, chunkdim(1:7)), rank)
        case default
            write(*,"(a,i0,a)") "=== Error: Image <",rank,"> report that Rank of Images must be 1~7 ==="
            stop
        end select

        deallocate(images)

    end subroutine calculate_domain_location

    subroutine calculate_domain_offset(domain_calc, domainrank, dimsfout)
        integer, intent(in)  :: domainrank
        class(DomainDump), intent(inout) :: domain_calc
        class(DomainDump), allocatable :: domaintemp
        integer(HSIZE_T), dimension(domainrank), intent(out) :: dimsfout   !总尺寸
        integer, dimension(domainrank) :: dimsftemp
        integer, allocatable :: SumTemp(:, :)
        integer :: i, rank, size_image, ierr
        integer :: win, win_count
        integer(kind=MPI_ADDRESS_KIND) :: win_size, disp_aint
        integer, dimension(domainrank) :: location_recv, location_send, count_recv, count_send
        integer, dimension(domainrank) :: dimsfout_bcast

        call MPI_COMM_SIZE(MPI_COMM_WORLD, size_image, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        allocate (SumTemp(domainrank, size_image))
        allocate (domaintemp)
        call domaintemp%init(domainrank)
        dimsfout = 0
        domaintemp = domain_calc
        location_send = domaintemp%location
        count_send = domaintemp%count

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        win_size = 4 * domainrank
        call MPI_WIN_CREATE(location_send, win_size, 4, MPI_INFO_NULL, MPI_COMM_WORLD, win, ierr)
        call MPI_WIN_CREATE(count_send, win_size, 4, MPI_INFO_NULL, MPI_COMM_WORLD, win_count, ierr)

        do i = 1, size_image
            disp_aint = 0
            call MPI_WIN_LOCK(MPI_LOCK_SHARED, i-1, 0, win, ierr)
            call MPI_GET(location_recv, domainrank, MPI_INTEGER, i-1, disp_aint, domainrank, MPI_INTEGER, win, ierr)
            call MPI_WIN_UNLOCK(i-1, win, ierr)

            if (all(location_recv <= domaintemp%location)) then
                disp_aint = 0
                call MPI_WIN_LOCK(MPI_LOCK_SHARED, i-1, 0, win_count, ierr)
                call MPI_GET(count_recv, domainrank, MPI_INTEGER, i-1, disp_aint, domainrank, MPI_INTEGER, win_count, ierr)
                call MPI_WIN_UNLOCK(i-1, win_count, ierr)

                SumTemp(:, i) = count_recv
            else
                SumTemp(:, i) = 0
            end if
        end do

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        dimsftemp = sum(SumTemp, dim=2)
        do i = 1, domainrank
            dimsftemp(i) = dimsftemp(i)/(PRODUCT(domaintemp%location)/domaintemp%location(i))
        end do
        domain_calc%offset = dimsftemp - domain_calc%count

        if (rank+1 == size_image) dimsfout = dimsftemp

        dimsfout_bcast = dimsfout
        call MPI_Bcast(dimsfout_bcast, domainrank, MPI_INTEGER, size_image-1, MPI_COMM_WORLD, ierr)
        dimsfout = dimsfout_bcast

        call MPI_WIN_FREE(win, ierr)
        call MPI_WIN_FREE(win_count, ierr)
        if (allocated(SumTemp)) deallocate (SumTemp)
        if (allocated(domaintemp)) deallocate (domaintemp)

        return
    end subroutine calculate_domain_offset

end module ModuleParallelDump
