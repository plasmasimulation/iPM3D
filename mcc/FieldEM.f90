Module ModuleFieldEM
    use mpi
	! use ModuleControlFlow
	! use ModuleFileName
    ! use ModuleParallelDump
    ! use ModuleDomain

	implicit none

	Type FieldEM			! Definition of electromagnetic field
		! Type(FileName) 	:: 	IOName
        ! Type(Domain), Pointer :: DM => Null()

		Real(8), Allocatable :: Ex(:,:,:),Ey(:,:,:),Ez(:,:,:)	! [Lz, Lr], Field in z, r, theta direction
		Real(8), Allocatable :: Bx(:,:,:),By(:,:,:),Bz(:,:,:)

	contains
	
	! 	procedure :: Init => InitFieldEM
	! 	! procedure :: Dump => DumpFieldEM
	! 	! procedure :: Load => LoadFieldEM
          procedure :: Destroy => DestroyFieldEM
    !     ! procedure :: Reset => ResetFieldEM
    !     procedure :: Zero => ZeroFieldEM
    !     ! procedure :: Diag => DiagFieldEM

    !     ! procedure, private :: LoadFieldEMHDF5
    !     ! procedure, private :: LoadFieldEMDAT

    !     ! procedure, private :: DumpFieldEMHDF5
         procedure ::Dump =>DumpFieldEMDAT

	 End Type FieldEM

    ! ! Type(HDF5_PDump),private :: hdf5FiledEMDump

	contains
    subroutine DestroyFieldEM(FG)
			Class(FieldEM), intent(inout) :: FG

            if (Allocated(FG%Ex)) Deallocate(FG%Ex)
            if (Allocated(FG%Ey)) Deallocate(FG%Ey)
            if (Allocated(FG%Ez)) Deallocate(FG%Ez)
            if (Allocated(FG%Bx)) Deallocate(FG%Bx)
            if (Allocated(FG%By)) Deallocate(FG%By)
            if (Allocated(FG%Bz)) Deallocate(FG%Bz)
           

        end subroutine DestroyFieldEM

        subroutine DumpFieldEMDAT(FG,xstart, xend, ystart, yend,zstart,zend)
            integer(4) :: xstart, xend, ystart, yend,zstart,zend
            		Class(FieldEM),intent(inout) :: FG
                    integer(4) :: size, rank, ierr,i
                    character(len=99) :: file_name
        
                    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
                    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        
                    if (0 == rank) write(*, '(a)') "The fieldEM is save as dat file."
                 write(file_name, '(i1)') rank
                    open(10, file="./field/field"//trim(file_name)//".txt")
                    ! open(10, file=FG%IOName%FullName%str)
                        write(10,*)"x方向电场"
                        do i=zstart,zend+1
                        write(10,*)""
                        write(10, "(3f9.3,/)") FG%Ex(xstart:xend+1,ystart:yend+1,i)
                        end do 
                        write(10,*)"y方向电场"
                        do i=zstart,zend+1
                            write(10,*)""
                            write(10, "(3f9.3,/)") FG%Ey(xstart:xend+1,ystart:yend+1,i)
                            end do 
                        write(10,*)"z方向电场"
                        do i=zstart,zend+1
                            write(10,*)""
                            write(10, "(3f9.3,/)") FG%Ez(xstart:xend+1,ystart:yend+1,i)
                            end do 
                        ! write(10, *) FG%Er
                        ! write(10, *) FG%Et
                        ! write(10, *) FG%Bz
                        ! write(10, *) FG%Br
                        ! write(10, *) FG%Bt
                    close(10)
        
                end subroutine DumpFieldEMDAT
	! 	subroutine InitFieldEM(FG, CF, fgname)
	! 		Class(FieldEM), intent(inout) :: FG
    !         Type(ControlFlow), intent(in) :: CF
    !         Character(*), optional, intent(in) :: fgname

    !         if (present(fgname)) then
    !             call FG%IOName%Init(fgname, RESTART_FILE_NAME)
    !         else
    !             call FG%IOName%Init("FieldEM", RESTART_FILE_NAME)
    !         end if

    !         call FG%Reset(CF%DM)
    !         call FG%Zero()

	! 	end subroutine InitFieldEM


    !     subroutine ResetFieldEM(FG, DM)
	! 		Class(FieldEM), intent(inout) :: FG
    !         Type(Domain), intent(in), target :: DM

    !         ! Setting
    !         FG%DM => DM

    !         ! load or init
    !         call FG%Destroy()

    !         Associate(zstart => FG%DM%CornerIndex(BOUNDARY_LEFT, 1), &
    !                   zend => FG%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
    !                   rstart => FG%DM%CornerIndex(BOUNDARY_LEFT, 2), &
    !                   rend => FG%DM%CornerIndex(BOUNDARY_RIGHT, 2))

    !             Allocate(FG%Ez(zstart-1 : zend, rstart   : rend))
    !             Allocate(FG%Er(zstart   : zend, rstart-1 : rend))
    !             Allocate(FG%Et(zstart   : zend, rstart   : rend))
    !             Allocate(FG%Bz(zstart   : zend, rstart   : rend))
    !             Allocate(FG%Br(zstart   : zend, rstart   : rend))
    !             Allocate(FG%Bt(zstart   : zend, rstart   : rend))

    !         end Associate

    !     end subroutine ResetFieldEM


	! 	subroutine DumpFieldEM(FG)
	! 		Class(FieldEM),intent(inout) :: FG

    !         if (FILE_EXTENSION_MODE_H5 == FG%IOName%ExtensionMode) then
    !             call FG%DumpFieldEMHDF5()

    !         else if (FILE_EXTENSION_MODE_DAT == FG%IOName%ExtensionMode) then
    !             call FG%DumpFieldEMDAT()
            
    !         else
    !             call FG%DumpFieldEMHDF5()

    !         end if

	! 	endsubroutine DumpFieldEM


    !     subroutine LoadFieldEM(FG)
	! 		Class(FieldEM),intent(inout) :: FG
    !         logical :: alive

    !         Inquire(file=FG%IOName%FullName%str, exist=alive)
    !         if (alive) then
    !             if (FILE_EXTENSION_MODE_H5 == FG%IOName%ExtensionMode) then
    !                 call FG%LoadFieldEMHDF5()

    !             else if (FILE_EXTENSION_MODE_DAT == FG%IOName%ExtensionMode) then
    !                 call FG%LoadFieldEMDAT()

    !             end if

    !         end if

	! 	end subroutine LoadFieldEM


    !     subroutine LoadFieldEMHDF5(FG)
    !         Class(FieldEM),intent(inout) :: FG
    !         real(8),ALLOCATABLE :: Temp2D(:,:)
    !         integer(4) :: size, rank, ierr

    !         call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    !         call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    !         if (0 == rank) write(*, '(a)') "The fieldEM is load from hdf5 file."
            
    !         Associate(Ez => FG%Ez, Er => FG%Er, Et => FG%Et, &
    !                   Bz => FG%Bz, Br => FG%Br, Bt => FG%Bt)

    !             call hdf5FiledEMDump%init(filename=FG%IOName%FullName%str, mode='read')
    !             call hdf5FiledEMDump%open()

    !             if (Allocated(Temp2D)) Deallocate(Temp2D)
    !             Allocate(Temp2D, mold=Ez)
    !             call hdf5FiledEMDump%read('Ez',Temp2D)
    !             Ez=Temp2D

    !             if (Allocated(Temp2D)) Deallocate(Temp2D)
    !             Allocate(Temp2D, mold=Er)
    !             call hdf5FiledEMDump%read('Er',Temp2D)
    !             Er=Temp2D

    !             if (Allocated(Temp2D)) Deallocate(Temp2D)
    !             Allocate(Temp2D, mold=Et)
    !             call hdf5FiledEMDump%read('Et',Temp2D)
    !             Et=Temp2D

    !             if (Allocated(Temp2D)) Deallocate(Temp2D)
    !             Allocate(Temp2D, mold=Bz)
    !             call hdf5FiledEMDump%read('Bz',Temp2D)
    !             Bz=Temp2D

    !             if (Allocated(Temp2D)) Deallocate(Temp2D)
    !             Allocate(Temp2D, mold=Br)
    !             call hdf5FiledEMDump%read('Br',Temp2D)
    !             Br=Temp2D

    !             if (Allocated(Temp2D)) Deallocate(Temp2D)
    !             Allocate(Temp2D, mold=Bt)
    !             call hdf5FiledEMDump%read('Bt',Temp2D)
    !             Bt=Temp2D
                
    !             Deallocate(Temp2D)
    !             call hdf5FiledEMDump%close()
    !         End Associate

    !     end subroutine LoadFieldEMHDF5


    !     subroutine LoadFieldEMDAT(FG)
	! 		Class(FieldEM),intent(inout) :: FG
	! 		Integer(4) :: i, j
    !         integer(4) :: size, rank, ierr

    !         call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    !         call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    !         if (0 == rank) write(*, '(a)') "The fieldEM is load from dat file."

    !         open(10, file=FG%IOName%FullName%str)
    !             read(10, *) FG%Ez
    !             read(10, *) FG%Er
    !             read(10, *) FG%Et
    !             read(10, *) FG%Bz
    !             read(10, *) FG%Br
    !             read(10, *) FG%Bt
    !         close(10)

    !     end subroutine LoadFieldEMDAT


    !     subroutine ZeroFieldEM(FG)
    !         Class(FieldEM), intent(inout) :: FG

    !         FG%Ez = 0.d0
    !         FG%Er = 0.d0
    !         FG%Et = 0.d0
    !         FG%Bz = 0.d0
    !         FG%Br = 0.d0
    !         FG%Bt = 0.d0

    !     end subroutine ZeroFieldEM


    !     subroutine DiagFieldEM(FG, name)
    !         Class(FieldEM),intent(inout) :: FG
    !         Character(*), intent(in) :: name
    !         Type(HDF5_PDump) :: hdf5Dump
    !         Type(FileName) 	:: 	IOName
    !         integer(4) :: zstart_on_bound, rstart_on_bound
    !         integer(4) :: zend_on_bound, rend_on_bound

    !         call IOName%Init(name, DIAG_FILE_NAME)

    !         call hdf5Dump%init(filename=IOName%FullName%str, mode='write')
    !         call hdf5Dump%open()

    !         Associate(zstart => FG%DM%CornerIndex(BOUNDARY_LEFT, 1), &
    !                   zend   => FG%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
    !                   rstart => FG%DM%CornerIndex(BOUNDARY_LEFT, 2), &
    !                   rend   => FG%DM%CornerIndex(BOUNDARY_RIGHT, 2), &
    !                   zRightImageType => FG%DM%NeighborType(BOUNDARY_RIGHT, 1), &
    !                   rRightImageType => FG%DM%NeighborType(BOUNDARY_RIGHT, 2), &
    !                   image_shape => FG%DM%ImageShape)

    !             call hdf5Dump%write('/Ezz', FG%Ez(zstart-1 : zend, rstart   : rend), chunkdim=image_shape)
    !             call hdf5Dump%write('/Err', FG%Er(zstart   : zend, rstart-1 : rend), chunkdim=image_shape)

    !             zstart_on_bound = zstart
    !             zend_on_bound = zend
    !             rstart_on_bound = rstart
    !             rend_on_bound = rend
    !             if (image_shape(1) > 1 .and. NEIGHBOR_TYPE_DOMAIN == zRightImageType) zend_on_bound = zend - 2
    !             if (image_shape(2) > 1 .and. NEIGHBOR_TYPE_DOMAIN == rRightImageType) rend_on_bound = rend - 1
    !             call hdf5Dump%write('/Ez', FG%Ez(zstart-1 : zend_on_bound, rstart   : rend_on_bound), chunkdim=image_shape)

    !             zstart_on_bound = zstart
    !             zend_on_bound = zend
    !             rstart_on_bound = rstart
    !             rend_on_bound = rend
    !             if (image_shape(1) > 1 .and. NEIGHBOR_TYPE_DOMAIN == zRightImageType) zend_on_bound = zend - 1
    !             if (image_shape(2) > 1 .and. NEIGHBOR_TYPE_DOMAIN == rRightImageType) rend_on_bound = rend - 2
    !             call hdf5Dump%write('/Er', FG%Er(zstart   : zend_on_bound, rstart-1 : rend_on_bound), chunkdim=image_shape)

    !         end Associate

    !         call hdf5Dump%close()

    !     end subroutine DiagFieldEM


    !     subroutine DumpFieldEMHDF5(FG)
    !         Class(FieldEM),intent(inout) :: FG
    !         integer(4) :: size, rank, ierr

    !         call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
    !         call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    !         if (0 == rank) write(*, '(a)') "The fieldEM is save as hdf5 file."

    !         Associate(dz => FG%DM%SpaceStep(1), dr => FG%DM%SpaceStep(2), &
    !                   Nz => FG%DM%GlobalShape(1), Nr => FG%DM%GlobalShape(2), &
    !                   Lz => FG%DM%LocalShape(1), Lr => FG%DM%LocalShape(2), &
    !                   Ez => FG%Ez, Er => FG%Er, Et => FG%Et, &
    !                   Bz => FG%Bz, Br => FG%Br, Bt => FG%Bt)

    !             call hdf5FiledEMDump%init(filename=FG%IOName%FullName%str, mode='write', serial=.True.)
    !             call hdf5FiledEMDump%open()

    !             call hdf5FiledEMDump%writeattr('/','Nz',Nz)
    !             call hdf5FiledEMDump%writeattr('/','Nr',Nr)
    !             call hdf5FiledEMDump%writeattr('/','Lz',Lz)
    !             call hdf5FiledEMDump%writeattr('/','Lr',Lr)
    !             call hdf5FiledEMDump%writeattr('/','dz',dz)
    !             call hdf5FiledEMDump%writeattr('/','dr',dr)
    !             call hdf5FiledEMDump%write('Ez',Ez)
    !             call hdf5FiledEMDump%write('Er',Er)
    !             call hdf5FiledEMDump%write('Et',Et)
    !             call hdf5FiledEMDump%write('Bz',Bz)
    !             call hdf5FiledEMDump%write('Br',Br)
    !             call hdf5FiledEMDump%write('Bt',Bt)

    !             call hdf5FiledEMDump%close()

    !         End Associate

    !     end subroutine DumpFieldEMHDF5


    


    

End Module ModuleFieldEM
