Module ModuleFieldSource
    use ModuleControlFlow
    use ModuleFieldOne
    use ModuleDomain

    implicit none
    
    Type FieldSource
        Type(Domain), Pointer :: DM => Null()

        real(8), Allocatable ::  Phi(:, :)
        real(8), Allocatable ::  Rho(:, :)
        real(8), Allocatable ::  Chi(:, :)
    
    contains

        procedure :: Init => InitFieldSource
        procedure :: Zero => ZeroFieldSource
        procedure :: Reset => ResetFieldSource
        procedure :: Destroy => DestroyFieldSource
        procedure :: Sum => SumFieldSource
        procedure :: SetE => SetElectrostaticFieldFromPotential
        procedure :: Diag => DiagFieldSource

    end Type FieldSource

    contains

        subroutine InitFieldSource(FT, CF)
            Class(FieldSource), intent(inout) :: FT
            Type(ControlFlow), intent(in) :: CF

            call FT%Reset(CF%DM)

        end subroutine InitFieldSource


        subroutine ResetFieldSource(FT, DM)
            Class(FieldSource), intent(inout) :: FT
            Type(Domain), intent(in), target :: DM

            FT%DM => DM
            call FT%Destroy()

            Associate(zstart => FT%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                      zend => FT%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                      rstart => FT%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                      rend => FT%DM%CornerIndex(BOUNDARY_RIGHT, 2))

                Allocate(FT%Phi(zstart-1 : zend+1, rstart-1 : rend+1))

                Allocate(FT%Rho(zstart : zend, rstart : rend))
                Allocate(FT%Chi(zstart-1 : zend+1, rstart-1 : rend+1))

            end Associate

            return
        end subroutine ResetFieldSource


        subroutine DestroyFieldSource(FT)
            Class(FieldSource), intent(inout) :: FT

            if (Allocated(FT%Rho)) Deallocate(FT%Rho)
            if (Allocated(FT%Phi)) Deallocate(FT%Phi)
            if (Allocated(FT%Chi)) Deallocate(FT%Chi)

            return
        end subroutine DestroyFieldSource


        subroutine SumFieldSource(FT, FO)
            Class(FieldSource), intent(inout) :: FT
            Type(FieldOne), intent(in) :: FO(:)
            Integer(4) :: i

            call FT%Zero()

            Associate(zstart => FT%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                      zend   => FT%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                      rstart => FT%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                      rend   => FT%DM%CornerIndex(BOUNDARY_RIGHT, 2))

                do i = 1, size(FO)
                    FT%Rho = FT%Rho + FO(i)%RhoOne
                    FT%Chi(zstart:zend, rstart:rend) = FT%Chi(zstart:zend, rstart:rend) + FO(i)%ChiOne(zstart:zend, rstart:rend)
                end do

            end Associate

            return
        end subroutine SumFieldSource


        subroutine SetElectrostaticFieldFromPotential(FT, FG)
            Class(FieldSource), intent(inout) :: FT
            Type(FieldEM), intent(inout) :: FG
            integer(4) :: i

            Associate(zstart => FT%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                      zend   => FT%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                      rstart => FT%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                      rend   => FT%DM%CornerIndex(BOUNDARY_RIGHT, 2), &
                      dz  => FG%DM%SpaceStep(1), dr => FG%DM%SpaceStep(2), &
                      Phi => FT%Phi, Ez => FG%Ez, Er => FG%Er)

                do i = zstart-1, zend
                    Ez(i, :) = Phi(i, rstart:rend) - Phi(i+1, rstart:rend)
                end do
                Ez = Ez / dz

                do i = rstart-1, rend
                    Er(:, i) = Phi(zstart:zend, i) - Phi(zstart:zend, i+1)
                end do
                Er = Er / dr

            end Associate

            return
        end subroutine SetElectrostaticFieldFromPotential


        subroutine DiagFieldSource(FT, name)
            Class(FieldSource),intent(inout) :: FT
            Character(*), intent(in) :: name
            Type(HDF5_PDump) :: hdf5Dump
            Type(FileName) 	:: 	IOName
            integer(4) :: zend_on_bound, rend_on_bound

            call IOName%Init(name, PathMode=FILE_PATH_MODE_DIAG, &
                                    ExtensionMode=FILE_EXTENSION_MODE_H5, &
                                    ParallelMode=FILE_PARALLEL_MODE_CLOSE, &
                                    DynamicIndex=FILE_DYNAMIC_MODE_CLOSE)

            call hdf5Dump%init(filename=IOName%FullName%str, mode='write')
            call hdf5Dump%open()

            Associate(zstart => FT%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                      zend   => FT%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                      rstart => FT%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                      rend   => FT%DM%CornerIndex(BOUNDARY_RIGHT, 2), &
                      zRightImageType => FT%DM%NeighborType(BOUNDARY_RIGHT, 1), &
                      rRightImageType => FT%DM%NeighborType(BOUNDARY_RIGHT, 2), &
                      image_shape => FT%DM%ImageShape)

                zend_on_bound = zend
                rend_on_bound = rend
                if (image_shape(1) > 1 .and. NEIGHBOR_TYPE_DOMAIN == zRightImageType) zend_on_bound = zend - 1
                if (image_shape(2) > 1 .and. NEIGHBOR_TYPE_DOMAIN == rRightImageType) rend_on_bound = rend - 1

                call hdf5Dump%write('/Phii', FT%Phi(zstart-1 : zend+1, rstart-1 : rend+1), chunkdim=image_shape)
                call hdf5Dump%write('/Phi', FT%Phi(zstart:zend_on_bound, rstart:rend_on_bound), chunkdim=image_shape)

                call hdf5Dump%write('/Rho', FT%Rho(zstart:zend_on_bound, rstart:rend_on_bound), chunkdim=image_shape)
                
                call hdf5Dump%write('/Chii', FT%Chi(zstart-1 : zend+1, rstart-1 : rend+1), chunkdim=image_shape)
                call hdf5Dump%write('/Chi', FT%Chi(zstart:zend_on_bound, rstart:rend_on_bound), chunkdim=image_shape)

            end Associate

            call hdf5Dump%close()

            return
        end subroutine DiagFieldSource


        subroutine ZeroFieldSource(FT)
            Class(FieldSource), intent(inout) :: FT

            FT%Rho = 0.d0
            FT%Phi = 0.d0
            FT%Chi = 0.d0
            
            return
        end subroutine ZeroFieldSource

EndModule ModuleFieldSource
