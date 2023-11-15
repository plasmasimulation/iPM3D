Module ModuleFieldOne
    use ModuleControlFlow
    use ModuleParticleBundle
    use ModuleBspline
    use ModuleDomain
    use ModuleFieldEM
    use ModuleFileName
    use ModuleGeometry

    implicit none

    Type FieldOne
        Type(Domain), Pointer :: DM => Null()
        real(8) :: RhoFactor, JonFactor, ChiFactor

        real(8), Allocatable :: RhoOne(:, :)
        real(8), Allocatable :: JonOne(:, :)
        real(8), Allocatable :: ChiOne(:, :)

    contains

        procedure :: Init => InitFieldOne
        procedure :: Zero => ZeroFieldOne
        procedure :: Rescale => RescaleFieldOne
        procedure :: Reset => ResetFieldOne
        procedure :: Destroy => DestroyFieldOne
        procedure :: Diag => DiagFieldOne

        procedure :: P2CES => WeightingFieldOneElectrostatic

    end Type FieldOne

    contains

        subroutine InitFieldOne(FO, PB, CF)
            Class(FieldOne),intent(inout) :: FO
            Type(ParticleBundle),intent(in) :: PB
            Type(ControlFlow), intent(in) :: CF

            call FO%Reset(PB, CF%DM)

            return
        end subroutine InitFieldOne


        subroutine ResetFieldOne(FO, PB, DM)
            Class(FieldOne),intent(inout) :: FO
            Type(ParticleBundle),intent(in) :: PB
            Type(Domain), intent(in), target :: DM

            FO%DM => DM
            FO%RhoFactor = PB%Charge
            FO%ChiFactor = 0.5d0 * PB%Charge / PB%Mass * PB%dt * PB%dt / Epsilon

            ! init
            call FO%Destroy()

            Associate(zstart => FO%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                      zend => FO%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                      rstart => FO%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                      rend => FO%DM%CornerIndex(BOUNDARY_RIGHT, 2))

                Allocate(FO%RhoOne(zstart : zend, rstart : rend))
                Allocate(FO%JonOne(zstart : zend, rstart : rend))
                Allocate(FO%ChiOne(zstart : zend, rstart : rend))

            end Associate

            return
        end subroutine ResetFieldOne


        subroutine ZeroFieldOne(FO)
            Class(FieldOne), intent(inout) :: FO

            FO%RhoOne=0.d0
            FO%JonOne=0.d0
            FO%ChiOne=0.d0

            return
        end subroutine ZeroFieldOne


        subroutine RescaleFieldOne(FO)
            Class(FieldOne),intent(inout) :: FO

            FO%RhoOne = FO%RhoOne * FO%RhoFactor
            FO%JonOne = FO%JonOne * FO%JonFactor
            FO%ChiOne = FO%ChiOne * FO%ChiFactor

            return
        end subroutine RescaleFieldOne


        subroutine DiagFieldOne(FO, name)
            Class(FieldOne),intent(inout) :: FO
            Character(*), intent(in) :: name
            Type(HDF5_PDump) :: hdf5Dump
            Type(FileName) 	:: 	IOName
            integer(4) :: zend_on_bound, rend_on_bound

            call IOName%Init(name, DIAG_FILE_NAME)

            call hdf5Dump%init(filename=IOName%FullName%str, mode='write')
            call hdf5Dump%open()

            Associate(zstart => FO%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                      zend   => FO%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                      rstart => FO%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                      rend   => FO%DM%CornerIndex(BOUNDARY_RIGHT, 2), &
                      zRightImageType => FO%DM%NeighborType(BOUNDARY_RIGHT, 1), &
                      rRightImageType => FO%DM%NeighborType(BOUNDARY_RIGHT, 2), &
                      image_shape => FO%DM%ImageShape)

                zend_on_bound = zend
                rend_on_bound = rend
                if (image_shape(1) > 1 .and. NEIGHBOR_TYPE_DOMAIN == zRightImageType) zend_on_bound = zend - 1
                if (image_shape(2) > 1 .and. NEIGHBOR_TYPE_DOMAIN == rRightImageType) rend_on_bound = rend - 1

                call hdf5Dump%write('/Rho', FO%RhoOne(zstart:zend_on_bound, rstart:rend_on_bound), chunkdim=image_shape)
                call hdf5Dump%write('/Jon', FO%JonOne(zstart:zend_on_bound, rstart:rend_on_bound), chunkdim=image_shape)
                call hdf5Dump%write('/Chi', FO%ChiOne(zstart:zend_on_bound, rstart:rend_on_bound), chunkdim=image_shape)

            end Associate

            call hdf5Dump%close()

            return
        end subroutine DiagFieldOne


        subroutine WeightingFieldOneElectrostatic(FO, PB, Geom)
            Class(FieldOne), intent(inout) :: FO
            Type(ParticleBundle), intent(in) :: PB
            type(Geometry), intent(inout) :: Geom
            integer(4) :: i
            integer(4) :: NzL, NzU, NrL, NrU
            real(8)    :: ZL, ZU, RL, RU
            real(8)    :: Volume
            
            call FO%Zero()
            
            Associate (zstart => FO%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                       zend   => FO%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                       rstart => FO%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                       rend   => FO%DM%CornerIndex(BOUNDARY_RIGHT, 2), &
                       Npar   => PB%Npar, &
                       Rho    => FO%RhoOne, Jon => FO%JonOne, Chi => FO%ChiOne,&
                       dx     => FO%DM%SpaceStep(1))

                do i = 1, Npar
                    NzL = Ceiling(PB%PO(i)%Z)
                    NzU = NzL + 1
                    ZL  = Dble(NzL) - PB%PO(i)%Z
                    ZU  = 1.d0 - ZL

                    NrL = Ceiling(PB%PO(i)%R)
                    NrU = NrL + 1
                    RL  = Dble(NrL) - PB%PO(i)%R
                    RU  = 1.d0 - RL

                    Rho(NzL, NrL) = Rho(NzL, NrL) + ZL * RL * PB%PO(i)%WQ
                    Rho(NzL, NrU) = Rho(NzL, NrU) + ZL * RU * PB%PO(i)%WQ
                    Rho(NzU, NrL) = Rho(NzU, NrL) + ZU * RL * PB%PO(i)%WQ
                    Rho(NzU, NrU) = Rho(NzU, NrU) + ZU * RU * PB%PO(i)%WQ
                end do

                Rho(zstart:zend, rstart:rend) = Rho(zstart:zend, rstart:rend) * Geom%point_inv_volume(zstart:zend, rstart:rend)
                Rho = Rho * FO%RhoFactor

                if (pic_type_explicit == pic_type) then
                    Chi = 0.d0
                else if (pic_type_implicit == pic_type) then
                    Chi = Rho * FO%Chifactor
                end if

            end Associate
            
            return
        end subroutine WeightingFieldOneElectrostatic


        subroutine DestroyFieldOne(FO)
            Class(FieldOne),intent(inout) :: FO

            if (Allocated(FO%RhoOne)) Deallocate(FO%RhoOne)
            if (Allocated(FO%JonOne)) Deallocate(FO%JonOne)
            if (Allocated(FO%ChiOne)) Deallocate(FO%ChiOne)

            return
        end subroutine DestroyFieldOne

end Module ModuleFieldOne
