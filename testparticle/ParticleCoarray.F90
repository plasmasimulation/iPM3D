Module ModuleParallelParticle
    use ModuleParticleBundle
    use ModuleDomain
    use ModuleParticleBoundary1D

    implicit none

contains

    subroutine CommunicationParticleCoarray(PBound, PB, Ns)
        integer(4), intent(in) :: Ns
        Type(ParticleBoundary), intent(inout) :: PBound(0:Ns)[*]
        Type(ParticleBundle), intent(inout) :: PB(0:Ns)
        integer(4) :: k, NCom, dim1D
        integer(4) :: leftImageId, rightImageId, leftImageType, rightImageType

        do k = 0, Ns
            call PBound(k)%AborpZ(PB(k))
        end do

        sync all

        ! Z
        dim1D = 1
        do k = 0, Ns
            leftImageId    = PBound(k)%DM%NeighborImageId(BOUNDARY_LEFT, dim1D)
            rightImageId   = PBound(k)%DM%NeighborImageId(BOUNDARY_RIGHT, dim1D)
            leftImageType  = PBound(k)%DM%NeighborType(BOUNDARY_LEFT, dim1D)
            rightImageType = PBound(k)%DM%NeighborType(BOUNDARY_RIGHT, dim1D)

            if (NEIGHBOR_TYPE_DOMAIN == leftImageType) then
                NCom = PBound(k)[leftImageId]%PBDump(BOUNDARY_RIGHT, dim1D)%NPar
                
                if (NCom > 0) then

#ifdef OPEN_CHECK_PARTICLE_OVERSTEP
                    call LogOverstep(PB(k)%NPar, NCom, PB(k)%NParNormal, k)
#endif

                    PB(k)%PO(PB(k)%NPar+1 : PB(k)%NPar+NCom) = PBound(k)[leftImageId]%PBDump(BOUNDARY_RIGHT, dim1D)%PO(1:NCom)
                    PB(k)%NPar = PB(k)%NPar+NCom
                end if
            end if

            if (NEIGHBOR_TYPE_DOMAIN == rightImageType) then
                NCom = PBound(k)[rightImageId]%PBDump(BOUNDARY_LEFT, dim1D)%NPar

                if (NCom > 0) then

#ifdef OPEN_CHECK_PARTICLE_OVERSTEP
                    call LogOverstep(PB(k)%NPar, NCom, PB(k)%NParNormal, k)
#endif

                    PB(k)%PO(PB(k)%NPar+1 : PB(k)%NPar+NCom) = PBound(k)[rightImageId]%PBDump(BOUNDARY_LEFT, dim1D)%PO(1:NCom)
                    PB(k)%NPar = PB(k)%NPar+NCom
                end if
            end if
        end do

        sync all

        do k = 0, Ns
            call PBound(k)%AborpR(PB(k))
        end do

        ! R
        dim1D = 2
        do k = 0, Ns
            leftImageId    = PBound(k)%DM%NeighborImageId(BOUNDARY_LEFT, dim1D)
            rightImageId   = PBound(k)%DM%NeighborImageId(BOUNDARY_RIGHT, dim1D)
            leftImageType  = PBound(k)%DM%NeighborType(BOUNDARY_LEFT, dim1D)
            rightImageType = PBound(k)%DM%NeighborType(BOUNDARY_RIGHT, dim1D)

            if (NEIGHBOR_TYPE_DOMAIN == leftImageType) then
                NCom = PBound(k)[leftImageId]%PBDump(BOUNDARY_RIGHT, dim1D)%NPar

#ifdef OPEN_CHECK_PARTICLE_OVERSTEP
                    call LogOverstep(PB(k)%NPar, NCom, PB(k)%NParNormal, k)
#endif

                if (NCom > 0) then
                    PB(k)%PO(PB(k)%NPar+1 : PB(k)%NPar+NCom) = PBound(k)[leftImageId]%PBDump(BOUNDARY_RIGHT, dim1D)%PO(1:NCom)
                    PB(k)%NPar = PB(k)%NPar+NCom
                end if
            end if

            if (NEIGHBOR_TYPE_DOMAIN == rightImageType) then
                NCom = PBound(k)[rightImageId]%PBDump(BOUNDARY_LEFT, dim1D)%NPar
                if (NCom > 0) then

#ifdef OPEN_CHECK_PARTICLE_OVERSTEP
                    call LogOverstep(PB(k)%NPar, NCom, PB(k)%NParNormal, k)
#endif

                    PB(k)%PO(PB(k)%NPar+1 : PB(k)%NPar+NCom) = PBound(k)[rightImageId]%PBDump(BOUNDARY_LEFT, dim1D)%PO(1:NCom)
                    PB(k)%NPar = PB(k)%NPar+NCom
                end if
            end if
        end do

#ifdef OPEN_CHECK_PARTICLE_OVERSTEP
        sync all
        call CheckParticlesAcrossDomains(PB, Ns, PBound(0)%DM)
#endif

        return
    end subroutine CommunicationParticleCoarray


    subroutine CommunicationParticleCoarray1DZ(PBound, PB)
        Type(ParticleBoundary), intent(inout) :: PBound[*]
        Type(ParticleBundle), intent(inout) :: PB
        integer(4) :: NCom
    
        Associate(leftImageId    => PBound%DM%NeighborImageId(BOUNDARY_LEFT, 1), &
                  rightImageId   => PBound%DM%NeighborImageId(BOUNDARY_RIGHT, 1), &
                  leftImageType  => PBound%DM%NeighborType(BOUNDARY_LEFT, 1), &
                  rightImageType => PBound%DM%NeighborType(BOUNDARY_RIGHT, 1))

            if (NEIGHBOR_TYPE_DOMAIN == leftImageType) then
                NCom = PBound[leftImageId]%PBDump(BOUNDARY_RIGHT, 1)%NPar
                if (NCom > 0) then
                    PB%PO(PB%NPar+1 : PB%NPar+NCom) = PBound[leftImageId]%PBDump(BOUNDARY_RIGHT, 1)%PO(1:NCom)
                    PB%NPar = PB%NPar+NCom
                end if
            end if

            if (NEIGHBOR_TYPE_DOMAIN == rightImageType) then
                NCom = PBound[rightImageId]%PBDump(BOUNDARY_LEFT, 1)%NPar
                if (NCom > 0) then
                    PB%PO(PB%NPar+1 : PB%NPar+NCom) = PBound[rightImageId]%PBDump(BOUNDARY_LEFT, 1)%PO(1:NCom)
                    PB%NPar = PB%NPar+NCom
                end if
            end if

        end Associate

        return
    end subroutine CommunicationParticleCoarray1DZ


    subroutine CommunicationParticleCoarray1DR(PBound, PB)
        Type(ParticleBoundary), intent(inout) :: PBound[*]
        Type(ParticleBundle), intent(inout) :: PB
        integer(4) :: NCom

        Associate(leftImageId    => PBound%DM%NeighborImageId(BOUNDARY_LEFT, 2), &
                  rightImageId   => PBound%DM%NeighborImageId(BOUNDARY_RIGHT, 2), &
                  leftImageType  => PBound%DM%NeighborType(BOUNDARY_LEFT, 2), &
                  rightImageType => PBound%DM%NeighborType(BOUNDARY_RIGHT, 2))

            if (NEIGHBOR_TYPE_DOMAIN == leftImageType) then
                NCom = PBound[leftImageId]%PBDump(BOUNDARY_RIGHT, 2)%NPar
                if (NCom > 0) then
                    PB%PO(PB%NPar+1 : PB%NPar+NCom) = PBound[leftImageId]%PBDump(BOUNDARY_RIGHT, 2)%PO(1:NCom)
                    PB%NPar = PB%NPar+NCom
                end if
            end if

            if (NEIGHBOR_TYPE_DOMAIN == rightImageType) then
                NCom = PBound[rightImageId]%PBDump(BOUNDARY_LEFT, 2)%NPar
                if (NCom > 0) then
                    PB%PO(PB%NPar+1 : PB%NPar+NCom) = PBound[rightImageId]%PBDump(BOUNDARY_LEFT, 2)%PO(1:NCom)
                    PB%NPar = PB%NPar+NCom
                end if
            end if

        end Associate

        return
    end subroutine CommunicationParticleCoarray1DR


    function getParticleNum(Ns, PB)
        integer(4), intent(in) :: Ns
        Type(ParticleBundle), intent(inout) :: PB(0:Ns)
        integer(4) :: getParticleNum(0:Ns)
        integer(4), save :: Npar[*]
        integer(4) :: k
        
        do k = 0, Ns
            Npar = PB(k)%NPar
            call co_sum(Npar)
            getParticleNum(k) = Npar
        end do

    end


    subroutine LogOverstep(NPar, NCom, NParNormal, k)
        integer(4), intent(in) :: NPar, NCom, NParNormal, k

        if (dble(NPar+NCom) > PARTICLE_NUMBER_MAX_RATIO * dble(NParNormal)) then
            write(*, '(a, i, 4x, a, i, 4x, a, i4, 2x, a, i4, 2x, a)') &
                    'The number ', NPar, ' + ', NCom, &
                    ' of particle ', k, ' in image ', this_image(), ' exceeds the limit.'
        end if

        return
    end subroutine LogOverstep


    subroutine CheckParticlesAcrossDomains(PB, Ns, dm)
        integer(4), intent(in) :: Ns
        Type(ParticleBundle), intent(inout) :: PB(0:Ns)
        Type(Domain), intent(in) :: dm
        Type(ParticleOne), save :: POtmp[*]
        integer(4) :: k, i, targetId

        do k = 0, Ns
            do i = 1, PB(k)%NPar
                targetId = getParticleImageId(PB(k)%PO(i), dm)
                if (targetId /= dm%MyId) then

                    write(*, '(i, i, i, f, f)') k, dm%MyId, targetId, PB(k)%PO(i)%Z, PB(k)%PO(i)%R

                end if
            end do
        end do

    end subroutine CheckParticlesAcrossDomains


    function getParticleImageId(PO, dm)
        Class(ParticleOne), intent(in) :: PO
        Type(Domain), intent(in) :: dm
        integer(4) :: getParticleImageId
        integer(4) :: dmindex(1:2), i
        real(8) :: tmp
        real(8) :: zmin, zmax, rmin, rmax

        zmin = dble(dm%CornerIndex(BOUNDARY_LEFT, 1)  - 1)
        zmax = dble(dm%CornerIndex(BOUNDARY_RIGHT, 1) - 1)
        rmin = dble(dm%CornerIndex(BOUNDARY_LEFT, 2)  - 1)
        rmax = dble(dm%CornerIndex(BOUNDARY_RIGHT, 2) - 1)

        if (PO%Z > zmin .and. PO%Z < zmax .and. PO%R > rmin .and. PO%R < rmax) then
            getParticleImageId = dm%MyId

        else if (PO%Z <= 0.d0 .or. PO%Z >= dble(dm%GlobalShape(1)) .or. &
                 PO%R <= 0.d0 .or. PO%R >= dble(dm%GlobalShape(2))) then
            getParticleImageId = 0

        else
            dmindex = 1

            ! z image index
            tmp = dble(dm%LocalShapeList(1, 1) - 1)
            do i = 2, dm%ImageShape(1)
                if (PO%Z < tmp) then
                    exit
                else
                    dmindex(1) = dmindex(1) + 1
                    tmp = tmp + dble(dm%LocalShapeList(1, i) - 1)
                end if
            end do

            ! r image index
            tmp = dble(dm%LocalShapeList(2, 1) - 1)
            do i = 2, dm%ImageShape(2)
                if (PO%R < tmp) then
                    exit
                else
                    dmindex(2) = dmindex(2) + 1
                    tmp = tmp + dble(dm%LocalShapeList(2, i) - 1)
                end if
            end do

            getParticleImageId = getImageIndex(dmindex, dm%ImageShape, 2)

        end if

    end

end Module ModuleParallelParticle
