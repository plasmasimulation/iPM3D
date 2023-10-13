Module ModuleParticleBundleIndex
    use ModuleParticleOneIndex

    implicit none

    type :: ParticleBundleIndex
        integer(4) :: NPar
        integer(4) :: size
        integer(4) :: npar_min, npar_max
        real(8) :: k_lb, k_dec, k_inc

        type(ParticleOneIndex), allocatable :: POI(:)

    contains

        procedure :: init => initParticleBundleIndex
        procedure :: realloc => reallocParticleBundleIndex
        procedure :: destroy => destroyParticleBundleIndex
        procedure :: addone => addParticleOneIndex
        procedure :: delone => delParticleOneIndex
        procedure :: adjust => adjustParticleBundleIndex

    end type ParticleBundleIndex

    contains

        subroutine initParticleBundleIndex(this, init_size, min_size, max_size, k_lb, k_dec, k_inc)
            class(ParticleBundleIndex), intent(inout) :: this
            integer(4), optional, intent(in) :: init_size, min_size, max_size
            real(8), optional, intent(in) :: k_lb, k_dec, k_inc

            this%NPar = 0
            this%size = 1000
            this%npar_min = 10
            this%npar_max = 10000
            this%k_lb = 0.25d0
            this%k_dec = 0.5d0
            this%k_inc = 1.5d0

            if (present(init_size) .and. init_size > 0) this%size = init_size
            if (present(min_size)  .and. min_size > 0)  this%npar_min = min_size
            if (present(max_size)  .and. max_size > 0)  this%npar_max = max_size

            if (present(k_lb)  .and. k_lb > 0.d0 .and. k_lb <= 1.d0) this%k_lb = k_lb
            if (present(k_dec) .and. k_dec > 0.d0 .and. k_dec <= 1.d0) this%k_dec = k_dec
            if (present(k_inc) .and. k_inc >= 1.d0) this%k_inc = k_inc

            call this%destroy()
            call this%realloc()

        end subroutine initParticleBundleIndex


        subroutine reallocParticleBundleIndex(this, newsize)
            class(ParticleBundleIndex), intent(inout) :: this
            integer(4), optional, intent(in) :: newsize
            integer(4) :: size_tmp
            type(ParticleOneIndex), allocatable :: temp(:)
            
            if (allocated(this%POI)) then
                if (present(newsize)) then
                    size_tmp = newsize

                    if (size_tmp < this%npar_min) size_tmp = this%npar_min
                    if (size_tmp > this%npar_max) size_tmp = this%npar_max

                    allocate(temp(size_tmp))
                    if (size_tmp <= this%size) then
                        temp(1:size_tmp) = this%POI(1:size_tmp)
                        deallocate(this%POI)
                        call move_alloc(temp, this%POI)

                    else if (size_tmp > this%size) then
                        temp(1:this%size) = this%POI(1:this%size)
                        deallocate(this%POI)
                        call move_alloc(temp, this%POI)

                    end if

                    if (allocated(temp)) deallocate(temp)
                    this%size = size_tmp

                else
                    deallocate(this%POI)
                    if (this%size > 0) then
                        allocate(this%POI(this%size))
                    else
                        write(*, *) "The particle bundle size needs to be greater than 0."
                        stop
                    end if
                end if

            else
                if (this%size > 0) then
                    allocate(this%POI(this%size))
                else
                    write(*, *) "The particle bundle size needs to be greater than 0."
                    stop
                end if
            end if

        end subroutine reallocParticleBundleIndex


        subroutine destroyParticleBundleIndex(this)
            class(ParticleBundleIndex), intent(inout) :: this

            if (allocated(this%POI)) deallocate(this%POI)
            this%NPar = 0

        end subroutine destroyParticleBundleIndex


        subroutine addParticleOneIndex(this, one)
            class(ParticleBundleIndex), intent(inout) :: this
            type(ParticleOneIndex), intent(in) :: one

            if (this%NPar+1 > this%size) call this%adjust(this%NPar+1)
            call this%POI(this%NPar+1)%Copy(one)
            this%NPar = this%NPar + 1

        end subroutine addParticleOneIndex


        subroutine delParticleOneIndex(this, index)
            class(ParticleBundleIndex), intent(inout) :: this
            integer(4), intent(in) :: index

            if (this%NPar > 0) then
                if (index >= 1 .and. index <= this%NPar) then
                    call this%POI(index)%Copy(this%POI(this%NPar))
                    this%NPar = this%NPar - 1                    
                    call this%adjust()

                else
                    write(*, *) "The index is out of the range."
                    stop
                end if

            else
                write(*, *) "No particle in particle bundle can be deleted."
                stop
            end if

        end subroutine delParticleOneIndex


        subroutine adjustParticleBundleIndex(this, newsize)
            class(ParticleBundleIndex), intent(inout) :: this
            integer(4), optional, intent(in) :: newsize
            integer(4) :: tmp_size

            if (present(newsize) .and. newsize > 0) then
                if (newsize >= this%size) then
                    tmp_size = max(newsize, min(this%npar_max, int(this%k_inc*dble(this%size))))
                    if (tmp_size > this%npar_max) then
                        write(*, *) "The size of particle bundle reaches the maximum."
                        stop
                    end if

                    call this%realloc(tmp_size)
                else
                    call this%realloc(newsize)
                end if

            else
                if (this%NPar == 0) then
                    call this%realloc(this%npar_min)
                else if (this%NPar < int(this%k_lb * dble(this%size))) then
                    tmp_size = int(max(this%k_lb, this%k_dec)*dble(this%size))
                    call this%realloc(tmp_size)
                end if
            end if

        end subroutine adjustParticleBundleIndex

end Module ModuleParticleBundleIndex
