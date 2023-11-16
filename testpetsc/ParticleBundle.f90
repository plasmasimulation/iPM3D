Module ModuleParticleBundle
    use ModuleParticleOne
    ! use ModuleFieldEM
    ! use ModuleBspline
     !use ModuleParallelDump

    implicit none

     real(8), save :: particle_number_max_ratio = 0.d0
     real(8), save :: particle_limit_max_velocity = 0.d0
    real(8), save :: particle_number_lower_boundary_ratio = 0.25d0
    real(8), save :: particle_number_decrease_ratio = 0.5d0
    real(8), save :: particle_number_increase_ratio = 1.5d0

    integer(4), save :: init_z_start = 0
    integer(4), save :: init_r_start = 0
    integer(4), save :: init_z_stop  = 0
    integer(4), save :: init_r_stop  = 0

    Type :: ParticleBundle
        ! Type(FileName) :: IOName
        integer(4) :: NPar
        integer(4) :: size
        integer(4) :: npar_min, npar_max
        real(8) :: k_lb, k_dec, k_inc

        real(8) :: XFactor, VFactor
        logical :: LXScaled=.False., LVScaled=.False.
        real(8) :: Charge, Mass, Weight
        real(8) :: dx, dt

        Type(ParticleOne), Allocatable :: PO(:)
        ! Type(Bspline1N2D), Allocatable :: BS(:)

    contains

       ! procedure :: AllInit => AllInitializationParticleBundle

        procedure :: Init => initParticleBundle
        procedure :: Adjust => adjustParticleBundle

        procedure :: Destroy => DestroyParticleBundle
        procedure :: Realloc => reallocParticleBundle

        procedure :: AddOne  => AddParticleOneParticleBundle
        procedure :: AddBun  => AddParticleBundleParticleBundle
        procedure :: DelOne  => DelParticleOneParticleBundle
        procedure :: DelOne2  => DelParticleOne2ParticleBundle
        procedure :: DelOneReady  => DelParticleReadyParticleBundle
        !procedure :: NumRes2=>ParticleBundleNumberRescale2

        !procedure :: PosRes  => PositionRescaleParticleBundle
        ! procedure :: VelRes  => VelocityRescaleParticleBundle
        ! procedure :: MoveES  => MoveElectrostaticParticleBundle

        ! procedure :: MoveESEx  => MoveElectrostaticParticleBundleExplicit
        ! procedure :: MoveESIm  => MoveElectrostaticParticleBundleImplicit

        ! procedure :: MoveEM  => MoveElectroMagneticParticleBundle
        ! !procedure :: WeightP2C=>WeightP2CParticleBundle

        ! procedure :: VelLimit => VelocityLimitParticleBundle

        ! procedure :: Dump => DumpParticleBundle
        ! procedure :: Load => LoadParticleBundle
        procedure :: Copy => CopyParticleBundle

        ! procedure :: Diag => DiagParticleBundle

        ! procedure, private :: DumpParticleBundleHDF5
        ! procedure, private :: DumpParticleBundleDat
        ! procedure, private :: LoadParticleBundleHDF5
        ! procedure, private :: LoadParticleBundleDat

        ! procedure :: Norm => ParticleBundleNormalization
        ! !procedure :: RescaleFieldOne=>RescaleFieldOneParticleBundle

    end Type ParticleBundle

    ! Type(HDF5_PDump), private :: hdf5ParticleBundleDump

    contains
    subroutine initParticleBundle(this, init_size, min_size, max_size, k_lb, k_dec, k_inc)
        class(ParticleBundle), intent(inout) :: this
        integer(4), optional, intent(in) :: init_size, min_size, max_size
        real(8), optional, intent(in) :: k_lb, k_dec, k_inc

        this%npar = 0
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

        call this%Destroy()
        call this%Realloc()

    end subroutine initParticleBundle


    subroutine CopyParticleBundle(this, src)
        class(ParticleBundle), intent(inout) :: this
        Type(ParticleBundle), intent(in) :: src

        this%NPar = src%NPar
        call this%Init(src%size, src%npar_min, src%npar_max, src%k_lb, src%k_dec, src%k_inc)

        this%XFactor = src%XFactor
        this%VFactor = src%VFactor
        this%LXScaled = src%LXScaled
        this%LVScaled = src%LVScaled

        this%Charge = src%Charge
        this%Mass = src%Mass
        this%Weight = src%Weight
        this%dx = src%dx
        this%dt = src%dt

    end subroutine CopyParticleBundle


    subroutine reallocParticleBundle(this, newsize)
        class(ParticleBundle), intent(inout) :: this
        integer(4), optional, intent(in) :: newsize
        integer(4) :: size_tmp
        type(ParticleOne), allocatable :: temp(:)
        
        if (allocated(this%PO)) then
            if (present(newsize)) then
                size_tmp = newsize

                if (size_tmp < this%npar_min) size_tmp = this%npar_min
                if (size_tmp > this%npar_max) size_tmp = this%npar_max

                allocate(temp(size_tmp))
                if (size_tmp <= this%size) then
                    temp(1:size_tmp) = this%PO(1:size_tmp)
                    deallocate(this%PO)
                    call move_alloc(temp, this%PO)

                else if (size_tmp > this%size) then
                    temp(1:this%size) = this%PO(1:this%size)
                    deallocate(this%PO)
                    call move_alloc(temp, this%PO)

                end if

                if (allocated(temp)) deallocate(temp)
                this%size = size_tmp

            else
                deallocate(this%PO)
                if (this%size > 0) then
                    allocate(this%PO(this%size))
                else
                    write(*, *) "The particle bundle size needs to be greater than 0."
                    stop
                end if
            end if

        else
            if (this%size > 0) then
                allocate(this%PO(this%size))
            else
                write(*, *) "The particle bundle size needs to be greater than 0."
                stop
            end if
        end if

    end subroutine reallocParticleBundle


    subroutine DestroyParticleBundle(this)
        class(ParticleBundle), intent(inout) :: this

        if(allocated(this%PO)) deallocate(this%PO)
        ! if(allocated(this%BS)) deallocate(this%BS)

        this%NPar = 0

    end subroutine DestroyParticleBundle


    subroutine AddParticleOneParticleBundle(this, PO)
        class(ParticleBundle),intent(inout) ::  this
        class(ParticleOne),intent(in) ::  PO

        if (this%NPar+1 > this%size) call this%adjust(this%NPar+1)
        this%PO(this%NPar+1) = PO
        this%NPar = this%NPar + 1

    end subroutine AddParticleOneParticleBundle


    subroutine AddParticleBundleParticleBundle(this, bun)
        class(ParticleBundle), intent(inout) :: this
        class(ParticleBundle), intent(in) :: bun

       if (this%NPar+bun%NPar > this%size) call this%adjust(this%NPar+bun%NPar)
        this%PO(this%NPar+1:this%NPar+bun%NPar) = bun%PO(1:bun%NPar)
        this%NPar = this%NPar + bun%NPar

    end subroutine AddParticleBundleParticleBundle


    subroutine DelParticleOneParticleBundle(this, NDel)
        class(ParticleBundle),intent(inout) :: this
        integer(4),intent(in) :: NDel

        if (this%NPar > 0) then
            if (NDel >= 1 .and. NDel <= this%NPar) then
                this%PO(NDel) = this%PO(this%NPar)
                this%NPar = this%NPar - 1                    
                call this%adjust()

            else
                write(*, *) "The NDel is out of the range."
                stop
            end if

        else
            write(*, *) "No particle in particle bundle can be deleted."
            stop
        end if

    end subroutine DelParticleOneParticleBundle
    subroutine  DelParticleReadyParticleBundle(this)
        class(ParticleBundle),intent(inout) :: this
       integer(4):: i
        if (this%NPar > 0) then
            i=1
            ! write(*,*) "this%NPar before",this%NPar
            do while(i<this%NPar-1) 
                if(this%PO(i)%tempdelflag==1) then
                       if(this%PO(this%NPar)%tempdelflag==1) then
                        do while(this%PO(this%NPar)%tempdelflag==1)
                             this%NPar = this%NPar - 1 
                        end do  
                        end if
                       
                    this%PO(i) = this%PO(this%NPar)
                    this%NPar = this%NPar - 1
                    
                end if
                i=i+1
            end do
            if(this%PO(this%NPar)%tempdelflag==1) then
                this%NPar = this%NPar - 1
            end if
            ! write(*,*) "this%NPar after",this%NPar
            call this%adjust()

        end if
    end subroutine DelParticleReadyParticleBundle

    subroutine DelParticleOne2ParticleBundle(this, NDel)
        class(ParticleBundle),intent(inout) :: this
        integer(4),intent(in) :: NDel

        if (this%NPar > 0) then
            if (NDel >= 1 .and. NDel <= this%NPar) then
                this%PO(NDel)%tempdelflag = 1
                !this%NPar = this%NPar - 1                    
                !call this%adjust()

            else
                write(*, *) "The NDel 2is out of the range.",NDel,this%NPar

                !stop
            end if

        else
            write(*, *) "No particle in particle bundle can be deleted."
            stop
        end if

    end subroutine DelParticleOne2ParticleBundle


    subroutine adjustParticleBundle(this, newsize)
        class(ParticleBundle), intent(inout) :: this
        integer(4), optional, intent(in) :: newsize
        integer(4) :: tmp_size

        if (present(newsize) .and. newsize > 0) then
            if (newsize >= this%size) then
                tmp_size = max(newsize, min(this%npar_max, int(this%k_inc*dble(this%size))))
                if (tmp_size > this%npar_max) then
                    write(*, *) "The size of particle bundle reaches the maximum."
                    stop
                end if

                call this%Realloc(tmp_size)
            else
                call this%Realloc(newsize)
            end if

        else
            if (this%NPar == 0) then
                call this%Realloc(this%npar_min)
            else if (this%NPar < int(this%k_lb * dble(this%size))) then
                tmp_size = int(max(this%k_lb, this%k_dec)*dble(this%size))
                call this%Realloc(tmp_size)
            end if
        end if

    end subroutine adjustParticleBundle

        ! subroutine DumpParticleBundle(this)
        !     class(ParticleBundle), intent(inout) :: this

        !     if (FILE_EXTENSION_MODE_H5 == this%IOName%ExtensionMode) then
        !         call this%DumpParticleBundleHDF5()

        !     else if (FILE_EXTENSION_MODE_DAT == this%IOName%ExtensionMode) then
        !         call this%DumpParticleBundleDat()
            
        !     else
        !         call this%DumpParticleBundleHDF5()

        !     end if

        ! end subroutine DumpParticleBundle


        ! subroutine LoadParticleBundle(this, Status)
		! 	class(ParticleBundle), intent(inout) :: this
        !     logical, intent(inout), optional :: Status
        !     logical :: alive

        !     Inquire(file=this%IOName%FullName%str, exist=alive)
        !     if (alive) then
        !         if (FILE_EXTENSION_MODE_H5 == this%IOName%ExtensionMode) then
        !             call this%LoadParticleBundleHDF5()

        !         else if (FILE_EXTENSION_MODE_DAT == this%IOName%ExtensionMode) then
        !             call this%LoadParticleBundleDat()
                
        !         end if
        !     end if

        !     if (present(Status)) Status = alive

		! end subroutine LoadParticleBundle


        ! subroutine DumpParticleBundleHDF5(this)
		! 	class(ParticleBundle), intent(inout) :: this
        !     integer(4) :: size, rank, ierr

        !     call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        !     call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        !     if (0 == rank) write(*, '(a)') "The ParticleBundle is save as hdf5 file."

        !     call hdf5ParticleBundleDump%init(filename=this%IOName%FullName%str, mode='write', serial=.True.)
        !     call hdf5ParticleBundleDump%open()

        !     call hdf5ParticleBundleDump%writeattr('/', 'NPar', this%NPar)
        !     call hdf5ParticleBundleDump%writeattr('/', 'size', this%size)
        !     call hdf5ParticleBundleDump%writeattr('/', 'npar_min', this%npar_min)
        !     call hdf5ParticleBundleDump%writeattr('/', 'npar_max', this%npar_max)
        !     call hdf5ParticleBundleDump%writeattr('/', 'k_lb', this%k_lb)
        !     call hdf5ParticleBundleDump%writeattr('/', 'k_dec', this%k_dec)
        !     call hdf5ParticleBundleDump%writeattr('/', 'k_inc', this%k_inc)
        !     call hdf5ParticleBundleDump%writeattr('/', 'Charge', this%Charge)
        !     call hdf5ParticleBundleDump%writeattr('/', 'Mass', this%Mass)
        !     call hdf5ParticleBundleDump%writeattr('/', 'Weight', this%Weight)
        !     call hdf5ParticleBundleDump%writeattr('/', 'XFactor', this%XFactor)
        !     call hdf5ParticleBundleDump%writeattr('/', 'VFactor', this%VFactor)
        !     call hdf5ParticleBundleDump%writeattr('/', 'LXScaled', this%LXScaled)
        !     call hdf5ParticleBundleDump%writeattr('/', 'LVScaled', this%LVScaled)

        !     call hdf5ParticleBundleDump%write('X', this%PO%X)
        !     call hdf5ParticleBundleDump%write('Y', this%PO%Y)
        !     call hdf5ParticleBundleDump%write('Z', this%PO%Z)
        !     call hdf5ParticleBundleDump%write('R', this%PO%R)
        !     call hdf5ParticleBundleDump%write('Vx', this%PO%Vx)
        !     call hdf5ParticleBundleDump%write('Vy', this%PO%Vy)
        !     call hdf5ParticleBundleDump%write('Vz', this%PO%Vz)
        !     call hdf5ParticleBundleDump%write('Ax', this%PO%Ax)
        !     call hdf5ParticleBundleDump%write('Ay', this%PO%Ay)
        !     call hdf5ParticleBundleDump%write('Az', this%PO%Az)
        !     call hdf5ParticleBundleDump%write('WQ', this%PO%WQ)

        !     call hdf5ParticleBundleDump%close()

        ! end subroutine DumpParticleBundleHDF5


        ! subroutine DumpParticleBundleDat(this)
        !     class(ParticleBundle), intent(inout) :: this
        !     integer(4) :: size, rank, ierr

        !     call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        !     call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        !     if (0 == rank) write(*, '(a)') "The ParticleBundle is save as hdf5 file."

        !     open(10, file=this%IOName%FullName%str)
        !         write(10, *) this%NPar, this%size
        !         write(10, *) this%npar_min, this%npar_max
        !         write(10, *) this%k_lb, this%k_dec, this%k_inc
        !         write(10, *) this%Charge, this%Mass, this%Weight
        !         write(10, *) this%XFactor, this%VFactor
        !         write(10, *) this%LXScaled, this%LVScaled
        !         write(10, *) this%PO%X
        !         write(10, *) this%PO%Y
        !         write(10, *) this%PO%Z
        !         write(10, *) this%PO%R
        !         write(10, *) this%PO%Vx
        !         write(10, *) this%PO%Vy
        !         write(10, *) this%PO%Vz
        !         write(10, *) this%PO%Ax
        !         write(10, *) this%PO%Ay
        !         write(10, *) this%PO%Az
        !         write(10, *) this%PO%WQ
        !     close(10)

        ! end subroutine DumpParticleBundleDat


        ! subroutine LoadParticleBundleHDF5(this)
        !     class(ParticleBundle), intent(inout) :: this
        !     real(8),ALLOCATABLE :: Temp1D(:)
        !     integer(4) :: size, rank, ierr, npar

        !     call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        !     call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        !     if (0 == rank) write(*, '(a)') "The ParticleBundle is load from hdf5 file."

        !     call hdf5ParticleBundleDump%init(filename=this%IOName%FullName%str, mode='read')
        !     call hdf5ParticleBundleDump%open()

        !     call hdf5ParticleBundleDump%readattr('/', 'NPar', npar)
        !     call hdf5ParticleBundleDump%readattr('/', 'size', this%size)
        !     ! call hdf5ParticleBundleDump%readattr('/', 'npar_min', this%npar_min)
        !     ! call hdf5ParticleBundleDump%readattr('/', 'npar_max', this%npar_max)
        !     ! call hdf5ParticleBundleDump%readattr('/', 'k_lb', this%k_lb)
        !     ! call hdf5ParticleBundleDump%readattr('/', 'k_dec', this%k_dec)
        !     ! call hdf5ParticleBundleDump%readattr('/', 'k_inc', this%k_inc)
        !     call hdf5ParticleBundleDump%readattr('/', 'Charge', this%Charge)
        !     call hdf5ParticleBundleDump%readattr('/', 'Mass', this%Mass)
        !     call hdf5ParticleBundleDump%readattr('/', 'Weight', this%Weight)
        !     call hdf5ParticleBundleDump%readattr('/', 'XFactor', this%XFactor)
        !     call hdf5ParticleBundleDump%readattr('/', 'VFactor', this%VFactor)
        !     call hdf5ParticleBundleDump%readattr('/', 'LXScaled', this%LXScaled)
        !     call hdf5ParticleBundleDump%readattr('/', 'LVScaled', this%LVScaled)

        !     call this%Init(this%size, this%npar_min, this%npar_max, &
        !                    this%k_lb, this%k_dec, this%k_inc)

        !     this%NPar = npar

        !     Allocate(Temp1D(this%size))
        !     call hdf5ParticleBundleDump%read('X', Temp1D)
        !     this%PO%X = Temp1D

        !     call hdf5ParticleBundleDump%read('Y', Temp1D)
        !     this%PO%Y = Temp1D

        !     call hdf5ParticleBundleDump%read('Z', Temp1D)
        !     this%PO%Z = Temp1D

        !     call hdf5ParticleBundleDump%read('R', Temp1D)
        !     this%PO%R = Temp1D

        !     call hdf5ParticleBundleDump%read('Vx', Temp1D)
        !     this%PO%Vx = Temp1D

        !     call hdf5ParticleBundleDump%read('Vy', Temp1D)
        !     this%PO%Vy = Temp1D

        !     call hdf5ParticleBundleDump%read('Vz', Temp1D)
        !     this%PO%Vz = Temp1D

        !     call hdf5ParticleBundleDump%read('Ax', Temp1D)
        !     this%PO%Ax = Temp1D

        !     call hdf5ParticleBundleDump%read('Ay', Temp1D)
        !     this%PO%Ay = Temp1D

        !     call hdf5ParticleBundleDump%read('Az', Temp1D)
        !     this%PO%Az = Temp1D

        !     call hdf5ParticleBundleDump%read('WQ', Temp1D)
        !     this%PO%WQ = Temp1D

        !     Deallocate(Temp1D)
        !     call hdf5ParticleBundleDump%close()

        ! end subroutine LoadParticleBundleHDF5


        ! subroutine LoadParticleBundleDat(this)
        !     class(ParticleBundle), intent(inout) :: this
        !     integer(4) :: size, rank, ierr, npar
        !     integer(4) :: tmpint1, tmpint2
        !     real(8) :: tmpreal1, tmpreal2, tmpreal3

        !     call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        !     call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        !     if (0 == rank) write(*, '(a)') "The ParticleBundle is load from dat file."

        !     open(10, file=this%IOName%FullName%str)
        !         read(10, *) npar, this%size
        !         read(10, *) tmpint1, tmpint2
        !         read(10, *) tmpreal1, tmpreal2, tmpreal3
        !         read(10, *) this%Charge, this%Mass, this%Weight
        !         read(10, *) this%XFactor, this%VFactor
        !         read(10, *) this%LXScaled, this%LVScaled
        !         call this%Init(this%size, this%npar_min, this%npar_max, &
        !                        this%k_lb, this%k_dec, this%k_inc)

        !         this%NPar = npar
        !         read(10, *) this%PO%X
        !         read(10, *) this%PO%Y
        !         read(10, *) this%PO%Z
        !         read(10, *) this%PO%R
        !         read(10, *) this%PO%Vx
        !         read(10, *) this%PO%Vy
        !         read(10, *) this%PO%Vz
        !         read(10, *) this%PO%Ax
        !         read(10, *) this%PO%Ay
        !         read(10, *) this%PO%Az
        !         read(10, *) this%PO%WQ
        !     close(10)

        ! end subroutine LoadParticleBundleDat


        ! subroutine DiagParticleBundle(this, name)
		! 	class(ParticleBundle), intent(inout) :: this
        !     Character(*), intent(in) :: name
        !     Type(FileName) :: tmpName

        !     call tmpName%Init(name, PathMode=FILE_PATH_MODE_DIAG, &
        !                             ExtensionMode=FILE_EXTENSION_MODE_H5, &
        !                             ParallelMode=FILE_PARALLEL_MODE_OPEN, &
        !                             DynamicIndex=FILE_DYNAMIC_MODE_CLOSE)

        !     call hdf5ParticleBundleDump%init(filename=tmpName%FullName%str, mode='write', serial=.True.)
        !     call hdf5ParticleBundleDump%open()

        !     call hdf5ParticleBundleDump%writeattr('/', 'NPar', this%NPar)
        !     call hdf5ParticleBundleDump%writeattr('/', 'size', this%size)
        !     call hdf5ParticleBundleDump%writeattr('/', 'npar_min', this%npar_min)
        !     call hdf5ParticleBundleDump%writeattr('/', 'npar_max', this%npar_max)
        !     call hdf5ParticleBundleDump%writeattr('/', 'k_lb', this%k_lb)
        !     call hdf5ParticleBundleDump%writeattr('/', 'k_dec', this%k_dec)
        !     call hdf5ParticleBundleDump%writeattr('/', 'k_inc', this%k_inc)
        !     call hdf5ParticleBundleDump%writeattr('/', 'Charge', this%Charge)
        !     call hdf5ParticleBundleDump%writeattr('/', 'Mass', this%Mass)
        !     call hdf5ParticleBundleDump%writeattr('/', 'Weight', this%Weight)
        !     call hdf5ParticleBundleDump%writeattr('/', 'XFactor', this%XFactor)
        !     call hdf5ParticleBundleDump%writeattr('/', 'VFactor', this%VFactor)
        !     call hdf5ParticleBundleDump%writeattr('/', 'LXScaled', this%LXScaled)
        !     call hdf5ParticleBundleDump%writeattr('/', 'LVScaled', this%LVScaled)

        !     call hdf5ParticleBundleDump%write('X', this%PO(1:this%NPar)%X)
        !     call hdf5ParticleBundleDump%write('Y', this%PO(1:this%NPar)%Y)
        !     call hdf5ParticleBundleDump%write('Z', this%PO(1:this%NPar)%Z)
        !     call hdf5ParticleBundleDump%write('R', this%PO(1:this%NPar)%R)
        !     call hdf5ParticleBundleDump%write('Vx', this%PO(1:this%NPar)%Vx)
        !     call hdf5ParticleBundleDump%write('Vy', this%PO(1:this%NPar)%Vy)
        !     call hdf5ParticleBundleDump%write('Vz', this%PO(1:this%NPar)%Vz)
        !     call hdf5ParticleBundleDump%write('Ax', this%PO(1:this%NPar)%Ax)
        !     call hdf5ParticleBundleDump%write('Ay', this%PO(1:this%NPar)%Ay)
        !     call hdf5ParticleBundleDump%write('Az', this%PO(1:this%NPar)%Az)
        !     call hdf5ParticleBundleDump%write('WQ', this%PO(1:this%NPar)%WQ)

        !     call hdf5ParticleBundleDump%close()

        ! end subroutine DiagParticleBundle


        ! subroutine AllInitializationParticleBundle(this, CF, charge, mass, temperature, density, ppg, pbname)
        !     class(ParticleBundle), intent(inout) :: this
        !     class(ControlFlow), intent(in) :: CF
        !     real(8), intent(in) :: charge, mass, temperature, density
        !     integer(4), intent(in) :: ppg
        !     character(*), optional, intent(in) :: pbname
        !     integer(4) :: i, j, k, NPArMax
        !     real(8) :: XFactor, VFactor, Wq
        !     logical ::  Status

        !     ! if (present(pbname)) then
        !     !     call this%IOName%Init(pbname, RESTART_FILE_NAME)
        !     ! else
        !     !     call this%IOName%Init("ParticleBundle", RESTART_FILE_NAME)
        !     ! end if

        !     associate (Nz => CF%DM%GlobalShape(1), &
        !                Nr => CF%DM%GlobalShape(2), &
        !                Lz => CF%DM%LocalShape(1), &
        !                Lr => CF%DM%LocalShape(2), &
        !                dz => CF%DM%SpaceStep(1), &
        !                dr => CF%DM%SpaceStep(2), &
        !                NzL => CF%DM%CornerIndex(BOUNDARY_LEFT, 1), &
        !                NzU => CF%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
        !                NrL => CF%DM%CornerIndex(BOUNDARY_LEFT, 2), &
        !                NrU => CF%DM%CornerIndex(BOUNDARY_RIGHT, 2))

        !         this%dx = dz
        !         this%dt = CF%dt
        !         this%XFactor = this%dx
        !         this%VFactor = this%dx / this%dt

        !         this%Charge = charge
        !         this%Mass   = mass
        !         this%Weight = density / Dble(ppg)

        !         call this%Init(ppg * (Lz-1) * (Lr-1), &
        !                        10, &
        !                        Ceiling(particle_number_max_ratio * ppg * (Lz-1) * (Lr-1)), &
        !                        particle_number_lower_boundary_ratio, &
        !                        particle_number_decrease_ratio, &
        !                        particle_number_increase_ratio)

        !         this%NPar = 0
        !         if (init_z_stop >= init_z_start .and. init_r_stop >= init_r_start) then
        !             do i = init_z_start, init_z_stop
        !                 do j = init_r_start, init_r_stop
        !                     do k = 1, ppg
        !                         this%NPar = this%NPar + 1

        !                         this%LXScaled = .True.
        !                         call this%PO(this%NPar)%PosInit(dble(i-1), dble(i), dble(j-1), dble(j))

        !                         this%LVScaled = .True.
        !                         call this%PO(this%NPar)%VelMaxInit(mass, temperature)
        !                         VFactor = 1.d0 / this%VFactor
        !                         call this%PO(this%NPar)%VelRes(VFactor)
        !                         call this%PO(this%NPar)%AccInpInit()

        !                         Wq = this%PO(this%NPar)%R
        !                         if(Wq <= 1.d0) then
        !                             this%PO(this%NPar)%Wq = 1.d0 / 4.d0 * this%Weight * 2.d0 * PI * this%dx * this%dx * this%dx
        !                         else
        !                             this%PO(this%NPar)%Wq = Wq * this%Weight * 2.d0 * PI * this%dx * this%dx * this%dx
        !                         end if

        !                     end do
        !                 end do
        !             end do
        !         end if

        !         call this%Adjust()

        !     end associate

        ! end subroutine AllInitializationParticleBundle



       

        ! subroutine PositionRescaleParticleBundle(this)
        !     class(ParticleBundle),intent(inout) :: this
        !     integer(4) :: i
        !     real(8) :: XFactor

        !     if(this%LXScaled) then
        !         XFactor=this%XFactor
        !         this%LXScaled=.False.
        !     else
        !         XFactor=1.d0/this%XFactor
        !         this%LXScaled=.True.
        !     end If
            
        !     do i=1,this%NPar
        !         call this%PO(i)%PosRes(XFactor)
        !     end do

        ! end subroutine PositionRescaleParticleBundle


        ! subroutine VelocityRescaleParticleBundle(this)
        !     class(ParticleBundle),intent(inout) :: this
        !     integer(4) :: i
        !     real(8) :: VFactor

        !     if(this%LVScaled) then
        !         VFactor=this%VFactor
        !         this%LVScaled=.False.
        !     else
        !         VFactor=1.d0/this%VFactor
        !         this%LVScaled=.True.
        !     end if

        !     do i=1,this%NPar
        !         call this%PO(i)%VelRes(VFactor)
        !     end do

        ! end subroutine VelocityRescaleParticleBundle


        ! subroutine MoveElectrostaticParticleBundle(this, FG)
        !     class(ParticleBundle),intent(inout) :: this
        !     Type(FieldEM),intent(inout) :: FG

        !     if (pic_type_explicit == pic_type) then
        !         call this%MoveESEx(FG)

        !     else if (pic_type_implicit == pic_type) then
        !         call this%MoveESIm(FG)

        !     end if

        ! end subroutine MoveElectrostaticParticleBundle


        ! subroutine MoveElectrostaticParticleBundleExplicit(this, FG)
        !     Class(ParticleBundle),intent(inout) :: this
        !     Type(FieldEM),intent(inout) :: FG
        !     Integer(4) :: i
        !     Integer(4) :: N1, N2, N1H, N2H
        !     Real(8)    :: Z1, Z2, R1, R2
        !     Real(8)    :: Z1H, Z2H, R1H, R2H
        !     Real(8)    :: E1, E2, E3
        !     Real(8)    :: Ep(3)
        !     Real(8)    :: CosA, SinA
        !     Real(8)    :: EFactor

        !     EFactor = this%Charge / this%Mass * this%dt / (this%dx / this%dt)

        !     Associate(NPar => this%NPar, Ez => FG%Ez, Er => FG%Er)
        !         do i = 1, NPar
        !             N1  = Ceiling(this%PO(i)%Z)
        !             Z1  = Dble(N1) - this%PO(i)%Z
        !             Z2  = 1.d0 - Z1
        !             N1H = Int(this%PO(i)%Z + 0.5d0)
        !             Z1H = Dble(N1H) + 0.5d0 - this%PO(i)%Z
        !             Z2H = 1.d0 - Z1H

        !             N2  = Ceiling(this%PO(i)%R)
        !             R1  = Dble(N2) - this%PO(i)%R
        !             R2  = 1.d0 - R1
        !             N2H = Int(this%PO(i)%R + 0.5d0)
        !             R1H = Dble(N2H) + 0.5d0 - this%PO(i)%R
        !             R2H = 1.d0 - R1H

        !             E1 = Er(N1,N2H)*Z1*R1H + Er(N1+1,N2H)*Z2*R1H + Er(N1,N2H+1)*Z1*R2H + Er(N1+1,N2H+1)*Z2*R2H
        !             E2 = Ez(N1H,N2)*Z1H*R1 + Ez(N1H+1,N2)*Z2H*R1 + Ez(N1H,N2+1)*Z1H*R2 + Ez(N1H+1,N2+1)*Z2H*R2
        !             E3 = 0.d0
                    
        !             E1 = E1 * EFactor
        !             E2 = E2 * EFactor
        !             E3 = E3 * EFactor
                    
        !             If(this%PO(i)%R < MinReal) then
        !                 CosA = this%PO(i)%X / (Abs(this%PO(i)%X) + MinReal)
        !                 SinA = 0.d0
        !             Else
        !                 CosA = this%PO(i)%X / this%PO(i)%R
        !                 SinA = this%PO(i)%Y / this%PO(i)%R
        !             Endif
        
        !             Ep(1) = E1*CosA - E3*SinA
        !             Ep(2) = E1*SinA + E3*CosA
        !             Ep(3) = E2
                    
        !             this%PO(i)%Vx = this%PO(i)%Vx + Ep(1)
        !             this%PO(i)%X  = this%PO(i)%X  + this%PO(i)%Vx
   
        !             this%PO(i)%Vy = this%PO(i)%Vy +Ep(2)
        !             this%PO(i)%Y  = this%PO(i)%Y  + this%PO(i)%Vy
 
        !             this%PO(i)%Vz = this%PO(i)%Vz + Ep(3)
        !             this%PO(i)%Z  = this%PO(i)%Z  + this%PO(i)%Vz

        !             this%PO(i)%R  = DSQRT(this%PO(i)%X*this%PO(i)%X + this%PO(i)%Y*this%PO(i)%Y)
        !         end do
        !     end Associate
            
        ! end subroutine MoveElectrostaticParticleBundleExplicit


        ! subroutine MoveElectrostaticParticleBundleImplicit(this, FG)
        !     class(ParticleBundle),intent(inout) :: this
        !     Type(FieldEM),intent(inout) :: FG
        !     Integer(4) :: i
        !     Integer(4) :: N1, N2, N1H, N2H
        !     Real(8)    :: Z1, Z2, R1, R2
        !     Real(8)    :: Z1H, Z2H, R1H, R2H
        !     Real(8)    :: E1, E2, E3
        !     Real(8)    :: Ep(3)
        !     Real(8)    :: CosA, SinA
        !     Real(8)    :: EFactor

        !     EFactor = this%Charge / this%Mass * this%dt / (this%dx / this%dt)

        !     associate(NPar => this%NPar, Ez => FG%Ez, Er => FG%Er)
        !         do i = 1, NPar
        !             N1  = Ceiling(this%PO(i)%Z)
        !             Z1  = Dble(N1) - this%PO(i)%Z
        !             Z2  = 1.d0 - Z1
        !             N1H = Int(this%PO(i)%Z + 0.5d0)
        !             Z1H = Dble(N1H) + 0.5d0 - this%PO(i)%Z
        !             Z2H = 1.d0 - Z1H

        !             N2  = Ceiling(this%PO(i)%R)
        !             R1  = Dble(N2) - this%PO(i)%R
        !             R2  = 1.d0 - R1
        !             N2H = Int(this%PO(i)%R + 0.5d0)
        !             R1H = Dble(N2H) + 0.5d0 - this%PO(i)%R
        !             R2H = 1.d0 - R1H

        !             E1 = Er(N1,N2H)*Z1*R1H + Er(N1+1,N2H)*Z2*R1H + Er(N1,N2H+1)*Z1*R2H + Er(N1+1,N2H+1)*Z2*R2H
        !             E2 = Ez(N1H,N2)*Z1H*R1 + Ez(N1H+1,N2)*Z2H*R1 + Ez(N1H,N2+1)*Z1H*R2 + Ez(N1H+1,N2+1)*Z2H*R2
        !             E3 = 0.d0
                    
        !             E1 = E1 * EFactor
        !             E2 = E2 * EFactor
        !             E3 = E3 * EFactor
                    
        !             If(this%PO(i)%R < MinReal) then
        !                 CosA = this%PO(i)%X / (Abs(this%PO(i)%X) + MinReal)
        !                 SinA = 0.d0
        !             Else
        !                 CosA = this%PO(i)%X / this%PO(i)%R
        !                 SinA = this%PO(i)%Y / this%PO(i)%R
        !             Endif
        
        !             Ep(1) = E1*CosA - E3*SinA
        !             Ep(2) = E1*SinA + E3*CosA
        !             Ep(3) = E2
                    
        !             this%PO(i)%Vx = this%PO(i)%Vx + Ep(1)
        !             this%PO(i)%X  = this%PO(i)%X  + Ep(1)
        !             this%PO(i)%Ax = 0.5d0 * (this%PO(i)%Ax + Ep(1))

        !             this%PO(i)%Vy = this%PO(i)%Vy +Ep(2)
        !             this%PO(i)%Y  = this%PO(i)%Y  +Ep(2)
        !             this%PO(i)%Ay = 0.5d0 * (this%PO(i)%Ay + Ep(2))
 
        !             this%PO(i)%Vz = this%PO(i)%Vz + Ep(3)
        !             this%PO(i)%Z  = this%PO(i)%Z  + Ep(3)
        !             this%PO(i)%Az = 0.5d0 * (this%PO(i)%Az + Ep(3))
                    
        !             this%PO(i)%Vx = this%PO(i)%Vx + this%PO(i)%Ax
        !             this%PO(i)%X  = this%PO(i)%X  + this%PO(i)%Vx
        !             this%PO(i)%Vy = this%PO(i)%Vy + this%PO(i)%Ay
        !             this%PO(i)%Y  = this%PO(i)%Y  + this%PO(i)%Vy
        !             this%PO(i)%Vz = this%PO(i)%Vz + this%PO(i)%Az
        !             this%PO(i)%Z  = this%PO(i)%Z  + this%PO(i)%Vz

        !             this%PO(i)%R  = DSQRT(this%PO(i)%X*this%PO(i)%X + this%PO(i)%Y*this%PO(i)%Y)
        !         end do
        !     end associate

        ! end subroutine MoveElectrostaticParticleBundleImplicit


        ! subroutine MoveElectroMagneticParticleBundle(this, FG)
        !     class(ParticleBundle),intent(inout) :: this
        !     Type(FieldEM),intent(inout) :: FG
        !     integer(4) :: i,N1,N2,N1H,N2H
        !     real(8) :: Z1,Z2,R1,R2,R
        !     real(8) :: Z1H,Z2H,R1H,R2H
        !     real(8) :: E1,E2,E3,B1,B2,B3
        !     real(8) :: X,Y,CosA,SinA
        !     real(8) :: SVelocity(1:3)=0.d0,Ep(1:3)=0.d0,Bp(1:3)=0.d0
        !     real(8) :: Omega(1:3)=0.d0,OmegaT,TransB(3,3)=0.d0
        !     real(8) :: EFactor,BFactor
        !     real(8) :: RhoFactor,ChiFactor

        !     EFactor=this%Charge/this%Mass*this%dt/(this%dx/this%dt)
        !     BFactor=this%Charge/this%Mass*this%dt

        !     associate(Ez=>FG%Ez, Er=>FG%Er, Bz=>FG%Bz, Br=>FG%Br) !,this%PO(i)%Z=>Z,this%PO(i)%R=>R
        !         do i=1,this%NPar
        !             X=this%PO(i)%X
        !             Y=this%PO(i)%Y
        !             R=this%PO(i)%R

        !             N1=Ceiling(this%PO(i)%Z)
        !             Z1=Dble(N1)-this%PO(i)%Z
        !             Z2=1.d0-Z1
        !             N1H=Int(this%PO(i)%Z+0.5d0)
        !             Z1H=Dble(N1H)+0.5d0-this%PO(i)%Z
        !             Z2H=1.d0-Z1H

        !             N2=Ceiling(this%PO(i)%R)
        !             R1=Dble(N2)-this%PO(i)%R
        !             R2=1.d0-R1
        !             N2H=Int(this%PO(i)%R+0.5d0)
        !             R1H=Dble(N2H)+0.5d0-this%PO(i)%R
        !             R2H=1.d0-R1H

        !             E1=Er(N1,N2H)*Z1*R1H+Er(N1+1,N2H)*Z2*R1H+Er(N1,N2H+1)*Z1*R2H+Er(N1+1,N2H+1)*Z2*R2H
        !             E2=Ez(N1H,N2)*Z1H*R1+Ez(N1H+1,N2)*Z2H*R1+Ez(N1H,N2+1)*Z1H*R2+Ez(N1H+1,N2+1)*Z2H*R2
        !             E3=0.d0
        !             E1=E1*EFactor
        !             E2=E2*EFactor
        !             E3=E3*EFactor

        !             !E3=Et(N1,N2)*Z1*R1+Et(N1+1,N2)*Z2*R1+Et(N1,N2+1)*Z1*R2+Et(N1+1,N2+1)*Z2*R2

        !             B1=Br(N1,N2H)*Z1*R1H+Br(N1+1,N2H)*Z2*R1H+Br(N1,N2H+1)*Z1*R2H+Br(N1+1,N2H+1)*Z2*R2H
        !             B2=Bz(N1H,N2)*Z1H*R1+Bz(N1H+1,N2)*Z2H*R1+Bz(N1H,N2+1)*Z1H*R2+Bz(N1H+1,N2+1)*Z2H*R2
        !             B1=B1*BFactor
        !             B2=B2*BFactor

        !             If(this%PO(i)%R<MinReal) then
        !                 CosA=X/(Abs(X)+MinReal)
        !                 SinA=0.d0
        !             Else
        !                 CosA=X/this%PO(i)%R
        !                 SinA=Y/this%PO(i)%R
        !             Endif

        !             Ep(1)=E1*CosA-E3*SinA
        !             Ep(2)=E1*SinA+E3*CosA
        !             Ep(3)=E2

        !             Bp(1)=B1*CosA
        !             Bp(2)=B1*SinA
        !             Bp(3)=B2

        !             Omega(1)=Bp(1)*0.5d0
        !             Omega(2)=Bp(2)*0.5d0
        !             Omega(3)=Bp(3)*0.5d0
        !             OmegaT=1.d0+Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3)

        !             TransB(1,1)=1.d0+Omega(1)*Omega(1)
        !             TransB(2,1)=Omega(1)*Omega(2)+Omega(3)
        !             TransB(3,1)=Omega(1)*Omega(3)-Omega(2)
        !             TransB(1,2)=Omega(1)*Omega(2)-Omega(3)
        !             TransB(2,2)=1.d0+Omega(2)*Omega(2)
        !             TransB(3,2)=Omega(2)*Omega(3)+Omega(1)
        !             TransB(1,3)=Omega(1)*Omega(3)+Omega(2)
        !             TransB(2,3)=Omega(2)*Omega(3)-Omega(1)
        !             TransB(3,3)=1.d0+Omega(3)*Omega(3)
        !             TransB(1:3,1:3)=TransB(1:3,1:3)/OmegaT

        !             SVelocity(1)=(TransB(1,1)*Ep(1)+TransB(2,1)*Ep(2)+TransB(3,1)*Ep(3))*0.5d0
        !             SVelocity(2)=(TransB(1,2)*Ep(1)+TransB(2,2)*Ep(2)+TransB(3,2)*Ep(3))*0.5d0
        !             SVelocity(3)=(TransB(1,3)*Ep(1)+TransB(2,3)*Ep(2)+TransB(3,3)*Ep(3))*0.5d0

        !             this%PO(i)%Vx=this%PO(i)%Vx+SVelocity(1)
        !             this%PO(i)%Vy=this%PO(i)%Vy+SVelocity(2)
        !             this%PO(i)%Vz=this%PO(i)%Vz+SVelocity(3)

        !             this%PO(i)%X=this%PO(i)%X+SVelocity(1)
        !             this%PO(i)%Y=this%PO(i)%Y+SVelocity(2)
        !             this%PO(i)%Z=this%PO(i)%Z+SVelocity(3)

        !             this%PO(i)%Ax=0.5d0*(this%PO(i)%Ax+Ep(1))
        !             this%PO(i)%Ay=0.5d0*(this%PO(i)%Ay+Ep(2))
        !             this%PO(i)%Az=0.5d0*(this%PO(i)%Az+Ep(3))

        !             SVelocity(1)=this%PO(i)%Vx+0.5d0*this%PO(i)%Ax+ &
        !                             this%PO(i)%Vy*Omega(3)-this%PO(i)%Vz*Omega(2)
        !             SVelocity(2)=this%PO(i)%Vy+0.5d0*this%PO(i)%Ay+ &
        !                             this%PO(i)%Vz*Omega(1)-this%PO(i)%Vx*Omega(3)
        !             SVelocity(3)=this%PO(i)%Vz+0.5d0*this%PO(i)%Az+ &
        !                             this%PO(i)%Vx*Omega(2)-this%PO(i)%Vy*Omega(1)

        !             this%PO(i)%Vx=TransB(1,1)*SVelocity(1)+TransB(2,1)*SVelocity(2)+TransB(3,1)*SVelocity(3)
        !             this%PO(i)%Vy=TransB(1,2)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(3,2)*SVelocity(3)
        !             this%PO(i)%Vz=TransB(1,3)*SVelocity(1)+TransB(2,3)*SVelocity(2)+TransB(3,3)*SVelocity(3)

        !             this%PO(i)%X=this%PO(i)%X+this%PO(i)%Vx
        !             this%PO(i)%Y=this%PO(i)%Y+this%PO(i)%Vy
        !             this%PO(i)%Z=this%PO(i)%Z+this%PO(i)%Vz

        !             this%PO(i)%R=Dsqrt(this%PO(i)%X*this%PO(i)%X+this%PO(i)%Y*this%PO(i)%Y)
        !             !=R
        !         End do
        !     End associate

        ! end subroutine MoveElectroMagneticParticleBundle


        ! subroutine VelocityLimitParticleBundle(this, limitfactor, scalefactor)
        !     class(ParticleBundle),intent(inout) :: this
        !     real(8), intent(in) :: limitfactor, scalefactor
        !     integer(4) :: i
        !     real(8) :: VelocityMax, limitfactorNorm, scalefactorNorm

        !     limitfactorNorm = limitfactor
        !     if (limitfactorNorm <= 0.d0 .or. limitfactorNorm > 1.d0) limitfactorNorm = 0.75d0

        !     scalefactorNorm = scalefactor
        !     if (scalefactorNorm <= 0.d0 .or. scalefactorNorm > 1.d0) scalefactorNorm = 0.25d0

        !     VelocityMax = limitfactorNorm * particle_limit_max_velocity
        !     do i = 1, this%NPar
        !         if (sqrt(this%PO(i)%Vx**2 + this%PO(i)%Vy**2 + this%PO(i)%Vz**2) > VelocityMax) then
        !             this%PO(i)%Vx = this%PO(i)%Vx * scalefactorNorm
        !             this%PO(i)%Vy = this%PO(i)%Vy * scalefactorNorm
        !             this%PO(i)%Vz = this%PO(i)%Vz * scalefactorNorm
        !         end if
        !     end do

        ! end subroutine VelocityLimitParticleBundle


        ! subroutine MoveParticleOneElectrostaticImplicit(this, PB, FG)
        !     class(ParticleOne), intent(inout) :: this
        !     type(ParticleBundle), intent(in) :: PB
        !     type(FieldEM), intent(in) :: FG
        !     integer(4) :: N1, N2, N1H, N2H
        !     real(8)    :: Z1, Z2, R1, R2
        !     real(8)    :: Z1H, Z2H, R1H, R2H
        !     real(8)    :: E1, E2, E3
        !     real(8)    :: Ep(3)
        !     real(8)    :: CosA, SinA
        !     real(8)    :: EFactor

        !     EFactor = PB%Charge / PB%Mass * PB%dt / (PB%dx / PB%dt)

        !     associate(Ez => FG%Ez, Er => FG%Er)
        !         N1  = Ceiling(this%Z)
        !         Z1  = Dble(N1) - this%Z
        !         Z2  = 1.d0 - Z1
        !         N1H = Int(this%Z + 0.5d0)
        !         Z1H = Dble(N1H) + 0.5d0 - this%Z
        !         Z2H = 1.d0 - Z1H

        !         N2  = Ceiling(this%R)
        !         R1  = Dble(N2) - this%R
        !         R2  = 1.d0 - R1
        !         N2H = Int(this%R + 0.5d0)
        !         R1H = Dble(N2H) + 0.5d0 - this%R
        !         R2H = 1.d0 - R1H

        !         E1 = Er(N1,N2H)*Z1*R1H + Er(N1+1,N2H)*Z2*R1H + Er(N1,N2H+1)*Z1*R2H + Er(N1+1,N2H+1)*Z2*R2H
        !         E2 = Ez(N1H,N2)*Z1H*R1 + Ez(N1H+1,N2)*Z2H*R1 + Ez(N1H,N2+1)*Z1H*R2 + Ez(N1H+1,N2+1)*Z2H*R2
        !         E3 = 0.d0
                
        !         E1 = E1 * EFactor
        !         E2 = E2 * EFactor
        !         E3 = E3 * EFactor
                
        !         If(this%R < MinReal) then
        !             CosA = this%X / (Abs(this%X) + MinReal)
        !             SinA = 0.d0
        !         Else
        !             CosA = this%X / this%R
        !             SinA = this%Y / this%R
        !         Endif
    
        !         Ep(1) = E1*CosA - E3*SinA
        !         Ep(2) = E1*SinA + E3*CosA
        !         Ep(3) = E2
                
        !         this%Vx = this%Vx + Ep(1)
        !         this%X  = this%X  + Ep(1)
        !         this%Ax = 0.5d0 * (this%Ax + Ep(1))

        !         this%Vy = this%Vy +Ep(2)
        !         this%Y  = this%Y  +Ep(2)
        !         this%Ay = 0.5d0 * (this%Ay + Ep(2))

        !         this%Vz = this%Vz + Ep(3)
        !         this%Z  = this%Z  + Ep(3)
        !         this%Az = 0.5d0 * (this%Az + Ep(3))
                
        !         this%Vx = this%Vx + this%Ax
        !         this%X  = this%X  + this%Vx
        !         this%Vy = this%Vy + this%Ay
        !         this%Y  = this%Y  + this%Vy
        !         this%Vz = this%Vz + this%Az
        !         this%Z  = this%Z  + this%Vz

        !         this%R  = DSQRT(this%X*this%X + this%Y*this%Y)
        !     end associate

        ! end subroutine MoveParticleOneElectrostaticImplicit

End Module ModuleParticleBundle
