Module ModuleSpecyOne
    use mpi
    use Constants
    ! use ModuleControlFlow

    implicit none

    integer(4), parameter :: NgMax = NG_MAX         ! the max gas number
    integer(4), parameter :: NsMax = NS_MAX         ! the max specy number per gas

    integer(4), parameter :: MC_MODEL_NULL = 1      ! the density and temperature will not change over time
    integer(4), parameter :: MC_MODEL_NEAUTRAL = 2  ! depleting is considered, but the energy balance or heating is not
    integer(4), parameter :: MC_MODEL_FLUID = 3     ! plasma+continous Fluid Gas
    integer(4), parameter :: MC_MODEL_DSMC = 4      ! plasma + DSMC Gas

    integer(4), parameter :: CS_MODEL_LXCAT = 1     ! crosssection data source
    integer(4), parameter :: CS_MODEL_PEGASUS = 2

    Type SpecyOnePhysical                           ! define a physical particle
        character(99)   :: Name  = "A+"
        integer(4)      :: NAtom = 1, Charge = 0
        real(8)         :: Mass  = 0.d0, Radius = 1.d0, polarizability = 0.d0
    end Type SpecyOnePhysical

    Type SpecyOnePegasus                            ! Define a pegasus particle.
        integer(4)      :: Index=0                  !粒子编号
        real(8)         :: Mass=0.d0                !粒子质量(AMU单位)
        integer(4)      :: Charge=0                 !粒子单位电荷
        real(8)         :: Radius=1.d0              !粒子半径（pm单位）
        real(8)         :: Polarizability=0.d0      !极化率(化学反应性碰撞，单位为Bohr半径的倍数)
        character(99)   :: Name="Ar+"               !粒子符号
        
    contains

        procedure :: Load => LoadSpecyOnePegasus
        procedure :: Convert2Specy => ConvertSpecyOnePegasus2SpecyOnePhysical
        procedure :: Convert2Gas   => ConvertSpecyOnePegasus2GasPhysical

        Generic :: Convert => Convert2Specy, Convert2Gas

    end Type SpecyOnePegasus 

    Type GasPhysical                                ! define a physical gas
        character(99)   :: Name = "Cu"
        integer(4)      :: MCModel = MC_MODEL_NULL, NAtom, Ns
        integer(4)      :: CSModel = CS_MODEL_PEGASUS
        real(8)         :: Mass, Radius = 0.d0, BetaMax = 0.d0, polarizability = 0.d0

        Type(SpecyOnePhysical) :: SP(NsMax)
        Type(SpecyOnePegasus)  :: SOP(0:NsMax)
    end Type GasPhysical

    Type SpecyOne                                   ! define a numerical particle.
        character(99)   :: Name
        integer(4)      :: SpecyIndex = -1, GasIndex = -1
        real(8)         :: Charge, Mass, Radius, NAtom
        real(8)         :: InitDensity, Density, InitTemperature, Temperature
    end Type SpecyOne

    Type GasOne                                     ! define a numerical particle.
        character(99)   :: Name
        integer(4)      :: MCModel = MC_MODEL_NULL, Ns, IndexStart
        integer(4)      :: CSModel = CS_MODEL_PEGASUS
        real(8)         :: Mass, Radius, NAtom, BetaMax, polarizability
        real(8)         :: InitDensity, Density, InitTemperature, Temperature
        
        Type(SpecyOnePegasus) :: SOP(0:NsMax)
    end Type GasOne

    Type(SpecyOnePhysical), parameter :: ElectronPhysical = &
         SpecyOnePhysical('Electron', 1, -1, ElectronMass/AtomicMass, 0.d0)

    Type(SpecyOne), Allocatable,save :: SpecyGlobal(:)
    Type(GasOne), Allocatable,save   :: GasGlobal(:)


    contains
    
        ! subroutine GasInit(Ng,Ns ,pressure, gas_temperature, ele_temperature, gas_name, gas_ratio)
        !     ! class(ControlFlow), intent(inout) :: CF
        !     real(8), intent(in) ::Ng,Ns , pressure, gas_temperature, ele_temperature
        !     character(len=99), dimension(CF%Ng) :: gas_name
        !     real(8), intent(in) :: gas_ratio(CF%Ng)
        !     type(GasPhysical) :: GP
        !     type(GasPhysical), Allocatable :: TempGP(:)
        !     integer(4) :: Nss = 0
        !     integer(4) :: i, j, k, Index = 0
        !     logical    :: alive
        !     character(len=99) :: tmpFilePath
        !     NAMELIST /GasPhysicsfile/ GP
        !     integer(4) :: size, rank, ierr

        !     call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        !     call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        !     call DestroyGas()

        !     ! load gas
        !     Allocate(GasGlobal(Ng))
        !     Allocate(TempGP(Ng))

        !     do i = 1, Ng
        !         tmpFilePath = "./input/gas/lxcat/"//Trim(gas_name(i))//"/"//Trim(gas_name(i))//".txt"
        !         inquire(file=tmpFilePath, exist=alive)
        !         if(alive) then
        !             open(20, FILE=tmpFilePath)
        !                 read(20, NML=GasPhysicsfile)
        !             close(20)
        !         else
        !             if (0 == rank) write(*, '(a, a, a)') 'The gas file ', tmpFilePath, ' is not find.'
        !             stop 1

        !         end if

        !         TempGP(i) = GP
        !         Ns = Ns + GP%Ns
        !     end do

        !     ! load specy
        !     ! CF%Ns=Ns
        !     Allocate(SpecyGlobal(0:Ns))

        !     call SpecyOnePhysics2SpecyOne(SpecyGlobal(0), ElectronPhysical)
        !     SpecyGlobal(0)%SpecyIndex       = 0
        !     SpecyGlobal(0)%GasIndex         = 0
        !     SpecyGlobal(0)%InitDensity      = 1.d0
        !     SpecyGlobal(0)%Density          = SpecyGlobal(0)%InitDensity
        !     SpecyGlobal(0)%InitTemperature  = ele_temperature
        !     SpecyGlobal(0)%Temperature      = SpecyGlobal(0)%InitTemperature

        !     do i = 1, Ng
        !         call GasPhysics2Gas(GasGlobal(i), TempGP(i))

        !         GasGlobal(i)%InitTemperature = gas_temperature
        !         GasGlobal(i)%Temperature     = GasGlobal(i)%InitTemperature
        !         GasGlobal(i)%InitDensity     = (pressure * mTorrtoPa) / (kB * GasGlobal(i)%InitTemperature) * gas_ratio(i)
        !         GasGlobal(i)%Density         = GasGlobal(i)%InitDensity
        !         GasGlobal(i)%IndexStart      = Index

        !         do j = 1, GasGlobal(i)%Ns
        !             Index = Index + 1
        !             SpecyGlobal(Index)%SpecyIndex       = Index
        !             SpecyGlobal(Index)%GasIndex         = i
        !             SpecyGlobal(Index)%InitDensity      = gas_ratio(i)
        !             SpecyGlobal(Index)%Density          = SpecyGlobal(Index)%InitDensity
        !             SpecyGlobal(Index)%InitTemperature  = gas_temperature
        !             SpecyGlobal(Index)%Temperature      = SpecyGlobal(Index)%InitTemperature
        !             call SpecyOnePhysics2SpecyOne(SpecyGlobal(Index), TempGP(i)%SP(j))
        !         end do
        !     end do

        !     if (Allocated(TempGP)) Deallocate(TempGP)
        !     if (0 == rank) write(*, *) "Init lxcat gas."
            
        ! end subroutine GasInit

        Subroutine  GasInitPegasus(Ng,Ns, pressure, gas_temperature, ele_temperature, gas_name, gas_ratio, model_type)
            ! class(ControlFlow), intent(inout) :: CF
            real(8), intent(in) :: pressure, gas_temperature, ele_temperature
            character(len=99), dimension(Ng) :: gas_name
            real(8), intent(in) :: gas_ratio(Ng)
            integer(4) :: model_type,Ng,Ns
            type(GasPhysical) :: GP
            type(GasPhysical), Allocatable :: TempGP(:)
            ! integer(4) :: Ns = 0
            integer(4) :: i, j, k, Index = 0
            logical    :: alive, INited
            character(len=99) :: tmpFilePath
            integer(4) :: size, rank, ierr

            call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

            call DestroyGas()

            Allocate(GasGlobal(CF%Ng))
            Allocate(TempGP(CF%Ng))

            do i=1,Ng
                if (model_type == CS_MODEL_PEGASUS) then
                    Call LoadGasPhysicalPegasus(TempGP(i),gas_name(i),INited)
                    if (0 == rank) write(*, *) "Init pegasus gas."

                else if (model_type == CS_MODEL_LXCAT) then
                    Call LoadGasPhysicalLxCat(TempGP(i),gas_name(i),INited)
                    if (0 == rank) write(*, *) "Init lxcat gas."

                else
                    if (0 == rank) write(*, *) "The cs model is invaild."

                end If

                Ns=Ns+TempGP(i)%Ns
            End Do

            ! CF%Ns=Ns
            Allocate(SpecyGlobal(0:Ns))

            Call SpecyOnePhysics2SpecyOne(SpecyGlobal(0),ElectronPhysical)
            SpecyGlobal(0)%SpecyIndex=0
            SpecyGlobal(0)%GasIndex=0

            SpecyGlobal(0)%InitDensity=1.d0
            SpecyGlobal(0)%Density=SpecyGlobal(0)%InitDensity
            SpecyGlobal(0)%InitTemperature = ele_temperature
            SpecyGlobal(0)%Temperature=SpecyGlobal(0)%InitTemperature

            Do i=1,CF%Ng
                Select Case (TempGP(i)%CSModel)
                    Case(1)
                        Call GasPhysics2Gas(GasGlobal(i),TempGP(i))
                    Case(2)
                        Call GasPhysics2GasPegasus(GasGlobal(i),TempGP(i))
                ENd Select

                GasGlobal(i)%InitTemperature=gas_temperature
                GasGlobal(i)%Temperature=GasGlobal(i)%InitTemperature
                GasGlobal(i)%InitDensity=(pressure*mTorrtoPa)/(kB*GasGlobal(i)%InitTemperature)*gas_ratio(i)
                GasGlobal(i)%Density=GasGlobal(i)%InitDensity
                GasGlobal(i)%IndexStart=Index

                Do j=1,GasGlobal(i)%Ns
                    Index=Index+1
                    SpecyGlobal(Index)%SpecyIndex=Index
                    SpecyGlobal(Index)%GasIndex=i
                    SpecyGlobal(Index)%InitDensity=gas_ratio(i)
                    SpecyGlobal(Index)%Density = SpecyGlobal(Index)%InitDensity
                    SpecyGlobal(Index)%InitTemperature = gas_temperature
                    SpecyGlobal(Index)%Temperature=SpecyGlobal(Index)%InitTemperature
                    Call SpecyOnePhysics2SpecyOne(SpecyGlobal(Index),TempGP(i)%SP(j))
                End Do
            ENd do

            if (Allocated(TempGP)) Deallocate(TempGP)

        End  subroutine  GasInitPegasus

        subroutine DestroyGas()

            if (Allocated(GasGlobal))   Deallocate(GasGlobal)
            if (Allocated(SpecyGlobal)) Deallocate(SpecyGlobal)

            return
        end subroutine DestroyGas


        subroutine SpecyOnePhysics2SpecyOne(SO, SP)
            Type(SpecyOne), intent(inout) :: SO
            Type(SpecyOnePhysical), intent(in) :: SP

            SO%Name   = SP%Name
            SO%Charge = Dble(SP%Charge) * ElectronCharge
            SO%Mass   = SP%Mass * AtomicMass
            SO%NAtom  = Dble(SP%NAtom)
            SO%Radius = SP%Radius

            return
        end subroutine SpecyOnePhysics2SpecyOne


        subroutine GasPhysics2Gas(GO, GP)
            Type(GasOne), intent(inout) :: GO
            Type(GasPhysical), intent(in) :: GP

            GO%Name    = GP%Name
            GO%MCModel = GP%MCModel
            GO%CSModel = GP%CSModel
            GO%Ns      = GP%Ns
            GO%Mass    = GP%Mass * AtomicMass
            GO%NAtom   = Dble(GP%NAtom)
            GO%Radius  = GP%Radius
            GO%BetaMax = GP%BetaMax
            GO%polarizability = 14.d0 * a0 * a0 * a0

            return
        end subroutine GasPhysics2Gas


        subroutine GasPhysics2GasPegasus(GO, GP)
            Type(GasOne), intent(inout) :: GO
            Type(GasPhysical), intent(in) :: GP
            integer(4) :: i

            GO%Name    = GP%Name
            GO%CSModel = GP%CSModel
            GO%MCModel = GP%MCModel
            GO%Ns      = GP%Ns
            GO%Mass    = GP%Mass * AtomicMass
            GO%NAtom   = Dble(GP%NAtom)
            GO%Radius  = GP%Radius
            GO%BetaMax = GP%BetaMax
            GO%SOP     = GP%SOP
            GO%polarizability = 14.d0 * a0 * a0 * a0

            return
        end subroutine  GasPhysics2GasPegasus


        ! subroutine ShowGas(Ng,Ns)
        !     ! Type(ControlFlow), intent(inout) :: CF
        !     integer(4) :: i, j,Ng,Ns
        !     integer(4) :: size, rank, ierr

        !     call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
        !     call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        !     if (0 == rank) write(*, '(a)') "The Gas Info: "

        !     if (0 == rank) write(*, '(4x, a8, i4, i4, *(es12.4))') SpecyGlobal(0)%Name, &
        !                                         SpecyGlobal(0)%SpecyIndex, SpecyGlobal(0)%GasIndex, &
        !                                         SpecyGlobal(0)%Charge, SpecyGlobal(0)%Mass, &
        !                                         SpecyGlobal(0)%Radius, SpecyGlobal(0)%NAtom, &
        !                                         SpecyGlobal(0)%InitDensity, SpecyGlobal(0)%Density, &
        !                                         SpecyGlobal(0)%InitTemperature, SpecyGlobal(0)%Temperature
        !     if (0 == rank) write(*, *) ""
            
        !     do i = 1, Ng
        !         if (0 == rank) write(*, '(a8, i4, i4, i4, *(es12.4))') GasGlobal(i)%Name, &
        !                                         GasGlobal(i)%MCModel, GasGlobal(i)%Ns, &
        !                                         GasGlobal(i)%IndexStart, GasGlobal(i)%Mass, GasGlobal(i)%Radius, &
        !                                         GasGlobal(i)%NAtom, GasGlobal(i)%BetaMax, &
        !                                         GasGlobal(i)%InitDensity, GasGlobal(i)%Density, &
        !                                         GasGlobal(i)%InitTemperature, GasGlobal(i)%Temperature
                
        !         do j = GasGlobal(i)%IndexStart+1, GasGlobal(i)%IndexStart + GasGlobal(i)%Ns
        !             if (0 == rank) write(*, '(4x, a6, i4, i4, *(es12.4))') SpecyGlobal(j)%Name, &
        !                                         SpecyGlobal(j)%SpecyIndex, SpecyGlobal(j)%GasIndex, &
        !                                         SpecyGlobal(j)%Charge, SpecyGlobal(j)%Mass, &
        !                                         SpecyGlobal(j)%Radius, SpecyGlobal(j)%NAtom, &
        !                                         SpecyGlobal(j)%InitDensity, SpecyGlobal(j)%Density, &
        !                                         SpecyGlobal(j)%InitTemperature, SpecyGlobal(j)%Temperature
        !         end do

        !         if (0 == rank) write(*, *) ""
        !     end do

        !     return
        ! end subroutine ShowGas


        subroutine ConvertSpecyOnePegasus2SpecyOnePhysical(SOP, SP)
            Class(SpecyOnePegasus), intent(inout) :: SOP
            Type(SpecyOnePhysical), Intent(inout) :: SP

            SP%Mass   = SOP%Mass
            SP%Charge = SOP%Charge
            SP%Radius = SOP%Radius*1.d-10
            SP%Polarizability = SOP%Polarizability
            SP%Name = SOP%Name

            return
        end subroutine ConvertSpecyOnePegasus2SpecyOnePhysical


        subroutine ConvertSpecyOnePegasus2GasPhysical(SOP, GP)
            Class(SpecyOnePegasus), intent(inout) :: SOP
            Type(GasPhysical) :: GP

            GP%Name   = SOP%Name
            GP%Mass   = SOP%Mass
            GP%Radius = SOP%Radius*1.d-10

            return
        end subroutine ConvertSpecyOnePegasus2GasPhysical


        subroutine LoadSpecyOnePegasus(SOP, IONumber)
            Class(SpecyOnePegasus), intent(inout) :: SOP
            Integer(4) :: IONumber

            read(IONumber,*) SOP%Index, SOP%Mass, SOP%Charge, SOP%Radius, SOP%Polarizability, SOP%Name
            
            return
        end subroutine LoadSpecyOnePegasus

end Module ModuleSpecyOne
