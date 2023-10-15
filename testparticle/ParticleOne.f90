Module ModuleParticleOne
    use Constants

    implicit none

    Type :: ParticleOne

        real(8) ::  X, Y, Z, R, Vx, Vy, Vz, Ax, Ay, Az, WQ
        integer(4):: tempdelflag=0

    contains

        procedure :: PosInit    => PositionRandomInitializationParticleOne
        procedure :: VelInpInit => VelocityInputInitializationParticleOne
        procedure :: VelMaxInit => VelocityMaxwellianInitializationParticleOne
        procedure :: VelRanInit => VelocityRandomInitializationParticleOne
        procedure :: AccInpInit => AccelerationInputInitializationParticleOne
        procedure :: WqInpInit  => WeightInputInitializationParticleOne

        procedure :: PosRes => PositionRescaleParticleOne
        procedure :: VelRes => VelocityRescaleParticleOne
        procedure :: AccRes => AccelerationRescaleParticleOne
        procedure :: Energy => EnergyParticleOne

        procedure :: Copy => CopyParticleOne
        procedure :: Swap => SwapParticleOne

    end Type ParticleOne

    contains

        subroutine PositionRandomInitializationParticleOne(PO,ZL,ZU,RL,RU)
            Class(ParticleOne), intent(inout) :: PO
            real(8), intent(in) :: ZL,ZU,RL,RU
            real(8) :: Theta

            call RANDOM_NUMBER(R)
            PO%R=RL+(RU-RL)*R
            call RANDOM_NUMBER(R)
            PO%Z=ZL+(ZU-ZL)*R

            call RANDOM_NUMBER(R)
            Theta=2.d0*PI*R
            PO%X=PO%R*DSin(Theta)
            PO%Y=PO%R*DCos(Theta)

            return
        end subroutine

        subroutine VelocityInputInitializationParticleOne(PO,Vx,Vy,Vz)
            Class(ParticleOne), intent(inout) :: PO
            real(8), intent(in), optional :: Vx,Vy,Vz

            if (present(Vx)) then
                PO%Vx=Vx
            else
                PO%Vx=0.d0
            end if

            if (present(Vy)) then
                PO%Vy=Vy
            else
                PO%Vy=0.d0
            end if

            if (present(Vz)) then
                PO%Vz=Vz
            else
                PO%Vz=0.d0
            end if

            return
        end subroutine VelocityInputInitializationParticleOne

        subroutine VelocityMaxwellianInitializationParticleOne(PO,Mass,Temperature)
            Class(ParticleOne), intent(inout) :: PO
            real(8), intent(in) :: Mass,Temperature
            real(8) :: V,Beta,FuncA,FuncB
            real(8) :: Theta,CosTheta,SinTheta,Phi

            Beta=1.d0/(kB*Temperature)
            FuncA=1.d0
            FuncB=0.d0
            do while(FuncA>FuncB)
                call RANDOM_NUMBER(R)
                FuncA=R*R
                call RANDOM_NUMBER(R)
                FuncB=-exp*R*Dlog(R)
            end do
            V=DSqrt(-3.d0*Dlog(R)/Beta/Mass)
            call VelocityRandomInitializationParticleOne(PO,V)

            return
        end subroutine VelocityMaxwellianInitializationParticleOne

        subroutine VelocityRandomInitializationParticleOne(PO,V)
            Class(ParticleOne), intent(inout) :: PO
            real(8), intent(in) ::  V
            real(8) :: Phi,CosPhi,SinPhi
            real(8) :: Theta,CosTheta,FcosTheta,SinTheta

            call RANDOM_NUMBER(R)
            Associate(Vx=>PO%Vx, Vy=>PO%Vy, Vz=>PO%Vz)
                CosTheta=1.d0-2.d0*R
                SinTheta=Dsqrt(1.d0-cosTheta*cosTheta)
                call RANDOM_NUMBER(R)
                Phi=2.d0*PI*R
                PO%Vx=V*CosTheta
                PO%Vy=V*SinTheta*DCos(Phi)
                PO%Vz=V*SinTheta*Dsin(Phi)
            end Associate

            return
        end subroutine VelocityRandomInitializationParticleOne

        subroutine AccelerationInputInitializationParticleOne(PO,Ax,Ay,Az)
            Class(ParticleOne), intent(inout) :: PO
            real(8), intent(in), optional :: Ax,Ay,Az

            if(present(Ax)) then
                PO%Ax=Ax
            else
                PO%Ax=0.d0
            end if

            if(present(Ay)) then
                PO%Ay=Ay
            else
                PO%Ay=0.d0
            end if
            
            if(present(Az)) then
                PO%Az=Az
            else
                PO%Az=0.d0
            end if

            return
        end subroutine AccelerationInputInitializationParticleOne

        subroutine WeightInputInitializationParticleOne(PO,Wq)
            Class(ParticleOne),intent(inout) :: PO
            real(8),intent(in) :: Wq

            PO%Wq=Wq

            return
        end subroutine WeightInputInitializationParticleOne

        subroutine PositionRescaleParticleOne(PO,XFactor)
            Class(ParticleOne), intent(inout) :: PO
            real(8), intent(in) :: XFactor

            Associate(X=>PO%X, Y=>PO%Y, Z=>PO%Z, R=>PO%R)
                X=X*XFactor
                Y=Y*XFactor
                Z=Z*XFactor
                R=R*XFactor
            end Associate

            return
        end subroutine PositionRescaleParticleOne

        subroutine VelocityRescaleParticleOne(PO,VFactor)
            Class(ParticleOne), intent(inout) :: PO
            real(8), intent(in) :: VFactor

            Associate(Vx=>PO%Vx,Vy=>PO%Vy,Vz=>PO%Vz)
                Vx=Vx*VFactor
                Vy=Vy*VFactor
                Vz=Vz*VFactor
            end Associate

            return
        end subroutine VelocityRescaleParticleOne

        subroutine AccelerationRescaleParticleOne(PO,AFactor)
            Class(ParticleOne), intent(inout) :: PO
            real(8), intent(in) :: AFactor
            Associate(Ax=>PO%Ax, Ay=>PO%Ay, Az=>PO%Az)
                Ax=Ax*AFactor
                Ay=Ay*AFactor
                Az=Az*AFactor
            end Associate

            return
        end subroutine AccelerationRescaleParticleOne

        subroutine CopyParticleOne(POD,POC)
            Class(ParticleOne), intent(inout) :: POD
            Class(ParticleOne), intent(in) :: POC

            select Type(POD)
                Type is(ParticleOne)
                    select Type(POC)
                        Type is(ParticleOne)
                            POD=POC
                    end select
            end select

            return
        end subroutine CopyParticleOne

        subroutine SwapParticleOne(POD,POC)
            Class(ParticleOne), intent(inout) :: POD
            Class(ParticleOne), intent(inout) :: POC
            Type(ParticleOne) :: POT

            select Type(POD)
                Type is(ParticleOne)
                    select Type(POC)
                        Type is(ParticleOne)
                            POT=POC
                            POC=POD
                            POD=POT
                    end select
            end select

            return
        end subroutine SwapParticleOne

        Elemental function EnergyParticleOne(PO,Mass,VFactor)
            Class(ParticleOne), intent(in) :: PO
            real(8), intent(in) :: Mass
            real(8), intent(in), optional :: VFactor
            real(8) :: EnergyParticleOne

            EnergyParticleOne=0.5d0*Mass*(PO%Vx*PO%Vx+PO%Vy*PO%Vy+PO%Vz*PO%Vz)
            if(present(VFactor)) EnergyParticleOne=EnergyParticleOne*VFactor*VFactor

            return
        end function EnergyParticleOne

end Module ModuleParticleOne
