module MCCInterface
use iso_c_binding
 Use ModuleTypeMCC
! use ModuleMCCInitialization
! Use ModuleParticleOne
 Use ModuleMCCElectron
 Use ModuleMCCIon
! Use ModuleMCCSigma
! Use ModuleParticleBundle
! use ModuleParticleBundleIndex
 use ModuleSpecyOne
 use ModuleReactionOnePegasus
 use ModuleMCCInitialization
 use ModuleParticleOne
implicit none
Type(MCCBundle),Allocatable :: MCCBundleGlobal(:,:)
real(8) :: gas_pressure, gas_temperature, ele_temperature,VFactor
character(len=99), dimension(1) :: gas_name
real(8) iongas_ratio,elegas_ratio
integer(4) model_type,Ng,Ns,collision_section_type
real(8) ,dimension(1)::gas_ratio
type(ParticleOne) aParticle

contains

subroutine  MCCBundleInit(coll_ratio,dx,dt)bind(C,name="MCCBundleInit")
   real(c_double) coll_ratio(2),dx(3),dt
 
gas_pressure=25
gas_temperature=300
ele_temperature=30000
gas_name(1)="Ar"
collision_section_type=2
model_type=2
Ns=0
Ng=1
gas_ratio(1)=1
! dt=1d-10
VFactor=dx(1)/dt



 
 call GasInitPegasus(Ng,Ns,gas_pressure, gas_temperature, ele_temperature, gas_name, gas_ratio, collision_section_type)
 Allocate(MCCBundleGlobal(0:Ns,Ng))
!  Ns=1
! Ng=1
 call MCCBundleElelctronInit(Ns,Ng,dt,SpecyGlobal(0),GasGlobal(1:Ng),MCCBundleGlobal(0,1:Ng))
! write(*,*)"Collisionratio 1 ",MCCBundleGlobal(0,1)%CollisionRatio
! Ns=1
! Ng=1
   call MCCBundleIonInit(Ns,Ng,dt,SpecyGlobal(1:Ns),GasGlobal(1:Ng),MCCBundleGlobal(1:Ns,1:Ng)) 
!    write(*,*)"Collisionratio 2 ",MCCBundleGlobal(1,1)%CollisionRatio
coll_ratio(1)=MCCBundleGlobal(0,1)%CollisionRatio
coll_ratio(2)=MCCBundleGlobal(1,1)%CollisionRatio
! write(*,*)"Ng,Ns",Ng,Ns,"  "
 end subroutine MCCBundleInit

subroutine MCC(x1,v1,x2,v2,x3,v3,ispecies)bind(C ,name="MCC")
    Implicit none
     real(c_double) :: x1(3),v1(3),x2(3),v2(3),x3(3),v3(3)
     Real(8) :: CollisionRatio,R
     integer(8):: flag
     integer(c_int) Index,ispecies(3)
     type(MCCParticleOne) :: One
     type(ParticleOne),target::  particle
     Call Random_Number(R)
   !   write(*,*)"output",x1(1),x1(2),x1(3),v1(1),v1(2),v1(3);
   !   x1(1)=2
   !   One%POT%X=x1(1)
     particle%X=x1(1)
     particle%Y=x1(2)
     particle%Z=x1(3)
     particle%Vx=v1(1)
     particle%Vy=v1(2)
     particle%Vz=v1(3)
      One%POI=>particle
       VFactor=1

     If (ispecies(1)==0) Then
      call  One%POI%VelRes(VFactor)
      call One%Updater(SpecyGlobal(0),GasGlobal(1))
      !   One%POT%X=0
        Call SelectProbility(One,MCCBundleGlobal(0,1))
     Index=One%ReactionIndex
     Call  SelectCollisionElectron(One,SpecyGlobal(0),GasGlobal(1),MCCBundleGlobal(0,1)%Reaction(Index))
   !   if(Index>1) then
   !   write(*,*)"Newnumber ispecies(1)",ispecies(1),Index,MCCBundleGlobal(0,1)%Reaction(Index)%ReactionType 
   !   end if
   else
      call  One%POI%VelRes(VFactor)
      call One%Updater(SpecyGlobal(1),GasGlobal(1))
      !   One%POT%X=0
        Call SelectProbility(One,MCCBundleGlobal(1,1))
        Index=One%ReactionIndex
        Call  SelectCollisionIon(One,SpecyGlobal(1),GasGlobal(1),MCCBundleGlobal(1,1)%Reaction(Index))
        call One%POI%VelRes(1.d0/VFactor)
     
 Index=One%ReactionIndex

     end if 



     if(One%NPONew>0)Then
      ! write(*,*)"newparticle",One%PON(1)%X,One%PON(2)%X
      if(One%NPONew>=1) then
         ispecies(1)=One%PON(1)%Index+1
         v2(1)=One%PON(1)%Vx
         v2(2)=One%PON(1)%Vy
         v2(2)=One%PON(1)%Vz
     end if
     if(One%NPONew==2) then
         ispecies(2)=One%PON(2)%Index+1
         v3(1)=One%PON(2)%Vx
         v3(2)=One%PON(2)%Vy
         v3(2)=One%PON(2)%Vz
     end if
     end if

     x1(1)=particle%X
     x1(2)=particle%Y
     x1(3)=particle%Z
     v1(1)=particle%Vx
     v1(2)=particle%Vy
     v1(3)=particle%Vz
   !   Do m=1,MCCParticleBundle(k)%NPONew
      ! Call PBIGlobal%AddOne(MCCParticleBundle(k)%PON(m))
   !   if(One%NPONew>0) then
     
   !   end if
   !   Vx=Vx*VFactor
   !   Vy=Vy*VFactor
   !   Vz=Vz*VFactor
    !  write(*,*)One%POT%X
   !   PB%VFactor = PB%dx / PB%dt
   !   VFactor=1
   !   write(*,*)"output",x1(1);
   !  call  One%POI%VelRes(VFactor)
   !   call One%Updater(SpecyGlobal(0),GasGlobal(1))

   !    Call SelectProbility(One,MCCBundleGlobal(0,1))
   !   Index=One%ReactionIndex
      ! If (Index>0) Then
   !   Index=13
   !  !  write(*,*)"11_4",11_4
   !  !  MCCBundleGlobal(0,1)%Reaction(Index)=1
         
   !    ! End If
   !        Index=113
   !        ! MCCBundleGlobal(1,1)%Reaction(Index)=1
   !        Call  SelectCollisionIon(One,SpecyGlobal(1),GasGlobal(1),MCCBundleGlobal(1,1)%Reaction(Index))
   !    call One%POI%VelRes(1.d0/VFactor)

  
    ! ! j is the Particle Index;
    ! Integer(4),intent(in) :: Ns,Ng
    ! Type(ParticleBundle),intent(inout) :: PB(0:Ns)
    ! Type(SpecyOne),intent(in) :: SO(0:Ns)
    ! Type(GasOne),intent(in) :: GO(Ng)
    ! Type(MCCBundle) :: MCCB(0:Ns,Ng)
    ! Integer(4) :: i,j,k,m,Index
    ! Integer(4) :: CounterParticleAnnihilationIndex,NParTemp


    ! Type(MCCBundle) :: MCCB
    ! x1=x2
    ! v1=v2
    ! x3=v3
    ! MCCB%SigmaMax=
    ! MCCB%CollisionRatio
    ! call Random_Number(R)
    ! if( R<CollisionRatio) then
    !     Call SelectProbility(MCCParticleBundle(k),MCCB(j,i))
    !        Index=MCCParticleBundle(k)%ReactionIndex
    !       If (Index>0) Then
    !      Call  SelectCollisionElectron(MCCParticleBundle(k),SO(j),GO(i),MCCB(j,i)%Reaction(Index)) 
    !       end if
    ! end if

end subroutine 



subroutine InitParticleOne(v1,ispecies)bind(C,name="InitParticleOne")
   real(c_double) :: v1(3)
   integer(c_int)::ispecies
   real(8) mass,temperature
if(ispecies==0) then
   mass=ElectronMass
   temperature=30000000
   call aParticle%VelMaxInit(mass,temperature)
   ! write(*,*)"Velosity",aParticle%Vx,aParticle%Vy,aParticle%Vz
   v1(1)=aParticle%Vx
   v1(2)=aParticle%Vy
   v1(3)=aParticle%Vz
end if
   


end subroutine

end module MCCInterface

 ! init MCCB  //colloctionrato
    ! Gas
!     call InitGasJson(json_para, ControlFlowGlobal)
!     CF%Ng = Ng
!     "gas": {
!     "ng": 1,
!     "pressure": 25,
!     "temperature": 300,
!     "name": [
!         "Ar"
!     ],
!     "ratio": [
!         1.0
!     ],

!     "collision_section_type": 1,
!     "@collision_section_type": "1=collision_section_type_lxcat | 2=collision_section_type_pegasus",
!     "ng_max": 3,
!     "ns_max": 9
! }
! "init_density": 1e15,
! "init_electron_temperature": 30000,