module MCCInterface
use iso_c_binding
 Use ModuleTypeMCC
! use ModuleMCCInitialization
! Use ModuleParticleOne
! Use ModuleMCCElectron
! Use ModuleMCCIon
! Use ModuleMCCSigma
! Use ModuleParticleBundle
! use ModuleParticleBundleIndex
 use ModuleSpecyOne
 use ModuleReactionOnePegasus
 use ModuleMCCInitialization
implicit none

contains

subroutine  MCCBundleInit()bind(C,name="MCCBundleInit")
   
 Type(MCCBundle),Allocatable :: MCCBundleGlobal(:,:)
real(8) :: gas_pressure, gas_temperature, ele_temperature
character(len=99), dimension(1) :: gas_name
real(8) iongas_ratio,elegas_ratio
integer(4) model_type,Ng,Ns,collision_section_type
real(8) ,dimension(1)::gas_ratio
gas_pressure=25
gas_temperature=300
ele_temperature=30000
gas_name="Ar"
collision_section_type=2
model_type=2
Ns=1
Ng=1
gas_ratio(1)=1


 Allocate(MCCBundleGlobal(0:Ns,Ng))
 call GasInitPegasus(Ng,Ns,gas_pressure, gas_temperature, ele_temperature, gas_name, gas_ratio, collision_section_type)
 call MCCBundleElelctronInit(Ns,Ng,SpecyGlobal(0),GasGlobal(1:Ng),MCCBundleGlobal(0,1:Ng))
write(*,*)"specy",SpecyGlobal(0)%Charge
 ! call MCCBundleIonInit(Ns,Ng,SO(1:Ns),GO(1:Ng),MCCBundleGlobal(1:Ns,1:Ng)) 





 





 end subroutine MCCBundleInit

subroutine MCC(x1,v1,x2,v2,x3,v3,flag)bind(C ,name="MCC")
    real(c_double) :: x1,v1,x2,v2,x3,v3,flag
    Real(8) :: CollisionRatio,R
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