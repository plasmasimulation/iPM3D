! Module ModuleRecombination
!     Use ModuleParticleBundle
!     Use ModuleFieldOne
!     Use ModuleSpecyOne
!     Use ModuleMCCPublic
!     Implicit none
     
!         Type ElectronIonRecombination
!             Integer(4) :: Ns=0,Nx=NxMax,Nt=1,Ng=1
!             Integer(4) :: Timer=0
!             Integer(4),allocatable :: NsIndex(:)
!             Logical,allocatable :: doRecombination(:)!=.false.
!             Real(8) :: Dx=Inputdx,Dt=Inputdt
!             Real(8),allocatable :: RecombinationRate(:)
!         contains
!                 Procedure :: Init=>InitElectronIonRecombination
!                 Procedure :: InitConstant=>InitConstantElectronIonRecombination
!                 Procedure :: Update=>UpdateElectronIonRecombination
!                 Procedure :: Onestep=>OnestepElectronIonRecombination
!             ENd Type ElectronIonRecombination
           
!         Type IonIonRecombination
!            Integer(4) :: Ns=0,Nx=NxMax,Nt=1,Ng=1
!            Integer(4) :: Timer=0
!         Integer(4),allocatable :: NsIndex(:)
!            Real(8) :: Dx=Inputdx,Dt=Inputdt
!            Logical,allocatable :: doRecombination(:,:)!=.false.
!            Real(8),allocatable :: RecombinationRate(:,:)
!            contains
!                 Procedure :: Init=>InitIonIonRecombination
!                 Procedure :: InitConstant=>InitConstantIonIonRecombination
!                 Procedure :: Update=>UpdateIonIonRecombination
!                 Procedure :: Onestep=>OnestepIonIonRecombination
!         ENd Type IonIonRecombination
     
!      !Type()
!      !            
!       !      Real(8) :: K0=1.d-13 !G文献中 都是1d-13
     
     
!     contains
!         subroutine  OnestepElectronIonRecombination(EIR,PB,FO)
!             Class(ElectronIonRecombination),intent(inout) :: EIR
!             Type(ParticleBundle),intent(inout) :: PB(0:)
!             Type(FieldOne),intent(inout) :: FO(0:EIR%Ns)
!             Integer(4) :: i
!             EIR%Timer=EIR%Timer+1
                        
!             If(EIR%Timer==EIR%Nt) Then
!                 do i=0,EIR%Ns 
!                     Call PB(i)%WeightP2C2(FO(i)) 
!                 ENd Do
!                 Call  EIR%Update(PB,FO)
!                 EIR%Timer=0
!             ENd If
!            return
!         End Subroutine OnestepElectronIonRecombination
        
!         subroutine  OnestepIonIonRecombination(IIR,PB,FO)
!             Class(IonIonRecombination),intent(inout) :: IIR
!             Type(ParticleBundle),intent(inout) :: PB(0:)
!             Type(FieldOne),intent(inout) :: FO(0:IIR%Ns)
!             Integer(4) :: i
!             IIR%Timer=IIR%Timer+1
                        
!             If(IIR%Timer==IIR%Nt) Then
!                 do i=0,IIR%Ns 
!                     Call PB(i)%WeightP2C2(FO(i)) 
!                 ENd Do
!                 Call  IIR%Update(PB,FO)
!                 IIR%Timer=0
!             ENd If
!            return
!         End Subroutine OnestepIonIonRecombination
    
!         subroutine  InitConstantElectronIonRecombination(EIR,PB)
!             Class(ElectronIonRecombination),intent(inout) :: EIR
!             Type(ParticleBundle),intent(in) :: PB(0:)
!             Integer(4) :: i   
!             Do i=1,EIR%Ns
!                 If (PB(i)%Charge>0.d0) Then
!                     EIR%doRecombination(i)=.True.
!                     EIR%RecombinationRate(i)=3.95d-15 !1.d-13
!                 ELse
!                     EIR%doRecombination(i)=.false.
!                     EIR%RecombinationRate(i)=0.d0
!                 ENd If
!             ENd Do

!            return
!         End Subroutine InitConstantElectronIonRecombination
        
!         subroutine  InitElectronIonRecombination(EIR,CF,GO)
!             Class(ElectronIonRecombination),intent(inout) :: EIR
!             Type(ControlFlow),intent(in) :: CF
!             Type(GasOne),intent(in) :: GO(:)
!             Integer(4) :: i
!             EIR%Dx=CF%Dx
!             EIR%Dt=CF%Dt
            
!             EIR%Ns=CF%Ns
!             EIR%Nx=CF%Nx
!             EIR%Nt=CF%Nt
!             EIR%Timer=0
            
!             EIR%Ng=CF%Ng
            
!             Allocate(EIR%NsIndex(EIR%Ng))
!             Do i=1,Size(GO)
!                EIR%NsIndex(i)=GO(i)%IndexStart
!             ENd Do
            
!             Allocate(EIR%doRecombination(EIR%Ns))
!             EIR%doRecombination=.false.
!             Allocate(EIR%RecombinationRate(EIR%Ns))
!             EIR%RecombinationRate=0.d0
!            return
!         End Subroutine InitElectronIonRecombination
        
!         subroutine  InitConstantIonIonRecombination(IIR,PB)
!             Class(IonIonRecombination),intent(inout) :: IIR
!             Type(ParticleBundle),intent(in) :: PB(0:)
!             Integer(4) :: i,j
!             Logical :: isPositive(IIR%Ns)
!             Do i=1,IIR%Ns
!                 If (PB(i)%Charge>0.d0) Then
!                     isPositive(i)=.true.
!                 ELse If (PB(i)%Charge<0.d0) Then
!                     isPositive(i)=.false.
!                 ENd If
!             ENd Do
!             IIR%doRecombination=.false.
!             IIR%RecombinationRate=0.d0 
!             Do i=1,IIR%Ns
!                 DO j=1,IIR%Ns
!                     If (j<=i) Then
!                         IIR%doRecombination(j,i)=.false.
!                     ELse
!                         If(isPositive(j)==.true..and.isPositive(i)==.false.) Then
!                             IIR%doRecombination(j,i)=.true.
!                             IIR%RecombinationRate(j,i)=1.d-13
!                         Else If(isPositive(j)==.false..and.isPositive(i)==.true.) Then
!                             IIR%doRecombination(j,i)=.true.
!                             IIR%RecombinationRate(j,i)=1.d-13
!                         ELse
!                             IIR%doRecombination(j,i)=.false.
!                             IIR%RecombinationRate(j,i)=0.d0
!                         End IF
!                     End IF    
!                 End DO
!             ENd do
            
!            return
!         End Subroutine InitConstantIonIonRecombination
    
!         subroutine  InitIonIonRecombination(IIR,CF,GO)
!             Class(IonIonRecombination),intent(inout) :: IIR
!             Type(ControlFlow),intent(in) :: CF
!             Type(GasOne),intent(in) :: GO(:)
!             Integer(4) :: i
!             IIR%Dx=CF%Dx
!             IIR%Dt=CF%Dt
            
!             IIR%Ns=CF%Ns
!             IIR%Nx=CF%Nx
!             IIR%Nt=CF%Nt
!             IIR%Timer=0
            
!             IIR%Ng=CF%Ng
!             Allocate(IIR%NsIndex(IIR%Ng))
!             Do i=1,Size(GO)
!                IIR%NsIndex(i)=GO(i)%IndexStart
!             ENd Do
            
!             Allocate(IIR%doRecombination(IIR%Ns,IIR%Ns))
!             IIR%doRecombination=.false.
!             Allocate(IIR%RecombinationRate(IIR%Ns,IIR%Ns))
!             IIR%RecombinationRate=0.d0
!            return
!         End Subroutine InitIonIonRecombination
        
!         !subroutine  InitConstantElectronIonRecombination(IIR,PB)
!         !    Class(IonIonRecombination),intent(inout) :: IIR
!         !    Type(ParticleBundle),intent(in) :: PB(:)
!         !    
!         !    
!         !      
!         !   return
!         !End Subroutine InitConstantElectronIonRecombination
    
!         subroutine  UpdateIonIonRecombination(IIR,PB,FO)
!             Class(IonIonRecombination),intent(inout) :: IIR
!             Type(ParticleBundle),intent(inout) :: PB(0:IIR%Ns)
!             Type(FieldOne),intent(in) :: FO(0:IIR%Ns)
!             Integer(4) :: i,j,Nx,N,Ns
!             Real(8) :: NumberPositive(IIR%Nx),NumberNegative(IIR%Nx),NumberRecombination(IIR%Nx),NumberRecombinationTemp(IIR%Nx)
! 			Real(8) :: Weighting,Dx
!             Weighting = PB(1)%Weight*IIR%Dx
!             Dx=IIR%Dx
!             Nx=IIR%Nx
!             Ns=IIR%Ns

!             Do i=1,Ns
!                DO j=1,Ns 
!                   If(IIR%doRecombination(j,i)) Then
!                         NumberPositive=0.d0
!                         NumberNegative=0.d0
!                         NumberRecombination=0.d0
!                         NumberRecombinationTemp=0.d0
!                         NumberPositive = abs(FO(j)%RhoOne/PB(j)%charge)*Dx
!                         NumberNegative = abs(FO(i)%RhoOne/PB(i)%charge)*Dx
!                         NumberRecombination(1:Nx) = NumberNegative(1:Nx)*IIR%RecombinationRate(j,i)*IIR%dt*Dble(IIR%Nt)*(NumberPositive(1:Nx)/Dx)/Weighting
!                         NumberRecombination(1) = 0.5*NumberRecombination(1)
!                         NumberRecombination(Nx) = 0.5*NumberRecombination(Nx)
!                         NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx)
!                         Call Recombination(PB(j),Nx,NumberRecombinationTemp) !处理复合的正离子-del
!                         NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx) !因为上面的调用函数已经改变了NumberRecombinationTemp值，需要重新赋值
!                         Call Recombination(PB(i),Nx,NumberRecombinationTemp)  !处理复合的负离子-del
!                     ENd IF
!                 End DO
!             ENd do      
!            return
!         End Subroutine UpdateIonIonRecombination
        
!         !Weighting = PB(0)%Weight*CF%Dx
! 		!		Nx = CF%Nx
! 		!		
!   !          NumberPositive=0.d0
!   !          NumberNegative=0.d0
!   !          NumberRecombination=0.d0
!   !          NumberRecombinationTemp=0.d0
!   !          NumberPositive = abs(FO(TpyePositive)%RhoOne/PB(TpyePositive)%charge)*CF%Dx
!   !          NumberNegative = abs(FO(TpyeNegative)%RhoOne/PB(TpyeNegative)%charge)*CF%Dx
! 		!		NumberRecombination(1:Nx) = NumberNegative(1:Nx)*K0*dt*(NumberPositive(1:Nx)/CF%Dx)/Weighting
!   !          NumberRecombination(1) = 0.5*NumberRecombination(1)
! 		!		NumberRecombination(Nx) = 0.5*NumberRecombination(Nx)
!   !          NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx)
!   !          Call Recombination(PB(TpyePositive),Nx,NumberRecombinationTemp) !处理复合的正离子-del
!   !          
!   !          NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx) !因为上面的调用函数已经改变了NumberRecombinationTemp值，需要重新赋值
!   !          Call Recombination(PB(TpyeNegative),Nx,NumberRecombinationTemp)  !处理复合的负离子-del 
        
        
!         subroutine UpdateElectronIonRecombination(EIR,PB,FO)
!             Class(ElectronIonRecombination),intent(inout) :: EIR
!             Type(ParticleBundle),intent(inout) :: PB(0:EIR%Ns)
!             Type(FieldOne),intent(in) :: FO(0:EIR%Ns)
!             Integer(4) :: i,j,Nx,N
!             Real(8) :: K0(EIR%Nx)
!             Real(8) :: NumberPositive(EIR%Nx),NumberElectron(EIR%Nx),NumberRecombination(EIR%Nx),NumberRecombinationTemp(EIR%Nx)
!             Real(8) :: Weighting,Dx !Add By YuShimin
				
!             Weighting = PB(0)%Weight*EIR%Dx
! 			Nx = EIR%Nx
!             Dx=EIR%Dx
!             Do i=1,EIR%Ns
!                 If(EIR%doRecombination(i)) Then
!                 NumberElectron=0.d0
!                 NumberPositive=0.d0
!                 NumberRecombination=0.d0
!                 NumberRecombinationTemp=0.d0
!                 NumberElectron = abs(FO(0)%RhoOne/PB(0)%charge)*Dx
!                 NumberPositive = abs(FO(i)%RhoOne/PB(i)%charge)*Dx
!                 Do j=1,Nx
! 					if((DSQRT(FO(0)%TOne(j))*FO(i)%TOne(j)) > 0.d0) then
! 					   !K0(i) = a_ElectronIonRecombinationCoefficient * (1.d0/DSQRT((2.d0/3.d0)*FO%TOne(i)))
! 					   K0(j) = EIR%RecombinationRate(i)/(DSQRT(FO(0)%TOne(j))*FO(i)%TOne(j) )
! 				  else 
! 					   K0(j) =0.d0
! 				  end if
!                 End Do
! 				NumberRecombination(1:Nx) = NumberElectron(1:Nx)*K0(1:Nx)*EIR%dt*Dble(EIR%Nt)*(NumberPositive(1:Nx)/EIR%Dx)/Weighting
! 				NumberRecombination(1) = 0.5*NumberRecombination(1)
! 				NumberRecombination(Nx) = 0.5*NumberRecombination(Nx)
!                 NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx)
!                 Call Recombination(PB(i),Nx,NumberRecombinationTemp) !处理复合的正离子-del
!                 NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx) !因为上面的调用函数已经改变了NumberRecombinationTemp值，需要重新赋值
!                 Call Recombination(PB(0),Nx,NumberRecombinationTemp)  !处理复合的负离子-del
!                 ENd IF
!             End Do
!                 !Write (*,*)  ParticlePositive%Name,ParticlePositive%NPar,ParticleNegative%Name,ParticleNegative%NPar,'After'       
!         return
!        end subroutine UpdateElectronIonRecombination
        
    
!      subroutine  IonIonRecombinationRate(CF,PB,FO,TpyePositive,TpyeNegative,dt) 
!         implicit none
! 		      Type(ControlFlow),INTENT(IN) :: CF
!             Type(ParticleBundle),intent(inout) :: PB(0:CF%Ns)
! 				Type(FieldOne),intent(inout) :: FO(0:CF%Ns)
! 				Integer(4),intent(in) :: TpyeNegative,TpyePositive
!             Real(8),intent(in) :: dt
!             Real(8) :: K0=1.d-13 !G文献中 都是1d-13
!             Integer(4) :: i,Nx,N
!             Real(8) :: NumberPositive(CF%Nx),NumberNegative(CF%Nx),NumberRecombination(CF%Nx),NumberRecombinationTemp(CF%Nx)
! 				Real(8) :: Weighting
				
!             Weighting = PB(0)%Weight*CF%Dx
! 				Nx = CF%Nx
				
!             NumberPositive=0.d0
!             NumberNegative=0.d0
!             NumberRecombination=0.d0
!             NumberRecombinationTemp=0.d0
!             NumberPositive = abs(FO(TpyePositive)%RhoOne/PB(TpyePositive)%charge)*CF%Dx
!             NumberNegative = abs(FO(TpyeNegative)%RhoOne/PB(TpyeNegative)%charge)*CF%Dx
! 				NumberRecombination(1:Nx) = NumberNegative(1:Nx)*K0*dt*(NumberPositive(1:Nx)/CF%Dx)/Weighting
!             NumberRecombination(1) = 0.5*NumberRecombination(1)
! 				NumberRecombination(Nx) = 0.5*NumberRecombination(Nx)
!             NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx)
!             Call Recombination(PB(TpyePositive),Nx,NumberRecombinationTemp) !处理复合的正离子-del
            
!             NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx) !因为上面的调用函数已经改变了NumberRecombinationTemp值，需要重新赋值
!             Call Recombination(PB(TpyeNegative),Nx,NumberRecombinationTemp)  !处理复合的负离子-del 
!            return
!        End Subroutine  IonIonRecombinationRate
!   !     
!          subroutine ElectronIonRecombinationRate(CF,PB,FO,TpyeIon,dt) 
!         implicit none
! 		    Type(ControlFlow),INTENT(IN) :: CF
!             Type(ParticleBundle),intent(inout) :: PB(0:CF%Ns)
! 				Type(FieldOne),intent(inout) :: FO(0:CF%Ns)
! 				Integer(4),intent(in) :: TpyeIon
!             Real(8),intent(in) :: dt 
!             Integer(4) :: i,Nx,N
!             Real(8) :: K0(CF%Nx),NumberPositive(CF%Nx),NumberElectron(CF%Nx),NumberRecombination(CF%Nx),NumberRecombinationTemp(CF%Nx)
!             Real(8) :: Weighting,a_ElectronIonRecombinationCoefficient !Add By YuShimin
				
!             Weighting = PB(0)%Weight*CF%Dx
! 				Nx = CF%Nx
!             K0=1.d-13
!             NumberPositive=0.d0
!             NumberElectron=0.d0
!             NumberRecombination=0.d0
!             NumberRecombinationTemp=0.d0
! 				NumberElectron = abs(FO(0)%RhoOne/PB(0)%charge)*CF%Dx
  
! 				NumberPositive = abs(FO(TpyeIon)%RhoOne / PB(TpyeIon)%charge) *CF%Dx
!             a_ElectronIonRecombinationCoefficient = 3.95d-15 !e-Ar+
!             Do i=1,Nx
! 					if((DSQRT(FO(0)%TOne(i))*FO(TpyeIon)%TOne(i)) > 0.d0) then
! 					   !K0(i) = a_ElectronIonRecombinationCoefficient * (1.d0/DSQRT((2.d0/3.d0)*FO%TOne(i)))
! 					   K0(i) = a_ElectronIonRecombinationCoefficient/(DSQRT(FO(0)%TOne(i))*FO(TpyeIon)%TOne(i) )
! 				  else 
! 					   K0(i) =0.d0
! 				  end if
!             End Do
! 				NumberRecombination(1:Nx) = NumberElectron(1:Nx)*K0*dt*(NumberPositive(1:Nx)/CF%Dx)/Weighting
! 				NumberRecombination(1) = 0.5*NumberRecombination(1)
! 				NumberRecombination(Nx) = 0.5*NumberRecombination(Nx)
!             NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx)
				
!             Call Recombination(PB(TpyeIon),Nx,NumberRecombinationTemp) !处理复合的正离子-del
            
!             NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx) !因为上面的调用函数已经改变了NumberRecombinationTemp值，需要重新赋值
!             Call Recombination(PB(0),Nx,NumberRecombinationTemp)  !处理复合的负离子-del
!             !Write (*,*)  ParticlePositive%Name,ParticlePositive%NPar,ParticleNegative%Name,ParticleNegative%NPar,'After'       
!         return
!   end subroutine ElectronIonRecombinationRate
!   !     
!        Subroutine  Recombination(PB,Nx,NumberRecombination)
!                Implicit none
!                Integer(4),intent(in) ::  Nx
!                Type(ParticleBundle),intent(inout) :: PB
!                Real(8),intent(inout) :: NumberRecombination(Nx)
!                Integer(4) :: i,N,IndexBegin,j
!                Real(8) :: Ratio,S1,S2
!                    CALL RANDOM_NUMBER(R)
!                    IndexBegin=Ceiling(R*dble(PB%Npar))  !随机一个遍历起始点
!                    do i=IndexBegin,1,-1   !这里是分成两部分，将数组都遍历了一下，---改进：可以做一个总数，附加判断条件一旦总数达到 就不进行进一步的判断了
!                        N=Ceiling(PB%PO(i)%X) 
! 						     S1=Dble(N)-PB%PO(i)%X
! 						     !S2 = 1.D0-S1
! 							  CALL RANDOM_NUMBER(R)
! 						     IF (R<S1) then
!                          If(NumberRecombination(N+1)>MinReal) Then
!                            If(NumberRecombination(N+1)>=1.d0)  Then
!                                NumberRecombination(N+1)=NumberRecombination(N+1)-1.d0
!                                Call PB%DelOne(i)
!                                !DelParticle(i,InputParticle)
!                             Else
!                                  CALL RANDOM_NUMBER(R)
!                                  Ratio=NumberRecombination(N+1)
!                                  NumberRecombination(N+1)=NumberRecombination(N+1)-1.d0
!                                  If(R<Ratio)  Then
!                                     Call PB%DelOne(i)
!                                     !DelParticle(i,InputParticle)
!                                  End If 
!                              End if
! 								 End if
! 							  else
! 								  If(NumberRecombination(N)>MinReal) Then
!                            If(NumberRecombination(N)>=1.d0)  Then
!                                NumberRecombination(N)=NumberRecombination(N)-1.d0
!                                Call PB%DelOne(i)
!                                !DelParticle(i,InputParticle)
!                             Else
!                               CALL RANDOM_NUMBER(R)
!                               Ratio=NumberRecombination(N)
!                               NumberRecombination(N)=NumberRecombination(N)-1.d0
!                               If(R<Ratio)  Then
!                                     Call PB%DelOne(i)
!                                     !DelParticle(i,InputParticle)
!                               End If 
!                             End if
! 								 End if
! 							  end if
!                  End do
!                 i=IndexBegin+1
!                 do while (i<=PB%NPar)
!                        N=Ceiling(PB%PO(i)%X) 
!                        S1=Dble(N)-PB%PO(i)%X
! 						     !S2 = 1.D0-S1
! 							  CALL RANDOM_NUMBER(R)
! 						     IF (R<S1) then
!                          If(NumberRecombination(N+1)>MinReal) Then
!                            If(NumberRecombination(N+1)>=1.d0)  Then
!                                NumberRecombination(N+1)=NumberRecombination(N+1)-1.d0
!                                Call PB%DelOne(i)
!                                !DelParticle(i,InputParticle)
!                             Else
!                                  CALL RANDOM_NUMBER(R)
!                                  Ratio=NumberRecombination(N+1)
!                                  NumberRecombination(N+1)=NumberRecombination(N+1)-1.d0
!                                  If(R<Ratio)  Then
!                                     Call PB%DelOne(i)
!                                     !DelParticle(i,InputParticle)
!                                  End If 
!                              End if
! 								 End if
! 							  else
! 								  If(NumberRecombination(N)>MinReal) Then
!                            If(NumberRecombination(N)>=1.d0)  Then
!                                NumberRecombination(N)=NumberRecombination(N)-1.d0
!                                Call PB%DelOne(i)
!                                !DelParticle(i,InputParticle)
!                             Else
!                               CALL RANDOM_NUMBER(R)
!                               Ratio=NumberRecombination(N)
!                               NumberRecombination(N)=NumberRecombination(N)-1.d0
!                               If(R<Ratio)  Then
!                                     Call PB%DelOne(i)
!                                     !DelParticle(i,InputParticle)
!                               End If 
!                             End if
! 								 End if
! 							  end if
!                     i=i+1
!                  End do      
!                return
! 		 End Subroutine  Recombination
!   !
! 		! 
! 	 !  Subroutine  IonIonRecombinationRateOld(Nx,ParticlePositive,ParticleNegative,Weighting,dt)
!   !          Implicit none
!   !          Type(ParticleBundle),intent(inout) :: ParticlePositive,ParticleNegative
!   !          Real(8),intent(in) :: Weighting,dt
!   !          Real(8),parameter :: K0=1.d-13 !G文献中 都是1d-13
!   !          Integer(4) :: i,Nx,N
!   !          Real(8) :: NumberPositive(Nx), NumberNegative(Nx),NumberRecombination(Nx),NumberRecombinationTemp(Nx)
!   !          
!   !          NumberPositive=0.d0
!   !          NumberNegative=0.d0
!   !          NumberRecombination=0.d0
!   !          NumberRecombinationTemp=0.d0
!   !          !Write (*,*)  ParticlePositive%Name,ParticlePositive%NPar,ParticleNegative%Name,ParticleNegative%NPar,'Before'
!   !          do i=1,ParticlePositive%NPar
!   !                     N=Ceiling(ParticlePositive%PO(i)%X)
!   !                     NumberPositive(N)=NumberPositive(N)+1.d0   !求出每个网格内粒子数 N
!   !          End do
!   !         
!   !          do i=1,ParticleNegative%NPar
!   !                     N=Ceiling(ParticleNegative%PO(i)%X) !同上
!   !                     NumberNegative(N)=NumberNegative(N)+1.d0 
!   !          End do
!   !          NumberRecombination(1:Nx)=Weighting*K0*dt*NumberPositive(1:Nx)*NumberNegative(1:Nx) !Na*w*k*dt/V*Nb, Na*w*k*dt/V=n*v*sigame 注意(k=v*sigma)
!   !          NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx)
!   !          Call RecombinationOld(ParticlePositive,Nx,NumberRecombinationTemp) !处理复合的正离子-del
!   !          NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx) !因为上面的调用函数已经改变了NumberRecombinationTemp值，需要重新赋值
!   !          Call RecombinationOld(ParticleNegative,Nx,NumberRecombinationTemp)  !处理复合的负离子-del
!   !          !Write (*,*)  ParticlePositive%Name,ParticlePositive%NPar,ParticleNegative%Name,ParticleNegative%NPar,'After'       
!   !         return
!   !     End Subroutine  IonIonRecombinationRateOld
!   !     
!   !     subroutine ElectronIonRecombinationRateOld(Nx,FO,ParticlePositive,Electron,Weighting,dt) 
!   !      implicit none
!   !      Type(ParticleBundle),intent(inout) :: ParticlePositive,Electron
!   !          Real(8),intent(in) :: Weighting,dt 
!   !          Integer(4) :: i,Nx,N
!   !          Real(8) :: K0(Nx),NumberPositive(Nx), NumberElectron(Nx),NumberRecombination(Nx),NumberRecombinationTemp(Nx)
!   !          
!   !          Real(8) :: a_ElectronIonRecombinationCoefficient !Add By YuShimin
!   !          !Integer(4),intent(in) :: Ns
!   !          Type(FieldOne),intent(inout) :: FO
!   !          K0=1.d-13
!   !          NumberPositive=0.d0
!   !          NumberElectron=0.d0
!   !          NumberRecombination=0.d0
!   !          NumberRecombinationTemp=0.d0
!   !          !Write (*,*)  ParticlePositive%Name,ParticlePositive%NPar,ParticleNegative%Name,ParticleNegative%NPar,'Before'
!   !          do i=1,ParticlePositive%NPar
!   !                     N=Ceiling(ParticlePositive%PO(i)%X)
!   !                     NumberPositive(N)=NumberPositive(N)+1.d0   !求出每个网格内粒子数 N
!   !          End do
!   !          do i=1,Electron%NPar
!   !                     N=Ceiling(Electron%PO(i)%X) !同上
!   !                     NumberElectron(N)=NumberElectron(N)+1.d0 
!   !          End do
!   !          a_ElectronIonRecombinationCoefficient = 0.5d-13 !e-Ar+
!   !          !do i=0,ControlFlowGlobal%Ns
!   !          !     !Call ParticleAborption(ParticleGlobal(i),ParticleBDOneGlobal(i))
!   !          !     Call ParticleGlobal(i)%WeightP2C(FieldOneGlobal(i))
!   !          !end do
!   !          Do i = 1,Nx
!   !              K0(i) = a_ElectronIonRecombinationCoefficient * (1.d0/DSQRT((2.d0/3.d0)*FO%TOne(i)))
!   !              !K0(i) = 1.d-13
!   !          End Do
!   !          NumberRecombination(1:Nx)=Weighting*K0(1:Nx)*dt*NumberPositive(1:Nx)*NumberElectron(1:Nx) !Na*w*k*dt/V*Nb, Na*w*k*dt/V=n*v*sigame 注意(k=v*sigma)
!   !          NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx)
!   !          Call RecombinationOld(ParticlePositive,Nx,NumberRecombinationTemp) !处理复合的正离子-del
!   !          NumberRecombinationTemp(1:Nx)=NumberRecombination(1:Nx) !因为上面的调用函数已经改变了NumberRecombinationTemp值，需要重新赋值
!   !          Call RecombinationOld(Electron,Nx,NumberRecombinationTemp)  !处理复合的负离子-del
!   !          !Write (*,*)  ParticlePositive%Name,ParticlePositive%NPar,ParticleNegative%Name,ParticleNegative%NPar,'After'       
!   !      return
!   !end subroutine ElectronIonRecombinationRateOld
!   !     
!   !     Subroutine  RecombinationOld(InputParticle,Nx,NumberRecombination)
!   !             Implicit none
!   !             Integer(4),intent(in) ::  Nx
!   !             Type(ParticleBundle),intent(inout) :: InputParticle
!   !             Real(8),intent(inout) :: NumberRecombination(Nx)
!   !             Integer(4) :: i,N,IndexBegin
!   !             Real(8) :: Ratio
!   !             
!   !                 CALL RANDOM_NUMBER(R)
!   !                 IndexBegin=Ceiling(R*dble(InputParticle%Npar))  !随机一个遍历起始点
!   !                   
!   !                 do i=IndexBegin,1,-1   !这里是分成两部分，将数组都遍历了一下，---改进：可以做一个总数，附加判断条件一旦总数达到 就不进行进一步的判断了
!   !                     N=Ceiling(InputParticle%PO(i)%X) 
!   !                     If(NumberRecombination(N)>MinReal) Then
!   !                         If(NumberRecombination(N)>=1.d0)  Then
!   !                             NumberRecombination(N)=NumberRecombination(N)-1.d0
!   !                             Call InputParticle%DelOne(i)
!   !                             !DelParticle(i,InputParticle)
!   !                          Else
!   !                               CALL RANDOM_NUMBER(R)
!   !                               Ratio=NumberRecombination(N)
!   !                               NumberRecombination(N)=NumberRecombination(N)-1.d0
!   !                               If(R<Ratio)  Then
!   !                                  Call InputParticle%DelOne(i)
!   !                                  !DelParticle(i,InputParticle)
!   !                               End If 
!   !                           End if
!   !                      End if
!   !               End do
!   !               
!   !              i=IndexBegin+1
!   !              do while (i<=InputParticle%NPar)
!   !                     N=Ceiling(InputParticle%PO(i)%X) 
!   !                     If(NumberRecombination(N)>MinReal) Then
!   !                         If(NumberRecombination(N)>=1.d0)  Then
!   !                             NumberRecombination(N)=NumberRecombination(N)-1.d0
!   !                             Call InputParticle%DelOne(i)
!   !                             !DelParticle(i,InputParticle)
!   !                          Else
!   !                               CALL RANDOM_NUMBER(R)
!   !                               Ratio=NumberRecombination(N)
!   !                               NumberRecombination(N)=NumberRecombination(N)-1.d0
!   !                             If(R<Ratio)  Then
!   !                               Call InputParticle%DelOne(i)
!   !                               !DelParticle(i,InputParticle)
!   !                              End If 
!   !                           End if
!   !                      End if
!   !                        i=i+1
!   !               End do      
!   !             return
!   !        End Subroutine  RecombinationOld	 
		 
		 
		 
    
! End Module ModuleRecombination