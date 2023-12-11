Module ModuleReactionOnePegasus
    use Constants
    use Numrical
    use ModuleMCCSigma
    
    implicit none

    Integer(4),parameter, private :: NspecyMax=9
    Integer(4),parameter, private :: NcrosssectionMax=500

    
    Type  ReactionOneBase    ! Define a collosion process.
          Real(8) :: Threshold !最小反应能量
          Integer(4) :: NSpecy !参加粒子种类数量
          Integer(4) :: SpecyIndex(NspecyMax) !参加反应的粒子编号
          Integer(4) :: NumberChange(NspecyMax) !反应过程中消失一个-1，生成一个+1，不变为0
          Integer(4) :: ReactionType !反应类型
          contains
          Procedure :: LoadReaction=>LoadReactionOneBase
          Procedure :: Convert=>ConvertReactionOneBase2ReactionOne
          Procedure :: FindIndex=>FindResultantIndex
    End Type ReactionOneBase

    Type,extends(ReactionOneBase) :: ReactionOnePegasus
          Integer(4) :: Index !编号
          Character(CharLenthMax) :: Name="Ar" !反应方程式文字表达式
          Integer(4) :: Ncrosssection !碰撞截面数据行数
          Real(8) :: Crosssection(2,NcrosssectionMax) !能量   碰撞面积（cm2）
          contains
          Procedure :: Load=>LoadReactionOnePegasus
          Procedure :: Convert2SR=>ConvertReactionOnePegasus2SigmaRaw
    End Type ReactionOnePegasus

    contains
    
       subroutine  SigmaNormalizationPegasus(SN,SO,GO,Inited)
                    Implicit none
                    Type(SigmaNormalized),Intent(out) :: SN !, intent(inout)
                    Type(SpecyOne),Intent(in)  :: SO
                    Type(GasOne),Intent(in) :: GO
                    Logical,Intent(inout) :: Inited
                    Type(SigmaRaw),Allocatable :: SR(:)
                    Type(ReactionOnePegasus) :: ROP
                    Integer(4) :: Model,NReaction,i
                    !Integer(4) :: CSIndex
                    Character(len=99) :: Filename,CSname,CSIndex,ReactionType,Threshold
                    !Logical :: Alive
                    !NAMELIST /MRList/ Model,NReaction
                    !NAMELIST /SRList/ SR
                    
                    !GO%Name="Ar"
                        
                    Filename="./input/gas/pegasus/"//Trim(GO%Name)//"/"//Trim(SO%Name)//Trim(GO%Name)//".txt"

                    Inquire(file=Filename,exist=Inited)
                    If(Inited) then
                         OPEN(17,FILE=Filename)
                         Read(17,*) Model,NReaction
                         Call SigmaNormalizedInit(SN,Model,NReaction)
                         Allocate(SR(NReaction))
                         
                        Do i=1,NReaction
                            Read(17,*) CSIndex,ReactionType,Threshold
                            CSname="./input/gas/pegasus/"//Trim(GO%Name)//"/"//Trim(CSIndex)
                            !OPEN(18,FILE=CSname)
                            Call ROP%Load(CSname)
                            Call ROP%Convert2SR(GO%SOP(1:),SR(i))
                            !close(18)
                        End Do
                           close(17)
                            Select case(Model)
                              Case(1)
                                Call SigmaNormalizedUpdate(SN,SR)
                              Case(2)
                                Call SigmaNormalizedUpdate(SN,SR)
                              Case(3)
                                Call SigmaNormalizedUpdateReactive(SN,SR,SO,GO)
                            End Select
                            Call UpdateReactionIndex(SN,GO)
                            Call SigmaNormalizedFinalization(SN)
                            
                          DeAllocate(SR)
                    Else
                        Call SigmaNormalizedInit(SN,2,1)
                          !Inited=.true.
                  End If
                  return
      End  subroutine SigmaNormalizationPegasus
       
    
        
    
        subroutine LoadReactionOneBase(ROB,IONumber)
                Class(ReactionOneBase),intent(inout) :: ROB
                Integer(4) :: IONumber
                    Call skip_comments(IONumber)
                    Read(IONumber,*) ROB%Threshold
                    Read(IONumber,*) ROB%NSpecy
                    Read(IONumber,*) ROB%SpecyIndex(1:ROB%NSpecy)
                    Read(IONumber,*) ROB%NumberChange(1:ROB%NSpecy)
                    Read(IONumber,*) ROB%ReactionType
                return
        end subroutine LoadReactionOneBase
    
         subroutine LoadReactionOnePegasus(ROP,Filename)
                Class(ReactionOnePegasus),intent(inout) :: ROP
                Character(len=99) :: Filename
                Logical :: alive,Status
                Integer(4) :: i
                integer(4) :: size, rank, ierr

                call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
                call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

                !Logical :: alive
                !Character(len=99) :: Filename
                !Character(len=99) :: TempChar(10)
                !Integer(4) :: i,NameIndex,NPArMax
                !NameIndex=DefaultNameIndex
                !Write(filename,*) NameIndex,trim(ROP%SO%Name),".dat"
                !Write(filename,*) "180000.dat"
                !Open (10,file=filename)
                !Write(10,*) 11
                !close(10)
                !!
                Inquire(file=Filename,exist=alive)
                If(alive) then
                       Open (10,file=filename)
                       Call skip_comments(10)
                       Read(10,*) ROP%Index
                       Read(10,*) ROP%Name
                       !Read(10,*) ROP%Threshold
                       !Read(10,*) ROP%NSpecy
                       !Read(10,*) ROP%SpecyIndex(1:ROP%NSpecy)
                       !Read(10,*) ROP%NumberChange(1:ROP%NSpecy)
                       !Read(10,*) ROP%ReactionType
                       Call ROP%ReactionOneBase%LoadReaction(10)
                       Read(10,*) ROP%Ncrosssection
                       do i=1,ROP%Ncrosssection
                             Read(10,*)  ROP%Crosssection(:,i)
                       end do
                       !Status=.False.
                       if (0 == rank) Write(*,*) 'Loading ',Trim(ROP%Name),' Complete!' 
                       Close(10)
                 else
                       if (0 == rank) Write(*,*) 'Can not find the file, the particle will be randomly initilalized.'    
                       Status=.True. 
                 End if    
                 
                return
         end subroutine  LoadReactionOnePegasus
         
         !subroutine ConvertReactionType(ROB,SOP,ReactionType)
         !       Class(ReactionOneBase),intent(inout) :: ROB
         !       Type(SpecyOnePegasus),intent(in) :: SOP(1:)
         !       Integer(4),intent(inout) :: ReactionType
         !       Logical :: isElectron=.false.
         !       Integer(4) :: i
         !       Integer(4) :: NResultant=0,IndexResultant(2)=-1
         !       DO i=1,2
         !           If(ROB%SpecyIndex(i)==0) Then
         !              isElectron=.True.
         !           Else
         !              isElectron=.false. 
         !           End If    
         !       End Do
         !                       
         !       !DO i=3,NSpecy
         !       !    
         !       !NReaction=0    
         !       !
         !       If(isElectron==.True.) Then
         !           Select case(ROB%ReactionType)
         !               Case(1)
         !                   ReactionType=11
         !               Case(100)
         !                   ReactionType=12
         !               Case(2)
         !                   ReactionType=13
         !               Case(3)
         !                   ReactionType=14
         !               Case(4)
         !                   ReactionType=15
         !               End Select
         !       Else
         !           Select case(ROB%ReactionType)
         !               Case(1)
         !                   ReactionType=101
         !               Case(6)
         !                   ReactionType=103
         !               End Select
         !       ENd If
         !       return
         !end subroutine ConvertReactionType
         
         subroutine ConvertReactant(ROB,SOP,Reactant)
                Class(ReactionOneBase),intent(inout) :: ROB
                Type(SpecyOnePegasus),intent(in) :: SOP(0:)
                Integer(4),intent(inout) :: Reactant
                Logical :: isElectron=.false.
                Integer(4) :: i
                DO i=1,2
                   If(ROB%SpecyIndex(i)==0) Then
                      isElectron=.True.
                      Reactant=0
                   Else
                      isElectron=.false. 
                      Reactant=-1
                   End If    
                End Do

                If(isElectron==.false. ) Then
                    DO i=1,2
                        If (ROB%SpecyIndex(i)>0) Then
                           Reactant=ROB%SpecyIndex(i)
                        ENd If
                    End Do
                ENd If
                return
         end subroutine ConvertReactant
        ! 
        !subroutine ConvertResultant(ROB,Resultant)
        !        Class(ReactionOneBase),intent(inout) :: ROB
        !        Integer(4),intent(inout) :: Resultant
        !        Logical :: isElectron=.false.
        !        Integer(4) :: i
        !        DO i=1,2
        !           If(ROB%SpecyIndex(i)==0) Then
        !              isElectron=.True.
        !           Else
        !              isElectron=.false. 
        !           End If    
        !        End Do
        !        
        !        If(isElectron==.True.) Then
        !            Select case(ROB%ReactionType)
        !                Case(1)
        !                    Resultant=0
        !                Case(100)
        !                    Resultant=0
        !                Case(2)
        !                    Resultant=ROB%SpecyIndex(3)
        !                Case(3)
        !                    Resultant=ROB%SpecyIndex(3)
        !                Case(4)
        !                    Resultant=10*ROB%SpecyIndex(3)+ROB%SpecyIndex(4)
        !                End Select
        !            
        !        Else
        !            Resultant=ROB%SpecyIndex(3)
        !        ENd If
        !        return
        ! end subroutine ConvertResultant
         
        Subroutine ConvertReactionOneBase2ReactionOne(ROB,SOP,RO)
                implicit none
                Class(ReactionOneBase),intent(inout) :: ROB
                Type(SpecyOnePegasus),intent(in) :: SOP(1:)
                Type(ReactionOne),intent(inout) :: RO
                Logical :: isElectron=.false.
                Integer(4) :: i,NResultant=0,IndexResultant(2)=-1
                
                Call ConvertReactionOneIndex(ROB,SOP)
                Associate (Reactant=>RO%Reactant,ReactionType=>RO%ReactionType,Resultant=>RO%Resultant,&
                    Threshold=>RO%Threshold)

                    Threshold=ROB%Threshold

                    DO i=1,2
                        If(ROB%SpecyIndex(i)==0) Then
                           isElectron=.True.
                        Else
                           isElectron=.false. 
                        End If    
                    End Do
                        
                    NResultant=0    
                    IndexResultant(2)=-1
                    
                    If(isElectron==.True.) Then
                        Select case(ROB%ReactionType)
                            Case(1)
                                ReactionType=11
                                Resultant=0
                            Case(100)
                                ReactionType=12
                                Resultant=0
                            Case(2)
                                ReactionType=13
                                Call ROB%FindIndex(NResultant,IndexResultant)
                                Resultant=IndexResultant(1)
                            Case(3)
                                Call ROB%FindIndex(NResultant,IndexResultant)
                                Select Case(NResultant)
                                Case(0)
                                    ReactionType=12
                                    Resultant=0
                                Case(1)
                                    ReactionType=14
                                    Resultant=IndexResultant(1)
                                Case(2)
                                    ReactionType=15
                                    Resultant=IndexResultant(1)*10+IndexResultant(2)
                                CASE DEFAULT
                                    Write(*,*) "Reaction Conversiong Error!!!"
                                    pause
                                ENd Select
                             Case(4)
                                Call ROB%FindIndex(NResultant,IndexResultant)
                                Select Case(NResultant)
                                Case(0)
                                    ReactionType=12
                                    Resultant=0
                                Case(1)
                                    DO i=1,2
                                        If(ROB%SpecyIndex(i)==0) Then
                                            Select Case(ROB%NumberChange(i))
                                            Case(0)
                                                ReactionType=12
                                                Resultant=IndexResultant(1)
                                            Case(-1)
                                                ReactionType=14
                                                Resultant=IndexResultant(1) 
                                            CASE DEFAULT
                                            Write(*,*) "Reaction Conversiong Error!!!"   
                                            pause
                                            ENd Select
                                        ENd IF
                                    ENd DO
                                Case(2)
                                    ReactionType=15
                                    Resultant=IndexResultant(1)*10+IndexResultant(2)
                                CASE DEFAULT
                                    Write(*,*) "Reaction Conversiong Error!!!"
                                    pause
                                ENd Select
                             Case(255)
                                Call ROB%FindIndex(NResultant,IndexResultant)
                                Select Case(NResultant)
                                Case(0)
                                    ReactionType=12
                                    Resultant=0
                                Case(1)
                                    DO i=1,2
                                        If(ROB%SpecyIndex(i)==0) Then
                                            Select Case(ROB%NumberChange(i))
                                            Case(1)
                                                ReactionType=13
                                                Resultant=IndexResultant(1)
                                            Case(0)
                                                ReactionType=12
                                                Resultant=0
                                            Case(-1)
                                                ReactionType=14
                                                Resultant=IndexResultant(1) 
                                            ENd Select
                                        ENd IF
                                    ENd DO
                                Case(2)
                                    ReactionType=15
                                    Resultant=IndexResultant(1)*10+IndexResultant(2)
                                CASE DEFAULT
                                    Write(*,*) "Reaction Conversiong Error!!!"
                                ENd Select 
                            End Select
                    Else
                        Select case(ROB%ReactionType)
                            Case(1)
                                ReactionType=101
                            Case(6)
                                ReactionType=103
                            Case(255)
                                Call ROB%FindIndex(NResultant,IndexResultant)
                                Select Case(NResultant)
                                Case(0)
                                    DO i=1,2
                                        If(ROB%SpecyIndex(i)==Reactant) Then
                                            Select Case(ROB%NumberChange(i))
                                            Case(0)
                                                ReactionType=101_4
                                                Resultant=ROB%SpecyIndex(i) 
                                            Case(-1)
                                                ReactionType=100_4
                                                Resultant=-1 
                                            ENd Select
                                        ENd IF
                                    ENd DO
                                Case(1)
                                    ReactionType=112
                                    Resultant=IndexResultant(1)
                                CASE DEFAULT
                                    Write(*,*) "Reaction Conversiong Error!!!"
                                ENd Select 
                            End Select
                    ENd If
   
                    End Associate
                    
                    
                    
                return
        end subroutine ConvertReactionOneBase2ReactionOne
        
        Subroutine FindResultantIndex(ROB,NResultant,IndexResultant)
                Class(ReactionOneBase),intent(inout) :: ROB
                !Type(SpecyOnePegasus),intent(in) :: SOP(1:)
                Integer(4),intent(inout) :: NResultant,IndexResultant(2)
                Integer(4) :: i,j

                NResultant=0
                IndexResultant(2)=-1

                Do i=3,ROB%NSpecy
                    If(ROB%SpecyIndex(i)>0) Then
                        NResultant=NResultant+1
                        IndexResultant(NResultant)=ROB%SpecyIndex(i)
                    ENd If
                ENd Do
                return
        end subroutine FindResultantIndex           
        
        Subroutine ConvertReactionOneIndex(ROB,SOP)
                implicit none
                Class(ReactionOneBase),intent(inout) :: ROB
                Type(SpecyOnePegasus),intent(in) :: SOP(1:)
                Integer(4) :: i,j,TempIndex(1:ROB%NSpecy)

                TempIndex=ROB%SpecyIndex(1:ROB%NSpecy)
                ROB%SpecyIndex=-1
                Do i=1,ROB%NSpecy
                    If(TempIndex(i)==1) Then
                       ROB%SpecyIndex(i)=0  
                    ELse   
                       DO j=1,Size(SOP)
                            If(TempIndex(i)==SOP(j)%Index) Then
                                ROB%SpecyIndex(i)=j
                                Exit                          
                            End If  
                        ENd Do
                    End If
                ENd Do    
                return
        end subroutine ConvertReactionOneIndex
        
        
        
        Subroutine ConvertReactionOnePegasus2SigmaRaw(ROG,SOP,SR)
                Class(ReactionOnePegasus),intent(inout) :: ROG
                Type(SpecyOnePegasus),intent(in) :: SOP(1:)
                Type(SigmaRaw),intent(inout) :: SR
                Call ROG%ReactionOneBase%Convert(SOP,SR%Reaction)
                
                SR%EnergySigma(1:2,1:ROG%NCrosssection)=ROG%Crosssection(1:2,1:ROG%NCrosssection)
                SR%EnergySigma(2,1:ROG%NCrosssection)=SR%EnergySigma(2,1:ROG%NCrosssection)/1.d4
                return
        end subroutine ConvertReactionOnePegasus2SigmaRaw
        
        Subroutine LoadGasPhysicalPegasus(GP,GasName,Inited)
                Type(GasPhysical),intent(inout) :: GP
                Character(len=99),intent(in) :: GasName
                logical,intent(inout) :: Inited

                Character(len=99) :: Filename
                Integer(4) ::i,j,NSOP
                logical :: Alive
                
                Filename="./input/gas/pegasus/"//Trim(GasName)//"/"//Trim(GasName)//".txt"
                Inquire(file=Filename,exist=alive)
                If(Alive) then
                    Inited=.false.
                    GP%CSModel=2
                    OPEN(20,FILE=Filename)
                        !Call skip_comments(10)
                        Read(20,*) GP%Name
                        Read(20,*) GP%MCModel
                        Read(20,*) GP%Natom
                        Read(20,*) NSOP
                        GP%Ns=NSOP-1
                        !Allocate(GP%SOP(0:NSOP-1))
                        Do j=0,NSOP-1
                              Call GP%SOP(j)%Load(20)
                        ENd Do
                    Close (20)
                Else
                   Inited=.True. 
                ENd IF
                
                If(Alive) then
                    Call GP%SOP(0)%Convert(GP)
                    Do j=1,NSOP-1
                        Call GP%SOP(j)%Convert(GP%SP(j))
                    ENd Do
                End IF
                return
        end subroutine LoadGasPhysicalPegasus
        
        Subroutine LoadGasPhysicalLxCat(GP,GasName,Inited)
                Type(GasPhysical),intent(inout) :: GP
                Character(len=99),intent(in) :: GasName
                logical,intent(inout) :: Inited
                NAMELIST /GasPhysicsfile/ GP
                Character(len=99) :: SOPName,Filename
                Integer(4) ::i,j,NSOP
                logical :: Alive
                
                Filename="./input/gas/lxcat/"//Trim(GasName)//"/"//Trim(GasName)//".txt"
                          Inquire(file=Filename,exist=alive)
                          If(alive) then
                               Inited=.false.
                               GP%CSModel=1
                               OPEN(20,FILE=Filename)
                               Read(20,NML=GasPhysicsfile)
                               Close (20)
                          ELse
                               Inited=.True.
                          End If
                return
        end subroutine LoadGasPhysicalLxCat      

        
        Subroutine  GasInitPegasus(CF, pressure, gas_temperature, ele_temperature, gas_name, gas_ratio, model_type)
            class(ControlFlow), intent(inout) :: CF
            real(8), intent(in) :: pressure, gas_temperature, ele_temperature
            character(len=99), dimension(CF%Ng) :: gas_name
            real(8), intent(in) :: gas_ratio(CF%Ng)
            integer(4), intent(in) :: model_type
            type(GasPhysical) :: GP
            type(GasPhysical), Allocatable :: TempGP(:)
            integer(4) :: Ns = 0
            integer(4) :: i, j, k, Index = 0
            logical    :: alive, INited
            character(len=99) :: tmpFilePath
            integer(4) :: size, rank, ierr

            call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

            call DestroyGas()

            Allocate(GasGlobal(CF%Ng))
            Allocate(TempGP(CF%Ng))

            do i=1,CF%Ng
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

            CF%Ns=Ns
            Allocate(SpecyGlobal(0:CF%Ns))

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

    End Module ModuleReactionOnePegasus
