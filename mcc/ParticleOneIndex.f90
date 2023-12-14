Module ModuleParticleOneIndex
    Use ModuleParticleOne
    Implicit none
    Type,extends(ParticleOne) :: ParticleOneIndex
        Integer(4) ::  Index
    contains
        procedure ::  IndexInit=>IndexInitializationParticleOneIndex
        procedure :: Copy=>CopyParticleOneIndex
        procedure :: Swap=>SwapParticleOneIndex
    EndType ParticleOneIndex

contains
    subroutine IndexInitializationParticleOneIndex(POI,Index)
        Class(ParticleOneIndex),intent(inout) :: POI
        Integer(4),intent(in) :: Index
        POI%Index=Index
    Endsubroutine IndexInitializationParticleOneIndex

    subroutine CopyParticleOneIndex(POD,POC)
        Class(ParticleOneIndex),intent(inout) :: POD
        Class(ParticleOne),intent(in) :: POC
        Select Type(POD)
        Type is(ParticleOneIndex)
            Select Type(POC)
            Type is(ParticleOneIndex)
                POD=POC
            ENdSelect
        ENdSelect
    Endsubroutine CopyParticleOneIndex

    Subroutine SwapParticleOneIndex(POD,POC)
        Class(ParticleOneIndex),intent(inout) :: POD
        Class(ParticleOne),intent(inout) :: POC
        Type(ParticleOneIndex) :: POT
        Select Type(POD)
        Type is(ParticleOneIndex)
            Select Type(POC)
            Type is(ParticleOneIndex)
                POT=POC
                POC=POD
                POD=POT
            ENdSelect
        ENdSelect
    Endsubroutine SwapParticleOneIndex
ENDModule ModuleParticleOneIndex
