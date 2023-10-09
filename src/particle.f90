module particles

type Species
    character(8):: id             !species id
    real(8)     :: mass           !molecular mass
    real(8)     :: charge         !multiple of electron charge
end type Species

type OnePart
    integer(8)  :: id                  !particle ID
    integer(8)  :: ispecies            !particle species index
    real(8)     :: x(3),v(3),a(3)      !position,velocity,accelation
    real(8)     :: weight              !particle or cell weight, if weighting enabled
end type OnePart

type particles
type(OnePart),allocatable:: onepart(:) 
end type partices

end module particles