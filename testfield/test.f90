program main 
    implicit none
    integer(4)::c(3,3,3),b(3,3,3)
    c=1
    b=2
    c(1,1:3,1:3)=b(1:3,1,1:3)
    write(*,*)c
end program main 