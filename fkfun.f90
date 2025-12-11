
subroutine fkfun(x,f,ier)

    use fcnpointer 

    implicit none

    ! input arguments 
    integer*4, intent(inout) :: ier
    real*8, intent(in) :: x(*)
    real*8, intent(inout) :: f(*)

    ! local arguments

    call fcnptr(x,f,ier)

    return

end subroutine fkfun


