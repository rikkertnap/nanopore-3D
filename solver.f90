subroutine solver(xvec, xvecguess, ier)

    use system, only : neqs, methodflag
    use anderson
    use kinsol, only : neq 


    implicit none

    real*8, intent(inout) :: xvec(neqs), xvecguess(neqs)
    integer*4, intent(inout) :: ier ! error flag type given  by kinsol

    ! local : put into system  
    real*8 :: accuracy, residual
    logical :: isSolution
    integer :: maxfkfunevals
  
    
    if(methodflag/=1) then 
        maxfkfunevals = 1000 ! == defined in kinsol.f90 
        accuracy = 1.0d-6    ! == defined fnormtol defined in kinsol.f90
    endif

    if(methodflag==1) then

        call call_kinsol(xvec, xvecguess, ier)

    else if(methodflag==2) then

        call anderson_min_loop(xvecguess, xvec, accuracy, residual, isSolution, maxfkfunevals, neq)

    else if(methodflag==3) then
    
        call simple_min_loop(xvecguess, xvec, accuracy, residual, isSolution, maxfkfunevals, neq)

    else 
        print*,"Solver method incorrect"
        stop
    endif

    if(methodflag/=1) then 
        if(iSsolution) then 
            ier = 1
        else
            ier = 0
        endif    
    endif

end subroutine solver



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subrutina que llama a kinsol
     

!subroutine call_fkfun(x1_old)
!    use system
!    use MPI
!
!    integer i

!    real*8 :: x1_old(neqs) 

    ! == local arguments 
!    real*8 :: x1(neqs)     ! == this make a local copy call x1 independetly of the one defien din solve in 3D.f90
!    real*8 :: f(neqs)

    ! MPI

!    integer tag
!    parameter(tag = 0)
!    integer err

!    x1 = 0.0
!    do i = 1,neqs
!        x1(i) = x1_old(i)
!    enddo

!    CALL MPI_BCAST(x1, neqs , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)

!    call fkfun(x1,f, ier) ! todavia no hay solucion => fkfun 
!end
