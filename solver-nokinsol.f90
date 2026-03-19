subroutine solver(xvec, xvecguess, ier)

    use system, only : neqs, methodflag 
    use kinsol, only : neq 
    use anderson

    implicit none

    real*8, intent(inout) :: xvec(neqs), xvecguess(neqs)
    integer*4, intent(inout) :: ier ! error flag type given  by kinsol

    ! local : put into system  
    real*8 :: accuracy, residual
    logical :: isSolution
    integer :: maxfkfunevals
    
    print*,"neq=",neq,"neqs=",neqs

    maxfkfunevals = 1000 ! == defined in kinsol.f90 need to connected 
    accuracy = 1.0d-6    ! == defined fnormtol defined in kinsol.f90
    ier = 0              ! == retrun value solver  

    if(methodflag==1) then

        ! call call_kinsol(xvec, xvecguess, ier)

    else if(methodflag==2) then

        call anderson_min_loop(xvecguess, xvec, accuracy, residual, isSolution, maxfkfunevals, neq)
        
    
    else if(methodflag==3) then
    
        call simple_min_loop(xvecguess, xvec, accuracy, residual, isSolution, maxfkfunevals, neq)

    else 
        print*,"Solver method incorrect"
        stop
    endif

    if(isSolution) ier=1 
      

end subroutine solver



!  == Subroutine that calls Kinsol   

subroutine call_fkfun(x)

    use system, only : neqs
    use kinsol, only : ier
    use MPI ! need for mpi definitions 

    implicit none 

    real*8 :: x(neqs)

    ! == local arguments 
    real*8 :: f(neqs)
    integer :: err  ! == mpi 
    
    
    CALL MPI_BCAST(x, neqs , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)

    call fkfun(x, f, ier) 

end subroutine call_fkfun

