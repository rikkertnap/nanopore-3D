subroutine solver(xvec, xvecguess, ier)

    use system, only : eqs, dimx, dimy, dimz, method 

    implicit none

    real*8, intent(inout) :: xvec(eqs*dimx*dimy*dimz), xvecguess(eqs*dimx*dimy*dimz)
    integer*4, intent(inout) :: ier ! error flag type given  by kinsol

    ! local : put into system  
    real*8 :: accuracy, residual
    logical :: isSolution
    integer :: maxfkfunevals, neq
    
    method=1

    neq= eqs * dimx * dimy * dimz  
    maxfkfunevals = 1000 ! == defined in kinsol.f90 need to connected 
    accuracy = 1.0d-6    ! == defined fnormtol defined in kinsol.f90
   
    if(method==1) then

        call call_kinsol(xvec, xvecguess, ier)

    else if(method==2) then

        call anderson_min_loop(xvecguess, xvec, accuracy, residual, isSolution, maxfkfunevals, neq)

    else if(method==3) then
    
        call simple_min_loop(xvecguess, xvec, accuracy, residual, isSolution, maxfkfunevals, neq)

    else 
        print*,"Solver method incorrect"
        stop
    endif

end subroutine solver



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subrutina que llama a kinsol
     

subroutine call_fkfun(x1_old)
    use system
    use MPI

    integer i

    real*8 :: x1_old(eqs*dimx*dimy*dimz) 

    ! == local arguments 
    real*8 :: x1(eqs*dimx*dimy*dimz)     ! == this make a local copy call x1 independetly of the one defien din solve in 3D.f90
    real*8 :: f(eqs*dimx*dimy*dimz)

    ! MPI

    integer tag
    parameter(tag = 0)
    integer err

    x1 = 0.0
    do i = 1,eqs*dimx*dimy*dimz
        x1(i) = x1_old(i)
    enddo

    CALL MPI_BCAST(x1, eqs*dimx*dimy*dimz , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)

    call fkfun(x1,f, ier) ! todavia no hay solucion => fkfun 
end
