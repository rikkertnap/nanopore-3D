!  Anderson.c module 
!  routine to solve SCF equations by means of Anderson mixing
!  an update of the iterations vector is a linear combination 
!  of  previous  iteration vectors  
!  see V. Eyert J. Comp. Phys.1996,124, 271-285 */

module anderson 

    !use precision_definition
    !use mathconst
    !use myutils

    implicit none

    real*8, parameter ::  BETA_A = 1.0 ! 0.9         ! mixing parameter for anderson mixing 
    real*8, parameter ::  BETA_S = -0.2  !  -0.2    ! mixing parameter for simple mixing 
    real*8, parameter ::  TOL_DELTA = 0.0001 ! 0.01     ! treshold for Anderson to start 
    integer, parameter  ::  NUMBER_SOL= 6       ! number of previous solutions used in Anderson mixing */

    integer ::  step   ! counter 

contains

! computes L2norm of vector f 

    function l2norm(f,n)result(norm)
      
        implicit none

        integer, intent(in)  :: n 
        real*8, intent(in) :: f(n)
        real*8             :: norm ! output

        integer :: i ! dummy index

        norm=0.0
        do i=1,n
            norm = norm + f(i)**2
        enddo
        norm=dsqrt(norm)
        
    end function l2norm


! shift columns of matrix  
subroutine cycle_matrix_ptr(Mat,MM,NN)
    
    real*8, intent(inout), dimension(:,:) :: Mat
    integer, intent(in) :: MM, NN
    
    
    integer  ::i
    real*8:: tmp_ptr(NN)
    
    tmp_ptr = Mat(1,:) 
    do i = 1, MM-1
        Mat(i,:)= Mat(i+1,:)
    enddo    
    Mat(MM,:) = tmp_ptr

end subroutine 

subroutine  cycle_matrix(Mat, MM, NN)

    real*8, intent(inout), dimension(:,:) :: Mat
    integer, intent(in) :: MM, NN

    integer :: i,j
	
    do i = 1, MM-1
        do j = 1, NN
            Mat(i,j) = Mat(i+1,j)
        enddo
    enddo        
    
end subroutine 

! matrix of  residual vectors 
! U[n][m] = < F[L] - F[L-n]|F[L] - F[L-m]>  

subroutine deviation_matrix(U,d,MM,NN)
 

    real*8, intent(inout), dimension(:,:) :: U
    real*8, intent(in), dimension(:,:) :: d
    integer, intent(in) :: MM, NN
    
    integer :: n, m, i
    real*8:: tmp
    
    do n = 1, MM
        do  m = 1, MM
            tmp = 0.0
            do i = 1, NN
                tmp = tmp+ (d(MM+1,i) - d(MM+1-n,i))*(d(MM+1,i) - d(MM+1-m,i))
            enddo
            U(n,m) = tmp
        enddo    
    enddo
    
end subroutine

! vector  of  residual vectors 
!  V[n] =  < F[L] - F[L-n]|F[L]> 

subroutine deviation_vector(V,d,MM,NN)

    real*8, intent(inout), dimension(:) :: V
    real*8, intent(in), dimension(:,:) :: d
    integer, intent(in) :: MM, NN

    integer :: n, i
    real*8::  tmp
    
    do n = 1, MM
        tmp = 0.0
        do i = 1, NN
            tmp = tmp + (d(MM+1,i) - d(MM+1-n,i))*(d(MM +1,i)) 
        enddo
        V(n) = tmp
    enddo
    
end subroutine 

! calculated the mixing coefficient by solving  U theta = V 
! uses LAPACK routines 

subroutine  mixing_coefficients(theta, U, V,  M)

    real*8, dimension(:) :: theta
    real*8, dimension(:,:) :: U 
    real*8, dimension(:) :: V
    integer, intent(in) ::  M
    
    integer ::  NRHS, N, INFO ! long int
    integer,  dimension(:), allocatable ::  IPIV
    character :: TRANS
    
    N = M
    TRANS ='N'
    NRHS = 1
    allocate(IPIV(M))   ! pivot matrix 
    

   call dgetrf(N,N, U, N, IPIV, INFO) 

    if(INFO/=0) then
        print*,"LU decomposition of UF failed ."
        print*,"INFO=",INFO
        stop
    endif 
    
    call dgetrs(TRANS,N,NRHS,U,N,IPIV,V,N,INFO )
   
    if(INFO/=0)then
        print*,"Solving linear equation  U theta = V  failed .\n"
        print*,"INFO=",INFO
        stop 
    endif
    ! on succesfull solving V contains solution theta 
   
    theta = V
    deallocate(IPIV)

end subroutine

! Iterative solver for scf eq 
! using  Anderson method 
subroutine anderson_min_loop(xguess, x, TOL, fnorm, isSolution, MAX_INT, N)

   ! use vectornorm
   ! use fcnpointer
    
    real*8, intent(inout), dimension(:) :: xguess
    real*8, intent(inout), dimension(:) :: x
    real*8, intent(in) :: TOL
    real*8, intent(out) :: fnorm 
    logical,  intent(out) :: isSolution
    integer, intent(in) :: MAX_INT
    integer(8), intent(in) :: N

    integer ::  i, j, r,  NN, M
    logical :: conv
    real*8 :: delta
    real*8, allocatable, dimension(:)  :: fvec, Ftmp, xoutM
    real*8::  tmpx, tmpF, xbar, Fbar
    real*8, allocatable, dimension(:,:) ::   xin, xout, F, U
    real*8, allocatable, dimension(:)  :: V, theta
    character(len=256) :: istr, rstr, text

    integer*4 :: ier
    
    NN = int(N,kind(NN))
    M = NUMBER_SOL         ! number of previous solution 

    ! memory allocation 
            
    allocate(fvec(NN))     ! tmp F(x) 
    allocate(xin(M+1,NN) ) 
    allocate(xout(M+1,NN)) 
    allocate(xoutM(NN)) 
    allocate(F(M+1,NN))    !  F(x) 
    allocate(Ftmp(NN))
    allocate(U(M,M))
    allocate(V(M))
    allocate(theta(M))
    
    !  set fields and initial  starting solution  
    do  i = 1, NN
        x(i) =  xguess(i)
    enddo
    
    ! Solver for F(X) = 0 
    
    conv = .false.                     
    step = 1   
    j = 1         
    
    ! the first M iterations must be simple mixing   

    do while( (j <= M+1) .and. (.not.conv .and. (step < MAX_INT) ))
    
        !call fcnptr(x,fvec,N)
        call fkfun(x,fvec,ier)

        do i = 1, NN
            xin(j,i) = x(i)
            F(j,i) = fvec(i)
            xout(j,i) = fvec(i) + x(i) 
            x(i) = x(i) + BETA_S * fvec(i) ! update 
        enddo

        fnorm = L2norm(fvec,NN)
        conv = (fnorm <= TOL)  ! 1=true if found solution 

        print*,"step = ",step,"norm = ",fnorm
        step = step + 1
        j = j+1

    enddo
    
    do while( (.not.(conv) .and. (step < MAX_INT) ))
    
        !call fcnptr(x,fvec,N)
        call fkfun(x,fvec,ier)
        
        call cycle_matrix(xin,M+1,NN)
        call cycle_matrix(xout,M+1,NN)
        call cycle_matrix(F,M+1,NN)

        do i = 1, NN
            xin(M+1,i) = x(i)
            F(M+1,i) = fvec(i)
            xout(M+1,i) = fvec(i) +x(i) 
            xoutM(i) = xout(M+1,i)
        enddo   
     
        fnorm = L2norm(fvec,NN)
        delta = fnorm/L2norm(xoutM,NN)
        conv = (fnorm <= TOL)       ! 1=true if found solution
        
        print*,"step = ",step,"norm = ",fnorm,"delta = ",delta
        
        if( delta < TOL)  then ! Anderson mixing
        
            call deviation_matrix(U, F, M, NN)
            call deviation_vector(V, F, M, NN)

            call mixing_coefficients(theta, U, V, M)
        
            do i =1, NN
                tmpx = 0.0
                tmpF = 0.0
                do r = 1, M ! linear combination of M previous vectors
                    tmpx = tmpx + theta(r) * (xin(M+1-r,i) - xin(M+1,i))
                    tmpF = tmpF + theta(r) * (F(M+1-r,i) - F(M+1,i))
                enddo
                xbar = xin(M+1,i) + tmpx
                Fbar = F(M+1,i)   + tmpF
                x(i) =  xbar  + BETA_A *  Fbar  ! update  
            enddo
            
        else ! simple mixing
            do i =1, NN
                x(i) = x(i) + BETA_S * fvec(i) ! update
            enddo
        endif
        
        step = step + 1
        
    enddo ! end second while loop i.e. end of iteration */ 
    
    isSolution = conv

    if(conv) then
        write(rstr,'(E25.16)')fnorm
        text="final L2 norm of the residuals = = "//trim(adjustl(rstr))
        !call print_to_log(LogUnit,text)
        write(istr,'(I8)')step
        text="number of iterations  = "//trim(adjustl(istr))
        !call print_to_log(LogUnit,text)
    endif    

    if (step>=MAX_INT) then
        text="program exceeded maximuum number of interations, iteration aborted"
        !call print_to_log(LogUnit,text)
        !call write_last_fcn_eval(x)
    endif     
    
    ! free memory 
    deallocate(fvec)
    deallocate(xin)  
    deallocate(xout)   
    deallocate(xoutM) 
    deallocate(F)    !  F(x) 
    deallocate(Ftmp)
    deallocate(U)
    deallocate(V)
    deallocate(theta)
    
end subroutine 

! Iterative solver for scf eq
! simple mixing method 
subroutine  simple_min_loop(xguess, x, TOL, fnorm, isSolution, MAX_INT, N)

    !use vectornorm
    !use fcnpointer
    
    real*8, intent(inout), dimension(:) :: xguess
    real*8, intent(inout), dimension(:) :: x
    real*8, intent(in) :: TOL
    real*8, intent(out) :: fnorm
    logical, intent(out) :: isSolution
    integer, intent(in) :: MAX_INT
    integer(8), intent(in) :: N 

    integer :: i
    integer :: NN
    logical :: conv
    real*8, allocatable, dimension(:)  :: fvec
    character(len=256) :: rstr, istr, text
    integer*4 :: ier
        
    NN = int(N,kind(NN))
          
    !memory allocation 
    allocate(fvec(NN))!     ! tmp F(x) 
         
    ! set fields and initial  starting solution  
    do i = 1,  NN
        x(i) = xguess(i)
    enddo
        
    ! Solver for F(X) = 0  
    
    conv = .false.                      
    step = 1  
    
    !    simple mixing   

    do while( (.not.(conv) .and. (step < MAX_INT)))
        ! call fcnptr(x,fvec,N) 
        call fkfun(x,fvec,ier)
        do i = 1, NN
            x(i) = x(i) + BETA_S * fvec(i) ! update 
        enddo
   
        fnorm = L2norm(fvec,NN)
        conv = (fnorm <= TOL)  ! 1=true if found solution */

!        print*,"step = ",step,"L2norm = ",fnorm

        step = step+1
    enddo    
    
    isSolution = conv

    if(conv) then
        write(rstr,'(E25.16)')fnorm
        text="final L2 norm of the residuals = = "//trim(adjustl(rstr))
        !call print_to_log(LogUnit,text)
        write(istr,'(I8)')step
        text="number of iterations  = "//trim(adjustl(istr))
        !call print_to_log(LogUnit,text)
    endif    
    
    if (step>=MAX_INT) then
        text="program exceeded maximuum number of interations, iteration aborted"
       ! call print_to_log(LogUnit,text)
        call write_last_fcn_eval(x)
    endif    

    ! free memory 
    deallocate(fvec)
  




  
end subroutine  simple_min_loop


subroutine  write_last_fcn_eval(x)

    real*8, intent(in), dimension(:) :: x

    integer :: un_out 
    character(len=5) :: outfilename
    integer :: ios, i


    outfilename = "x.out"
    open(newunit=un_out,file=outfilename, iostat=ios, action="write")
    
    if(ios > 0 ) then
        print*, 'Error opening file : iostat =', ios
    endif

    do i=1,size(x)
        write(un_out,*)x(i)
    enddo

    close(un_out)

end subroutine write_last_fcn_eval

end module 













