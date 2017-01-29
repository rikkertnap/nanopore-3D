!     The routine kpsol is the preconditioner solve routine. It must have
!     that specific name be used in order that the c code can find and link
!     to it.  The argument list must also be as illustrated below:

subroutine fkpsol(udata, uscale, fdata, fscale, vv, ftem, ier)
use system
use mkinsol
implicit none

integer ier
integer*8 neq, i
double precision udata(*), uscale(*), fdata(*), fscale(*)
double precision vv(*), ftem(*)

common /psize/ neq

do  i = 1, neq
   vv(i) = vv(i) * pp(i)
enddo
ier = 0

return
end

!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!* * *
!c     The routine kpreco is the preconditioner setup routine. It must have
!c     that specific name be used in order that the c code can find and link
!c     to it.  The argument list must also be as illustrated below:

subroutine fkpset(udata, uscale, fdata, fscale,vtemp1,vtemp2, ier)
use system
use mkinsol
use mparameters_monomer
implicit none

integer ier
integer*8 neq, i
double precision udata(*), uscale(*), fdata(*), fscale(*)
double precision vtemp1(*), vtemp2(*)
integer ncells

common /psize/ neq

ncells =  dimx*dimy*dimz

do i = 1, ncells
   pp(i) = 0.1 / (1.0+exp(1.0-udata(i)))
! pp(i) = 1
enddo

do i = ncells+1, (N_poorsol+1)*ncells
!   pp(i) = 0.1 / (1.0+exp(1.0-udata(i)))
pp(i) = 0.1 / (1.0+exp(1.0-udata(i)))
enddo

if(electroflag.eq.1) then
do i = ncells*(N_poorsol+1), ncells*(N_poorsol+2)
   pp(i) = 1.0
enddo
endif

   ier = 0
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subrutina que llama a kinsol
     

subroutine call_fkfun(x1_old)
use system
use MPI

integer i

real*8 x1_old(eqs*dimx*dimy*dimz)
real*8 x1(eqs*dimx*dimy*dimz)
real*8 f(eqs*dimx*dimy*dimz)

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

subroutine call_kinsol(x1_old, xg1_old, ier)
use system
use const
use mparameters_monomer
implicit none
integer i
real*8 x1(eqs*dimx*dimy*dimz), xg1(eqs*dimx*dimy*dimz)
real*8 x1_old(eqs*dimx*dimy*dimz), xg1_old(eqs*dimx*dimy*dimz)
integer*8 iout(15) ! Kinsol additional output information
real*8 rout(2) ! Kinsol additional out information
integer*8 msbpre
real*8 fnormtol, scsteptol
real*8 scale(eqs*dimx*dimy*dimz)
real*8 constr(eqs*dimx*dimy*dimz)
integer*4  globalstrat, maxl, maxlrst
integer*4 ier ! Kinsol error flag
integer*8 neq ! Kinsol number of equations
integer*4 max_niter
common /psize/ neq ! Kinsol
integer ierr
integer ncells


! INICIA KINSOL
ncells = dimx*dimy*dimz
neq = eqs*dimx*dimy*dimz
msbpre  = 10 ! maximum number of iterations without prec. setup (?)
fnormtol = 1.0d-6 ! Function-norm stopping tolerance
scsteptol = 1.0d-6 ! Function-norm stopping tolerance

maxl = 500 ! maximum Krylov subspace dimesion (?!?!?!) ! Esto se usa para el preconditioner
maxlrst = 5 ! maximum number of restarts
max_niter = 1000
globalstrat = 0

call fnvinits(3, neq, ier) ! fnvinits inits NVECTOR module
if (ier .ne. 0) then       ! 3 for Kinsol, neq ecuantion number, ier error flag (0 is OK)
  write(stdout,*) 'call_kinsol: SUNDIALS_ERROR: FNVINITS returned IER = ', ier
  call MPI_FINALIZE(ierr) ! finaliza MPI
  stop
endif

call fkinmalloc(iout, rout, ier)    ! Allocates memory and output additional information
if (ier .ne. 0) then
   write(stdout,*) 'call_kinsol: SUNDIALS_ERROR: FKINMALLOC returned IER = ', ier
   call MPI_FINALIZE(ierr) ! finaliza MPI
   stop
 endif

call fkinsetiin('MAX_SETUPS', msbpre, ier)  ! Additional input information
call fkinsetrin('FNORM_TOL', fnormtol, ier)
call fkinsetrin('SSTEP_TOL', scsteptol, ier)
call fkinsetiin('MAX_NITER', max_niter, ier)

do i = 1, ncells  !constraint vector
   constr(i) = 2.0 ! xh > 0
enddo

do i = ncells+1, (N_poorsol+1)*ncells
   constr(i) = 1.0 ! xtotal >= 0
enddo

if(electroflag.eq.1) then
do i = ncells*(N_poorsol+1), ncells*(N_poorsol+2)  !constraint vector
   constr(i) = 0.0 ! no contraint for psi
enddo
endif

call fkinsetvin('CONSTR_VEC', constr, ier) ! constraint vector
! CALL FKINSPTFQMR (MAXL, IER)
call fkinspgmr(maxl, maxlrst, ier) !  Scale Preconditioned GMRES solution of linear system (???)
!call fkinspbcg(maxl, ier) !  Scale Preconditioned BCG

if (ier .ne. 0) then
  write(stdout,*) 'call_kinsol: SUNDIALS_ERROR: FKINSPGMR returned IER = ', ier
  call fkinfree ! libera memoria
  call MPI_FINALIZE(ierr) ! finaliza MPI
  stop
endif
call fkinspilssetprec(1, ier) ! preconditiones

do i = 1, neq ! scaling vector
  scale(i) = 1.0
enddo

do i = 1, neq ! Initial guess
      x1(i) = x1_old(i)
      xg1(i) = x1(i)  
enddo


call fkinsol(x1, globalstrat, scale, scale, ier)         ! Llama a kinsol
print*, 'ier', ier

if (ier .lt. 0) then
      write(stdout,*) 'call_kinsol: SUNDIALS_ERROR: FKINSOL returned IER = ', ier
      write(stdout,*) 'call_kinsol: Linear Solver returned IER = ', iout(9)
!      call fkinfree
endif

do i = 1, neq ! output
  x1_old(i) = x1(i)
  xg1_old(i) = x1(i)
enddo

return
end



