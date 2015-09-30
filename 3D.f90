subroutine solve

use system
use const
use kai
use chainsdat
use molecules
use results
use kinsol
use bulk
use MPI
use ellipsoid
use ematrix
implicit none
external fcn
integer i, ix, iy, iz

!-----  varables de la resolucion -----------

real*8 x1(2*dimx*dimy*dimz),xg1(2*dimx*dimy*dimz)
real*8 f(2*dimx*dimy*dimz)
       
integer n

! Volumen fraction
real*8 xh(dimx, dimy, dimz)
real*8 psi(dimx, dimy, dimz) ! potencial

! MPI
integer tag, source
parameter(tag = 0)
integer err
integer ier_tosend
double  precision norma_tosend

! number of equations

n = dimx*dimy*dimz

! Initial guess

if((infile.eq.2).or.(infile.eq.-1)) then
  do i = 1, 2*n  
      xg1(i) = xflag(i)     
      x1(i) = xflag(i)
  enddo
endif

if(infile.eq.0) then
  do i=1,n
    xg1(i)=xsolbulk
    x1(i)=xsolbulk
  enddo
  do i=n+1, n*2
    xg1(i)=0.0d0
    x1(i)=0.0d0
  enddo
endif

!--------------------------------------------------------------
! Solve               
!--------------------------------------------------------------

! JEFE
if(rank.eq.0) then ! solo el jefe llama al solver
   iter = 0
   print*, 'solve: Enter solver ', 2*n, ' eqs'

   if(infile.ge.0) then
    call call_kinsol(x1, xg1, ier)
   endif
   if(infile.eq.-1) then
    call fkfun(x1, f, ier)
   endif
   flagsolver = 0
   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
endif
  
! Subordinados

if(rank.ne.0) then
  do
     flagsolver = 0
     source = 0
     CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
     if(flagsolver.eq.1) then
        call call_fkfun(x1) ! todavia no hay solucion => fkfun 
     endif ! flagsolver
     if(flagsolver.eq.0) exit ! Detiene el programa para este nodo
   enddo
endif

! Recupero el valor de ier y de la norma
! Asi los subordinados se enteran si el solver convergio o si hay que
! cambiar la   estrategia...
! Jefe

if (rank.eq.0) then
   norma_tosend = norma
   ier_tosend = ier ! distinto tipo de integer
   CALL MPI_BCAST(norma_tosend, 1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
   CALL MPI_BCAST(ier_tosend,1, MPI_INTEGER,0,MPI_COMM_WORLD,err)
endif

! Subordinados

if (rank.ne.0) then
   CALL MPI_BCAST(norma_tosend, 1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
   CALL MPI_BCAST(ier_tosend, 1, MPI_INTEGER,0,MPI_COMM_WORLD,err)
   norma = norma_tosend
   ier = ier_tosend
endif

! Recupera xh y psi (NO SON COMMON!)
do ix=1,dimx
   do iy=1,dimy
      do iz=1,dimz
       xh(ix,iy,iz)=x1(ix+dimx*(iy-1)+dimx*dimy*(iz-1))
       psi(ix,iy,iz)=x1(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+n)
      enddo
   enddo  
enddo

! Chequea si exploto... => Sistema anti-crash

if(infile.ne.5) then
  if((ier.lt.0).or.(.not.((norma.gt.0).or.(norma.lt.0))).or.(norma.gt.error)) then ! exploto...
    if(rank.eq.0)print*, 'solve: Error in solver: ', ier
    if(rank.eq.0)print*, 'solve: norma ', norma
    call MPI_FINALIZE(ierr) ! finaliza MPI
    stop
  endif
endif    

! No exploto, guardo xflag
do i = 1, 2*n
  xflag(i) = x1(i) ! xflag sirve como input para la proxima iteracion
enddo
infile = 2 ! no vuelve a leer infile

!----------------------------------------------------------
!  OUTPUT
!----------------------------------------------------------

!if(rank.eq.0) then ! solo el jefe escribe a disco....
  ! Guarda infile
!  write(filename,'(A4, I3.3, A4)')'out.', cccc, '.dat'
!  open(unit=45, file=filename)
!   do i = 1, 2*n
!    write(45, *)x1(i)
!   enddo
!  close(45)
! endif

end subroutine


