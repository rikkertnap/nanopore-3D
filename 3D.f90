subroutine solve(flagcrash)

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
use mparameters_monomer
implicit none
external fcn
integer i, ix, iy, iz, ip
integer flagcrash

!-----  varables de la resolucion -----------

real*8 x1(eqs*dimx*dimy*dimz),xg1(eqs*dimx*dimy*dimz)
real*8 f(eqs*dimx*dimy*dimz)
       
integer ncells

! Volumen fraction
real*8 xh(dimx, dimy, dimz), xtotal(dimx,dimy,dimz,N_poorsol)
real*8 psi(dimx, dimy, dimz) ! potencial
real*8 temp

! MPI
integer tag, source
parameter(tag = 0)
integer err
integer ier_tosend
double  precision norma_tosend

! number of equations

ncells = dimx*dimy*dimz

! Initial guess

if((infile.eq.2).or.(infile.eq.-1).or.(infile.eq.3)) then
  do i = 1, eqs*ncells  
      xg1(i) = xflag(i)     
      x1(i) = xflag(i)
  enddo
endif

if(infile.eq.0) then
  do i=1,ncells
    xg1(i)=xsolbulk
    x1(i)=xsolbulk
  enddo

  do i = ncells+1,(N_poorsol+1)*ncells
    xg1(i)=0.0
    x1(i)=0.0
  enddo

if(electroflag.eq.1) then

  do i=(N_poorsol+1)*ncells+1, (N_poorsol+2)*ncells
    xg1(i)=0.0d0
    x1(i)=0.0d0
  enddo

endif ! electroflag
endif ! infile

!--------------------------------------------------------------
! Solve               
!--------------------------------------------------------------

! JEFE
if(rank.eq.0) then ! solo el jefe llama al solver
   iter = 0
   write(stdout,*) 'solve: Enter solver ', eqs*ncells, ' eqs'

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

       do ip=1, N_poorsol
          xtotal(ix,iy,iz,ip)=x1(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells)
       enddo

       if(electroflag.eq.1)psi(ix,iy,iz)=x1(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsol+1)*ncells)
      enddo
   enddo  
enddo

! Chequea si exploto... => Sistema anti-crash

if(infile.ne.-1) then
  if((ier.lt.0).or.(.not.((norma.gt.0).or.(norma.lt.0))).or.(norma.gt.error)) then ! exploto...
    if(rank.eq.0)write(stdout,*) 'solve: Error in solver: ', ier
    if(rank.eq.0)write(stdout,*) 'solve: norma ', norma
    flagcrash = 1
    return
  endif
endif    

! No exploto, guardo xflag
do i = 1, eqs*ncells
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

flagcrash = 0
return

end subroutine


