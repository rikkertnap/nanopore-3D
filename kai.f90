
!#####################################################################
!
! Este programa calcula los kai para poor-solvent en 3D (geometria 
! esferica) usando un metodo de MC
!
!#####################################################################


subroutine kais
use system
use const
use kai
use molecules
use MPI
use chainsdat
use transform
implicit none

real*8 l ! medio lseg, radio del segmento
 
integer MCsteps ! numero de steps de MC
integer ix, iy , iz, jx, jy, jz
real*8 x,y,z, radio
integer limit
real*8,allocatable :: matriz(:,:,:)
integer i
real*8 rands
real*8 suma
real*8 or
real*8 volume
real*8 xx(3), vv(3)
integer temp

Xulimit = int(cutoff/delta)+1
limit = Xulimit + 5 ! make it much larger, will check at the end.

ALLOCATE (matriz(-limit:limit, -limit:limit, -limit:limit)) ! matriz de kai
if(rank.eq.0)print*,'kais: Kai calculation'
suma = 0.0
matriz = 0.0
MCsteps = 200*Xulimit
volume = dfloat(MCsteps)**3

!if(rank.eq.0)print*, 'kais: CORREGIR MCSTEPS!!!!'

l = lseg 

do jx = 1, MCsteps
do jy = 1, MCsteps
do jz = 1, MCsteps

x = 2.0*cutoff*(((float(jx)-0.5)/float(MCsteps))-0.5) ! number between -cutoff and +cutoff
y = 2.0*cutoff*(((float(jy)-0.5)/float(MCsteps))-0.5)
z = 2.0*cutoff*(((float(jz)-0.5)/float(MCsteps))-0.5)

 radio = sqrt(x**2 + y**2 + z**2) ! espacio real

 if(radio.gt.cutoff) cycle ! No esta dentro de la esfera del cut-off   
 if(radio.lt.l) cycle ! esta dentro de la esfera del segmento

 ! celda 

 vv(1) = x ! real space
 vv(2) = y
 vv(3) = z

 xx = MATMUL(MAT,vv) ! xx in transformed space

 ix = int(anint(xx(1)/delta))   ! espacio de la grilla
 iy = int(anint(xx(2)/delta))
 iz = int(anint(xx(3)/delta))

 matriz(ix, iy, iz) = matriz(ix, iy, iz) + (l/radio)**6

enddo
enddo
enddo

! Find minimum Xulimit

Xulimit = 0
do ix = -limit, limit
do iy = -limit, limit
do iz = -limit, limit
if(matriz(ix,iy,iz).ne.0.0) then
  if(Xulimit.lt.abs(ix))Xulimit=abs(ix)
  if(Xulimit.lt.abs(iy))Xulimit=abs(iy)
  if(Xulimit.lt.abs(iz))Xulimit=abs(iz)
endif
enddo
enddo
enddo

if(Xulimit.eq.limit) then
print*, 'kais: error Xulimit = limit'
call MPI_FINALIZE(ierr)
stop
endif

if(rank.eq.0)print*,'kais: New Xulimit', Xulimit

ALLOCATE (Xu(-Xulimit:Xulimit,-Xulimit:Xulimit,-Xulimit:Xulimit))
Xu = 0.0
sumXu = 0.0
do ix = -Xulimit, Xulimit
do iy = -Xulimit, Xulimit
do iz = -Xulimit, Xulimit
 Xu(ix, iy, iz) = matriz(ix, iy, iz)/volume*((2.0*cutoff)**3)
 suma = suma +  matriz(ix, iy, iz)/volume*((2.0*cutoff)**3)
 sumXu = sumXu + Xu(ix, iy, iz)
enddo
enddo
enddo

if(rank.eq.0)print*, 'kais: Sum Xulimit', suma

suma = 0.0
do ix = -limit, limit
do iy = -limit, limit
do iz = -limit, limit
 suma = suma +  matriz(ix, iy, iz)/volume*((2.0*cutoff)**3)
enddo
enddo
enddo

if(rank.eq.0)print*, 'kais: Total Sum', suma

do ix = -Xulimit, Xulimit
do iy = -Xulimit, Xulimit
do iz = -Xulimit, Xulimit
if(rank.eq.0)write(999,*)ix,iy,iz,Xu(ix, iy, iz)
enddo
enddo
enddo

 
end



