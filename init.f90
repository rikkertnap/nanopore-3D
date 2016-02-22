
subroutine initmpi
use MPI
use chainsdat
implicit none

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

end subroutine

subroutine initconst
use const
use molecules
use ellipsoid
implicit none
pi = acos(-1.0)
seed = 938121
lb = 0.714 ! bjerrum lenght in nm
zpos = 1.0
zneg = -1.0
vsol = vsol0
vsalt=((4.0/3.0)*pi*(0.2)**3)/vsol  ! volume salt in units of vsol 0.2=radius salt  
vpol= vpol0/vsol ! ((4.0/3.0)*pi*(0.2)**3)/vsol  ! volume polymer segment in units of vsol 
constq=delta*delta*4.0*pi*lb/vsol   ! multiplicative factor in poisson eq  
pKw = 14
error = 1e-4 ! para comparar con la norma...
errel=1d-6
itmax=200

if(electroflag.eq.0)eqs=1
if(electroflag.eq.1)eqs=2

end subroutine

subroutine initall
use molecules
use const
use bulk
use MPI
use ellipsoid
use chainsdat
use inputtemp
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open common files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

if(rank.eq.0) then
       open(unit=301, file='F_tot.dat', access='APPEND')
       open(unit=302, file='F_mixs.dat',  access='APPEND')
       open(unit=303, file='F_mixpos.dat',  access='APPEND')
       open(unit=304, file='F_mixneg.dat',  access='APPEND')
       open(unit=305, file='F_mixH.dat',  access='APPEND')
       open(unit=306, file='F_mixOH.dat',  access='APPEND')
       open(unit=307, file='F_conf.dat',  access='APPEND')
       open(unit=3071, file='F_gauche.dat',  access='APPEND')
       open(unit=308, file='F_eq.dat',  access='APPEND')
       open(unit=309, file='F_vdW.dat',  access='APPEND')
       open(unit=410, file='F_eps.dat',  access='APPEND')
       open(unit=311, file='F_electro.dat',  access='APPEND')
       open(unit=312, file='F_tot2.dat',  access='APPEND')
       open(unit=314, file='F_mixpos2.dat',  access='APPEND')
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Input-dependent variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!

constqE = vpol/(2.0d0*constq)
dielW = 78.54
dielPr = dielP/dielW
dielSr = dielS/dielW

Ka=10**(-pKa)
cHplus = 10**(-pHbulk)    ! concentration H+ in bulk
xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
pOHbulk= pKw -pHbulk
cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol  
xsalt=(csalt*Na/(1.0d24))*(vsalt*vsol)   ! volume fraction salt,csalt in mol/l 
if(pHbulk.le.7) then  ! pH<= 7
  xposbulk=xsalt/zpos
  xnegbulk=-xsalt/zneg+(xHplusbulk -xOHminbulk) *vsalt ! NaCl+ HCl  
else                  ! pH >7 
  xposbulk=xsalt/zpos +(xOHminbulk -xHplusbulk) *vsalt ! NaCl+ NaOH   
  xnegbulk= -xsalt/zneg 
endif

xsolbulk=1.0 -xHplusbulk -xOHminbulk -xnegbulk -xposbulk 
K0 = (Ka*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant 
expmupos = xposbulk /xsolbulk**vsalt
expmuneg = xnegbulk /xsolbulk**vsalt
expmuHplus = xHplusbulk /xsolbulk   ! vsol = vHplus 
expmuOHmin = xOHminbulk /xsolbulk   ! vsol = vOHmin 

end subroutine

subroutine endall
use MPI
implicit none

!!!!!!!!!!!!!!!!!!!!!!
! Close common files
!!!!!!!!!!!!!!!!!!!!!!

close(301)
close(302)
close(303)
close(304)
close(305)
close(306)
close(307)
close(3071)
close(308)
close(309)
close(310)
close(311)
close(312)
close(313)

call MPI_FINALIZE(ierr) ! finaliza MPI    
stop

end subroutine



subroutine savedata(cccc)
use system
use results
use const
use molecules
use chainsdat
use kai
use ematrix
use fields_fkfun
use MPI
use kinsol
use kaist
implicit none
integer cccc
character*20 filename
character*5  title
real*8 temp(dimx,dimy,dimz)
real*8 sumpol
integer ix,iy,iz
!----------------------------------------------------------
!  OUTPUT
!----------------------------------------------------------

if(rank.eq.0) then ! solo el jefe escribe a disco....
  ! Guarda infile
!  write(filename,'(A4, I3.3, A4)')'out.', cccc, '.dat'
!  open(unit=45, file=filename)
!   do i = 1, 2*n
!    write(45, *)x1(i)
!   enddo
!  close(45)

!!!!!!!!!!!!!!!!!!! Guarda archivos !!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Polimero

  temp(:,:,:) = avpol(:,:,:)*(1.0 - volprot(:,:,:))

  title = 'avpol'
  call savetodisk(temp, title, cccc)

! Solvente
  temp(:,:,:) = xh(:,:,:)*(1.0 - volprot(:,:,:))

  title = 'avsol'
  call savetodisk(temp, title, cccc)
! Cationes
!  title = 'avpos'
!  call savetodisk(xpos, title, cccc)
! Aniones
!  title = 'avneg'
!  call savetodisk(xneg, title, cccc)
! H+
!  title = 'avHpl'
!  call savetodisk(xHplus, title, cccc)
! OH-
!  title = 'avOHm'
!  call savetodisk(xOHmin, title, cccc)
! fdis
!  title = 'frdis'
!  call savetodisk(fdis, title, cccc)
! Potencial electrostatico

  temp(1:dimx,1:dimy, 1:dimz) = psi(1:dimx,1:dimy, 1:dimz)

  title = 'poten'
  call savetodisk(temp, title, cccc)
! Particle
  title = 'avpar'
  call savetodisk(volprot, title, cccc)

! save volprot for supercell

if(rank.eq.0) then
open (unit=8, file='out.par', form='unformatted')
do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
  xpar(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) = volprot(ix,iy,iz)
  enddo
 enddo
enddo
write(8)xpar
close(8)
endif

! system

  write(filename,'(A7, I3.3, A4)')'system.', cccc, '.dat'
  open (unit=310, file=filename)
  write(310,*)'st          = ',st ! residual size of iteration vector
  write(310,*)'fnorm       = ',norma ! residual size of iteration vector
  write(310,*)'length seg  = ',0.35 ! value see subroutine cadenas
  write(310,*)'delta       = ',delta
  write(310,*)'vsol        = ',vsol
  write(310,*)'vsalt       = ',vsalt*vsol
  write(310,*)'vpol       = ',vpol*vsol
  write(310,*)'pKw         = ',pKw
  write(310,*)'zpos        = ',zpos
  write(310,*)'zneg        = ',zneg
  write(310,*)'long        = ',long
  write(310,*)'iterations  = ',iter
  write(310,*)'sigma cad/nm2 = ',ncha/(dimx*dimy*delta*delta)
  write(310,*)'kai =          ', Xu
  write(310,*)'GIT version = ', _VERSION

  sumpol = 0.0
  do ix = 1, dimx
  do iy = 1, dimy
  do iz = 1, dimz
  sumpol = sumpol + avpol(ix,iy,iz)*(delta**3)*(1.0-volprot(ix,iy,iz))/vpol/vsol
  enddo
  enddo
  enddo

  write(310,*)'Number of segments =          ', sumpol
  close(310)

endif

end subroutine

subroutine store2disk(counter) ! saves state to disk
use ellipsoid
use kinsol
use montecarlo
use ematrix
use results
use MPI
use const
implicit none
integer counter
character*20 filename

if(rank.eq.0) then
open (unit=8, file='out.out', form='unformatted')
write(8)counter
write(8)free_energy
write(8)xflag
close(8)
endif

if(rank.eq.0) then
write(filename,'(A4, I3.3, A4)')'out.', counter, '.dat'
open(unit=8, file=filename, form='unformatted')
write(8)counter
write(8)free_energy
write(8)xflag
close(8)
endif
end subroutine

subroutine retrivefromdisk(counter) ! saves state to disk
use ellipsoid
use kinsol
use montecarlo
use ematrix
use results
use const
implicit none
integer counter

open (unit=8, file='in.in', form='unformatted')
read(8)counter
read(8)free_energy
read(8)xflag
close(8)
end subroutine


