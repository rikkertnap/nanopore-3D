!## Montecarlo - Molecular theory for the adsorption of a particle on a brush

use system
use MPI
use ellipsoid
use kinsol
use const
use montecarlo
use ematrix
use old
use kaist

implicit none
integer counter, counterr
integer MCpoints
integer saveevery 
real*8 maxmove
real*8 maxrot
integer moves,rots
real*8 rv(3)
real*8 temp
real*8 theta
real*8, external :: rands
logical flag
character*10 filename
integer j, i


!!!!!!!!!!!! global parameters ... !!!!!!
saveevery = 1

counter = 0
counterr = 1

call initmpi
if(rank.eq.0)print*, 'Program Crystal'
if(rank.eq.0)print*, 'GIT Version: ', _VERSION
if(rank.eq.0)print*, 'MPI OK'
call readinput

call initconst
call inittransf ! Create transformation matrixes
call initall
call allocation

!!! General files

do j = 1, NNN
write(filename,'(A3,I3.3, A4)')'pos',j,'.dat'
open(file=filename, unit=5000+j)
write(filename,'(A3,I3.3, A4)')'orn',j,'.dat'
open(file=filename, unit=6000+j)
write(filename,'(A3,I3.3, A4)')'rot',j,'.dat'
open(file=filename, unit=7000+j)
enddo

open(file='free_energy.dat', unit=9000)
open(file='acceptance.dat', unit=9001)

call kais
if(rank.eq.0)print*, 'Kai OK'

call update_matrix(flag) ! updates 'the matrix'

  if(flag.eqv..true.) then
    print*, 'Initial position of particle does not fit in z'
    print*, 'or particles collide'
    stop
  else
    if(rank.eq.0)print*, 'Particle OK'
  endif

call  graftpoints
if(rank.eq.0)print*, 'Graftpoints OK'

call creador ! Genera cadenas
if(rank.eq.0)print*, 'Creador OK'

if(infile.ne.0) then
   call retrivefromdisk(counter)
   if(rank.eq.0)print*, 'Load input from file'
   if(rank.eq.0)print*, 'Free energy', free_energy
   infile = 2
   call update_matrix(flag)
   if(flag.eqv..true.) then
    print*, 'Initial position of particle does not fit in z'
    print*, 'or particles collide'
    stop
   endif

endif

do i = 1, nst
 st = sts(i)
 if(rank.eq.0)print*,'Switch to st = ', st
 call solve
 counterr = counter + i
 call Free_Energy_Calc(counterr)
 if(rank.eq.0)print*, 'Free energy after solving', free_energy
 call savedata(counterr)
 if(rank.eq.0)print*, 'Save OK'
 call store2disk(counterr)
enddo

call endall
end


