!## Montecarlo - Molecular theory for the adsorption of a particle on a brush

use system
use MPI
use ellipsoid
use kinsol
use const
use montecarlo
use ematrix
use old

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
integer j


!!!!!!!!!!!! global parameters ... !!!!!!
verbose = 5
PBC = 1 ! Peridic boundary conditions 
        ! 0 = only x and y
        ! 1 = xyz

saveevery = 1

counter = 0
counterr = 1

call initmpi
call readinput
if(rank.eq.0)print*, 'MPI OK'

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

if(infile.eq.0) then
   call solve
   call Free_Energy_Calc(counter)
   call savedata(counter/saveevery)
else
   call retrivefromdisk(counter)
   counterr = counter
   if(rank.eq.0)print*, 'Load input from file'
   if(rank.eq.0)print*, 'Free energy', free_energy
   infile = 2
   call update_matrix(flag)
   if(flag.eqv..true.) then
    print*, 'Initial position of particle does not fit in z'
    print*, 'or particles collide'
    stop
   endif

   call solve
   call Free_Energy_Calc(counter)
   if(rank.eq.0)print*, 'Free energy after solving', free_energy
endif

counter = 1
call savedata(counter)
if(rank.eq.0)print*, 'Save OK'
call store2disk(counter)
call endall
end


