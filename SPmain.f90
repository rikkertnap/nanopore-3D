!## Montecarlo - Molecular theory for the adsorption of a particle on a brush

use system
use MPI
use ellipsoid
use kinsol
use const
use montecarlo
use ematrix
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
integer j, i, ii, iii
integer flagcrash
real*8 stOK

stdout = 6

!!!!!!!!!!!! global parameters ... !!!!!!
saveevery = 1

counter = 0
counterr = 1

call initmpi
if(rank.eq.0)write(stdout,*) 'Program Crystal'
if(rank.eq.0)write(stdout,*) 'GIT Version: ', _VERSION
if(rank.eq.0)write(stdout,*) 'MPI OK'

if(rank.eq.0)write(10,*) 'Program Crystal'
if(rank.eq.0)write(10,*) 'GIT Version: ', _VERSION
if(rank.eq.0)write(10,*) 'MPI OK'
call readinput

call monomer_definitions
call chains_definitions

call initconst
call inittransf ! Create transformation matrixes
call initellpos ! calculate real positions for ellipsoid centers
call initall
call allocation

!!! General files

if(systemtype.eq.1) then
 do j = 1, NNN
 write(filename,'(A3,I3.3, A4)')'pos',j,'.dat'
 open(file=filename, unit=5000+j)
 write(filename,'(A3,I3.3, A4)')'orn',j,'.dat'
 open(file=filename, unit=6000+j)
 write(filename,'(A3,I3.3, A4)')'rot',j,'.dat'
 open(file=filename, unit=7000+j)
 enddo
endif

open(file='free_energy.dat', unit=9000)

call kais
if(rank.eq.0)write(stdout,*) 'Kai OK'

if (systemtype.eq.1) then
call update_matrix(flag) ! updates 'the matrix'
elseif (systemtype.eq.2) then
call update_matrix_channel(flag) ! updates 'the matrix'
elseif (systemtype.eq.3) then
call update_matrix_channel_3(flag) ! updates 'the matrix'
elseif (systemtype.eq.4) then
call update_matrix_channel_4(flag) ! updates 'the matrix'
elseif (systemtype.eq.41) then
call update_matrix_channel_4(flag) ! updates 'the matrix'
endif

  if(flag.eqv..true.) then
    write(stdout,*) 'Initial position of particle does not fit in z'
    write(stdout,*) 'or particles collide'
    stop
  else
    if(rank.eq.0)write(stdout,*) 'Particle OK'
  endif

call  graftpoints
if(rank.eq.0)write(stdout,*) 'Graftpoints OK'

call creador ! Genera cadenas
if(rank.eq.0)write(stdout,*) 'Creador OK'

if(infile.ne.0) then
   call retrivefromdisk(counter)
   if(rank.eq.0)write(stdout,*) 'Load input from file'
   if(rank.eq.0)write(stdout,*) 'Free energy', free_energy
   if(infile.eq.3)call mirror
   if(infile.ne.-1)infile = 2
!   call update_matrix(flag)
   if(flag.eqv..true.) then
    write(stdout,*) 'Initial position of particle does not fit in z'
    write(stdout,*) 'or particles collide'
    stop
   endif

endif
 ii = 1
 sc = scs(ii)

st = 1.0d10+sts(1)
do i = 1, nst
do while (st.ne.sts(i))
 st = sts(i)
 if(rank.eq.0)write(stdout,*)'Switch to st = ', st
 flagcrash = 1
do while(flagcrash.eq.1)
 flagcrash = 0
 call solve(flagcrash)
 if(flagcrash.eq.1) then
    if(i.eq.1)stop
    st = (st + stOK)/2.0
 if(rank.eq.0)write(stdout,*)'Error, switch to st = ', st
 endif
enddo

 stOK = st ! last st solved OK
 if(rank.eq.0)write(stdout,*) 'Solved OK, st: ', stOK
 

enddo



 counterr = counter + i + ii  - 1
 call Free_Energy_Calc(counterr)
 if(rank.eq.0)write(stdout,*) 'Free energy after solving', free_energy
 call savedata(counterr)
 if(rank.eq.0)write(stdout,*) 'Save OK'
 call store2disk(counterr)

enddo

call endall
end


