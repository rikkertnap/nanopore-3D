subroutine chains_definitions

use chainsdat
use MPI
use branches
implicit none
integer i, ii
ALLOCATE (segtype(long)) 
!l=0
open(file="sequence.in",unit=333)
do i=1,long
 read(333,*) segtype(i)
! l=l+1
enddo
!if(l.ne.long) then
!    write(6,*) "error in sequence.in"
!    call MPI_FINALIZE(ierr) ! finaliza MPI
!    stop
!endif
close(333)
!segtype = 1


!if(branched.eq.1) then
!   segtype(1:longbb+longb(1)) = 2 ! backbone and first branch is hydrophobic
!endif

!if(branched.ne.1) segtype = 2

end




