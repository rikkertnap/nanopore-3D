subroutine chains_definitions

use chainsdat
use MPI
use branches
implicit none
integer i, ii
ALLOCATE (segtype(long)) 

segtype = 1


if(branched.eq.1) then
   segtype(longbb+1:longbb+longb(1)) = 2 ! first branch is hydrophobic
endif

end




