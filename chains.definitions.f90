subroutine chains_definitions

use chainsdat
use MPI
use branches
implicit none
integer i, ii
ALLOCATE (segtype(long)) 

segtype = 1


if(branched.eq.1) then
   segtype(1:longbb+longb(1)) = 2 ! backbone and first branch is hydrophobic
endif

if(branched.ne.1) segtype = 2

end




