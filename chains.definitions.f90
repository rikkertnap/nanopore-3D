subroutine chains_definitions

use chainsdat

implicit none
include 'mpif.h'
include 'MPI.h'
integer i, ii

ALLOCATE (segtype(long)) 
segtype(:) = 3

end




