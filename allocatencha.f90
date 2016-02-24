subroutine allocatencha

use fields_fkfun
use chainsdat
use conformations
use rotchain
implicit none

! fields_fkfun
ALLOCATE(q(ncha))
ALLOCATE(sumgauche(ncha))
ALLOCATE(ngauche(cuantas,ncha))

! chainsdat
allocate(posicion(ncha,3))
allocate(ngpol(ncha))
allocate(newcuantas(ncha))

end subroutine
