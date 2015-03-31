subroutine allocatencha

use fields_fkfun
use chainsdat
implicit none

! fields_fkfun
ALLOCATE(q(ncha))


! chainsdat
allocate(posicion(ncha,3))
allocate(ngpol(ncha))
allocate(newcuantas(ncha))

end subroutine
