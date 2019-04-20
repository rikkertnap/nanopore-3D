subroutine allocatecpp
use fields_fkfun
use conformations
use chainsdat
implicit none

ALLOCATE(px(cuantas, long, maxcpp))
ALLOCATE(py(cuantas, long, maxcpp))
ALLOCATE(pz(cuantas, long, maxcpp))
ALLOCATE(pro(cuantas, maxcpp))
ALLOCATE(zfinal(cuantas, maxcpp))
end subroutine
