!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Esta subrutina se encarga de poner a todas los segmentos dentro del slab

subroutine pxs

use system
use MPI
use chainsdat
use conformations
use const
use transform
implicit none
    
integer j, ii, jj,i
real*8 pxtemp(3,long)
real*8 xx(3)
real*8 x(3)
real*8 v(3)
integer testsystem
real*8 maxx(3)
integer flag
integer a1,a2,a3

integer, external :: PBCREFI, PBCSYMI

maxx(1) = float(dimx)*delta
maxx(2) = float(dimy)*delta
maxx(3) = float(dimz)*delta

do jj = 1, cpp(rank+1)
  ii = cppini(rank+1)+jj
  flag = 0

    do j=1,long
       x(1) = in1(j ,2)
       x(2) = in1(j, 3)
       x(3) = in1(j, 1)

       x = x + posicion(ii,:)
 
       if(testsystem(x).eq.-1) then ! if testsystem = -1,  there is a collision with all or particle 
         flag = -1
         exit
       endif

       if(testsystem(x).eq.-2) then ! if testsystem = -2, the polymer goes out-of-system
         print*, 'pxs: out-of-system'
         stop
       endif


    enddo ! j

    if(flag.eq.0) then
    newcuantas(ii) = newcuantas(ii)+1
    px(newcuantas(ii), :, jj) = int(pxtemp(1,:)/delta) + 1
            if(PBC(1).eq.1)px(newcuantas(ii), :, jj) = PBCSYMI(px(newcuantas(ii),:,jj),dimx)
            if(PBC(1).eq.3)px(newcuantas(ii), :, jj) = PBCREFI(px(newcuantas(ii),:,jj),dimx)
 
    py(newcuantas(ii), :, jj) = int(pxtemp(2,:)/delta) + 1
            if(PBC(3).eq.1)py(newcuantas(ii),:,jj) = PBCSYMI(py(newcuantas(ii),:,jj),dimy)
            if(PBC(3).eq.3)py(newcuantas(ii),:,jj) = PBCREFI(py(newcuantas(ii),:,jj),dimy)

    pz(newcuantas(ii), :, jj) = int(pxtemp(3,:)/delta) + 1
            if(PBC(5).eq.1)pz(newcuantas(ii),:,jj) = PBCSYMI(pz(newcuantas(ii),:,jj),dimz)
            if(PBC(5).eq.3)pz(newcuantas(ii),:,jj) = PBCREFI(pz(newcuantas(ii),:,jj),dimz)

    endif

enddo ! jj
return
end
      



