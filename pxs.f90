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
 
       v = MATMUL(MAT,x)

       do i = 1, 3
  
            pxtemp(i,j) = v(i)             ! new coordinate system
            do while (pxtemp(i,j).gt.maxx(i))
             pxtemp(i,j) = pxtemp(i,j) - maxx(i)
            enddo
            do while (pxtemp(i,j).lt.1.0d-20)
              pxtemp(i,j) = pxtemp(i,j) + maxx(i)
            enddo     
       enddo

       xx(:) = MATMUL(IMAT,pxtemp(:,j))
 
       if(testsystem(xx).eq.-1) then
         flag = -1
         exit
       endif

    enddo ! j

    if(flag.eq.0) then
    newcuantas(ii) = newcuantas(ii)+1
    px(newcuantas(ii), :, jj) = int(pxtemp(1,:)/delta) + 1
    py(newcuantas(ii), :, jj) = int(pxtemp(2,:)/delta) + 1
    pz(newcuantas(ii), :, jj) = int(pxtemp(3,:)/delta) + 1
    endif

enddo ! jj
return
end
      



