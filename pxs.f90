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

real, external :: PBCSYMR, PBCREFR
real*8, parameter :: erd =0.99999

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
 
       v = MATMUL(MAT,x) ! to transformed space

       do i = 1, 3
            pxtemp(i,j) = v(i)             ! new coordinate system
            if(PBC(2*i-1).eq.1)pxtemp(i,j) = PBCSYMR(pxtemp(i,j),maxx(i))
            if(PBC(2*i-1).eq.3)pxtemp(i,j) = PBCREFR(pxtemp(i,j),maxx(i))
       enddo

       xx(:) = MATMUL(IMAT,pxtemp(:,j)) ! to real space

       if(testsystem(xx).eq.-1) then ! if testsystem = -1,  there is a collision with all or particle 
         flag = -1
         exit
       endif

       if(testsystem(xx).eq.-2) then ! if testsystem = -2, the polymer goes out-of-system
         print*, 'pxs: out-of-system'
         stop
       endif


    enddo ! j

    if(flag.eq.0) then
    newcuantas(ii) = newcuantas(ii)+1
    px(newcuantas(ii), :, jj) = int(pxtemp(1,:)*erd/delta) + 1 ! erd compresses coordinates very slighly to prevent numerical errors
    py(newcuantas(ii), :, jj) = int(pxtemp(2,:)*erd/delta) + 1
    pz(newcuantas(ii), :, jj) = int(pxtemp(3,:)*erd/delta) + 1
    endif

enddo ! jj
return
end
      



