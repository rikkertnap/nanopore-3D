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
integer aa

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
       v = MATMUL(MAT,x)
       pxtemp(:,j) = v(:)

 
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
            do j = 1, long
            aa = floor(pxtemp(1,j)/delta) + 1
            if(PBC(1).eq.1)px(newcuantas(ii),j,jj) = PBCSYMI(aa,dimx)
            if(PBC(1).eq.3)px(newcuantas(ii),j,jj) = PBCREFI(aa,dimx)
            if((PBC(1).eq.0).or.(PBC(1).eq.2))px(newcuantas(ii),j,jj) = aa
            enddo
            do j = 1, long
            aa = floor(pxtemp(2,j)/delta) + 1
            if(PBC(3).eq.1)py(newcuantas(ii),j,jj) = PBCSYMI(aa,dimy)
            if(PBC(3).eq.3)py(newcuantas(ii),j,jj) = PBCREFI(aa,dimy)
            if((PBC(3).eq.0).or.(PBC(3).eq.2))py(newcuantas(ii),j,jj) = aa
            enddo
            do j = 1, long
            aa = floor(pxtemp(3,j)/delta) + 1
            if(PBC(5).eq.1)pz(newcuantas(ii),j,jj) = PBCSYMI(aa,dimz)
            if(PBC(5).eq.3)pz(newcuantas(ii),j,jj) = PBCREFI(aa,dimz)
            if((PBC(5).eq.0).or.(PBC(5).eq.2))pz(newcuantas(ii),j,jj) = aa
            enddo
    endif

enddo ! jj
return
end
      



