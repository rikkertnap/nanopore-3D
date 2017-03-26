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
integer testsystemc
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

       if((systemtype.eq.2).or.(systemtype.eq.3).or.(systemtype.eq.4).or.(systemtype.eq.41)   &
       .or.(systemtype.eq.42))call rot_chain_cyl(x,ii)

       x = x + posicion(ii,:)
       v = MATMUL(MAT,x)
       pxtemp(:,j) = v(:)

select case (systemtype)

case (1)
 
       if(testsystem(x).eq.-1) then ! if testsystem = -1,  there is a collision with all or particle 
         flag = -1
         exit
       endif

       if(testsystem(x).eq.-2) then ! if testsystem = -2, the polymer goes out-of-system
         write(stdout,*) 'pxs: out-of-system'
         stop
       endif

case (2, 3, 4, 41, 42)


       if(testsystemc(x).eq.-1) then ! if testsystem = -1,  there is a collision with all or particle 
         flag = -1
         exit
       endif

       if(testsystemc(x).eq.-2) then ! if testsystem = -2, the polymer goes out-of-system
         write(stdout,*) 'pxs: out-of-system'
         stop
       endif

endselect




    enddo ! j

    if(flag.eq.0) then

    newcuantas(ii) = newcuantas(ii)+1
    ngauche(newcuantas(ii),ii) = ing

            do j = 1, long
            aa = floor(pxtemp(1,j)/delta) + 1
            px(newcuantas(ii),j,jj) = aa
            if(aa.lt.1) then
              if(PBC(1).eq.1)px(newcuantas(ii),j,jj) = PBCSYMI(aa,dimx)
              if(PBC(1).eq.3)px(newcuantas(ii),j,jj) = PBCREFI(aa,dimx)
            endif
            if(aa.gt.dimx) then
              if(PBC(2).eq.1)px(newcuantas(ii),j,jj) = PBCSYMI(aa,dimx)
              if(PBC(2).eq.3)px(newcuantas(ii),j,jj) = PBCREFI(aa,dimx)
            endif

            aa = floor(pxtemp(2,j)/delta) + 1
            py(newcuantas(ii),j,jj) = aa
            if(aa.lt.1) then
              if(PBC(3).eq.1)py(newcuantas(ii),j,jj) = PBCSYMI(aa,dimy)
              if(PBC(3).eq.3)py(newcuantas(ii),j,jj) = PBCREFI(aa,dimy)
            endif
            if(aa.gt.dimy) then
              if(PBC(4).eq.1)py(newcuantas(ii),j,jj) = PBCSYMI(aa,dimy)
              if(PBC(4).eq.3)py(newcuantas(ii),j,jj) = PBCREFI(aa,dimy)
            endif

            aa = floor(pxtemp(3,j)/delta) + 1
            pz(newcuantas(ii),j,jj) = aa
            if(aa.lt.1) then
              if(PBC(5).eq.1)pz(newcuantas(ii),j,jj) = PBCSYMI(aa,dimz)
              if(PBC(5).eq.3)pz(newcuantas(ii),j,jj) = PBCREFI(aa,dimz)
            endif
            if(aa.gt.dimz) then
              if(PBC(6).eq.1)pz(newcuantas(ii),j,jj) = PBCSYMI(aa,dimz)
              if(PBC(6).eq.3)pz(newcuantas(ii),j,jj) = PBCREFI(aa,dimz)
            endif
 
            enddo
    endif

enddo ! jj
return
end
      


subroutine rot_chain_cyl(x,ii)
use rotchain
implicit none
integer ii
real*8 x(3)
real*8 y(3)
real*8 t

t = rotangle(ii)
y = x
x(1) = cos(t)*y(1)-sin(t)*y(2)
x(2) = sin(t)*y(1)+cos(t)*y(2)
x(3) = y(3)
end
