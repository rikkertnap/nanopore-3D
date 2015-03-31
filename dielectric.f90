subroutine dielectfcn(phi,prot,epsfcn,Depsfcn)

! determines the dielectric function using an average mixing rule

use const
use system
implicit none
integer ix,iy,iz

real*8 phi(dimx,dimy,dimz)

real*8 epsfcn(0:dimx+1,0:dimy+1,0:dimz+1)

real*8 prot(dimx,dimy,dimz)

real*8 Depsfcn(0:dimx+1,0:dimy+1,0:dimz+1)

do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimz
epsfcn(ix,iy,iz) = prot(ix,iy,iz)*dielSr + (1.0-prot(ix,iy,iz))*(phi(ix,iy,iz)*dielPr + (1.0-phi(ix,iy,iz))) 
Depsfcn(ix,iy,iz) = (1.0-prot(ix,iy,iz))*(dielPr -1.0)
enddo
enddo
enddo

do ix = 1, dimx
do iy = 1, dimy
epsfcn(ix,iy,0) = epsfcn(ix,iy,1)
epsfcn(ix,iy,dimz+1) = epsfcn(ix,iy,dimz)
Depsfcn(ix,iy,0) = Depsfcn(ix,iy,1)
Depsfcn(ix,iy,dimz+1) = Depsfcn(ix,iy,dimz)
enddo
enddo

do ix = 1, dimx
do iz = 1, dimz
epsfcn(ix,0,iz) = epsfcn(ix,dimy,iz)
epsfcn(ix,dimy+1,iz) = epsfcn(ix,1,iz)
Depsfcn(ix,0,iz) = Depsfcn(ix,dimy,iz)
Depsfcn(ix,dimy+1,iz) = Depsfcn(ix,1,iz)
enddo
enddo

do iy = 1, dimy
do iz = 1, dimz
epsfcn(0,iy,iz) = epsfcn(dimx,iy,iz)
epsfcn(dimx+1,iy,iz) = epsfcn(1,iy,iz)
Depsfcn(0,iy,iz) = Depsfcn(dimx,iy,iz)
Depsfcn(dimx+1,iy,iz) = Depsfcn(1,iy,iz)
enddo
enddo

end subroutine
