subroutine dielectfcn(phi,prot,epsfcn,Depsfcn)

! determines the dielectric function using an average mixing rule

use const
use system
implicit none
integer ix,iy,iz, jx,jy,jz
integer, external :: PBCSYMI, PBCREFI
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! Boundary conditions for dielectric function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reflection or PBC, (PBC = 1 or 3)

do jx = 0, dimx+1
do jy = 0, dimy+1
do jz = 0, dimz+1

ix=jx
iy=jy
iz=jz ! these lines are necessary for PBC = 0 or 2

if (PBC(1).eq.1)ix = PBCSYMI(jx,dimx) ! STANDARD PBC
if (PBC(3).eq.1)iy = PBCSYMI(jy,dimy)
if (PBC(5).eq.1)iz = PBCSYMI(jz,dimz)

if (PBC(1).eq.3)ix = PBCREFI(jx,dimx) ! REFLECTING PBC
if (PBC(3).eq.3)iy = PBCREFI(jy,dimy)
if (PBC(5).eq.3)iz = PBCREFI(jz,dimz)

   epsfcn(jx, jy, jz) = epsfcn(ix, iy, iz)
   Depsfcn(jx, jy, jz) = Depsfcn(ix, iy, iz)
enddo
enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bulk or Wall, PBC = 0 or 2

select case (PBC(1)) ! x = 0
case(0) ! set bulk 
   epsfcn(0,:,:) = 1.0 ! water dielectric
   Depsfcn(0,:,:) = 0.0 ! irrevant, shouldn't be polymer in bulk
case(2)
   epsfcn(0,:,:) = 0.0 ! dielectric is zero at wall
   Depsfcn(0,:,:) = 0.0 ! irrevant, shouldn't be polymer in wall
endselect

select case (PBC(2)) ! x = dimx
case(0) ! set bulk 
   epsfcn(dimx+1,:,:) = 1.0 ! water dielectric
   Depsfcn(dimx+1,:,:) = 0.0 ! irrevant, shouldn't be polymer in bulk
case(2)
   epsfcn(dimx+1,:,:) = 0.0 ! dielectric is zero at wall
   Depsfcn(dimx+1,:,:) = 0.0 ! irrevant, shouldn't be polymer in wall
endselect

select case (PBC(3)) ! y = 0
case(0) ! set bulk 
   epsfcn(:,0,:) = 1.0 ! water dielectric
   Depsfcn(:,0,:) = 0.0 ! irrevant, shouldn't be polymer in bulk
case(2)
   epsfcn(:,0,:) = 0.0 ! dielectric is zero at wall
   Depsfcn(:,0,:) = 0.0 ! irrevant, shouldn't be polymer in wall
endselect

select case (PBC(4)) ! y = dimy
case(0) ! set bulk 
   epsfcn(:,dimy+1,:) = 1.0 ! water dielectric
   Depsfcn(:,dimy+1,:) = 0.0 ! irrevant, shouldn't be polymer in bulk
case(2)
   epsfcn(:,dimy+1,:) = 0.0 ! dielectric is zero at wall
   Depsfcn(:,dimy+1,:) = 0.0 ! irrevant, shouldn't be polymer in wall
endselect

select case (PBC(5)) ! z = 0
case(0)
   epsfcn(:,:,0) = 1.0 ! water dielectric
   Depsfcn(:,:,0) = 0.0 ! irrevant, shouldn't be polymer in bulk
case(2)
   epsfcn(:,:,0) = 0.0 ! dielectric is zero at wall
   Depsfcn(:,:,0) = 0.0 ! irrevant, shouldn't be polymer in wall
endselect

select case (PBC(6)) ! z = dimz
case(0) ! set bulk 
   epsfcn(:,:,dimz+1) = 1.0 ! water dielectric
   Depsfcn(:,:,dimz+1) = 0.0 ! irrevant, shouldn't be polymer in bulk
case(2)
   epsfcn(:,:,dimz+1) = 0.0 ! dielectric is zero at wall
   Depsfcn(:,:,dimz+1) = 0.0 ! irrevant, shouldn't be polymer in wall
endselect
end subroutine
