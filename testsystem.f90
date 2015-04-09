integer function testsystem(x)
use system
use ellipsoid
implicit none
real*8 x(3), xx(3)
integer j
real*8 vect
real*8 mmmult

! collision with walls and out of system

if (x(1).le.0.0) then
 if (PBC(1).eq.2)testsystem = -1
 if (PBC(1).eq.0)testsystem = -2
endif

if (x(1).gt.dimx*delta) then
 if (PBC(2).eq.2)testsystem = -1
 if (PBC(2).eq.0)testsystem = -2
endif

if (x(2).le.0.0) then
 if (PBC(3).eq.2)testsystem = -1
 if (PBC(3).eq.0)testsystem = -2
endif

if (x(2).gt.dimy*delta) then
 if (PBC(4).eq.2)testsystem = -1
 if (PBC(4).eq.0)testsystem = -2
endif

if (x(3).le.0.0) then
 if (PBC(5).eq.2)testsystem = -1
 if (PBC(5).eq.0)testsystem = -2
endif

if (x(3).gt.dimz*delta) then
 if (PBC(6).eq.2)testsystem = -1
 if (PBC(6).eq.0)testsystem = -2
endif

! collision with particles
do j = 1, NNN
xx(:) = x(:) - Rell(:,j)
vect = mmmult(xx,AAA(:,:,j))

!print*, x
!print*, xx
!print*, Rell(:,j)
!print*, vect

 if(vect.le.1.0) then 
  testsystem = -1
  return
 endif
enddo

testsystem = 0
return

end function

