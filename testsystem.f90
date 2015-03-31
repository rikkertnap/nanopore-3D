integer function testsystem(x)
use system
use ellipsoid
implicit none
real*8 x(3), xx(3)
integer j
real*8 vect
real*8 mmmult

! collision with walls and out of system

if (PBC(1).ne.0) then
 if (x(1).le.0.0) then
  testsystem = -1
  print*, 'wall', 1
  stop
  return
 endif
endif

if (PBC(2).ne.0) then
 if (x(1).ge.dimx*delta) then
  testsystem = -1
  print*, 'wall', 2
  stop
  return
 endif
endif

if (PBC(3).ne.0) then
 if (x(2).le.0.0) then
  testsystem = -1
  print*, 'wall', 3
  stop
  return
 endif
endif

if (PBC(4).ne.0) then
 if (x(2).ge.dimy*delta) then
  testsystem = -1
  print*, 'wall', 4
  stop
  return
 endif
endif

if (PBC(5).ne.0) then
 if (x(3).le.0.0) then
  testsystem = -1
  print*, 'wall', 5
  stop
  return
 endif
endif

if (PBC(6).ne.0) then
 if (x(3).ge.dimz*delta) then
  testsystem = -1
  print*, 'wall', 6
  stop
  return
 endif
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

