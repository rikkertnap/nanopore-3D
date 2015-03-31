integer function testsystem(x)
use system
use ellipsoid
implicit none
real*8 x(3), xx(3)
integer j
real*8 vect
real*8 mmmult

! collision with walls and out of system

if (PBC.eq.0) then
 if ((x(3).ge.dimz*delta).or.(x(3).le.0.0)) then
 testsystem = -1
 print*, 'wall'
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

