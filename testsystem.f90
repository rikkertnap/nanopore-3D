integer function testsystem(x)
use system
use ellipsoid
use transform
implicit none
real*8 x(3), xx(3), v(3)
integer j
real*8 vect
real*8 mmmult

! collision with walls and out of system

testsystem = 0

if (x(1).le.0.0) then
 if (PBC(1).eq.2)testsystem = -1
 if (PBC(1).eq.0)testsystem = -2
endif

if (x(1).gt.(float(dimx)*delta)) then
 if (PBC(2).eq.2)testsystem = -1
 if (PBC(2).eq.0)testsystem = -2
endif

if (x(2).le.0.0) then
 if (PBC(3).eq.2)testsystem = -1
 if (PBC(3).eq.0)testsystem = -2
endif

if (x(2).gt.(float(dimy)*delta)) then
 if (PBC(4).eq.2)testsystem = -1
 if (PBC(4).eq.0)testsystem = -2
endif

if (x(3).le.0.0) then
 if (PBC(5).eq.2)testsystem = -1
 if (PBC(5).eq.0)testsystem = -2
endif

if (x(3).gt.(float(dimz)*delta)) then
 if (PBC(6).eq.2)testsystem = -1
 if (PBC(6).eq.0)testsystem = -2
endif

! collision with particles
do j = 1, NNN
xx(:) = x(:) - Rell(:,j)


v = xx
xx = MATMUL(MAT,v)
xx(1) = xx(1) - nint(xx(1) / (float(dimx)*delta)) * float(dimx)*delta
xx(2) = xx(2) - nint(xx(2) / (float(dimy)*delta)) * float(dimy)*delta
xx(3) = xx(3) - nint(xx(3) / (float(dimz)*delta)) * float(dimz)*delta
v(:) = MATMUL(IMAT,xx)
xx = v

vect = mmmult(xx,AAA(:,:,j))

 if(vect.le.1.0) then 
  testsystem = -1
  return
 endif
enddo

return

end function

