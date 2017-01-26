integer function testsystemc(x)
use system
use transform
use channel

implicit none
real*8 x(3), xx(3), v(3), maxx(3)
integer j, i
real*8 vect
real*8 mmmult
real, external :: PBCSYMR, PBCREFR
real*8 dims(3)

dims(1) = delta*dimx
dims(2) = delta*dimy
dims(3) = delta*dimz

maxx(1) = float(dimx)*delta
maxx(2) = float(dimy)*delta
maxx(3) = float(dimz)*delta

! collision with walls and out of system

testsystemc = 0

v = MATMUL(MAT,x) ! to transformed space

! out-of-boundaries check is performed in transformed space

if (v(1).le.0.0) then
 if (PBC(1).eq.2)testsystemc = -1
 if (PBC(1).eq.0)testsystemc = -2
endif

if (v(1).gt.(float(dimx)*delta)) then
 if (PBC(2).eq.2)testsystemc = -1
 if (PBC(2).eq.0)testsystemc = -2
endif

if (v(2).le.0.0) then
 if (PBC(3).eq.2)testsystemc = -1
 if (PBC(3).eq.0)testsystemc = -2
endif

if (v(2).gt.(float(dimy)*delta)) then
 if (PBC(4).eq.2)testsystemc = -1
 if (PBC(4).eq.0)testsystemc = -2
endif

if (v(3).le.0.0) then
 if (PBC(5).eq.2)testsystemc = -1
 if (PBC(5).eq.0)testsystemc = -2
endif

if (v(3).gt.(float(dimz)*delta)) then
 if (PBC(6).eq.2)testsystemc = -1
 if (PBC(6).eq.0)testsystemc = -2
endif

if (testsystemc.eq.0) then ! saves some time

do i = 1, 3 ! put into cell
if(v(i).lt.0.0) then
   if(PBC(2*i-1).eq.1)v(i) = PBCSYMR(v(i),maxx(i))
   if(PBC(2*i-1).eq.3)v(i) = PBCREFR(v(i),maxx(i))
endif
if(v(i).gt.dims(i)) then
   if(PBC(2*i).eq.1)v(i) = PBCSYMR(v(i),maxx(i))
   if(PBC(2*i).eq.3)v(i) = PBCREFR(v(i),maxx(i))
endif
enddo

! collision with the channel (real coordinates)

xx(1) = x(1) - originc(1) ! distance to the center of the channel
xx(2) = x(2) - originc(2) 

vect = xx(1)**2+xx(2)**2

 if(vect.ge.rchannel**2) then
  if((v(3).gt.float(RdimZ)*delta).and.(v(3).lt.(float(dimz-RdimZ)*delta))) then  ! reservoirs in transformed coordinates
  testsystemc = -1
  return
 endif
endif

endif ! testsystem = 0

return

end function

