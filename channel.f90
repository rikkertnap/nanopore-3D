subroutine update_matrix_channel_4(flag)
use system
use channel
use ematrix
use MPI
use const
use chainsdat
use molecules
use channel
use transform, only : MAT, IMAT
use rotchain

implicit none

real*8 rchannel2, rchannelL2, rchannelS2
real*8, external :: rands
integer npoints ! points per cell for numerical integration 
integer counter
character*5 title
logical flag
integer j,ix,iy,iz
real pnumber
real*8 area
real*8 sumpolseg 
real*8 sstemp,vvtemp, maxss
real*8 cutarea
real*8 temp
real*8 temp2
real*8 sumvoleps1, sumvolprot1, sumvolq1, sumvolx1
integer ncha1
real*8 volx1(maxvolx)
real*8 com1(maxvolx,3)
integer p1(maxvolx,3)
integer i
real*8 volxx1(dimx,dimy,dimz)
real*8 volxx(dimx,dimy,dimz)
real*8 x(3), v(3), hcyl
integer nbands

cutarea = 0.0 ! throw away cells that have less area than cutarea x area of the cell with largest area  
sumpolseg = 0.0

rchannel2 = rchannel**2
rchannelL2 = (rchannel - 3*delta)**2
rchannelS2 = (rchannel + delta)**2

! clear all
voleps = 0.0
volprot = 0.0
volq = 0.0
volx = 0.0
volxx = 0.0
com = 0.0
ncha = 0

! channel center in x, y plane

 originc(1) = float(dimx)*delta/2.0 
 originc(2) = float(dimy)*delta/2.0 

 npoints = 50

 flag = .false.

 call integrate_c(rchannelL2, RdimZ,originc ,npoints, voleps1 , sumvoleps1, flag)

 flag = .false. ! not a problem if eps lays outside boundaries

 call integrate_c(rchannel2,RdimZ, originc,npoints, volprot1, sumvolprot1, flag)

 call integrate_c(rchannelS2,RdimZ,originc,npoints, volq1, sumvolq1, flag)

 call newintegrateg_c_4(rchannel2,RdimZ,originc,npoints,volx1,sumvolx1, com1, p1, ncha1, volxx1, NBRUSH)

!! eps
 voleps1 = voleps1-volprot1
 voleps1 = voleps1*eepsc

! epstype


select case (epstype)
case (1)
nbands = dimz/8
do iz = 1, dimz
if (mod(int((iz-1)/nbands),2).eq.1) then
voleps1(:,:,iz) = 0.0
endif
enddo
endselect


!! charge
 volq1 = volprot1-volq1
 temp = sumvolprot1-sumvolq1
 volq1 = volq1/temp*echargec/(delta**3) ! sum(volq) is echarge

!! grafting

v(1) = 0.0
v(2) = 0.0
v(3) = float(dimz-2*RdimZ)*delta

! v in transformed space, x in real space
! only work for gam = 90, cdiva any value

x = MATMUL(IMAT,v)

hcyl = x(3) ! height of the cylinder

area = 2.0*pi*rchannel*hcyl

!! volume  
 volprot1 = volprot1 * 0.9999
 volprot = volprot+volprot1

! CHECK COLLISION HERE...
 if(maxval(volprot).gt.1.0) then ! collision
   flag=.true. 
 endif
 
 voleps = voleps + voleps1
 volq = volq + volq1 

! add com1 and volx to list

 volxx = volxx1

ncha = ncha1
do i = 1, ncha
volx(i)=volx1(i)
com(i,:)=com1(i,:)
p0(i,:)=p1(i,:)
rotangle(i) = atan2(com1(i,1)-originc(1), com1(i,2)-originc(2))
enddo

title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)

sumpolseg = ncha

if (verbose.ge.2) then
temp = pi*rchannel2*float(dimz)*delta
!do j = 1, NNN
!temp = temp + 4.0/3.0*pi*Aell(1,j)*Aell(2,j)*Aell(3,j)
!enddo
if (rank.eq.0) then
write(stdout,*) 'channel:', 'update_matrix: Total nanochannel volumen real space= ', temp
write(stdout,*) 'channel:', 'update_matrix: Total discretized volumen =', (dimx*dimy*dimz-sum(volprot))*delta**3
write(stdout,*) 'channel:', 'number of polymers in system =', sumpolseg 
write(stdout,*) 'channel:', 'surface area =', area
write(stdout,*) 'channel:', 'surface density =', sumpolseg/area
write(stdout,*) 'channel:', 'surface density expected from input =', &
 float(NBRUSH)/(2.0*pi*rchannel)/cos(30.0/180.0*pi)/(2.0*pi*rchannel/float(NBRUSH))
endif
endif

title = 'aveps'
counter = 1
call savetodisk(voleps, title, counter)

title = 'avcha'
counter = 1
call savetodisk(volq, title, counter)

title = 'avgrf'
counter = 1
call savetodisk(volxx, title, counter)
end subroutine





subroutine update_matrix_channel_3(flag)
use system
use channel
use ematrix
use MPI
use const
use chainsdat
use molecules
use channel
use transform, only : MAT, IMAT
use rotchain

implicit none

real*8 rchannel2, rchannelL2, rchannelS2
real*8, external :: rands
integer npoints ! points per cell for numerical integration 
integer counter
character*5 title
logical flag
integer j,ix,iy,iz
real pnumber
real*8 area
real*8 sumpolseg 
real*8 sstemp,vvtemp, maxss
real*8 cutarea
real*8 temp
real*8 temp2
real*8 sumvoleps1, sumvolprot1, sumvolq1, sumvolx1
integer ncha1
real*8 volx1(maxvolx)
real*8 com1(maxvolx,3)
integer p1(maxvolx,3)
integer i
real*8 volxx1(dimx,dimy,dimz)
real*8 volxx(dimx,dimy,dimz)
real*8 x(3), v(3), hcyl
integer nbands

cutarea = 0.0 ! throw away cells that have less area than cutarea x area of the cell with largest area  
sumpolseg = 0.0

rchannel2 = rchannel**2
rchannelL2 = (rchannel - 3*delta)**2
rchannelS2 = (rchannel + delta)**2

! clear all
voleps = 0.0
volprot = 0.0
volq = 0.0
volx = 0.0
volxx = 0.0
com = 0.0
ncha = 0

! channel center in x, y plane

 originc(1) = float(dimx)*delta/2.0 
 originc(2) = float(dimy)*delta/2.0 

 npoints = 50

 flag = .false.

 call integrate_c(rchannelL2,RdimZ, originc ,npoints, voleps1 , sumvoleps1, flag)

 flag = .false. ! not a problem if eps lays outside boundaries

 call integrate_c(rchannel2, RdimZ, originc,npoints, volprot1, sumvolprot1, flag)

 call integrate_c(rchannelS2,RdimZ, originc,npoints, volq1, sumvolq1, flag)

 call newintegrateg_c_3(rchannel2,RdimZ, originc,npoints,volx1,sumvolx1, com1, p1, ncha1, volxx1, NBRUSH)

!! eps
 voleps1 = voleps1-volprot1
 voleps1 = voleps1*eepsc

!! charge
 volq1 = volprot1-volq1
 temp = sumvolprot1-sumvolq1
 volq1 = volq1/temp*echargec/(delta**3) ! sum(volq) is echarge

!! grafting

v(1) = 0.0
v(2) = 0.0
v(3) = float(dimz-2*RdimZ)*delta

! v in transformed space, x in real space
! only work for gam = 90, cdiva any value

x = MATMUL(IMAT,v)

hcyl = x(3) ! height of the cylinder

area = 2.0*pi*rchannel*hcyl

!! volume  
 volprot1 = volprot1 * 0.9999
 volprot = volprot+volprot1

! CHECK COLLISION HERE...
 if(maxval(volprot).gt.1.0) then ! collision
   flag=.true. 
 endif
 
 voleps = voleps + voleps1
 volq = volq + volq1 

! add com1 and volx to list

 volxx = volxx1

ncha = ncha1
do i = 1, ncha
volx(i)=volx1(i)
com(i,:)=com1(i,:)
p0(i,:)=p1(i,:)
rotangle(i) = atan2(com1(i,1)-originc(1), com1(i,2)-originc(2))
enddo

title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)

sumpolseg = ncha

if (verbose.ge.2) then
temp = pi*rchannel2*float(dimz)*delta
!do j = 1, NNN
!temp = temp + 4.0/3.0*pi*Aell(1,j)*Aell(2,j)*Aell(3,j)
!enddo
if (rank.eq.0) then
write(stdout,*) 'channel:', 'update_matrix: Total nanochannel volumen real space= ', temp
write(stdout,*) 'channel:', 'update_matrix: Total discretized volumen =', (dimx*dimy*dimz-sum(volprot))*delta**3
write(stdout,*) 'channel:', 'number of polymers in system =', sumpolseg 
write(stdout,*) 'channel:', 'surface area =', area
write(stdout,*) 'channel:', 'surface density =', sumpolseg/area
write(stdout,*) 'channel:', 'surface density expected from input =', 1.0/(((2.0*pi*rchannel)/float(NBRUSH))**2)
endif
endif

title = 'aveps'
counter = 1
call savetodisk(voleps, title, counter)

title = 'avcha'
counter = 1
call savetodisk(volq, title, counter)

title = 'avgrf'
counter = 1
call savetodisk(volxx, title, counter)
end subroutine




subroutine update_matrix_channel(flag)
use system
use channel
use ematrix
use MPI
use const
use chainsdat
use molecules
use rotchain
implicit none

real*8 rchannel2, rchannelL2, rchannelS2
real*8, external :: rands
integer npoints ! points per cell for numerical integration 
integer counter
character*5 title
logical flag
integer j,ix,iy,iz
real pnumber
real*8 area
real*8 sumpolseg 
real*8 sstemp,vvtemp, maxss
real*8 cutarea
real*8 temp
real*8 temp2
real*8 sumvoleps1, sumvolprot1, sumvolq1, sumvolx1
integer ncha1
real*8 volx1(maxvolx)
real*8 com1(maxvolx,3)
integer p1(maxvolx,3)
integer i
real*8 volxx1(dimx,dimy,dimz)
real*8 volxx(dimx,dimy,dimz)
integer nbands


cutarea = 0.0 ! throw away cells that have less area than cutarea x area of the cell with largest area  
sumpolseg = 0.0

rchannel2 = rchannel**2
rchannelL2 = (rchannel - 3*delta)**2
rchannelS2 = (rchannel + delta)**2

! clear all
voleps = 0.0
volprot = 0.0
volq = 0.0
volx = 0.0
volxx = 0.0
com = 0.0
ncha = 0

! channel center in x, y plane

 originc(1) = float(dimx)*delta/2.0 
 originc(2) = float(dimy)*delta/2.0 

 npoints = 50

 flag = .false.


 call integrate_c(rchannelL2,RdimZ, originc ,npoints, voleps1 , sumvoleps1, flag)

 flag = .false. ! not a problem if eps lays outside boundaries

 call integrate_c(rchannel2, RdimZ, originc,npoints, volprot1, sumvolprot1, flag)

 call integrate_c(rchannelS2,RdimZ, originc,npoints, volq1, sumvolq1, flag)

 call newintegrateg_c(rchannel2,RdimZ,originc,npoints,volx1,sumvolx1, com1, p1, ncha1, volxx1)

!! eps
 voleps1 = voleps1-volprot1
 voleps1 = voleps1*eepsc

! epstype

select case (epstype)

case (1)
nbands = dimz/8
do iz = 1, dimz
if (mod(int((iz-1)/nbands),2).eq.1) then
voleps1(:,:,iz) = 0.0
endif
enddo

endselect


!! charge
 volq1 = volprot1-volq1
 temp = sumvolprot1-sumvolq1
 volq1 = volq1/temp*echargec/(delta**3) ! sum(volq) is echarge

!! grafting

 area = 2.0*pi*rchannel*float(dimz)*delta

 temp2 = maxval(volx1)

 where(volx1<temp2*cutarea) ! remove cells with very little area
 volx1 = 0.0
 end where 

 do i = 1, ncha1
 volx1(i) = volx1(i)/sumvolx1*area*(sigmac+sigmar*(rands(seed)-0.5))
 volxx1(p1(i,1),p1(i,2),p1(i,3)) = & 
     volxx1(p1(i,1),p1(i,2),p1(i,3))/sumvolx1*area*(sigmac+sigmar*(rands(seed)-0.5))
 enddo

 maxss = 1.0d100
 sumpolseg = sumpolseg + area*sigmac*long

!! volume  
 volprot1 = volprot1 * 0.9999
 volprot = volprot+volprot1

! CHECK COLLISION HERE...
 if(maxval(volprot).gt.1.0) then ! collision
   flag=.true. 
 endif
 
 voleps = voleps + voleps1
 volq = volq + volq1 

! add com1 and volx to list

 volxx = volxx + volxx1

ncha = ncha1
do i = 1, ncha
volx(i)=volx1(i)
com(i,:)=com1(i,:)
p0(i,:)=p1(i,:)
rotangle(i) = atan2(com1(i,1)-originc(1), com1(i,2)-originc(2))
enddo

title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)

if (verbose.ge.2) then
temp = pi*rchannel2*float(dimz)*delta
!do j = 1, NNN
!temp = temp + 4.0/3.0*pi*Aell(1,j)*Aell(2,j)*Aell(3,j)
!enddo
if (rank.eq.0) then
write(stdout,*) 'channel:', 'update_matrix: Total nanochannel volumen real space= ', temp
write(stdout,*) 'channel:', 'update_matrix: Total discretized volumen =', (dimx*dimy*dimz-sum(volprot))*delta**3
write(stdout,*) 'channel:', 'number of monomers in system =', sumpolseg 
endif
endif

title = 'aveps'
counter = 1
call savetodisk(voleps, title, counter)

title = 'avcha'
counter = 1
call savetodisk(volq, title, counter)

title = 'avgrf'
counter = 1
call savetodisk(volxx, title, counter)

end subroutine

subroutine integrate_c(rchannel2,RdimZ, originc, npoints,volprot,sumvolprot, flag)
use system
use transform

implicit none
real*8 sumvolprot
integer npoints
real*8 rchannel2, originc(2)
real*8 volprot(dimx,dimy,dimz)
real*8 dr(3), dxr(3)
integer ix,iy,iz,ax,ay,az
real*8 vect
logical flagin, flagout
real*8 intcell_c
real*8 mmmult
integer jx,jy, jz
logical flag
integer RdimZ

real*8 box(4)
real*8 x(3), v(3)
integer xmin,xmax,ymin,ymax,zmin,zmax
integer i,j

logical flagsym
real*8 voltemp

volprot = 0.0
sumvolprot = 0.0 ! total volumen, including that outside the system


! scan over all cells

do ix = 1, dimx
do iy = 1, dimy
do iz = RdimZ+1, dimz-RdimZ

flagin = .false.
flagout = .false.

do ax = 0,1
do ay = 0,1
do az = 0,1

! v in transformed space
v(1) = float(ax+ix-1)*delta
v(2) = float(ay+iy-1)*delta
v(3) = 0.0

! x in real space, v in transformed space
    x = MATMUL(IMAT,v)

x(1) = x(1) - originc(1)
x(2) = x(2) - originc(2)

if((x(1)**2+x(2)**2).lt.rchannel2)flagin=.true. ! inside the channel
if((x(1)**2+x(2)**2).gt.rchannel2)flagout=.true. ! outside the channel

enddo
enddo
enddo

if((flagin.eqv..true.).and.(flagout.eqv..false.)) then ! cell all inside channel
    voltemp = 0.0
endif
if((flagin.eqv..false.).and.(flagout.eqv..true.)) then ! cell all outside channel
    voltemp = 1.0
endif
if((flagin.eqv..true.).and.(flagout.eqv..true.)) then ! cell part inside annd outside channel
    voltemp = intcell_c(rchannel2, originc,ix,iy,iz, npoints)
endif

sumvolprot = sumvolprot + voltemp
volprot(ix,iy,iz) = voltemp

enddo ! ix
enddo ! iy
enddo ! iz


end subroutine

double precision function intcell_c(rchannel2,originc,ix,iy,iz,n)
use system
use transform

implicit none
real*8 rchannel2
real*8 originc(2)
integer ix,iy,iz,ax,ay,az
integer cc
real*8 vect
integer n
real*8 mmmult
real*8 dr(3), dxr(3)

cc = 0
do ax = 1, n
do ay = 1, n
do az = 1, n

dr(1) = ix*delta-(ax)*delta/float(n) 
dr(2) = iy*delta-(ay)*delta/float(n) 
dr(3) = iz*delta-(az)*delta/float(n) 

! dr in transformed space
dxr = MATMUL(IMAT, dr)

dxr(1) = dxr(1)-originc(1)
dxr(2) = dxr(2)-originc(2)

vect = dxr(1)**2+dxr(2)**2
if(vect.gt.rchannel2)cc=cc+1 ! outside channel, integrate

enddo
enddo
enddo

intcell_c = float(cc)/(float(n)**3)
end function

subroutine newintegrateg_c_4(rchannel2,RdimZ,originc, npoints,volx1,sumvolx1,com1,p1,ncha1,volxx1, NBRUSH)
use system
use transform
use chainsdat
use ematrix
use const
use channel, only : sigmar
implicit none
real*8 rtetha, rz
integer NBRUSH
real*8 sumvolx1
integer npoints
integer indexvolx(dimx,dimy,dimz)
integer listvolx(ncha,3)
real*8 radio
real*8 rchannel, rchannel2, originc(2)
real*8 phi, dphi, tetha,dtetha, as, ds
integer mphi, mtetha
integer ix,iy,iz,jx,jy,jz
real*8 x(3), v(3)
integer i,j
integer ncount
real*8 comshift ! how far from the surface of the sphere the grafting point is
integer ncha1 ! count for current sphere
real*8 volx1(maxvolx)
real*8 com1(maxvolx,3)
integer p1(maxvolx,3)
real*8 volxx1(dimx,dimy,dimz)
integer flagin
integer dims(3), is(3), js(3)
integer jjjz, jjjt, npointz, npointt
real*8 hcyl
real*8 hcyl0
real*8, external :: rands
real*8 tethaadd, disp
integer RdimZ

disp = delta

pi=acos(-1.0)

dims(1) = dimx
dims(2) = dimy
dims(3) = dimz

rchannel = sqrt(rchannel2)

indexvolx = 0
ncha1 = 0
volx1 = 0.0
sumvolx1 = 0.0 ! total volumen, including that outside system
com1 = 0.0
p1 = 0
volxx1 = 0.0

!
! This routine determines the surface coverage and grafting positions only for cylinder
!
!

v(1) = 0.0
v(2) = 0.0
v(3) = float(dimz-2*RdimZ)*delta

! v in transformed space, x in real space
! only work for gam = 90, cdiva any value

x = MATMUL(IMAT,v)

hcyl = x(3) ! height of the cylinder


npointt = NBRUSH    ! number of sites along the tetha coordinate



select case (systemtype)
case (4) ! uniformed coated cylinder
   npointz = nint(float(NBRUSH)*hcyl/(2.0*pi*rchannel)/cos(30.0/180.0*pi))   ! number of sites along the z - coordinate, first site is at bdist
   v(1) = 0.0
   v(2) = 0.0
   v(3) = float(RdimZ+1)*delta
   x = MATMUL(IMAT,v)
   hcyl0 = x(3) ! position of the base of the cylinder
case (41) ! only one row
   npointz = 1
   v(1) = 0.0
   v(2) = 0.0
   v(3) = float(dimz)*delta/2.0
   x = MATMUL(IMAT,v)
   hcyl0 = x(3) ! position of the base of the cylinder
endselect

do jjjz = 1, npointz
do jjjt = 1, npointt

select case (randominput)
 
 case(0)
 rtetha = 0.0
 rz = 0.0

 case(1)
 
 rtetha = disp/2.0/(2.0*pi*rchannel)*(rands(seed)-1.0)
 rz = disp/delta*hcyl/float(dimz-2*RdimZ)/2.0*(rands(seed)-1.0)

 case(2)
 rtetha = 0.0
 rz = disp*(float(mod(jjjz,2))-0.5)

 case(20)
 rtetha = 0.0
 rz = disp*(float(mod(jjjz,4))-0.5)


 case(3)
 rtetha = 2.0*disp*(float(mod(jjjt,2))-0.5)
 rz = 0.0

 case(4)
 rtetha = 0.0
 rz = (float(jjjt)/float(npointt)-0.5)*4.0*disp

 case(40)
 disp = hcyl/8.0
 rtetha = 0.0
 rz = (((float(jjjt)-1)/float(npointt))-0.5)*disp*2.0
 rz = rz + disp/2.0*(float(mod(jjjz,2))-0.5)

 case(50)
 disp = hcyl/8.0
 rtetha = 0.0

 if(jjjt.lt.npointt/2) then
 rz = (((float(jjjt)-1)/float(npointt/2))-0.5)*disp*4.0
 rz = rz + disp/2.0*(float(mod(jjjz,2))-0.5)
 else
 rz = (((float(jjjt)-1-npointt/2)/float(npointt/2))-0.5)*disp*4.0
 rz = rz + disp/2.0*(float(mod(jjjz,2))-0.5)
 endif

 case(51)
 disp = hcyl/8.0
 rtetha = 0.0
 rz = (((float(jjjt)-1)/float(npointt))-0.5)*disp*4.0
 rz = rz + disp/2.0*(float(mod(jjjz,2))-0.5)

 case(52)
 disp = hcyl/8.0
 rtetha = 0.0
 rz = (((float(jjjt)-1)/float(npointt))-0.5)*disp*3.0
 rz = rz + disp/2.0*(float(mod(jjjz,2))-0.5)


 case(41)
 rtetha = 2.0*disp*(float(mod(jjjt,2))-0.5)
 rz = (float(jjjt)/float(npointt)-0.5)*4.0*disp

endselect

tethaadd = mod(jjjz,2)*2.0*pi/float(npointt)*0.5

x(1) = cos(float(jjjt-1)/float(npointt)*2.0*pi+rtetha+tethaadd)*rchannel + originc(1)
x(2) = sin(float(jjjt-1)/float(npointt)*2.0*pi+rtetha+tethaadd)*rchannel + originc(2)
x(3) = float(jjjz-1)/float(npointz)*hcyl+rz+hcyl0

!x in  real space
v = MATMUL(MAT,x)
v(3) = v(3) + float(dimz/npointz)/2.0*delta ! centers the first row of polymers at the middle of the layer, useful to avoid numerical rounding errors.
x = MATMUL(IMAT,v) ! and recalculates x due to the change in v

do j = 1,3
    js(j) = floor(v(j)/delta)+1
enddo

js(3)=mod(js(3)+dimz-1,dimz)+1

jx = js(1)
jy = js(2)
jz = js(3)

do i = 1, 3
if((js(i).le.0).or.(js(i).gt.dims(i))) then
   write(stdout,*) 'error in channel', i, js(i), dims(i)
    write(stdout,*) v(1), v(2), v(3)
   stop
endif
enddo

! increase counter

 if(ncha1.eq.maxvolx) then
   write(stdout,*) 'channel: increase maxvolx'
   stop
 endif

 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz

 volxx1(jx,jy,jz) =  1.0
 volx1(indexvolx(jx,jy,jz)) = 1.0
 com1(indexvolx(jx,jy,jz),:) = x(:)
 sumvolx1 = sumvolx1 + 1.0

enddo ! jjjt
enddo ! jjjz

do i = 1, ncha1
! Moves the position of the first segment lseg/2 away from the surface to prevent collision due to round errors.
com1(i,1) = com1(i,1) - lseg*((com1(i,1)-originc(1)))/rchannel 
com1(i,2) = com1(i,2) - lseg*((com1(i,2)-originc(2)))/rchannel 
enddo
end



subroutine newintegrateg_c_3(rchannel2,RdimZ,originc, npoints,volx1,sumvolx1,com1,p1,ncha1,volxx1, NBRUSH)
use system
use transform
use chainsdat
use ematrix
use const
use channel, only : sigmar
implicit none
real*8 rtetha, rz
integer NBRUSH
real*8 sumvolx1
integer npoints
integer indexvolx(dimx,dimy,dimz)
integer listvolx(ncha,3)
real*8 radio
real*8 rchannel, rchannel2, originc(2)
real*8 phi, dphi, tetha,dtetha, as, ds
integer mphi, mtetha
integer ix,iy,iz,jx,jy,jz
real*8 x(3), v(3)
integer i,j
integer ncount
real*8 comshift ! how far from the surface of the sphere the grafting point is
integer ncha1 ! count for current sphere
real*8 volx1(maxvolx)
real*8 com1(maxvolx,3)
integer p1(maxvolx,3)
real*8 volxx1(dimx,dimy,dimz)
integer flagin
integer dims(3), is(3), js(3)
integer jjjz, jjjt, npointz, npointt
real*8 hcyl, hcyl0
real*8, external :: rands
integer RdimZ

pi=acos(-1.0)

dims(1) = dimx
dims(2) = dimy
dims(3) = dimz

rchannel = sqrt(rchannel2)

indexvolx = 0
ncha1 = 0
volx1 = 0.0
sumvolx1 = 0.0 ! total volumen, including that outside system
com1 = 0.0
p1 = 0
volxx1 = 0.0

!
! This routine determines the surface coverage and grafting positions only for cylinder
!
!

v(1) = 0.0
v(2) = 0.0
v(3) = float(dimz-2*RdimZ)*delta

! v in transformed space, x in real space
! only work for gam = 90, cdiva any value

x = MATMUL(IMAT,v)

hcyl = x(3) ! height of the cylinder

v(1) = 0.0
v(2) = 0.0
v(3) = float(RdimZ+1)*delta

! v in transformed space, x in real space
! only work for gam = 90, cdiva any value

x = MATMUL(IMAT,v)

hcyl0 = x(3) ! height of the cylinder

npointz = nint(float(NBRUSH)*hcyl/(2.0*pi*rchannel))   ! number of sites along the z - coordinate, first site is at bdist
                             ! rounds to nearest number to avoid rounding errors later

npointt = NBRUSH    ! number of sites along the tetha coordinate



do jjjz = 1, npointz
do jjjt = 1, npointt

if(sigmar.ne.0.0) then
 rtetha = delta/2.0/(2.0*pi*rchannel)*(rands(seed)-1.0)
 rz = hcyl/float(dimz-2*RdimZ)/2.0*(rands(seed)-1.0)
endif

x(1) = cos(float(jjjt-1)/float(npointt)*2.0*pi+rtetha)*rchannel + originc(1)
x(2) = sin(float(jjjt-1)/float(npointt)*2.0*pi+rtetha)*rchannel + originc(2)
x(3) = float(jjjz-1)/float(npointz)*hcyl+rz+hcyl0


!x in  real space
v = MATMUL(MAT,x)

do j = 1,3
    js(j) = floor(v(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

v(3) = v(3) + 0.5*delta ! centers the first row of polymers at the middle of the layer, useful to avoid numerical rounding errors.
x = MATMUL(IMAT,v) ! and recalculates x due to the change in v

do j = 1,3
    js(j) = floor(v(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

do i = 1, 3
if((js(i).le.0).or.(js(i).gt.dims(i))) then
   write(stdout,*) 'error in channel', i, js(i), dims(i)
    write(stdout,*) v(1), v(2), v(3)
   stop
endif
enddo

! increase counter

 if(ncha1.eq.maxvolx) then
   write(stdout,*) 'channel: increase maxvolx'
   stop
 endif

 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz

 volxx1(jx,jy,jz) =  1.0
 volx1(indexvolx(jx,jy,jz)) = 1.0
 com1(indexvolx(jx,jy,jz),:) = x(:)
 sumvolx1 = sumvolx1 + 1.0

enddo ! jjjt
enddo ! jjjz

do i = 1, ncha1
! Moves the position of the first segment lseg/2 away from the surface to prevent collision due to round errors.
com1(i,1) = com1(i,1) - lseg*((com1(i,1)-originc(1)))/rchannel 
com1(i,2) = com1(i,2) - lseg*((com1(i,2)-originc(2)))/rchannel 
enddo
end


subroutine newintegrateg_c(rchannel2,RdimZ,originc, npoints,volx1,sumvolx1,com1,p1,ncha1,volxx1)
use system
use transform
use chainsdat
use ematrix
use const
implicit none
real*8 sumvolx1
integer npoints
integer indexvolx(dimx,dimy,dimz)
integer listvolx(ncha,3)
real*8 radio
real*8 rchannel, rchannel2, originc(2)
real*8 phi, dphi, tetha,dtetha, as, ds
integer mphi, mtetha
integer ix,iy,iz,jx,jy,jz
real*8 x(3), v(3)
integer i,j
integer ncount
real*8 comshift ! how far from the surface of the sphere the grafting point is
integer ncha1 ! count for current sphere
real*8 volx1(maxvolx)
real*8 com1(maxvolx,3)
integer p1(maxvolx,3)
real*8 volxx1(dimx,dimy,dimz)
integer flagin
integer dims(3), is(3), js(3)
integer jjjz, jjjt, npointz, npointt
integer RdimZ

pi=acos(-1.0)

dims(1) = dimx
dims(2) = dimy
dims(3) = dimz

rchannel = sqrt(rchannel2)

indexvolx = 0
ncha1 = 0
volx1 = 0.0
sumvolx1 = 0.0 ! total volumen, including that outside system
com1 = 0.0
p1 = 0
volxx1 = 0.0

! This routine determines the surface coverage and grafting positions only for cylinder
!

npointz = npoints*dimz
npointt = int(2.0*pi*rchannel/delta)*npoints

do jjjz = 1, npointz-1
do jjjt = 1, npointt

x(1) = cos(float(jjjt)/float(npointt)*2.0*pi)*rchannel
x(2) = sin(float(jjjt)/float(npointt)*2.0*pi)*rchannel
x(3) = float(jjjz)/float(npointz)*float(dimz)*delta

x(1) = x(1) + originc(1)
x(2) = x(2) + originc(2)

do j = 1,3
    js(j) = floor(x(j)/delta)+1
enddo
jx = js(1)
jy = js(2)
jz = js(3)

! increase counter
if(indexvolx(jx,jy,jz).eq.0) then

 if(ncha1.eq.maxvolx) then
   write(stdout,*) 'channel: increase maxvolx'
   stop
 endif

 ncha1 = ncha1 + 1

 indexvolx(jx,jy,jz) = ncha1
 p1(ncha1,1)=jx
 p1(ncha1,2)=jy
 p1(ncha1,3)=jz
endif

volxx1(jx,jy,jz) =  volxx1(jx,jy,jz) + 1.0
volx1(indexvolx(jx,jy,jz)) = volx1(indexvolx(jx,jy,jz)) + 1.0
com1(indexvolx(jx,jy,jz),:) = com1(indexvolx(jx,jy,jz),:) + x(:)
sumvolx1 = sumvolx1 + 1.0

enddo ! jjjt
enddo ! jjjz

do i = 1, ncha1
com1(i,:) = com1(i,:)/volx1(i)

! Moves the position of the first segment lseg/2 away from the surface to prevent collision due to round errors.

com1(i,1) = com1(i,1) - 0.5*lseg*((com1(i,1)-originc(1)))/rchannel 
com1(i,2) = com1(i,2) - 0.5*lseg*((com1(i,2)-originc(2)))/rchannel 
enddo
end

