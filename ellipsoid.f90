subroutine randomvect(V)
use const
implicit none

real*8, external :: rands
real*8 u, w
real*8 V(3)
real*8 theta, phi

u = rands(seed)
w = rands(seed)
theta = 2*pi*u
phi = acos(2.0*w-1.0)
V(1) = cos(theta)*sin(phi)
V(2) = sin(theta)*sin(phi)
V(3) = cos(phi)
end subroutine

subroutine make_ellipsoid
use system
use ellipsoid
implicit none
integer i
integer j
real deltaX

! clear all
AAA = 0.0
AAAS = 0.0
AAAL = 0.0
AAAX = 0.0

deltaX = delta/2.0

! LOOP over particle
do j = 1, NNN

 ! orientation vector
 orient(:,j) = 0
 orient(1,j) = 1.0

 AellL(1,j) = Aell(1,j)+delta
 AellL(2,j) = Aell(2,j)+delta
 AellL(3,j) = Aell(3,j)+delta

 AellS(1,j) = Aell(1,j)-delta
 AellS(2,j) = Aell(2,j)-delta
 AellS(3,j) = Aell(3,j)-delta

 AellX(1,j) = Aell(1,j)+deltaX
 AellX(2,j) = Aell(2,j)+deltaX
 AellX(3,j) = Aell(3,j)+deltaX

 do i = 1,3
 AAA(i,i,j) = 1.0/(Aell(i,j)**2)
 AAAS(i,i,j) = 1.0/(AellS(i,j)**2)
 AAAL(i,i,j) = 1.0/(AellL(i,j)**2)
 AAAX(i,i,j) = 1.0/(AellX(i,j)**2)
 enddo

enddo
end subroutine

subroutine update_matrix(flag)
use system
use ellipsoid
use ematrix
use MPI
use const
use chainsdat
use molecules
implicit none
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


call make_ellipsoid ! update matrixes for all particles

!cutarea = 0.01 ! throw away cells that have less area than cutarea x area of the cell with largest area  
cutarea = 0.0 ! throw away cells that have less area than cutarea x area of the cell with largest area  

sumpolseg = 0.0

! clear all
voleps = 0.0
volprot = 0.0
volq = 0.0
volx = 0.0
volxx = 0.0
com = 0.0
ncha = 0

do j = 1, NNN

! rotate ellipsoid matrixes according to current rotation matrix

 call rotv(AAA(:,:,j), rotmatrix(:,:,j))
 call rotv(AAAS(:,:,j), rotmatrix(:,:,j))
 call rotv(AAAL(:,:,j), rotmatrix(:,:,j))
 call rotv(AAAX(:,:,j), rotmatrix(:,:,j))

 call rotvo(orient(:,j), rotmatrix(:,:,j))

 npoints = 50

 flag = .false.

 call integrate(AAAL(:,:,j),AellL(:,j), Rell(:,j),npoints, voleps1 , sumvoleps1, flag)
 flag = .false. ! not a problem if eps lays outside boundaries
 call integrate(AAA(:,:,j),Aell(:,j), Rell(:,j),npoints, volprot1, sumvolprot1, flag)

 call integrate(AAAS(:,:,j),AellS(:,j), Rell(:,j),npoints, volq1, sumvolq1, flag)

 npoints = 100000000
 call newintegrateg(Aell(:,j),Rell(:,j),npoints,volx1,sumvolx1, com1, p1, ncha1, volxx1)

!! volume
 temp = 4.0/3.0*pi*Aell(1,j)*Aell(2,j)*Aell(3,j)/(sumvolprot1*delta**3) ! rescales volume
 volprot1 = volprot1*temp                                                 ! OJO: transformation should mantain cell volumen
 sumvolprot1 = sumvolprot1*temp

!! eps
 voleps1 = voleps1-volprot1
 voleps1 = voleps1*eeps(j)

!! charge
 volq1 = volprot1-volq1
 temp = sumvolprot1-sumvolq1
 volq1 = volq1/temp*echarge(j)/(delta**3) ! sum(volq) is echarge

!! grafting
 pnumber = 1.6075

 area = (Aell(1,j)*Aell(2,j))**pnumber 
 area = area+(Aell(1,j)*Aell(3,j))**pnumber 
 area = area+(Aell(2,j)*Aell(3,j))**pnumber 
 area = 4.0*pi*(area/3.0)**(1.0/pnumber) ! approximate (< 1% error) area of elipsoid, see wikipedia

 temp2 = maxval(volx1)

 where(volx1<temp2*cutarea) ! remove cells with very little area
 volx1 = 0.0
 end where 

 volx1 = volx1/sumvolx1*area*sigma(j)
 volxx1 = volxx1/sumvolx1*area*sigma(j)

 maxss = 1.0d100

! do i = 1, ncha1
! vvtemp = (1.0-volprot1(ix,iy,iz))*(delta**3)/vpol/vsol
! sstemp = volx1(i)/sigma(j)
! if(sstemp.eq.0.0) then
! vvtemp = 1.0d100
! else
! write(stdout,*) 'ellipsoid:', ix,iy,iz,(1.0-volprot1(ix,iy,iz)),volx1(ix,iy,iz)
! vvtemp = vvtemp/sstemp
! endif
! if(vvtemp.lt.maxss)maxss=vvtemp
! enddo
! if(rank.eq.0)write(stdout,*) 'ellipsoid:', 'Maxsigma for ', j,'is ', maxss

 sumpolseg = sumpolseg + area*sigma(j)*long

!! volume  
 volprot1 = volprot1 * 0.99
 volprot = volprot+volprot1

! CHECK COLLISION HERE...
 if(maxval(volprot).gt.1.0) then ! collision
   flag=.true. 
   exit
 endif
 
 voleps = voleps + voleps1
 volq = volq + volq1 

! add com1 and volx to list

 volxx = volxx + volxx1


do i = 1, ncha1
ncha = ncha+1
volx(ncha)=volx1(i)
com(ncha,:)=com1(i,:)
p0(ncha,:)=p1(i,:)
enddo

enddo
title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)

if (verbose.ge.2) then
temp = 0
do j = 1, NNN
temp = temp + 4.0/3.0*pi*Aell(1,j)*Aell(2,j)*Aell(3,j)
enddo
if (rank.eq.0) then
write(stdout,*) 'ellipsoid:', 'update_matrix: Total volumen real space= ', temp
write(stdout,*) 'ellipsoid:', 'update_matrix: Total discretized volumen =', sum(volprot)*delta**3
write(stdout,*) 'ellipsoid:', 'number of monomers in system =', sumpolseg 
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

subroutine integrate(AAA,Aell, Rell, npoints,volprot,sumvolprot, flag)
use system
use transform
use const
implicit none
real*8 sumvolprot
integer npoints
real*8 AAA(3,3)
real*8 volprot(dimx,dimy,dimz)
real*8 Rell(3), Aell(3)
real*8 dr(3), dxr(3)
integer ix,iy,iz,ax,ay,az
real*8 vect
logical flagin, flagout
real*8 intcell
real*8 mmmult
integer jx,jy, jz
real*8 Rpos(3)
real*8 maxAell
logical flag

real*8 box(4)
real*8 x(3), v(3)
integer xmin,xmax,ymin,ymax,zmin,zmax
integer i,j

logical flagsym
real*8 voltemp

volprot = 0.0
sumvolprot = 0.0 ! total volumen, including that outside the system

maxAell = max(Aell(1),Aell(2),Aell(3)) ! maximum lenght/2.0 of a box enclosing the ellipsoid in Cartesian Coordinates

Rpos(1) = Rell(1)
Rpos(2) = Rell(2)
Rpos(3) = Rell(3)

! create a box in transformed coordinate enclosing the ellipsoid
!!!!!!!!!!!!!!!! xmin !!!!!!!!!!!!!!!
x(1) = Rpos(1)-maxAell
do j = 1, 2
do i = 1, 2
x(2) = Rpos(2)+maxAell*(-1.0)**j
x(3) = Rpos(3)+maxAell*(-1.0)**i
v = MATMUL(MAT,x)
box(i+2*(j-1)) = v(1)
enddo
enddo
xmin = int(minval(box)/delta)-2
!!!!!!!!!!!!!! xmax !!!!!!!!!!!!!!!!!!!!!1
x(1) = Rpos(1)+maxAell
do j = 1, 2
do i = 1, 2
x(2) = Rpos(2)+maxAell*(-1.0)**j
x(3) = Rpos(3)+maxAell*(-1.0)**i
v = MATMUL(MAT,x)
box(i+2*(j-1)) = v(1)
enddo
enddo
xmax = int(maxval(box)/delta)+2
!!!!!!!!!!!!!!!! ymin !!!!!!!!!!!!!!!
x(2) = Rpos(2)-maxAell
do j = 1, 2
do i = 1, 2
x(1) = Rpos(1)+maxAell*(-1.0)**j
x(3) = Rpos(3)+maxAell*(-1.0)**i
v = MATMUL(MAT,x)
box(i+2*(j-1)) = v(2)
enddo
enddo
ymin = int(minval(box)/delta)-2
!!!!!!!!!!!!!!!! ymax !!!!!!!!!!!!!!!
x(2) = Rpos(2)+maxAell
do j = 1, 2
do i = 1, 2
x(1) = Rpos(1)+maxAell*(-1.0)**j
x(3) = Rpos(3)+maxAell*(-1.0)**i
v = MATMUL(MAT,x)
box(i+2*(j-1)) = v(2)
enddo
enddo
ymax = int(maxval(box)/delta)+2
!!!!!!!!!!!!!!!! zmin !!!!!!!!!!!!!!!
x(3) = Rpos(3)-maxAell
do j = 1, 2
do i = 1, 2
x(1) = Rpos(1)+maxAell*(-1.0)**j
x(2) = Rpos(2)+maxAell*(-1.0)**i
v = MATMUL(MAT,x)
box(i+2*(j-1)) = v(3)
enddo
enddo
zmin = int(minval(box)/delta)-2
!!!!!!!!!!!!!!!! zmax !!!!!!!!!!!!!!!
x(3) = Rpos(3)+maxAell
do j = 1, 2
do i = 1, 2
x(1) = Rpos(1)+maxAell*(-1.0)**j
x(2) = Rpos(2)+maxAell*(-1.0)**i
v = MATMUL(MAT,x)
box(i+2*(j-1)) = v(3)
enddo
enddo
zmax = int(maxval(box)/delta)+2

! Make a list of the cells that have no ellipsoid, those that have part ellipsoid and those that have full ellipsoid
! Consider boundary conditions 

do ix = xmin, xmax
do iy = ymin, ymax
do iz = zmin, zmax

jx=ix
jy=iy
jz=iz

if(PBC(1).eq.1) then
 jx=mod(ix+dimx-1,dimx)+1
endif
if(PBC(3).eq.1) then
 jy=mod(iy+dimy-1,dimy)+1
endif
if(PBC(5).eq.1) then
 jz=mod(iz+dimz-1,dimz)+1
endif

flagin=.false.
flagout=.false.

do ax = -1,0
do ay = -1,0
do az = -1,0

dr(1) = (ix+ax)*delta 
dr(2) = (iy+ay)*delta 
dr(3) = (iz+az)*delta 

! dr is in transformed space, change to cartesian space

dxr = MATMUL(IMAT,dr) - Rell
vect = mmmult(dxr,AAA)

if(vect.le.1.0) then           ! inside the ellipsoid
  flagin=.true.
  if(flagout.eqv..true.) then

  flagsym = .false.

    if (jx.lt.1) then
       if(PBC(1).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: ix', ix
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jy.lt.1) then
       if(PBC(3).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: iy', iy
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jz.lt.1) then
       if(PBC(5).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: iz', iz
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jx.gt.dimx) then
       if(PBC(2).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: ix', ix
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jy.gt.dimy) then
       if(PBC(4).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: iy', iy
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jz.gt.dimz) then
       if(PBC(6).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: iz', iz
         stop
       else
         flagsym = .true.
       endif
    endif

    voltemp =  intcell(AAA, Rell, ix,iy,iz, npoints)
    sumvolprot = sumvolprot + voltemp

    if(flagsym.eqv..false.) then ! cell is not out of system due to reflection symmetry
         volprot(jx,jy,jz) = voltemp
    endif

      goto 999 ! one in and one out, break the cycle
  endif
else 
  flagout=.true.
  if(flagin.eqv..true.) then

  flagsym = .false.

    if (jx.lt.1) then
       if(PBC(1).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: ix', ix
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jy.lt.1) then
       if(PBC(3).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: iy', iy
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jz.lt.1) then
       if(PBC(5).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: iz', iz
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jx.gt.dimx) then
       if(PBC(2).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: ix', ix
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jy.gt.dimy) then
       if(PBC(4).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: iy', iy
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jz.gt.dimz) then
       if(PBC(6).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: iz', iz
         stop
       else
         flagsym = .true.
       endif
    endif

    voltemp = intcell(AAA, Rell, ix,iy,iz, npoints)
    sumvolprot = sumvolprot + voltemp

    if(flagsym.eqv..false.) then ! cell is not out of system due to reflection symmetry
         volprot(jx,jy,jz) = voltemp
    endif

    goto 999 ! one in and one out, break the cycle
  endif
endif

enddo
enddo
enddo

if((flagin.eqv..true.).and.(flagout.eqv..false.)) then 

  flagsym = .false.

    if (jx.lt.1) then
       if(PBC(1).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: ix', ix
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jy.lt.1) then
       if(PBC(3).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: iy', iy
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jz.lt.1) then
       if(PBC(5).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: iz', iz
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jx.gt.dimx) then
       if(PBC(2).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: ix', ix
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jy.gt.dimy) then
       if(PBC(4).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: iy', iy
         stop
       else
         flagsym = .true.
       endif
    endif
    if (jz.gt.dimz) then
       if(PBC(6).ne.3) then
         write(stdout,*) 'ellipsoid:','update_matrix: iz', iz
         stop
       else
         flagsym = .true.
       endif
    endif

         sumvolprot = sumvolprot + 1.0
    if(flagsym.eqv..false.) then ! cell is not out of system due to reflection symmetry
         volprot(jx,jy,jz)=1.0 ! all inside
    endif
endif
999 continue

enddo
enddo
enddo
end subroutine

double precision function intcell(AAA,Rell,ix,iy,iz,n)
use system
use transform

implicit none
real*8 AAA(3,3)
real*8 Rell(3)
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
dxr = MATMUL(IMAT, dr) - Rell
vect = mmmult(dxr,AAA)

if(vect.le.1.0)cc=cc+1

enddo
enddo
enddo

intcell = float(cc)/(float(n)**3)
end function

double precision function mmmult(V,A)
implicit none
real*8 V(3)
real*8 A(3,3)
real*8 C(3)
C(1) = A(1,1)*V(1)+A(1,2)*V(2)+A(1,3)*V(3)
C(2) = A(2,1)*V(1)+A(2,2)*V(2)+A(2,3)*V(3)
C(3) = A(3,1)*V(1)+A(3,2)*V(2)+A(3,3)*V(3)
mmmult = V(1)*C(1) + V(2)*C(2) + V(3)*C(3)
endfunction

subroutine rotv(A, B) ! applies rotation matrix B to ellipsoid matrix A
implicit none
real*8 A(3,3)
real*8 B(3,3)
real*8 BT(3,3)
BT = TRANSPOSE(B)
A = MATMUL(A, B)
A = MATMUL(BT, A)
end subroutine

subroutine rotvo(orient, B) ! applies rotation matrix B to vector orient
implicit none
real*8 B(3,3)
real*8 orient(3)
orient = MATMUL(B, orient)
end subroutine

subroutine rotvm(A, theta, V) ! rotates the rotation matrix A theta degress round V
implicit none
real*8 A(3,3)
real*8 theta
real*8 B(3,3)
real*8 V(3)
B(1,1)=cos(theta)+V(1)*V(1)*(1.0-cos(theta))
B(1,2)=V(1)*V(2)*(1.0-cos(theta))-V(3)*sin(theta)
B(1,3)=V(1)*V(3)*(1.0-cos(theta))+V(2)*sin(theta)
B(2,1)=V(2)*V(1)*(1.0-cos(theta))+V(3)*sin(theta)
B(2,2)=cos(theta)+V(2)*V(2)*(1.0-cos(theta))
B(2,3)=V(2)*V(3)*(1.0-cos(theta))-V(1)*sin(theta)
B(3,1)=V(3)*V(1)*(1.0-cos(theta))-V(2)*sin(theta)
B(3,2)=V(3)*V(2)*(1.0-cos(theta))+V(1)*sin(theta)
B(3,3)=cos(theta)+V(3)*V(3)*(1.0-cos(theta))
A = MATMUL(B, A)
end subroutine


subroutine newintegrateg(Aell,Rell, npoints,volx1,sumvolx1,com1,p1,ncha1,volxx1)
use system
use transform
use chainsdat
use ematrix
use const
implicit none
real*8 sumvolx1
integer npoints
!real*8 AAA(3,3), AAAX(3,3)
integer indexvolx(dimx,dimy,dimz)
integer listvolx(ncha,3)
real*8 Rell(3), Aell(3)
real*8 radio
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

pi=acos(-1.0)

dims(1) = dimx
dims(2) = dimy
dims(3) = dimz

indexvolx = 0
ncha1 = 0
volx1 = 0.0
sumvolx1 = 0.0 ! total volumen, including that outside system
com1 = 0.0
p1 = 0
volxx1 = 0.0

! This routine determines the surface coverage and grafting positions only for spheres
!

! check we have a sphere

radio=Aell(1)
if((radio.ne.Aell(2)).or.(radio.ne.Aell(3))) then 
write(stdout,*) 'newintegrateg: needs spherical particle... stop'
stop
endif

do i = 1, npoints
call randomvect(x) ! uniform points on a sphere

x(:) = x(:)*radio + Rell(:)

! x in real space, v in transformed space
    v = MATMUL(MAT,x)


! PBC 

flagin = 1

do j = 1,3

    is(j) = floor(v(j)/delta)+1
    js(j) = is(j)

select case (PBC((j-1)*2+1))
  case (0 , 2)
    if(is(j).lt.1) then
    write(stdout,*) 'Error in newintegrateg: out of boundary'
    endif
  case (1)
    js(j)=mod(is(j)+dims(j)-1,dims(j))+1
  case (3)
    if(v(j).lt.0.0)flagin=0
endselect

select case (PBC((j-1)*2+2))
  case (0 , 2)
    if(is(j).gt.dims(j)) then
    write(stdout,*) 'Error in newintegrateg: out of boundary'
    endif
  case (1)
    js(j)=mod(is(j)+dims(j)-1,dims(j))+1
  case (3)
    if(v(j).gt.float(dims(j))*delta)flagin=0
endselect
enddo
jx = js(1)
jy = js(2)
jz = js(3)

if(flagin.eq.1) then

! increase counter
if(indexvolx(jx,jy,jz).eq.0) then

 if(ncha1.eq.maxvolx) then
   write(stdout,*) 'ellipsoid: increase maxvolx'
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
endif
sumvolx1 = sumvolx1 + 1.0

enddo ! npoints

do i = 1, ncha1
com1(i,:) = com1(i,:)/volx1(i)
! Moves the position of the first segment lseg/2 away from the surface to prevent collision due to round errors.
com1(i,:) = com1(i,:) + 1.5*lseg*((com1(i,:)-Rell(:)))/Aell(:)
enddo

end





