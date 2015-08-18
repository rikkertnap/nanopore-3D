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

call make_ellipsoid ! update matrixes for all particles

!cutarea = 0.01 ! throw away cells that have less area than cutarea x area of the cell with largest area  
cutarea = 0.0 ! throw away cells that have less area than cutarea x area of the cell with largest area  

sumpolseg = 0.0

! clear all
voleps = 0.0
volprot = 0.0
volq = 0.0
volx = 0.0
com = 0.0

do j = 1, NNN

! rotate ellipsoid matrixes according to current rotation matrix

 call rotv(AAA(:,:,j), rotmatrix(:,:,j))
 call rotv(AAAS(:,:,j), rotmatrix(:,:,j))
 call rotv(AAAL(:,:,j), rotmatrix(:,:,j))
 call rotv(AAAX(:,:,j), rotmatrix(:,:,j))

 call rotvo(orient(:,j), rotmatrix(:,:,j))

 npoints = 50

 flag = .false.

 call integrate(AAAL(:,:,j),AellL(:,j), Rell(:,j),npoints, voleps1 ,flag)
 flag = .false. ! not a problem if eps lays outside boundaries
 call integrate(AAA(:,:,j),Aell(:,j), Rell(:,j),npoints, volprot1, flag)

 call integrate(AAAS(:,:,j),AellS(:,j), Rell(:,j),npoints, volq1, flag)

 call integrateg(AAA(:,:,j),Aell(:,j),AAAX(:,:,j),AellX(:,j), Rell(:,j),npoints, volx1, com)

! do ix = 1, dimx
! do iy = 1, dimy
! do iz = 1, dimz
! if(volx(ix,iy,iz).ne.0.0)print*, 'ellipsoid:', ix,iy,iz, volx(ix,iy,iz),com(ix,iy,iz,:)
! enddo
! enddo
! enddo
! stop

!! volume
 temp = 4.0/3.0*pi*Aell(1,j)*Aell(2,j)*Aell(3,j)/(sum(volprot1)*delta**3) ! rescales volume
 volprot1 = volprot1*temp                                                 ! OJO: transformation should mantain cell volumen

!! eps
 voleps1 = voleps1-volprot1
 voleps1 = voleps1*eeps(j)

!! charge
 volq1 = volprot1-volq1
 temp = sum(volq1)
 volq1 = volq1/temp*echarge(j)/(delta**3) ! sum(volq) is echarge

!! grafting
 temp = sum(volx1)

 pnumber = 1.6075

 area = (Aell(1,j)*Aell(2,j))**pnumber 
 area = area+(Aell(1,j)*Aell(3,j))**pnumber 
 area = area+(Aell(2,j)*Aell(3,j))**pnumber 
 area = 4.0*pi*(area/3.0)**(1.0/pnumber) ! approximate (< 1% error) area of elipsoid, see wikipedia

 temp2 = maxval(volx1)

 where(volx1<temp2*cutarea) ! remove cells with very little area
 volx1 = 0.0
 end where 

 volx1 = volx1/temp*area*sigma(j)

 maxss = 1.0d100

 do ix = 1, dimx
 do iy = 1, dimy
 do iz = 1, dimz

 vvtemp = (1.0-volprot1(ix,iy,iz))*(delta**3)/vpol/vsol
 sstemp = volx1(ix,iy,iz)/sigma(j)
 if(sstemp.eq.0.0) then
 vvtemp = 1.0d100
 else
! print*, 'ellipsoid:', ix,iy,iz,(1.0-volprot1(ix,iy,iz)),volx1(ix,iy,iz)
 vvtemp = vvtemp/sstemp
 endif

 if(vvtemp.lt.maxss)maxss=vvtemp
 enddo
 enddo
 enddo

 if(rank.eq.0)print*, 'ellipsoid:', 'Maxsigma for ', j,'is ', maxss

! print*, 'ellipsoid:', j, area,sigma(j)
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
 volx = volx + volx1

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
print*, 'ellipsoid:', 'update_matrix: Total volumen real space= ', temp
print*, 'ellipsoid:', 'update_matrix: Total discretized volumen =', sum(volprot)*delta**3
print*, 'ellipsoid:', 'number of monomers in system =', sumpolseg 
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
call savetodisk(volx, title, counter)
end subroutine

subroutine integrateg(AAA,Aell,AAAX, AellX, Rell, npoints,volprot, com)
! call integrateg(AAA(:,:,j),Aell(:,j),AAAX(:,:,j),AellX(:,j), Rell(:,j),npoints, volq1, com)
use system
use transform

implicit none
integer npoints
real*8 AAA(3,3), AAAX(3,3)
real*8 volprot(dimx,dimy,dimz)
real*8 com(dimx,dimy,dimz,3)
real*8 Rell(3), Aell(3), AellX(3)
real*8 dr(3), dxr(3)
integer ix,iy,iz,ax,ay,az
real*8 vect, vectx
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

real*8 com1(3), volprot1
volprot = 0.0

maxAell = max(AellX(1),AellX(2),AellX(3)) ! maximum lenght/2.0 of a box enclosing the ellipsoid in Cartesian Coordinates

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
vectx = mmmult(dxr,AAAX)

if((vect.ge.1.0).and.(vectx.le.1.0)) then           ! between ellipsoid 1 and 2
  flagin=.true.
  if(flagout.eqv..true.) then

    if ((jx.lt.1).or.(jx.gt.dimx)) then
         print*, 'ellipsoid:','update_matrix: ix', ix
         stop
    endif
    if ((jy.lt.1).or.(jy.gt.dimy)) then
         print*, 'ellipsoid:','update_matrix: iy', iy
         stop
    endif
    if ((jz.lt.1).or.(jz.gt.dimz)) then
         print*, 'ellipsoid:','update_matrix: iz', iz
         stop
    endif

      call intcellg(AAA,AAAX,Rell,ix,iy,iz,npoints,volprot1,com1)
      volprot(jx,jy,jz) = volprot1
      com(jx,jy,jz,:) = com1(:)

      goto 999 ! one in and one out, break the cycle
  endif
else 
  flagout=.true.
  if(flagin.eqv..true.) then

    if ((jx.lt.1).or.(jx.gt.dimx)) then
         print*, 'ellipsoid:','update_matrix: ix', ix
         stop
    endif
    if ((jy.lt.1).or.(jy.gt.dimy)) then
         print*, 'ellipsoid:','update_matrix: iy', iy
         stop
    endif
    if ((jz.lt.1).or.(jz.gt.dimz)) then
         print*, 'ellipsoid:','update_matrix: iz', iz
         stop
    endif

    call intcellg(AAA,AAAX,Rell,ix,iy,iz,npoints,volprot1,com1)
    volprot(jx,jy,jz) = volprot1
    com(jx,jy,jz,:) = com1(:)

    goto 999 ! one in and one out, break the cycle
  endif
endif

enddo
enddo
enddo

if((flagin.eqv..true.).and.(flagout.eqv..false.)) then 

    if ((jx.lt.1).or.(jx.gt.dimx)) then
         print*, 'ellipsoid:','update_matrix: ix', ix
         stop
    endif
    if ((jy.lt.1).or.(jy.gt.dimy)) then
         print*, 'ellipsoid:','update_matrix: iy', iy
         stop
    endif
    if ((jz.lt.1).or.(jz.gt.dimz)) then
         print*, 'ellipsoid:','update_matrix: iz', iz
         stop
    endif

    volprot(jx,jy,iz)=1.0 ! all inside
endif
999 continue

enddo
enddo
enddo

end subroutine



subroutine integrate(AAA,Aell, Rell, npoints,volprot, flag)
use system
use transform

implicit none
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



volprot = 0.0

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

    if ((jx.lt.1).or.(jx.gt.dimx)) then
         print*, 'ellipsoid:','update_matrix: ix', ix
         stop
    endif
    if ((jy.lt.1).or.(jy.gt.dimy)) then
         print*, 'ellipsoid:','update_matrix: iy', iy
         stop
    endif
    if ((jz.lt.1).or.(jz.gt.dimz)) then
         print*, 'ellipsoid:','update_matrix: iz', iz
         stop
    endif


         volprot(jx,jy,jz) = intcell(AAA, Rell, ix,iy,iz, npoints)

      goto 999 ! one in and one out, break the cycle
  endif
else 
  flagout=.true.
  if(flagin.eqv..true.) then

    if ((jx.lt.1).or.(jx.gt.dimx)) then
         print*, 'ellipsoid:','update_matrix: ix', ix
         stop
    endif
    if ((jy.lt.1).or.(jy.gt.dimy)) then
         print*, 'ellipsoid:','update_matrix: iy', iy
         stop
    endif
    if ((jz.lt.1).or.(jz.gt.dimz)) then
         print*, 'ellipsoid:','update_matrix: iz', iz
         stop
    endif

         volprot(jx,jy,jz) = intcell(AAA, Rell, ix,iy,iz, npoints)

    goto 999 ! one in and one out, break the cycle
  endif
endif

enddo
enddo
enddo

if((flagin.eqv..true.).and.(flagout.eqv..false.)) then 

    if ((jx.lt.1).or.(jx.gt.dimx)) then
         print*, 'ellipsoid:','update_matrix: ix', ix
         stop
    endif
    if ((jy.lt.1).or.(jy.gt.dimy)) then
         print*, 'ellipsoid:','update_matrix: iy', iy
         stop
    endif
    if ((jz.lt.1).or.(jz.gt.dimz)) then
         print*, 'ellipsoid:','update_matrix: iz', iz
         stop
    endif

         volprot(jx,jy,jz)=1.0 ! all inside
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


subroutine intcellg(AAA,AAAX,Rell,ix,iy,iz,n,vp,com)
use system
use transform

implicit none
real*8 AAA(3,3)
real*8 AAAX(3,3)
real*8 Rell(3)
integer ix,iy,iz,ax,ay,az
integer cc
real*8 vect, vectX
integer n
real*8 mmmult
real*8 dr(3), dxr(3)
real*8 vp
real*8 com(3)

com = 0.0
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
vectx = mmmult(dxr,AAAX)

if((vect.ge.1.0).and.(vectx.le.1.0)) then
cc=cc+1

dxr = MATMUL(IMAT, dr) 
com(:) = com(:) + dxr(:)
endif


enddo
enddo
enddo

vp = float(cc)/(float(n)**3)
com(:) =com(:)/float(cc)
end













