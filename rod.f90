subroutine update_matrix_rod_4(flag) 
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
rchannelL2 = (rchannel + 3*delta)**2 
rchannelS2 = (rchannel - delta)**2 
 
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
 call integrate_c(rchannelL2, RdimZ-3,originc ,npoints, voleps1 , sumvoleps1, flag) 
 
 do iz = 1, RdimZ-3 
  voleps1(:,:,iz) = 1.0 
 enddo 
 do iz = dimz-RdimZ+4, dimz 
  voleps1(:,:,iz) = 1.0 
 enddo 
 voleps1 = 1.0-voleps1
 flag = .false. ! not a problem if eps lays outside boundaries 
 
 call integrate_c(rchannel2,RdimZ, originc,npoints, volprot1, sumvolprot1, flag) 
 do iz = 1, RdimZ 
  volprot1(:,:,iz) = 1.0 
 enddo 
 do iz = dimz-RdimZ+1, dimz 
  volprot1(:,:,iz) = 1.0 
 enddo 
 volprot1 = 1.0-volprot1 
 
 call integrate_c(rchannelS2,RdimZ+1,originc,npoints, volq1, sumvolq1, flag) 
 do iz = 1, RdimZ+1 
  volq1(:,:,iz) = 1.0 
 enddo 
 do iz = dimz-RdimZ, dimz 
  volq1(:,:,iz) = 1.0 
 enddo 
 volq1 = 1.0-volq1 


!rchannel2 = (rchannel+delta)**2 
! grafting positions on channel surface are OK because first segment is not fixed to (0,0,0) 
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
 temp = sum(volq1) 
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
write(stdout,*) 'rod:', 'update_matrix: Total nanorod volumen real space= ', temp 
write(stdout,*) 'rod:', 'update_matrix: Total discretized volumen =', (dimx*dimy*dimz-sum(volprot))*delta**3 
write(stdout,*) 'rod:', 'number of polymers in system =', sumpolseg  
write(stdout,*) 'rod:', 'surface area =', area 
write(stdout,*) 'rod:', 'surface density =', sumpolseg/area 
write(stdout,*) 'rod:', 'surface density expected from input =', & 
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

