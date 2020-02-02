subroutine px2csr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine transforms the px,py,pz matrices into a csr format 
! to use the intel mkl libraries
!
! note that there is a csr matrix per grafting point
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


use mkl
use system
use chainsdat
use MPI
use mparameters_monomer
use maps
use conformations

implicit none
real*8 sparse(dimx*dimy*dimz) ! space matrix for a given conformation
real*8 nse(dimx*dimy*dimz) ! non-sparse element matrix for a given conformation
integer nsp(dimx*dimy*dimz) ! non-sparse position matrix for a given conformation
integer pos, posx, posy, posz
integer nonzero ! number of non-zeros in sparse matrix
integer i, j, k, idx
integer flag
integer jj, ii
integer im

if(rank.eq.0) print*, 'Create CSR matrixes for mkl'

! 0. Allocate array

ALLOCATE (csrm(ncha,N_monomer))
ALLOCATE (promkl(ncha)) 

! 1. Determine non-zeros

do im = 1, N_monomer
do jj = 1, cpp(rank+1)     ! local grafting point
ii = cppini(rank+1)+jj   ! grafting point

  nonzero = 0
  do i = 1, newcuantas(ii)
  sparse = 0.0
  do j = 1, long
  if(segtype(j).eq.im) then
     posx = px(i,j,jj)
     posy = py(i,j,jj)
     posz = pz(i,j,jj)
     pos = imap(posx,posy,posz)
   if(sparse(pos).eq.0.0) then
     nonzero = nonzero + 1 
   endif
    sparse(pos) = sparse(pos) + 1.0
   endif ! segtype(j) = im
   enddo ! j
   enddo ! i

   csrm(ii,im)%nonzeros = nonzero

!2. Allocate arrays for csr format

   ALLOCATE (csrm(ii,im)%inc_values(nonzero))
   ALLOCATE (csrm(ii,im)%inc_columns(nonzero))
   ALLOCATE (csrm(ii,im)%pntrb(newcuantas(ii)))
   ALLOCATE (csrm(ii,im)%pntre(newcuantas(ii)))

   ALLOCATE (promkl(ii)%pro(newcuantas(ii)))
   ALLOCATE (promkl(ii)%lnpro(newcuantas(ii)))

!3. Dump px,py,pz into csr arrays

  csrm(ii,im)%gidx = 0 ! global index
  do i = 1, newcuantas(ii)
  nse = 0.0 
  nsp = 0
  idx = 0

  do j = 1, long
  if(segtype(j).eq.im) then
   flag = 0

! get position on the grid
   posx = px(i,j,jj)
   posy = py(i,j,jj)
   posz = pz(i,j,jj)
   pos = imap(posx,posy,posz)

   do k = 1, idx 
     if(nsp(k).eq.pos) then ! position already saved
       nse(k) = nse(k) + 1.0 ! increment by one
       flag = 1 ! will not increase index
       exit ! exit k loop
     endif
   enddo ! k
   if(flag.eq.0) then ! increase index
    idx = idx + 1
    nse(idx) = 1.0
    nsp(idx) = pos
   endif 
 endif ! segtype(j) = im
 enddo ! j

! nse and nsp now contains the element and position of conformation i

 csrm(ii,im)%pntrb(i) = csrm(ii,im)%gidx + 1 ! position of first element in the row in global list  


 do j = 1, idx ! values to add
  csrm(ii,im)%gidx = csrm(ii,im)%gidx + 1 ! increase total counter by one
  csrm(ii,im)%inc_values(csrm(ii,im)%gidx) = nse(j) ! store value
  csrm(ii,im)%inc_columns(csrm(ii,im)%gidx) = nsp(j) ! store position (1...dimx*dimy)
 enddo ! j
 enddo ! i

csrm(ii,im)%pntre(1:newcuantas(ii)-1) = csrm(ii,im)%pntrb(2:newcuantas(ii))
csrm(ii,im)%pntre(newcuantas(ii)) =  csrm(ii,im)%gidx + 1

enddo ! im
enddo ! jj

deallocate (px, py, pz)  ! free mem space
end subroutine


subroutine calc_mkl(xpot)

!
! This routine calculates avpol from xpot using MKL's csr format
!
!
use mkl
use system
use chainsdat
use MPI
use fields_fkfun
use conformations
use molecules
use ematrix
use maps
use mparameters_monomer
use kaist
use results 

implicit none
integer ii,k,kx,ky,i, j, kk, jj
integer err
CHARACTER matdescra(6)

integer im,ix,iy,iz
real*8 avpol_tosend(dimx*dimy*dimz, N_monomer), avpol_tmp(dimx*dimy*dimz)
real*8 avpolmkl(dimx*dimy*dimz, N_monomer)
real*8 xpot(dimx, dimy, dimz,N_monomer), lnxpot(dimx*dimy*dimz, N_monomer)
real*8 q_tosend

matdescra(1) = "G"
matdescra(4) = "F"

avpol_tosend = 0.0

! 1. Reindex xpot and take log 

do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimz
do im = 1, N_monomer
lnxpot(imap(ix,iy,iz),im)=log(xpot(ix,iy,iz,im))
enddo
enddo
enddo
enddo

! 2. Get lnpro from xpotpos and CSR compressed matrix

do im = 1, N_monomer
do jj = 1, cpp(rank+1)     ! local grafting point
ii = cppini(rank+1)+jj   ! grafting point

 call mkl_dcsrmv('N', newcuantas(ii), dimx*dimy*dimz, 1.0d+0, matdescra, &
  csrm(ii,im)%inc_values(:), csrm(ii,im)%inc_columns(:), csrm(ii,im)%pntrb(:), &
  csrm(ii,im)%pntre(:),lnxpot(:,im), 0.0d+0, promkl(ii)%lnpro(:))
  
promkl(ii)%pro(:) = exp(promkl(ii)%lnpro(:))

! 3. Recover q 

q(ii) = sum(promkl(ii)%pro(:))

! 4.  Calculate avpol_tosend
 
 avpol_tmp = 0.0
 call mkl_dcsrmv('T', newcuantas(ii), dimx*dimy*dimz, 1.0d+0, matdescra, csrm(ii,im)%inc_values(:), &
 csrm(ii,im)%inc_columns(:), csrm(ii,im)%pntrb(:), csrm(ii,im)%pntre(:), & 
 promkl(ii)%pro(:), 0.0d+0, avpol_tmp(:)) ! already normalized

avpol_tosend(:,im) = avpol_tosend(:,im) + &
avpol_tmp(:) * vpol*vsol/(delta**3)/fvmkl(:)*sc/q(ii)*ngpol(ii)

enddo ! jj
enddo ! im

  call MPI_Barrier(MPI_COMM_WORLD, err)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5.  MPI 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jefe
if (rank.eq.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpolmkl, dimx*dimy*dimz*N_monomer, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
endif

! Subordinados
if(rank.ne.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpolmkl, dimx*dimy*dimz*N_monomer, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
endif

! 6. Calculate avpol

do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimz
avpol(ix,iy,iz,:) = avpolmkl(imap(ix,iy,iz),:)
enddo
enddo
enddo

end subroutine







