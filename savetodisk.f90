subroutine savetodisk(array, title, counter)

use system
use transform
use s2d
implicit none

integer ix, iy, iz, jx, jy, jz
real*8 array(dimx, dimy, dimz)
real*8 arraytemp(dimx, dimy, dimz)
real*8 arrayz(dimz)
integer counter
character*5 title
character*6 titlez
character*21 filename, tempc
real*4 singlepres
real*8 x(3), v(3)
integer, external :: PBCSYMI, PBCREFI

!-----  coordenadas -------------------------
! Variables

arraytemp = 0
do iz = 1, dimz
  arraytemp(:,:,iz) = array(:,:,iz)
enddo
 
array = -1000
do iz = 1, dimz
  array(:,:,iz) = arraytemp(:,:,iz)
enddo

! Graba material en crudo

!      write(filename,'(A5, A1, I3.3, A4)') title,'.', counter, '.raw' 
!      open(unit=45, file=filename)
!
!       do ix = 1, dimx
!        do iy = 1, dimy
!          do iz = 1, dimz
!            write(45,*)array(ix,iy,iz)
!          enddo
!        enddo
!      enddo
!      close (45)

! Integra y graba promedios en z

do iz=1,dimz
         arrayz(iz) = 0.0
      do ix=1,dimx
         do iy=1, dimy
           arrayz(iz)=arrayz(iz)+array(ix,iy,iz)
         enddo
      enddo
enddo

print*, title

titlez = title // 'z'
write(filename,'(A6,A1, I3.3, A4)') titlez,'.', counter, '.dat' 
open(unit=45, file=filename)
do iz=1,dimz
         write(45,*)(iz-0.5)*delta,arrayz(iz)/dimx/dimy
enddo
close(45)

! Archivo paraview

if(vtkflag.eq.1) then

      write(filename,'(A5, A1,I3.3, A4)')title,'.', counter, '.vtk'
      open (unit=45, file=filename)
      write(45,'(A)')'# vtk DataFile Version 2.0'
      write(45,'(A)')title
      write(45,'(A)')'ASCII'
      write(45,'(A)')'DATASET STRUCTURED_GRID '
      write(45,'(A, I5, A1, I5, A1, I5)')'DIMENSIONS', scz*dimz+1, ' ', scy*dimy+1, ' ',scx*dimx+1
      write(45,'(A, I8, A)')'POINTS ',(scx*dimx+1)*(scy*dimy+1)*(scz*dimz+1),' float'
      do ix = 0, scx*dimx
        do iy = 0, scy*dimy
          do iz = 0, scz*dimz

          v(1) = ix*delta ! transformed coordinates
          v(2) = iy*delta
          v(3) = iz*delta

          x = MATMUL(IMAT,v)
 
            write(45, *) x(1),'   ',  x(2), '   ', x(3) ! Cartesian coordinates 
          enddo
        enddo
      enddo

      write(45,'(A, I8)')'CELL_DATA ', scx*scy*scz*dimx*dimy*dimz
      tempc = 'SCALARS ' // title // ' float 1'
      write(45,'(A)')tempc
      write(45,'(A)')'LOOKUP_TABLE default'

       do ix = 1, scx*dimx
        do iy = 1, scy*dimy
          do iz = 1, scz*dimz

            if(PBC(2).ne.3)jx = PBCSYMI(ix,dimx)
            if(PBC(2).eq.3)jx = PBCREFI(ix,dimx)

            if(PBC(4).ne.3)jy = PBCSYMI(iy,dimy)
            if(PBC(4).eq.3)jy = PBCREFI(iy,dimy)

            if(PBC(6).ne.3)jz = PBCSYMI(iz,dimz)
            if(PBC(6).eq.3)jz = PBCREFI(iz,dimz)

            singlepres = array(jx, jy, jz) ! Lo necesito en single presicion
 
            write(45,*)singlepres
          enddo
        enddo
      enddo
      close (45)

endif
return
end


