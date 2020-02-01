subroutine makemaps
use maps
use system
use ematrix

implicit none
integer ix,iy,iz,i

do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimy
imap(ix,iy,iz) = ix+dimx*(iy-1)+dimx*dimy*(iz-1)
mapx(imap(ix,iy,iz))=ix
mapy(imap(ix,iy,iz))=iy
mapz(imap(ix,iy,iz))=iz
enddo
enddo
enddo
end

subroutine calcfv
use system
use ematrix
use maps

implicit none
integer ix,iy,iz,i

do ix = 1, dimx
do iy = 1, dimy
do iz = 1, dimy
fvstd(ix,iy,iz) = 1.0-volprot(ix,iy,iz)
fvmkl(imap(ix,iy,iz))=1.0-volprot(ix,iy,iz)
enddo
enddo
enddo

end

