
subroutine readinput
use molecules
use const
use bulk
use MPI
use ellipsoid
use chainsdat
use inputtemp
use transform
use system
implicit none
character basura
integer j


!!!!!!!!!!!!!!!!!!!!!!!!!
! Read input variables
!!!!!!!!!!!!!!!!!!!!!!!!!

read(8,*), basura
read(8,*), dimx,dimy,dimz

read(8,*), basura
read(8,*), gama0

read(8,*), basura
read(8, *), delta

read(8,*), basura
read(8,*), long

read(8,*), basura
read(8,*), lseg

read(8,*), basura
read(8,*), vpol0

read(8,*), basura
read(8,*), cuantas

read(8,*), basura
read(8,*), readchains

read(8,*), basura
read(8, *), dielP, dielS

read(8,*), basura
read(8, *), zpol

read(8, *), basura
read(8, *), csalt

read(8, *), basura
read(8, *), pKa

read(8, *), basura
read(8, *), pHbulk

read(8, *), basura
read(8, *), infile

read(8, *), basura
read(8, *), st

read(8, *), basura
read(8, *), randominput

read(8, *), basura
read(8, *), NNN

call allocateell

read(8, *), basura
do j = 1, NNN
read(8, *), Rell(1,j), Rell(2,j), Rell(3,j)
enddo

read(8, *), basura
do j = 1, NNN
read(8, *), Aell(1,j), Aell(2,j), Aell(3,j)
enddo

read(8, *), basura
do j = 1, NNN
read(8, *), rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
read(8, *), rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
read(8, *), rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
enddo

read(8, *), basura
do j = 1, NNN
read(8, *), sigma(j)
enddo

read(8, *), basura
do j = 1, NNN
read(8, *), echarge(j)
enddo

read(8, *), basura
do j = 1, NNN
read(8, *), eeps(j)
enddo

end subroutine



