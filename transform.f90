subroutine inittransf
use const
use transform 
use system
use MPI

implicit none
real*8 beta
real*8 temp(3,3)
integer dimn, j
real*8 x(3),xx(3), vect1(3),vect2(3),vect3(3), vol
real*8 vectc(3)
real*8 fix
real*8 cdivl
real*8, external :: NORMA

gama0 = gama0/180.0*pi

beta = (pi/2.0 - gama0)/2.0
cdivl = cdiva/(sqrt((cos(beta)**2) - (sin(beta)**2)))
fix = 1.0/(cdivl**(1.0/3.0))

MAT(1,1) = cos(beta)/fix
MAT(1,2) = -sin(beta)/fix
MAT(1,3) = 0.0
MAT(2,1) = -sin(beta)/fix
MAT(2,2) = cos(beta)/fix
MAT(2,3) = 0.0

MAT = MAT/(sqrt((cos(beta)**2) - (sin(beta)**2)))

MAT(3,1) = 0.0
MAT(3,2) = 0.0
MAT(3,3) = 1.0/cdivl/fix ! divide by cdivl to keep constant volume


TMAT = TRANSPOSE(MAT)

dimn = 3
temp = MAT

call inverse(temp,IMAT,dimn)

xx(1) = delta
xx(2) = 0.0
xx(3) = 0.0

x = MATMUL(IMAT,xx) 
vect1 = x

xx(2) = delta
xx(1) = 0.0
xx(3) = 0.0

x = MATMUL(IMAT,xx) !
vect2 = x

xx(3) = delta
xx(2) = 0.0
xx(1) = 0.0

x = MATMUL(IMAT,xx) ! to real space 
vect3 = x

call cross_product(vect2,vect3,vectc)

vol = DOT_PRODUCT(vect1,vectc)

if(rank.eq.0)write(stdout,*) 'transform:', 'Volume of a lattice cell in real space', vol, 'nm^3'
if(rank.eq.0)write(stdout,*) 'transform:', 'Volume of a lattice cell in transformed space', delta**3, 'nm^3'

xx(1) = delta*dfloat(dimx)
xx(2) = 0.0
xx(3) = 0.0

x = MATMUL(IMAT,xx) 
vect1 = x

xx(2) = delta*dfloat(dimy)
xx(1) = 0.0
xx(3) = 0.0

x = MATMUL(IMAT,xx) !
vect2 = x

xx(3) = delta*dfloat(dimz)
xx(2) = 0.0
xx(1) = 0.0

x = MATMUL(IMAT,xx) ! 
vect3 = x

call cross_product(vect2,vect3,vectc)

vol = DOT_PRODUCT(vect1,vectc)

if(rank.eq.0)write(stdout,*) 'transform:', 'a / nm ', NORMA(vect1)
if(rank.eq.0)write(stdout,*) 'transform:', 'b / nm ', NORMA(vect2)
if(rank.eq.0)write(stdout,*) 'transform:', 'c / nm ', NORMA(vect3)
if(rank.eq.0)write(stdout,*) 'transform:', 'gama ', gama0*180/3.14159 
if(rank.eq.0)write(stdout,*) 'transform:', 'c/a', NORMA(vect3)/NORMA(vect1)
if(rank.eq.0)write(stdout,*) 'transform:', 'c/b', NORMA(vect3)/NORMA(vect2)
if(rank.eq.0)write(stdout,*) 'transform:', 'cell volume ', vol, 'nm^3'

end subroutine

double precision function NORMA(x)
implicit none
real*8 x(3)
NORMA = (x(1)**2+x(2)**2+x(3)**2)**(0.5)
end function

subroutine initellpos ! transforms ellipsoid coordinates from fractional to real
use const
use transform
use system
use MPI
use ellipsoid
implicit none
integer j
real*8 vect(3), vect2(3)

do j = 1, NNN

vect(1) = Rellf(1,j) + dx*delta
vect(2) = Rellf(2,j) + dy*delta
vect(3) = Rellf(3,j) + dz*delta
vect2 = MATMUL(IMAT,vect)
Rell(:,j) = vect2(:)

if(rank.eq.0)write(stdout,*) 'transform:','Coordinates of particle',j,' in real space ', Rell(:,j)
enddo




end subroutine





  subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

subroutine cross_product(a, b, c)
    real*8, dimension(3) :: c
    real*8, dimension(3), intent(in) :: a, b
 
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - b(1)*a(2)
end subroutine cross_product
