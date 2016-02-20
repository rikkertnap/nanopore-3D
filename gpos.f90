implicit none
real*8 x1,x2,y1,y2,z1,z2
real*8 lx, ly
integer n
real*8 dd,dx,dy,dz
integer flag

lx = 10.0
ly = 10.0

open(file='pos001.dat', unit=10)
open(file='pos002.dat', unit=20)
open(file='dist12.dat', unit=999)

do while (.true.)

read(10, *, IOSTAT=flag)n,x1,y1,z1
read(20, *, IOSTAT=flag)n,x2,y2,z2

if(flag.lt.0)exit

dx = abs(x1-x2)
dx = mod(dx+lx, lx)

dy = abs(y1-y2)
dx = mod(dy+ly, ly)

dz = z1-z2

dd = sqrt(dx**2 + dy**2)
!dd = sqrt(dx**2 + dy**2 + dz**2)

write(stdout,*) dd

write(999,*)n,dd
flush(999)

enddo
close(30)
end
