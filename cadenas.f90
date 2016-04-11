subroutine creador

use system
use const
use system
use chainsdat
use MPI
use branches
implicit none
integer i,il,ll
integer j
real*8 indax, inday, indaz
real*8 chains(3,200,100), gauches(100)
real*8 altx,alty,altz,x(200),y(200),xp(200),yp(200)
real*8 theta,theta1
integer iglobal
integer nchas
integer ii, jj

indexncha = 1

newcuantas = 0

if((readchains.eq.-1).and.(rank.eq.0)) then
write(stdout,*) 'creador:', 'Saving conformations...'
open(file='cadenas.dat',unit=3113)
endif

if(readchains.eq.1) then
 if(rank.eq.0)write(stdout,*) 'creador:', 'Loading conformations'
 open(file='cadenas.dat',unit=3113)
 do i = 1, cuantas
 do j = 1, long
  read(3113,*)in1(j,1),in1(j,2),in1(j,3)
 enddo
 call pxs
! write(stdout,*) 'creador:', i, newcuantas(1)
 enddo

! do il = 1, ncha
! write(stdout,*) 'creador:', newcuantas(il)
! enddo
! stop

close(3113)
 return
endif

il=0
iglobal=1

do while (il.lt.cuantas)

 select case (branched)

 case (0)
  call cadenas72mr(chains,nchas,gauches)
 case (1)
  call cadenas_b(chains,nchas,gauches) ! branched chains
 case (2)
  call cadenas_b2(chains,nchas,gauches) ! branched chains type 2


endselect

  do i=1,nchas
      il=il+1
      if(il.gt.cuantas) goto 100
      ing = gauches(i)
      do j=1,long
         in1(j,2)=chains(2,j,i)
         in1(j,3)=chains(3,j,i)
         in1(j,1)=chains(1,j,i)
         if((readchains.eq.-1).and.(rank.eq.0))write(3113,*)in1(j,1),in1(j,2),in1(j,3)
      enddo
         call pxs
  enddo
enddo

! do il = 1, ncha
! write(stdout,*) 'creador:', newcuantas(il)
! enddo
! stop

 if((readchains.eq.-1).and.(rank.eq.0))close(3113)


 100 do jj = 1, cpp(rank+1)
  ii = cppini(rank+1)+jj
  write(9988,*)ii,newcuantas(ii)
enddo

return
end

subroutine cadenas72mr(chains,nchas,gauches)
use chainsdat
use const     
implicit none
real*8 chains(3,200,100), gauches(100)
integer i,state,ii,j,ive,jve
real*8 rn,state1,sitheta,cotheta,dista
real*8 siphip,cophip
character*1 test
real*8 m(3,3),mm(3,3),tt(3,3),tp(3,3),tm(3,3)
real*8 x(3),xend(3,200),xendr(3,200)
integer nchas
real*8 rands
integer ng

sitheta=sin(68.0*pi/180.0)
cotheta=cos(68.0*pi/180.0)
siphip=sin(120.0*pi/180.0)
cophip=cos(120.0*pi/180.0)
     
nchas=0

do while (nchas.eq.0) 
 x(1)=lseg
 x(2)=0.0
 x(3)=0.0
     
 xend(1,1)=lseg
 xend(2,1)=0.0
 xend(3,1)=0.0
      
 tt(1,1)=cotheta
 tt(1,2)=sitheta
 tt(1,3)=0.0
 tt(2,1)=sitheta
 tt(2,2)=-cotheta
 tt(2,3)=0.0
 tt(3,1)=0.0
 tt(3,2)=0.0
 tt(3,3)=-1.0
      
 tp(1,1)=cotheta
 tp(1,2)=sitheta
 tp(1,3)=0.0
 tp(2,1)=sitheta*cophip
 tp(2,2)=-cotheta*cophip
 tp(2,3)=siphip
 tp(3,1)=sitheta*siphip
 tp(3,2)=-cotheta*siphip
 tp(3,3)=-cophip
      
 tm(1,1)=cotheta
 tm(1,2)=sitheta
 tm(1,3)=0.0
 tm(2,1)=sitheta*cophip
 tm(2,2)=-cotheta*cophip
 tm(2,3)=-siphip
 tm(3,1)=-sitheta*siphip
 tm(3,2)=cotheta*siphip
 tm(3,3)=-cophip
      
 222  rn=rands(seed)
      
 state1=0.0
     
 m(1,1)=cotheta
 m(1,2)=sitheta
 m(1,3)=0.0
      
 m(2,1)=cos(state1)*sitheta
 m(2,2)=-cos(state1)*cotheta
 m(2,3)=sin(state1)
 m(3,1)=sin(state1)*sitheta
 m(3,2)=-sin(state1)*cotheta
 m(3,3)=-cos(state1)
      
 x(1)=m(1,1)*lseg
 x(2)=m(2,1)*lseg
 x(3)=m(3,1)*lseg
      
 xend(1,2)=lseg+x(1)
 xend(2,2)=x(2)
 xend(3,2)=x(3)
      
 ng = 0 ! number of trans-bonds

 do i=3,long
         rn=rands(seed)
         state=int(rn*3)
         if (state.eq.3) then 
            state=2
         endif

         if (state.eq.0) then ! trans
            call mrrrr(m,tt,mm)
          do ii=1,3
          do j=1,3
             m(ii,j)=mm(ii,j)
          enddo
          enddo
            ng = ng + 1
         elseif (state.eq.1) then
            call mrrrr(m,tp,mm)
          do ii=1,3
          do j=1,3
             m(ii,j)=mm(ii,j)
          enddo
          enddo
         elseif (state.eq.2) then
            call mrrrr(m,tm,mm)
          do ii=1,3
          do j=1,3
            m(ii,j)=mm(ii,j)
          enddo
          enddo
         endif
         
         x(1)=m(1,1)*lseg
         x(2)=m(2,1)*lseg
         x(3)=m(3,1)*lseg
         
         xend(1,i)=xend(1,i-1)+x(1)
         xend(2,i)=xend(2,i-1)+x(2)
         xend(3,i)=xend(3,i-1)+x(3)
 enddo       
      
 dista=0.0
 do ive=4,long
         do jve=1,ive-3
            dista=(xend(1,jve)-xend(1,ive))**(2.0)
            dista=dista+(xend(2,jve)-xend(2,ive))**(2.0)
            dista=dista+(xend(3,jve)-xend(3,ive))**(2.0)
            dista=sqrt(dista)
            if (dista.lt.lseg) then
               goto 222
            endif
         enddo
 enddo


 do i=1,300
         test='S'
         call rota36(xend,xendr,long,test)
         if (test.eq.'N')cycle
         nchas=nchas+1
         do j=1,long
            chains(1,j,nchas)=xendr(1,j)
            chains(2,j,nchas)=xendr(2,j)
            chains(3,j,nchas)=xendr(3,j)
         enddo
            gauches(nchas) = ng
         if (nchas.eq.25)exit
 enddo   
enddo
return
end

subroutine rota36(xend,xendr,n,test)
      
use system
use const
implicit none      
real*8 xend(3,200),rands,xendr(3,200)
character*1 test
integer n
real*8 fac,fac1,fac2,sbe,cbe,sal,cal,sga
real*8 a,b,c
real*8 alfa, cga
integer i
real*8 gama2

fac=rands(seed)
fac1=rands(seed)
fac2=rands(seed)
alfa=fac*2*pi
cbe=2.0d0*fac1-1.0d0
gama2=fac2*2*pi

sbe=(1-cbe**2)**0.5
cal=cos(alfa)
sal=sin(alfa)
cga=cos(gama2)
sga=sin(gama2)

do i=1,n
   a=xend(1,i)
   b=xend(2,i)
   c=xend(3,i)

   xendr(1,i)=a*(-cbe*sal*sga+cal*cga)-b*(cbe*sal*cga+cal*sga)+c*sbe*sal
   xendr(2,i)=a*(cbe*cal*sga+sal*cga)+b*(cbe*cal*cga-sal*sga)-c*sbe*cal
   xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe
enddo 

return
end
      
subroutine mrrrr(a,b,c)
implicit none
real*8 a(3,3),b(3,3),c(3,3)
integer i,j,k

do i=1,3
  do j=1,3
     c(i,j)=0
  enddo
enddo

do i=1,3
 do j=1,3
  do k=1,3
  c(i,j)=c(i,j)+a(i,k)*b(k,j)
  enddo
 enddo
enddo

return
end 
        
subroutine graftpoints
use system
use chainsdat
use const
use MPI
use ematrix
implicit none
integer ix,iy,iz,j
integer i

cpp = 0

call allocatencha

do i = 1, ncha
  ngpol(i) = volx(i)
  posicion(i, :) = com(i,:)
  cpp(mod(i,size)+1) = cpp(mod(i,size)+1) + 1
enddo

maxcpp = maxval(cpp)

cppini(1) = 0
do j = 2,size
cppini(j)=cppini(j-1)+cpp(j-1)
enddo

if(rank.eq.0) then
if(verbose.ge.1) then
 write(stdout,*) 'creador:', 'graftingpoints:'
 write(stdout,*) 'creador:', 'ncha = ', ncha
 do j = 1, size
 write(stdout,*) 'creador:', ' cpp    ', j, ' = ', cpp(j)
 write(stdout,*) 'creador:', ' cppini ', j, ' = ', cppini(j)+1
 write(stdout,*) 'creador:', ' cppfin ', j, ' = ', cppini(j)+cpp(j)
 write(stdout,*) 'creador:', '!!!!!!!!!!!!!!!!!!!!!!!!!'
 enddo
endif
endif

call allocatecpp

end


