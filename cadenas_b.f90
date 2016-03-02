subroutine cadenas_b(chains,nchas,gauches)
use chainsdat
use const     
use branches
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
         call print_ent(xendr)
         do j=1,long
            chains(1,j,nchas)=xendr(1,j)
            chains(2,j,nchas)=xendr(2,j)
            chains(3,j,nchas)=xendr(3,j)
         enddo

          

            gauches(nchas) = ng
         if (nchas.eq.25)exit
 enddo   
enddo

stop
return
end



subroutine print_ent(xend)

use const
use branches
use chainsdat
implicit none

real*8 xend(3,200)
integer i

indexncha = indexncha + 1

! Imprime cadenas en formato ENT

do i=1, long ! Imprime todo
WRITE((4400+indexncha),'(A6,I5,A3,I12,A4,F8.3,F8.3,F8.3)') &
"HETATM",i,"  C",i,"    ",xend(1, i)*10,  & 
xend(2, i)*10,xend(3, i)*10
end do

do i = 1, long ! Une segmentos
WRITE((4000+indexncha),'(A6,I5,I5)')"CONECT", i, i+1
end do

WRITE((4000+indexncha),*)"END"

close(4000+indexncha)
end subroutine

