subroutine cadenas_b2(chains,nchas,gauches)
use chainsdat
use const     
use branches
implicit none
real*8 chains(3,200,100), gauches(100)
integer i,state,statef,ii,j,ive,jve
real*8 rn,state1,sitheta,cotheta,dista
real*8 siphip,cophip
character*1 test
real*8 m(3,3),m1(3,3),mm(3,3),tt(3,3),tp(3,3),tm(3,3)
real*8 x(3),xend(3,200), xendr(3,200)
integer nchas
real*8 rands
integer ng
integer kkk, bbb, lll
real*8 tol

tol = 1.0d-4

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
      
 kkk = 2

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
      
 xend(1,kkk)=lseg+x(1)
 xend(2,kkk)=x(2)
 xend(3,kkk)=x(3)
      
 ng = 0 ! number of trans-bonds


!!! backbone 


 do i=3,longbb

 kkk = kkk + 1
 m1 = m ! store matrix to go back later...

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
         
         xend(1,kkk)=xend(1,kkk-1)+x(1)
         xend(2,kkk)=xend(2,kkk-1)+x(2)
         xend(3,kkk)=xend(3,kkk-1)+x(3)

 if((i.gt.longb(1)).and.(i.le.longb(2))) then ! start branching at long(1) and end at long(2)

 m = m1 ! go back one
 
 state=state+1  ! last chosen state
   if(state.eq.3)state=0 ! sets always the same order (chirality) but can start with anyone

   kkk = kkk + 1
   
   select case (state)
     case (0) 
        call mrrrr(m,tt,mm) ! trans
     case (1)
         call mrrrr(m,tp,mm)
     case (2)
         call mrrrr(m,tm,mm)
   end select

         m=mm

         x(1)=m(1,1)*lseg
         x(2)=m(2,1)*lseg
         x(3)=m(3,1)*lseg

         xend(1,kkk)=xend(1,kkk-2)+x(1)
         xend(2,kkk)=xend(2,kkk-2)+x(2)
         xend(3,kkk)=xend(3,kkk-2)+x(3)

 endif
 enddo       
 
 dista=0.0
 do ive=1,long ! check all possible pairs for self avoiding
         do jve=ive+1,long 
            dista=(xend(1,jve)-xend(1,ive))**(2.0)
            dista=dista+(xend(2,jve)-xend(2,ive))**(2.0)
            dista=dista+(xend(3,jve)-xend(3,ive))**(2.0)
            dista=sqrt(dista)+tol
            if (dista.lt.lseg) then
!               print*, jve, xend(:,jve)
!               print*, ive, xend(:,ive)
!               print*, dista
           
                goto 222
            endif
         enddo
 enddo

 do i=1,50
         test='S'
         call rota36(xend,xendr,long,test)
         if (test.eq.'N')cycle
         nchas=nchas+1
!         call print_ent2(xendr)
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



subroutine print_ent2(xend)

use const
use branches
use chainsdat
implicit none

real*8 xend(3,200)
integer i,j
character*21 filename



! Imprime cadenas en formato ENT

write(filename,'(A6,A1, I3.3, A4)') 'cadena','.', indexncha, '.ent'

open(unit=4400, file=filename)

do i=1, long ! Imprime todo
WRITE(4400,'(A6,I5,A3,I12,A4,F8.3,F8.3,F8.3)') &
"HETATM",i,"  C",i,"    ",xend(1, i)*10,  & 
xend(2, i)*10,xend(3, i)*10
end do

   WRITE((4400),'(A6,I5,I5)')"CONECT", 1, 2
   WRITE((4400),'(A6,I5,I5)')"CONECT", 2, 3

i = 3
do j = 3, longbb ! Une segmentos

 if((j.ge.longb(1)).and.(j.lt.longb(2))) then ! start branching at long(1) and end at long(2)
   WRITE((4400),'(A6,I5,I5)')"CONECT", i, i+1
   WRITE((4400),'(A6,I5,I5)')"CONECT", i, i+2
   i = i + 1
 else
   WRITE((4400),'(A6,I5,I5)')"CONECT", i, i+1
 endif

i = i + 1
end do

WRITE(4400,*)"END"

close(4400)
indexncha = indexncha + 1
if(indexncha.eq.1000)stop
end subroutine

