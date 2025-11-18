
subroutine creador

    call creadorA    
    call creadorB

end subroutine creador  

subroutine creadorA

    use const, only : stdout, seed
    use chainsdat, only : cuantasA, newcuantasA, nsegA, lsegA, inA1, ingA, readchains   
    use chainsdat, only : cpp, cppini, ncha
    use MPI, only : rank 
    use branches, only : indexncha, branched

    implicit none
    
    ! local arguments 
    integer :: i,il,ll, j, ii, jj
    real*8 :: chains(3,200,100), gauches(100)
    integer :: iglobal
    integer :: nchas

    indexncha = 1

    newcuantasA = 0

    if((readchains.eq.-1).and.(rank.eq.0)) then
        write(stdout,*) 'creador:', 'Saving conformations...'
        open(file='cadenasA.dat',unit=3113)
    endif

    if(readchains.eq.1) then
        if(rank.eq.0) write(stdout,*) 'creador:', 'Loading conformations'
        open(file='cadenasA.dat',unit=3113)
        do i = 1, cuantasA
            do j = 1, nsegA
                read(3113,*)inA1(j,1),inA1(j,2),inA1(j,3)
            enddo
            call pxsA
            !write(stdout,*) 'creador:', i, newcuantas(1)
        enddo

        ! do il = 1, ncha
        ! write(stdout,*) 'creador:', newcuantas(il)
        ! enddo
        ! stop

        close(3113)
        
        return ! == alternate return : chain generation by reading 
    
    endif

    il=0
    iglobal=1

    do while (il.lt.cuantasA)

        select case (branched)
        case (0)
            call cadenas72mr(chains,nchas,gauches,nsegA,lsegA)
        case (1)
            call cadenas_b(chains,nchas,gauches) ! branched chains
        case (2)
            call cadenas_b2(chains,nchas,gauches) ! branched chains type 2
        end select

        do i=1,nchas
            il=il+1
            if(il.gt.cuantasA) goto 100
            ingA = gauches(i)
            do j=1,nsegA
                inA1(j,2)=chains(2,j,i)
                inA1(j,3)=chains(3,j,i)
                inA1(j,1)=chains(1,j,i)
                if((readchains.eq.-1).and.(rank.eq.0)) write(3113,*)inA1(j,1),inA1(j,2),inA1(j,3)
            enddo
            !print*,"il=",il
            call pxsA 
        enddo
    enddo

    do il = 1, ncha
        write(stdout,*) 'creadorA:', rank, newcuantasA(il)
    enddo
    ! stop

    if((readchains.eq.-1).and.(rank.eq.0)) close(3113)


100    do jj = 1, cpp(rank+1)
            ii = cppini(rank+1)+jj
            write(9988,*)rank, ii, newcuantasA(ii)
        enddo

    return

end subroutine creadorA 


subroutine creadorB

    use const, only : stdout, seed
    use chainsdat, only : cuantasB, newcuantasB, nsegB, lsegB, inB1, ingB, readchains   
    use chainsdat, only : cpp, cppini, ncha
    use MPI, only : rank 
    use branches, only : indexncha, branched

    implicit none
    
    ! local arguments 
    integer :: i,il,ll, j, ii, jj
    real*8 :: chains(3,200,100), gauches(100)
    integer :: iglobal
    integer :: nchas

    indexncha = 1

    newcuantasB = 0

    if((readchains.eq.-1).and.(rank.eq.0)) then
        write(stdout,*) 'creador:', 'Saving conformations...'
        open(file='cadenasB.dat',unit=3113)
    endif

    if(readchains.eq.1) then
        if(rank.eq.0) write(stdout,*) 'creador:', 'Loading conformations'
        open(file='cadenasB.dat',unit=3113)
        do i = 1, cuantasB
            do j = 1, nsegB
                read(3113,*)inB1(j,1),inB1(j,2),inB1(j,3)
            enddo
            call pxsB
            !write(stdout,*) 'creador:', i, newcuantas(1)
        enddo

        ! do il = 1, ncha
        ! write(stdout,*) 'creador:', newcuantas(il)
        ! enddo
        ! stop

        close(3113)
        
        return ! == alternate return : chain generation by reading 
    
    endif

    il=0
    iglobal=1

    do while (il.lt.cuantasB)

        select case (branched)
        case (0)
            call cadenas72mr(chains,nchas,gauches,nsegB,lsegB)
        case (1)
            call cadenas_b(chains,nchas,gauches) ! branched chains
        case (2)
            call cadenas_b2(chains,nchas,gauches) ! branched chains type 2
        end select

        do i=1,nchas
            il=il+1
            if(il.gt.cuantasB) goto 111
            ingB = gauches(i)
            do j=1,nsegB
                inB1(j,2)=chains(2,j,i)
                inB1(j,3)=chains(3,j,i)
                inB1(j,1)=chains(1,j,i)
                if((readchains.eq.-1).and.(rank.eq.0)) write(3113,*)inB1(j,1),inB1(j,2),inB1(j,3)
            enddo
            call pxsB
        enddo
    enddo

    do il = 1, ncha
        write(stdout,*) 'creadorB:', newcuantasB(il)
    enddo
    ! stop

    if((readchains.eq.-1).and.(rank.eq.0)) close(3113)


 111    do jj = 1, cpp(rank+1)
            ii = cppini(rank+1)+jj
            write(9989,*)rank,ii,newcuantasB(ii)
        enddo

    return

end subroutine creadorB


subroutine cadenas72mr(chains,nchas,gauches,nseg,lseg)

    use const, only : pi , seed    

    implicit none

    real*8, intent(inout) :: chains(3,200,100) 
    integer, intent(inout) :: nchas
    real*8, intent(inout) :: gauches(100)
    integer, intent(in) :: nseg
    real*8, intent(in) :: lseg

    ! local variables 
    
    integer :: i,state,ii,j,ive,jve
    real*8 :: rn,state1,sitheta,cotheta,dista
    real*8 :: siphip,cophip
    character*1 :: test
    real*8 :: m(3,3),mm(3,3),tt(3,3),tp(3,3),tm(3,3)
    real*8 :: x(3),xend(3,200),xendr(3,200)
    real*8 :: rands
    integer :: ng

    sitheta=sin(68.0d0*pi/180.0d0)
    cotheta=cos(68.0d0*pi/180.0d0)
    siphip=sin(120.0d0*pi/180.0d0)
    cophip=cos(120.0d0*pi/180.0d0)
        
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

        do i=3,nseg
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
        do ive=4,nseg
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
            call rota36(xend,xendr,nseg,test)
            if (test.eq.'N')cycle
            nchas=nchas+1
            do j=1,nseg
                chains(1,j,nchas)=xendr(1,j)
                chains(2,j,nchas)=xendr(2,j)
                chains(3,j,nchas)=xendr(3,j)
            enddo
            gauches(nchas) = ng
            if (nchas.eq.25) exit
        enddo   
    enddo

end subroutine cadenas72mr

subroutine rota36(xend,xendr,n,test)
      
    use system
    use const

    implicit none      
    
    real*8, intent(in) :: xend(3,200)
    real*8, intent(inout) :: xendr(3,200)
    character*1, intent(inout) :: test
    integer, intent(in) :: n

    real*8 :: fac,fac1,fac2,sbe,cbe,sal,cal,sga
    real*8 :: a,b,c
    real*8 :: alfa, cga
    integer :: i
    real*8 :: gama2
    real*8 :: rands

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

end subroutine rota36 
      
subroutine mrrrr(a,b,c)

    implicit none
    real*8, intent(in) :: a(3,3), b(3,3) 
    real*8, intent(inout) :: c(3,3)
    integer :: i,j,k

    do i=1,3
        do j=1,3
            c(i,j)=0.0d0
        enddo
    enddo

    do i=1,3
        do j=1,3
            do k=1,3
                c(i,j)=c(i,j)+a(i,k)*b(k,j)
            enddo
        enddo
    enddo

end subroutine mrrrr
        
subroutine graftpoints

    use system
    use chainsdat, only : posicionA, posicionB, ngpol, cpp, cppini,  maxcpp, ncha 
    use chainsdat, only : hasGraftA, hasGraftB
    use const, only : stdout, verbose 
    use MPI, only : rank, size
    use ematrix, only : comA, comB, volx
    use graftpoint, only : read_pattern_grafts
    
    implicit none
    
    integer :: ix,iy,iz,j
    integer :: i
    integer :: sumpolAseg    ! total number graft points with A type polymer
    integer :: sumpolBseg    ! total number graft points with B type polymer
    integer :: info


    sumpolAseg = 0 
    sumpolBseg = 0 

    cpp = 0

    print*,"graftpoint ncha=",ncha
    call allocatencha

    do i = 1, ncha                                   ! ncha = number of graft point  
        ngpol(i) = volx(i)                           ! number of polymer 
        posicionA(i, :) = comA(i,:)                    ! position real in space  
        posicionB(i, :) = comB(i,:)                    ! position real in space  
        write(123,*) posicionA(i, :)
        write(124,*) posicionB(i, :) 
        cpp(mod(i,size)+1) = cpp(mod(i,size)+1) + 1  ! = distribution of graft point per processore 
    enddo

    maxcpp = maxval(cpp)
 
    cppini(1) = 0                        ! book keeping of where graft point 
    do j = 2,size                        ! cppini(2) first of the graft chain at second 
        cppini(j)=cppini(j-1)+cpp(j-1)
    enddo

    ! == init hasgraftA and hasgraftB 

    call read_pattern_grafts(ncha,hasGraftA,hasGraftB,info)

    do i=1,ncha
        if(hasGraftA(i)) sumpolAseg = sumpolAseg+1  !== total number of A polymr graft points
        if(hasGraftB(i)) sumpolBseg = sumpolBseg+1  !== total number of B polymr graft points
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
            write(stdout,*) 'creador: number of A polymers in system =', sumpolAseg 
            write(stdout,*) 'creador: number of B polymers in system =', sumpolBseg 
        endif
    endif

    call allocatecpp

end subroutine graftpoints


