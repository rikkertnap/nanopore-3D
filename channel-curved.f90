module channelcurved

    implicit none

    real*8 :: radiusL         ! largest radius == rchannelL  in module channel
    real*8 :: radiusS         ! smallest radius == rchannelS in module channel
    real*8 :: radiusC         ! curvature nanochannel 
    real*8 :: lengthchannel 
    real*8 :: Rzmax           ! for evaluation of radiusfz 

contains 

    subroutine set_radii_curved
         use system, only : dimx, dimy
        use channel, only : rchannelL, rchannelS

        Rzmax =  4.0d0*sqrt(1.0d8*(dimx**2+dimy**2)) ! == largest then larger radial radius possible
        radiusL = rchannelL
        radiusS = rchannelS
        lengthchannel = LengthChannelCurvature()
        radiusC = radiusCurvature(radiusL,radiusS,lengthchannel)

    end subroutine  

    function LengthChannelCurvature() result(lenChannel)

        use system, only : dimz, delta 
        use channel, only  : RdimZ      ! Rdimz is of in units of delta 
        
        ! return argument 
        real*8 :: lenChannel
        
        lenChannel = (dimz-2.0d0*Rdimz)*delta

    end function LengthChannelCurvature

    function radiusCurvature(radiusL,radiusS,lengthchannel) result(Rc)

        use const , only : stdout
        use MPI, only : ierr

        real*8, intent(in) :: radiusL, radiusS, lengthchannel
        ! return argument 
        real*8 :: Rc

        ! local argument 
        logical :: errorflag

        Rc = (radiusL-radiusS)/2.0d0 +lengthchannel**2/(8.0d0*(radiusL-radiusS))

        errorflag=.false.
        if(radiusL<radiusS) errorflag=.true.
        if(radiusS<0) errorflag=.true.
        if(radiusL<0) errorflag=.true. 

        if(errorflag) then 
            write(stdout,*) 'radiusCurvature: error in input : rL =',radiusL,' rS =',radiusS
            call MPI_FINALIZE(ierr) ! == end  MPI
            stop
        endif    

       !  print*,"RC=",RC," radiusL=",radiusL," radiusS=",radiusS," Lengthchannel=",lengthchannel

    end function radiusCurvature 

    ! Shape function of curved channel
    ! returns -axis dependent radius of curved channel
    ! input: real*8 z : argement betwween -L/2<=z<=L/2 qwith L lenghtchannel
    ! output : real*8 Rz : in nm
    ! pre radsiusC and radiusL  set with function radiusCurvature             

    function radiusfz(z,RL,RC) result(Rz)

        real*8, intent(in) :: z 
        real*8, intent(in) :: RL,RC
        ! return argument 
        real*8 :: Rz
        
        if(RC**2>z**2) then 
            Rz = sqrt(RC**2-z**2) -(RC-RL)    
        else 
            Rz = Rzmax 
            print*,"warning radiusfz: Rc<=z : z",z, "RC=",RC
        endif

    end function 

    function surface_area_channel(a,b,RL,RC,L) result(area)
        use const, only : pi
          
        real*8, intent(in) :: a,b ! integraton boundary along z-axis origin in middle
        real*8, intent(in) :: RL,RC,L
        ! return argument 
        real*8 :: area

        area= 2.0d0*pi*RC*(b-a)-2.0d0*pi*(RC-RL)*RC*(asin(b/RC) -asin(a/RC) )

    end function 

    function total_surface_area_curv(RL,RC,L) result(areatotal)
        use const, only : pi
        
        real*8, intent(in) :: RL,RC,L
        real*8 :: areatotal
       
        if(RC>RL) then 
            areatotal = 2.0d0*pi*RC*L-4.0d0*pi*(RC-RL)*RC*asin(L/(2.0d0*RC))
        else if(RC==RL) then 
            areatotal = 2.0d0*pi*RC*L
        else    
            areatotal =  0.0d0
            print*,"total_area_curv: arguments out of range."
        endif              

    end function total_surface_area_curv
    
    function total_volume_curv(RL,RC,L) result(voltotal)
        use const, only : pi
        
        real*8, intent(in) :: RL,RC,L
        real*8 :: voltotal

        voltotal = pi*(RC**2+(RC-RL)**2)*L/2.0d0 -(pi/3.0d0)*(L/2.0d0)**3
        voltotal = voltotal - pi*(RC-RL)*( (L/2.0d0)*sqrt(RC**2-L**2/4.0d0) +(RC**2)*asin(L/(2.0d0*RC)) )
        voltotal = 2.0d0*voltotal

       ! voltotal = pi*(RC**2+(RC-RL)**2)*L -2.0d0*(pi/3.0d0)*(L/2.0d0)**3
       ! voltotal = voltotal - pi*(RC-RL)*( (L)*sqrt(RC**2-L**2/4.0d0) +2.0d0*(RC**2)*asin(L/(2.0d0*RC)) )
       ! voltotal = voltotal

    end function total_volume_curv

    subroutine unit_test_area_channel(info)

        use system, only : delta 
        use const, only : pi, stdout
        use MPI, only : rank, ierr

        integer, intent(out) :: info

        ! local argument
        
        real*8  :: areatotalL, areatotalS, areachannel, sumarea, volchan, volchan_math
        real*8 :: a, b, x ,RC, RS, RL, L
        integer :: nsteps, i


        info = 0

        ! test 1 
        ! This tests addition and evaluation of surface_area_channel

        RL = 10.0d0        ! largest radius 
        RS =  5.0d0        ! smallest radius 
        L =  20.0d0         !  
        RC = radiusCurvature(RL,RS,L) 
        areachannel = total_surface_area_curv(RL,RC,L)

        a=-L/2.0d0
        b=-a
        nsteps=(b-a)/delta
        x=a
        sumarea=0.0d0

        do i=1,nsteps
            sumarea = sumarea + surface_area_channel(x,(x+delta),RL,RC,L)
            x = x + delta
        enddo 

        if(dabs(areachannel-sumarea)>0.000000001d0) info = 1
        if(info==1) then
            if(rank==0) then 
                write(stdout,*)"unit_test_area_channel: curved channel area=",areachannel,"sum area=",sumarea
            endif
        endif      

        ! test  2 cylinder 

        info = 0 ! reset 
        RC = 10.0d0
        RL = RC
        RS = RL
        L = 20.0d0

        areatotalL= 2.0d0*pi*RL*L
        areachannel= total_surface_area_curv(RL,RC,L)
        
        if(dabs(areachannel-areatotalL)>0.000000001d0) info = 2
        
        if(info==2) then 
            if(rank==0) write(stdout,*)"unit_test_area_channel: cylinder area=", &
                areachannel,"expected =",areatotalL
        endif

         ! test 3 : volumer of the channel

        info = 0
        L  =   3.0d0     
        RL =   1.25d0     
        rS =  0.5d0     
        rC =   1.8750d0   
        volchan= total_volume_curv(RL,RC,L)
        volchan_math= 10.31808099299808d0       ! from Mathematica 

        if(dabs(volchan-volchan_math)> 0.0000000001d0) info = 3
        if(info==3) then 
            if(rank==0) write(stdout,*)"unit_test_area_channel: volume channel=",volchan,&
                " expected =", volchan_math
        endif
        

        ! test 4 : for the channel

        info = 0

        call set_radii_curved() 

        RL = radiusL
        RS = radiusS
        RC = radiusC
        L = lengthchannel

        areachannel = total_surface_area_curv(RL,RC,L)
                                    
        areatotalL = 2.0d0*pi*RL*L 
        areatotalS = 2.0d0*pi*RS*L 

        if(.not.((areatotalL>areachannel).and.(areachannel>areatotalS))) info = 4


        ! area channel should be between area cylinder with radiusS and radiusL 
        
        if(info==4)  then 
            if(rank==0) write(stdout,*)"unit_test_area_channel: area=",areachannel,&
                " areaL =",areatotalL," areaS =",areatotalS
        endif


        if(info/=0) then  
            if(rank==0)write(stdout,*)"unit_test_area_channel: failure"
            call MPI_FINALIZE(ierr) ! ++end MPI
            stop
        endif

    end subroutine unit_test_area_channel

    ! == calculates volprot for systemtype =4 and curvedflag==1

    subroutine update_matrix_channel_4_curved(flag)

        use system
        use channel
        use ematrix
        use MPI
        use const
        use chainsdat
        use molecules
        use channel
        use transform, only : MAT, IMAT
        use rotchain

        implicit none

        logical, intent(inout) :: flag

        ! == local variables
        ! real*8 :: rchannel2, rchannelL2, rchannelS2
        real*8 :: radiusSeps, radiusLeps
        real*8 :: radiusSq, radiusLq
        real*8 :: lenC
        real*8, external :: rands
        integer :: npoints ! points per cell for numerical integration 
        integer :: counter
        character*5 :: title
        integer :: j,ix,iy,iz
        real :: pnumber
        real*8 :: areachannel, volumechannel
        real*8 :: sumpolseg 
        real*8 :: sstemp,vvtemp, maxss
        real*8 :: cutarea
        real*8 :: temp
       ! real*8 :: temp2
        real*8 :: sumvoleps1, sumvolprot1, sumvolq1, sumvolx1
        integer :: ncha1
        real*8 :: volx1(maxvolx)
        real*8 :: com1(maxvolx,3)
        integer :: p1(maxvolx,3)
        integer :: i
        real*8 :: volxx1(dimx,dimy,dimz)
        real*8 :: volxx(dimx,dimy,dimz)
        real*8 :: x(3), v(3), hcyl
        integer :: nbands

        real*8 :: originc_curv(3)


        cutarea = 0.0 ! throw away cells that have less area than cutarea x area of the cell with largest area  
        sumpolseg = 0.0

        !rchannel2 = rchannel**2
        !rchannelL2 = (rchannel - 3*delta)**2
        !rchannelS2 = (rchannel + delta)**2

        lenC = LengthChannelCurvature()

        ! radii for volprot given by radiusS & radiusL
    
        ! radii for voleps

        radiusSeps = radiusS-3*delta
        radiusLeps = radiusL-3*delta
         
        ! radii for volq

        radiusSq = radiusS + delta
        radiusLq = radiusL + delta


        ! clear all
        voleps = 0.0
        volprot = 0.0
        volq = 0.0
        volx = 0.0
        volxx = 0.0
        com = 0.0
        ncha = 0

        ! channel center in x, y plane

        originc(1) = float(dimx)*delta/2.0 
        originc(2) = float(dimy)*delta/2.0 
        
        ! channel center  in x, y AND z 
        originc_curv(1) = originc(1) 
        originc_curv(2) = originc(2) 
        originc_curv(3) = float(dimz)*delta/2.0  

        npoints = 50

        flag = .false.

        ! call integrate_channel_curved(rchannelL2, RdimZ,originc ,npoints, voleps1 , sumvoleps1, flag)
        
        call integrate_channel_curved(radiusSeps,radiusLeps,RdimZ,originc_curv ,npoints, voleps1 , sumvoleps1, flag)

        flag = .false. ! not a problem if eps lays outside boundaries

        call integrate_channel_curved(radiusS,radiusL,RdimZ, originc_curv,npoints, volprot1, sumvolprot1, flag)

        call integrate_channel_curved(radiusSq,radiusLq,RdimZ,originc_curv ,npoints, voleps1 , sumvoleps1, flag)
        

        ! call newintegrateg_c_4(radiusS,radiusL,RdimZ,originc_curv,npoints,volx1,sumvolx1, com1, p1, ncha1, volxx1, NBRUSH)
        call newintegrate_channel_curved(radiusS,radiusL,RdimZ,originc_curv, npoints,volx1,sumvolx1,com1,p1,&
            ncha1,volxx1, NBRUSH)
     
    
        !! eps
        voleps1 = voleps1-volprot1
        voleps1 = voleps1*eepsc

        ! epstype

        select case (epstype)
        case (1)
            nbands = dimz/8
        do iz = 1, dimz
            if (mod(int((iz-1)/nbands),2).eq.1) then
                voleps1(:,:,iz) = 0.0
            endif
        enddo
        end select


        !! charge
        volq1 = volprot1-volq1
        temp = sum(volq1)
        volq1 = volq1/temp*echargec/(delta**3) ! sum(volq) is echarge

        !! grafting

        v(1) = 0.0
        v(2) = 0.0
        v(3) = float(dimz-2*RdimZ)*delta

        ! v in transformed space, x in real space
        ! only work for gam = 90, cdiva any value

        x = MATMUL(IMAT,v)

        hcyl = x(3) ! height of the cylinder

        !area = 2.0*pi*rchannel*hcyl

        !! volume  
        volprot1 = volprot1 * 0.9999
        volprot = volprot+volprot1

        ! CHECK COLLISION HERE...
        if(maxval(volprot).gt.1.0) then ! collision
            print*,"collision"
            flag=.true. 
        endif
        
        voleps = voleps + voleps1
        volq = volq + volq1 

        ! add com1 and volx to list

        volxx = volxx1

        ncha = ncha1
        do i = 1, ncha
            volx(i)=volx1(i)
            com(i,:)=com1(i,:)
            p0(i,:)=p1(i,:)
            rotangle(i) = atan2(com1(i,1)-originc(1), com1(i,2)-originc(2))
        enddo

        title = 'avpro'
        counter = 1
        call savetodisk(volprot, title, counter)

        sumpolseg = ncha

        if (verbose.ge.2) then
           
            areachannel = total_surface_area_curv(radiusL,radiusC,LenC)
            volumechannel = total_volume_curv(radiusL,radiusC,LenC)
           
            if (rank.eq.0) then
                write(stdout,*) 'channel-curved: L  =',lenC
                write(stdout,*) 'channel-curved: rL =',radiusL
                write(stdout,*) 'channel-curved: rS =',radiusS
                write(stdout,*) 'channel-curved: rC =',radiusC
                write(stdout,*) 'channel-curved: size reservoir =',Rdimz*delta  
                write(stdout,*) 'channel-curved: update_matrix: Total volume =',dimx*dimy*dimz*delta**3
                write(stdout,*) 'channel-curved: update_matrix: Total free volume =',&
                    (dimx*dimy*dimz-sum(volprot))*delta**3
                write(stdout,*) 'channel-curved: update_matrix: volume channel   = ',volumechannel
                write(stdout,*) 'channel-curved: update_matrix: volume vprot     = ',sum(volprot)*delta**3
                write(stdout,*) 'channel-curved: update_matrix: volume sum vprot = ',sumvolprot1*delta**3
                write(stdout,*) 'channel-curved: update_matrix: volume channel via vprot = ',&
                    dimx*dimy*lenC*delta**2-sumvolprot1*delta**3
                write(stdout,*) 'channel-curved: number of polymers in system =', sumpolseg 
                write(stdout,*) 'channel-curved: surface area =', areachannel
                write(stdout,*) 'channel-curved: surface density =', sumpolseg/areachannel

       !         write(stdout,*) 'channel-curved:', 'surface density expected from input =', &
       !             float(NBRUSH)/(2.0*pi*rchannel)/cos(30.0/180.0*pi)/(2.0*pi*rchannel/float(NBRUSH))
            endif
        endif

        title = 'aveps'
        counter = 1
        call savetodisk(voleps, title, counter)

        title = 'avcha'
        counter = 1
        call savetodisk(volq, title, counter)

        title = 'avgrf'
        counter = 1
        call savetodisk(volxx, title, counter)

    end subroutine update_matrix_channel_4_curved


    subroutine integrate_channel_curved(radiusS,radiusL,RdimZ, origincurv, npoints,volprot,sumvolprot, flag)
    
        use system
        use transform

        implicit none

        real*8, intent(in)  :: radiusS,radiusL, origincurv(3)  
        integer, intent(in) :: RdimZ
        real*8, intent(inout) :: volprot(dimx,dimy,dimz)
        real*8, intent(inout) :: sumvolprot
        integer, intent(in) :: npoints
        logical, intent(inout) :: flag ! == in integrate_c flag has nor return !!1
        
        ! local variables
        integer ::  ix,iy,iz,ax,ay,az
        real*8 :: vect, Rz, Rz2, RC,RL,RS, LenC
        logical :: flagin, flagout
        real*8 :: x(3), v(3)
         ! logical :: flagsym
        real*8 :: voltemp

        LenC = LengthChannelCurvature()
        RL = radiusL
        RC = radiusCurvature(radiusL,radiusS,LenC)
        volprot = 0.0d0
        sumvolprot = 0.0d0 ! total volume, including that outside the system
      
        ! scan over all cells
            
        do ix = 1, dimx
            do iy = 1, dimy
                do iz = RdimZ+1, dimz-RdimZ                 ! == for iz< RdimZ+1 and iz> dimz-RdimZ outside : reservoir   

                flagin = .false.
                flagout = .false.

                do ax = 0,1
                    do ay = 0,1
                        do az = 0,1
                                                            ! == v in transformed space
                            v(1) = float(ax+ix-1)*delta     ! == scan over corners in cartesian coordiante aka transfomed coorindates
                            v(2) = float(ay+iy-1)*delta
                            v(3) = float(az+iz-1)*delta 

                                                            ! == x in real space, v in transformed space
                            x = MATMUL(IMAT,v)              ! not striclty  neccessary 

                            x(1) = x(1) - origincurv(1)
                            x(2) = x(2) - origincurv(2)
                            x(3) = x(3) - origincurv(3)

                            Rz = radiusfz(x(3),RL,RC)
                            Rz2 = Rz**2
                           ! print*,"iz=",iz,"z=",x(3),"Rz=",Rz

                            if((x(1)**2+x(2)**2).lt.Rz2) flagin=.true.  ! == inside the channel    .lt. == < 
                            if((x(1)**2+x(2)**2).gt.Rz2) flagout=.true. ! == outside the channel   .gt. == > 


                        enddo
                    enddo
                enddo

                if((flagin.eqv..true.).and.(flagout.eqv..false.)) then ! cell all inside channel
                    voltemp = 0.0
                endif
                if((flagin.eqv..false.).and.(flagout.eqv..true.)) then ! cell all outside channel
                    voltemp = 1.0
                endif
                if((flagin.eqv..true.).and.(flagout.eqv..true.)) then ! cell part inside and outside channel
                    voltemp = integration_cell(origincurv,ix,iy,iz,npoints,RC,RL)
                endif

                sumvolprot = sumvolprot + voltemp
                volprot(ix,iy,iz) = voltemp

                !write(567,*)ix,iy,iz,voltemp,flagin,flagout
                

            enddo ! iz
         !   write(567,*)""
        enddo ! iy
        !write(567,*)""
    enddo ! ix

end subroutine  integrate_channel_curved

! integration over volume cell ix,iy,iz inside curved channel; i.e not accesible by solvent polymer
! channel.f90 this function in called intcell_c

function integration_cell(origincurv,ix,iy,iz,n,RC,RL)result(intcell_c)

    use system, only : delta
    use transform, only : IMAT

    implicit none

    real*8, intent(in)  :: origincurv(3)  
    integer, intent(in) :: ix,iy,iz
    integer, intent(in) :: n
    real*8, intent(in)  :: RC, RL
   
    ! retrun argument 
    real*8 :: intcell_c

    ! local argument
    integer :: ax,ay,az
    integer :: cc
    real*8 :: vect
    real*8 :: dr(3), dxr(3), Rz, Rz2

    cc = 0

    do ax = 1, n
        do ay = 1, n
            do az = 1, n

                ! = points uniform and symetrically distrubuted over cell volume 

                dr(1) = ix*delta-(ax-0.5)*delta/float(n) 
                dr(2) = iy*delta-(ay-0.5)*delta/float(n) 
                dr(3) = iz*delta-(az-0.5)*delta/float(n) 

                ! dr in transformed space
                dxr = MATMUL(IMAT, dr)

                dxr(1) = dxr(1) - origincurv(1)
                dxr(2) = dxr(2) - origincurv(2)
                dxr(3) = dxr(3) - origincurv(3)

                Rz = radiusfz(dxr(3),RL,RC)
                Rz2 = Rz**2
                
                vect = dxr(1)**2+dxr(2)**2
              
                if(vect.gt.Rz2) cc = cc+1   ! outside channel, integrate

            enddo
        enddo
    enddo

    intcell_c = float(cc)/(float(n)**3)
    
end function integration_cell

! == In channel.f90 subroutine called newintegrateg_c_4(
! == calculates volprot for systemtype =4 and curvedflag==1
! == This routine determines the surface coverage and grafting positions only for cylinder

subroutine newintegrate_channel_curved(radiusS,radiusL,RdimZ,origincurv, npoints,volx1,sumvolx1,com1,p1,&
    ncha1,volxx1, NBRUSH)
    
    use system
    use transform
    use chainsdat
    use ematrix
    use const
    use channel, only : sigmar, Nrings, ringpos

    implicit none

    real*8, intent(inout) :: radiusS        ! == smallest radius channel  
    real*8, intent(inout) :: radiusL        ! == largest  radius channel   
    integer, intent(in) :: RdimZ            ! == size reservoir in delta  
    real*8, intent(in) :: origincurv(3)      ! == orgincurv center channel  
    integer, intent(in) :: npoints          ! == number of points
    real*8, intent(inout) :: volx1(maxvolx) 
    real*8, intent(inout) :: sumvolx1
    real*8, intent(inout) :: com1(maxvolx,3)
    integer, intent(inout)::  p1(maxvolx,3)
    real*8, intent(inout) :: volxx1(dimx,dimy,dimz)
    integer, intent(inout) :: ncha1            ! count for current sphere ?? == number of grafts ??
    integer, intent(in) :: NBRUSH           ! == number of graft in theta direction


    ! real*8 rchannel, rchannel2, originc(2),RdimZ,originc, npoints,volx1,sumvolx1,com1,p1,ncha1,volxx1, NBRUSH

    ! == local variables 

    real*8 :: rtetha, rz
    integer :: indexvolx(dimx,dimy,dimz)
    integer :: listvolx(ncha,3)
    real*8 :: phi, dphi, tetha,dtetha, as, ds
    integer :: mphi, mtetha
    integer :: ix,iy,iz,jx,jy,jz
    real*8 :: x(3), v(3)
    integer :: i,j
    integer :: ncount
    real*8 :: comshift ! how far from the surface of the sphere the grafting point i
    integer :: dims(3), is(3), js(3)
    integer :: jjjz, jjjt, npointz, npointt
    real*8 :: hcyl         ! == height cylinder/channel 
    real*8 :: hcyl0        ! == location base of cylinder/channel 
    real*8, external :: rands
    real*8 :: tethaadd     !  , disp

    real*8 :: zcoor       !== z coordinate relative to z origin
    real*8 :: Radiusz     !== z dependent radius of channel  
    real*8 :: RC,RL,RS,LenC,rchannelz

!    disp = delta

    dims(1) = dimx
    dims(2) = dimy
    dims(3) = dimz

   ! rchannel = sqrt(rchannel2)

    indexvolx = 0
    ncha1 = 0
    volx1 = 0.0
    sumvolx1 = 0.0 ! == total volume, including that outside system
    com1 = 0.0
    p1 = 0
    volxx1 = 0.0

    v(1) = 0.0
    v(2) = 0.0
    v(3) = float(dimz-2*RdimZ)*delta

    ! v in transformed space, x in real space
    ! only work for gam = 90, cdiva any value

    x = MATMUL(IMAT,v)

    hcyl = x(3)         ! == height of the cylinder

    npointt = NBRUSH    ! == number of sites along the tetha coordinate

    v(1) = 0.0
    v(2) = 0.0
    v(3) = float(RdimZ)*delta
    x = MATMUL(IMAT,v)
    hcyl0 = x(3)        ! == position of the base of the cylinder

    select case (systemtype)
    case (42)
        npointz = Nrings   
    case default
        write(stdout,*)"newintegrate_channel_curved: wrong systemtype"
        stop
    end select

    RS = radiusS
    RL = radiusL
    LenC = LengthChannelCurvature()
    RC = radiusCurvature(RL,RS,LenC)

    do jjjz = 1, npointz      ! == loop over number of graft point in z and  the
        do jjjt = 1, npointt

            select case (randominput) ! == random displacement of rz and rtheta see channel.f90 for complete list here no discplacement!!!
            case(0)
                rtetha = 0.0 
                rz = 0.0
            case default
                write(stdout,*)" newintegrate_channel_curved: wrong randominput value"
                stop
            end select


            !if ((systemtype.eq.42).or.(systemtype.eq.52).or.(systemtype.eq.60)) then
            
            tethaadd = 0.0
            
            ! == x and y coordiante for cylinder with no curvature update formula 

            ! x(1) = cos(float(jjjt-1)/float(npointt)*2.0*pi+rtetha+tethaadd)*rchannel + originc(1)
            ! x(2) = sin(float(jjjt-1)/float(npointt)*2.0*pi+rtetha+tethaadd)*rchannel + originc(2)
            ! x(3) = ringpos(jjjz)*hcyl+hcyl0 

            !zcoor =(ringpos(jjjz)*hcyl+hcyl0 -origincurv(3)) 
           
            !== zcoor for radiusfz needs to be centered around orgincurv channel 
            !== ringpos relative location graft point range [-0.5:0.5] 
            !== ringpos=-0.5 == base zrel=-hcyl/2: ringbase=0.5 to zrel=+hcyl/2 
           
            zcoor =((ringpos(jjjz)+0.5d0)*hcyl+hcyl0-origincurv(3)) 
            Radiusz = radiusfz(zcoor,RL,RC)

            x(1) = cos(float(jjjt-1)/float(npointt)*2.0*pi+rtetha+tethaadd)*Radiusz + origincurv(1)
            x(2) = sin(float(jjjt-1)/float(npointt)*2.0*pi+rtetha+tethaadd)*Radiusz + origincurv(2)

           ! x(3) = ringpos(jjjz)*hcyl+hcyl0  ! Coordinate system has unusally origin of middle of cylinder == hcyl0 == base cylinder
 
            x(3) = (ringpos(jjjz)+0.5d0)*hcyl+hcyl0 !  Coordinate system has orgin of lattice 


            ! == x in  real space

            v = MATMUL(MAT,x)


            ! == v(3) = v(3) + float((dimz-RdimZ*2))/2.0*delta  ! == RJN different for curvature
            ! == centers the first row of polymers at the middle of the layer, useful to avoid numerical rounding errors.
            ! == this translate v(3) by half height=(dimz-2Rdimz)delta/2 of cylinder !!!! 
            ! == if((systemtype.eq.42).or.(systemtype.eq.52).or.(systemtype.eq.60)) then
            ! ==   v(3) = v(3) + float((dimz-RdimZ*2))/2.0*delta ! centers the first row of polymers at the middle of the layer, useful to avoid numerical rounding errors.
            ! == else
            ! ==    v(3) = v(3) + float((dimz-RdimZ*2)/npointz)/2.0*delta ! centers the first row of polymers at the middle of the layer, useful to avoid numerical rounding errors.
            ! == endif
            ! == x = MATMUL(IMAT,v) ! and recalculates x due to the change in v ! == unneccsary

            do j = 1,3
               ! js(j) = floor(v(j)/delta)+1
                js(j) = int(v(j)/delta)+1
            enddo

            ! js(3)=mod(js(3)+dimz-1,dimz)+1    
            !== this is translate js integer (v)  by (dimz-1) and retrun remainder  of the division by dimz ??


            jx = js(1)
            jy = js(2)
            jz = js(3)

            do i = 1, 3
                if((js(i).le.0).or.(js(i).gt.dims(i))) then
                    write(stdout,*) 'neewintegrate_channel_curved: error in channel-curved', i, js(i), dims(i)
                    write(stdout,*) v(1), v(2), v(3)
                    stop
                endif
            enddo

             write(125,*) jjjz, jjjt,  jx,jy,jz, x
            ! increase counter

            if(ncha1.eq.maxvolx) then
                write(stdout,*) 'newintegrate_channel_curved:: increase maxvolx'
                stop
            endif

            ncha1 = ncha1 + 1

            indexvolx(jx,jy,jz) = ncha1
            p1(ncha1,1)=jx
            p1(ncha1,2)=jy
            p1(ncha1,3)=jz

            volxx1(jx,jy,jz) =  1.0
            volx1(indexvolx(jx,jy,jz)) = 1.0
            com1(indexvolx(jx,jy,jz),:) = x(:) ! == caredull use coordiante systmam assoicaite with vector x 
            sumvolx1 = sumvolx1 + 1.0

        enddo ! jjjt
    enddo ! jjjz

    do i = 1, ncha1
        ! Moves the position of the first segment lseg/2 away from the surface to prevent collision due to round errors.
        ! == lseg 
        
        rchannelz = dsqrt((com1(i,1)-origincurv(1))**2+(com1(i,2)-origincurv(2))**2)
        
       ! print*,"rchannelz=",rchannelz
        
        com1(i,1) = com1(i,1) - lseg*((com1(i,1)-origincurv(1)))/rchannelz ! replace rchannel with rchannelz : radius of channel is z dependent 
        com1(i,2) = com1(i,2) - lseg*((com1(i,2)-origincurv(2)))/rchannelz 
      
    enddo

end subroutine 


end module channelcurved