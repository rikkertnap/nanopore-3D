module channelcurved

    implicit none

    real*8 :: radiusL         ! largest radius 
    real*8 :: radiusS         ! smallest radius 
    real*8 :: lengthchannel    
    real*8 :: radiusC         ! curvature nanochannel 


contains 

    subroutine init_curved_var

        use system, only : dimz, delta 
        use channel, only  : RdimZ      ! Rdimz is of in units of delta 

        lengthchannel=dimz-2.0d0*Rdimz*delta
        radiusC=(radiusL-radiusS)/2.0d0 +lengthchannel**2/(8.0d0*(radiusL-radiusS))

    end subroutine init_curved_var


    function radiusf(z) result(Rz)

         
        real*8, intent(in) :: z 
        ! return argument 
        real*8 :: Rz
        if(Rz>=z) then 
             Rz = sqrt(radiusC*radiusC-z*z) -(radiusc-radiusl)
        else
           Rz=0.0d0 
        endif

    end function 

    function surface_area_channel(a,b) result(area)
        use const, only : pi
          
        real*8, intent(in) :: a,b ! integraton boundary along z-axis origin in middle
        ! return argument 
        real*8 :: area

        area= 2.0d0*pi*radiusC*(b-a)-2.0d0*pi*(radiusC-radiusL)*radiusC*(asin(b/radiusC) -asin(a/RadiusC) )

    end function 

    function total_surface_area_curv(RL,RC,RS,L) result(areatotal)
        use const, only : pi
        
        real*8, intent(in) :: RL,RC,RS,L
        real*8 :: areatotal

        areatotal= 2.0d0*pi*RC*L-4.0d0*pi*(RC-RL)*RC*asin(L/(2.0d0*RC)) 

    end function
    
    subroutine unit_test_area_channel(info)

        use system, only : delta 
        use const, only : pi, stdout

        integer, intent(out) :: info

        ! local argument
        
        real*8  :: areatotal
        real*8 :: sumarea, a, b, x ,RC,RS,RL,L
        integer :: nsteps, i

        info = 1

        ! need to change init 
        radiusL =10.0d0        ! largest radius 
        radiusS = 5.0d0        ! smallest radius 
        lengthchannel = 20.0d0 ! 
        radiusC=(radiusL-radiusS)/2.0d0 +lengthchannel**2/(8.0d0*(radiusL-radiusS))

        areatotal= total_surface_area_curv(radiusL,radiusC,radiusS,lengthchannel)

        a=-lengthchannel/2.0d0
        b=-a
        nsteps=(b-a)/delta
        x=a
        sumarea=0.0d0
        do i=1,nsteps
            sumarea=sumarea+surface_area_channel(x,(x+delta))
            x=x+delta
        enddo 
        if(dabs(areatotal-sumarea)<0.0000000001d0) info = 0

        if(info.eq.1)print*,"unit_test_area_channel: curved channel area=", sumarea,"expected area=",areatotal

        ! simple test case : cylinder 
        RC=10.0d0
        RL=RC
        RS=RL
        L=25.0d0

        areatotal= 4.0d0*pi*RL*L
        sumarea= total_surface_area_curv(RL,RC,RS,L)
        
        if(dabs(areatotal-sumarea)<0.0000000001d0) info = 0

        if(info.eq.1)write(stdout,*)"unit_test_area_channel: cylinder area=", sumarea,"sumarea=",areatotal
        

    end subroutine unit_test_area_channel



end module channelcurved