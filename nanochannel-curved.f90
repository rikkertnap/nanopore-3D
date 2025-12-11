module nanochannel-curved

    implicit none

    real*8 :: radiusL         ! largest radius 
    real*8 :: radiusS         ! smallest radius 
    real*8 :: lenghtchannel   ! total lenght channel 
    real*8 :: radiusC         ! curvature nanochannel 
    real*8 :: nsections       ! number of curved section 
    real*8 :: lenghtsection   ! lenght channel section  lenghtchannel/nsections


contains 

    ! not used init !!
    subroutine init_curved_var

        use system, only : dimz, delta 
        use channel :: RdimZ      ! in units of delta 

        lengthchannel = (dimz-2.0d0*Rdimz)*delta
        lenghtsection = lengthchannel/(1.0d0*nsections)

        !radiusC=(radiusL-radiusS)/2.0d0 +lenghtchannel**2/(8.0d0*(radiusL-radiusS))

        radiusC=(radiusL-radiusS)/2.0d0 +lenghtsection**2/(8.0d0*(radiusL-radiusS))

    end function  

    ! function return radius of channel 
    ! -L/2<= z <= ?/2 with L=lengthsection 

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

    ! area related to a ONE curved subsection of  nanochannel  
    function surface_area_channel(a,b) result(area)
        use const, only : pi
          
        real*8, intent(in) :: a,b ! integraton boundary along z-axis origin in middle
    
        ! return argument 
        real*8 :: area

        area = 2.0d0*pi*radiusC*(b-a)-2.0d0*pi*(radiusC-radiusL)*radiusC*(asin(b/radiusC) -asin(a/RadiusC) )

    end function 

    subroutine total_area_channel(areatotal,info)

        use system, only : delta 
        integer, intent(out) :: info
         real*8, intent(inout) :: areatotal

        ! local argument

        real*8 :: sumarea, a, b, x 
        integer :: nsteps, i

        ! need to change init 
        ! radiusL =10.0d0        ! largest radius 
        ! radiusS = 5.0d0        ! smallest radius 
        ! lenghtchannel = 20.0d0 ! 
        

        radiusC=(radiusL-radiusS)/2.0d0 +lenghtchannel**2/(8.0d0*(radiusL-radiusS))

        info = 1
        areatotal = 2.0d0*pi*radiusC*lenghtchannel-4.0d0*pi*(radiusC-radiusL)*radiusC*asin(lengthchannel/radiusC) 
        areatotal = nsection * areatotal 

        a=-lengthchannel
        b=a
        nsteps=(b-a)/delta
        x=a
        do i=1,nsteps
            sumarea=sumarea+surface_area_channel(x,(x+delta))
            x=x+delta
        enddo 
        if(dabs(areatotal-sumarea))<0.0000000001d0) info =0

        if(info.neq.0)print*,"check-total_area: areattotal+", areatotal,"sumarea=",sumarea

    end subroutine 



end module nanaochannel-curved