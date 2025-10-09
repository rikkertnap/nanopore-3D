! module computate divergence of diffusive (flux div_flux)
! and related quantities

module flux

    implicit none

    real*8, allocatable, dimension(:,:,:,:) :: divJ    ! divergence of flux range divJ(nsize,niontypes)
    real*8, allocatable, dimension(:,:,:)   :: mu      ! chemical potential range mu(nsize)
    real*8, allocatable, dimension(:,:,:,:,:) :: Jvec    ! vector flux in  range J(nsize,3,niontypes)
    real*8, allocatable, dimension(:,:,:,:) :: mu_ion  ! chemical potential range mu_ion(nsize,niontypes) 
    

    logical, parameter  :: DEBUG_ST=.true.
   
    character(len=5), parameter :: iontype(4)=(/"pos  ","neg  ","Hplus","OHmin"/) 
    integer, parameter :: niontypes=4 

    public :: divJ, mu, Jvec, mu_ion
    public :: div_flux, allocate_divJ, allocate_mu

contains

    subroutine allocate_divJ()

        use system, only : dimx, dimy, dimz
        
        allocate(divJ(dimx, dimy, dimz, niontypes))

    end subroutine allocate_divJ

    subroutine allocate_mu()

        use system, only : dimx, dimy, dimz

        allocate(mu(0:dimx+1, 0:dimy+1, 0:dimz+1))

    end subroutine allocate_mu

    subroutine allocate_Jvec()

        use system, only : dimx, dimy, dimz

        allocate(Jvec(dimx, dimy, dimz, 3, niontypes))

    end subroutine allocate_Jvec

    subroutine allocate_mu_ion()

        use system, only : dimx, dimy, dimz

        allocate(mu_ion(0:dimx+1, 0:dimy+1, 0:dimz+1, niontypes))

    end subroutine allocate_mu_ion


   ! position dependent chem pot  
   ! \beta \mu_i(r) = \ln(\rho_i(r) v_w) + \beta \pi(r) * v_i + \beta * q_i \psi(r)
   ! \beta \mu_i(r) = \ln(x_i(r) v_w/v_i) - \ln(x_w(r)) * v_i/v_w +  z_i *  e  * \beta \psi(r)
   ! input : real*8 xsol,xvol,psi,iontype)    

    subroutine chem_potential(mu,xsol,xvol,psi,iontype)

        use system, only : dimx, dimy, dimz
        use moleculeslist, only :  vol, zval, get_value_moleclist

        ! input arguments 
        real*8, intent(inout) :: mu(:,:,:)
        real*8, intent(in) :: xvol(:,:,:), xsol(:,:,:)
        real*8, intent(in) :: psi(:,:,:)
        character(len=*), intent(in)  :: iontype

       
        ! local variables

        integer :: ix, iy,iz
        character(len=5):: key
        real*8 :: volum, valence
        
        key = trim(iontype)
        volum   = get_value_moleclist(vol,key)
        valence = get_value_moleclist(zval,key)
        
        if(.true.) then
           print*,"key=",key  
           print*,"volum=",volum
           print*,"valence=",valence
        endif
        
        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx
                    mu(ix,iy,iz) = dlog(xvol(ix,iy,iz)/volum) -dlog(xsol(ix,iy,iz)) * volum +&
                         valence * psi(ix,iy,iz)
                enddo 
            enddo
        enddo   
 
    end subroutine chem_potential


    subroutine volumefraction(xvol,xsol,mu,psi,volum,valence )

        use system, only : dimx, dimy, dimz

        ! input arguments 
        real*8, intent(inout) :: xvol(:,:,:)
        real*8, intent(in) :: mu(:,:,:),xsol(:,:,:)
        real*8, intent(in) :: psi(:,:,:)
        real*8, intent(in) :: volum
        integer, intent(in) :: valence
    
        ! local variables
        integer ::  ix, iy, iz
        
        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx
                    xvol(ix,iy,iz)= volum * exp(  mu(ix,iy,iz) &
                        - valence * psi(ix,iy,iz)) *  ( xsol(ix,iy,iz) ** volum) 
                enddo 
            enddo
        enddo               
 
         
    end subroutine volumefraction

    ! boundary conditions chemical potential and xvol 

    subroutine bc_flux(mu,xvol,iontype)
        
        use system, only : dimx, dimy, dimz  
        use moleculeslist, only : mumin, mumax, xvolmin, xvolmax
        use moleculeslist, only : get_value_moleclist

        ! input arguments 
        real*8, intent(inout) :: xvol(0:dimx+1,0:dimy+1,0:dimz+1)
        real*8, intent(inout) ::   mu(0:dimx+1,0:dimy+1,0:dimz+1)
        character(len=*), intent(in) :: iontype
          
        ! local variables

        character(len=5):: key
        real*8 :: mu_zmin, mu_zmax, xvol_zmin, xvol_zmax

        key = trim(iontype)
        mu_zmin = get_value_moleclist(mumin,key)
        mu_zmax = get_value_moleclist(mumax,key)
        xvol_zmin = get_value_moleclist(xvolmin,key)
        xvol_zmax = get_value_moleclist(xvolmax,key)

        ! boundary conditions chemical potential and xvol 

        ! x = 0
        mu(0,:,:) = mu(1,:,:) 
        xvol(0,:,:) = xvol(1,:,:) 

        ! x = dimx
        mu(dimx+1,:,:) = mu(dimx,:,:)  
        xvol(dimx+1,:,:) = xvol(dimx,:,:)  

        ! y = 0
        mu(:,0,:) = mu(:,1,:) 
        xvol(:,0,:) = xvol(:,1,:) 

        ! y = dimy
        mu(:,dimy+1,:) = mu(:,dimy,:) 
        xvol(:,dimy+1,:) = xvol(:,dimy,:) 

        ! z = 0
        mu(:,:,0) = mu_zmin
        xvol(:,:,0) = xvol_zmin

        ! z = dimz
        mu(:,:,dimz+1) = mu_zmax
        xvol(:,:,dimz+1) = xvol_zmax

        ! internal membrane still to do ...
        print*,"Warning bc_flux subroutine : internal membrane symmetries not applied yet"

    end subroutine bc_flux



    subroutine div_flux(divJ,xsol,xion,psi,iontype)

        ! input arguments 
        real*8, intent(inout) :: divJ(:,:,:)
        real*8, intent(in) :: xion(:,:,:), xsol(:,:,:)
        real*8, intent(in) :: psi(:,:,:)
        character(len=*), intent(in)  :: iontype

        call div_flux_cubic(divJ,xsol,xion,psi,iontype)
                
    end subroutine div_flux
        
    ! Computes div.J with numerical scheme using Gauss's theorem 

    subroutine div_flux_cubic(divJ,xsol,xion,psi,iontype)

        use system, only : dimx, dimy, dimz
        use moleculeslist, only : vol, zval, mumin, mumax, xvolmin, xvolmax
        use moleculeslist, only : get_value_moleclist

        ! input arguments 

        real*8, intent(inout) :: divJ(:,:,:)
        real*8, intent(in) :: xion(:,:,:), xsol(:,:,:)
        real*8, intent(in) :: psi(:,:,:)
        character(len=*), intent(in) :: iontype
        
        ! local variables

        integer :: ix, iy, iz
        integer ::  id, idxpls, idxmin, idypls, idymin, idzpls, idzmin
        character(len=5):: key
        real*8 :: volum, valence, coeff_scaled
        real*8 :: mu_zmin, mu_zpls, xvol_zmin, xvol_zpls
        real*8 :: mu(0:dimx+1, 0:dimy+1, 0:dimz+1)
        real*8 :: xvol(0:dimx+1, 0:dimy+1, 0:dimz+1)
        real*8 :: Jdotxpls, Jdotxmin, Jdotypls, Jdotymin, Jdotzpls, Jdotzmin, divJtmp 

        key = trim(iontype)

        volum   = get_value_moleclist(vol,key)
        valence = get_value_moleclist(zval,key)
        mu_zmin = get_value_moleclist(mumin,key)
        mu_zpls = get_value_moleclist(mumax,key)
        xvol_zmin = get_value_moleclist(xvolmin,key)
        xvol_zpls = get_value_moleclist(xvolmax,key)
        
        !coeff_scaled = 1.0_dp/(volum*delta*2.0_dp)
        coeff_scaled =  1.0d0

        !  volum from conversion of volumefraction to density  
        !  in cylinder coordinates : 2 delta^2 from division by 2 \pi delta^2 {\bar r}_i 
        !  in cubic coordiantes :  2 delta  :  from interpolation of density

        if(DEBUG_ST) then
           print*,"key=",key  
           print*,"volum=",volum
           print*,"valence=",valence
           print*,"mu_zmin=",mu_zmin
           print*,"mu_zpls=",mu_zpls
           print*,"xvol_zpls=",xvol_zpls
           print*,"xvol_zmin=",xvol_zmin
           print*,"size(xvol)=",size(xvol)
        endif

        
        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx
                    xvol(ix,iy,iz) = xion(ix,iy,iz) ! xvol and xion not same range !!!
                    mu(ix,iy,iz) = log(xvol(ix,iy,iz)/volum) -log(xsol(ix,iy,iz))*volum + valence * psi(ix,iy,iz)
                enddo
            enddo
        enddo    
      
        ! apply boundary conditions to mu and xvol

        call bc_flux(mu,xvol,iontype)
        
    
        divJ=0.0d0
        if(DEBUG_ST) divJ= 123435600.0000d0 ! used to detect unassinged values of divJ 
        
        ! inside 
        
        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx
                    
                    Jdotxpls = (xvol(ix+1,iy,iz) + xvol(ix  ,iy,iz))*(mu(ix+1,iy,iz) - mu(ix,iy,iz)  ) 
                    Jdotxmin = (xvol(ix  ,iy,iz) + xvol(ix-1,iy,iz))*(mu(ix ,iy,iz)  - mu(ix-1,iy,iz))
                    
                    Jdotypls = (xvol(ix,iy+1,iz) + xvol(ix,iy  ,iz))*(mu(ix,iy+1,iz) - mu(ix,iy,iz)  ) 
                    Jdotymin = (xvol(ix,iy  ,iz) + xvol(ix,iy-1,iz))*(mu(ix,iy  ,iz) - mu(ix,iy-1,iz))
                   
                    Jdotzpls = (xvol(ix,iy,iz+1) + xvol(ix,iy,iz  ))*(mu(ix,iy,iz+1) - mu(ix,iy,iz)  ) 
                    Jdotzmin = (xvol(ix,iy,iz  ) + xvol(ix,iy,iz-1))*(mu(ix,iy,iz  ) - mu(ix,iy,iz-1))
               
                    divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
                    divJ(ix,iy,iz) = - coeff_scaled * divJtmp

                enddo
            enddo
        enddo    

    
        if(DEBUG_ST) then
            do iz=1,dimz
                do iy=1,dimy 
                    do ix=1,dimx    
                        if(divJ(ix,iy,iz)== 123435600.0000d0) print*,"divJ unassiged in ix=",ix,"iy=",iy,'iy=',iy
                    enddo
                enddo
            enddo             
        endif    

    end subroutine div_flux_cubic


    ! Computes flux J= - D rho * grad mu 

    subroutine fluxJ(Jvec,xsol,xvol,psi,iontype)

        use system, only : dimx, dimy, dimz
        use molecules, only : vsol
        use moleculeslist, only : vol, zval, mumin, mumax, xvolmin, xvolmax,  Diffcoeff
        use moleculeslist, only : get_value_moleclist

        ! input arguments 

        real*8, intent(inout) :: Jvec(:,:,:,:)
        real*8, intent(inout) :: xvol(:,:,:)
        real*8, intent(in) :: psi(:,:,:), xsol(:,:,:)
        character(len=*), intent(in)  :: iontype

        ! local arguments

        real*8 :: grad_mu(dimx,dimy,dimz,3)
        real*8 ::  volum, Diffconst, Jvec0
        integer :: i, ix, iy, iz, idx 
        character(len=5):: key
 
        call grad_chem_pot(grad_mu,xsol,xvol,psi,iontype)

        key = trim(iontype)
        volum = get_value_moleclist(vol,key)
        Diffconst = get_value_moleclist(Diffcoeff,key)

        ! pre factor in flux factor (1.0e-9)^2 arise from vsol in unit of nm and grad_mu in 1/nm
        ! J in 1/( nm^2) s)
        
        Jvec0 = - Diffconst / ( volum * vsol * (1.0d-9)**2) 
        
        do i=1,3 
            do iz=1,dimz
                do iy=1,dimy 
                    do ix=1,dimx 
                        Jvec(ix,iy,iz,i) = Jvec0 * xvol(ix,iy,iz) * grad_mu(ix,iy,iz,i) 
                    enddo
                enddo
            enddo
        enddo            

        if(DEBUG_ST) then
            if(key=="K") then
                do iz=1,dimz 
                    do iy=1,dimy
                        do ix=1,dimx
                            write(200,*)ix,iy,iz,Jvec(ix,iy,iz,1),Jvec(ix,iy,iz,2),Jvec(ix,iy,iz,3)
                            write(300,*)ix,iy,iz,grad_mu(ix,iy,iz,1),grad_mu(ix,iy,iz,2),grad_mu(ix,iy,iz,3)
                        enddo
                    enddo
                enddo            
            endif
        endif        

    end subroutine fluxJ 

    ! Computes  total current through plane z=nz/2 delta 

    function current_I() result(currI)

        use system, only : dimx, dimy, dimz
        use results, only : xpos, xneg, xHplus,xOHmin
        use fields_fkfun, only : xsol=>xh, psi
       
        ! return arguments 

        real*8 :: currI

        ! local arguments

        real*8 :: J(dimx, dimy, dimz,3)
        integer :: t
        real*8 :: current(niontypes)
        
        current = 0.0d0

        do t=1,niontypes
           ! if(isionselfconsistent(t)) then
                if(iontype(t)=="pos") then 
                    call fluxJ(J,xsol,xpos,psi,iontype(t))
                    current(t) = current_I_ion(J,iontype(t))
                endif    
                if(iontype(t)=="neg") then 
                    call fluxJ(J,xsol,xneg,psi,iontype(t))
                    current(t) = current_I_ion(J,iontype(t))
                endif 
                if(iontype(t)=="Hplus") then
                     call fluxJ(J,xsol,xHplus,psi,iontype(t))
                    current(t) = current_I_ion(J,iontype(t))
                endif    
                if(iontype(t)=="OHmin") then 
                    call fluxJ(J,xsol,xOHmin,psi,iontype(t))
                    current(t) = current_I_ion(J,iontype(t))
                endif
                    
            ! endif    
        enddo

        currI = sum(current)

    end function current_I 


    ! Computes  current through plane z=nz/2 delta  for given flux density J 
    ! for given charge velence zval 

    function current_I_ion(J,iontype) result(current)

        use system, only : dimx, dimy, dimz, delta
        ! use results, only : xpos,xneg, xHplus,xOHmin
        ! use fields_fkfun, only : xh,psi
        use moleculeslist, only :  zval, get_value_moleclist
        

        ! input arguments 

        real*8, intent(inout) :: J(:,:,:,:)
        character(len=*), intent(in)  :: iontype

        ! return arguments

        real*8 :: current

        ! local variables

        character(len=5):: key
        real*8 :: valence, sum_curJ
        integer ix, iy, izmidplane

        key = trim(iontype)
        valence = get_value_moleclist(zval,key)
        izmidplane = int(dimz/2.0d0)
        sum_curJ = 0.0d0

        do iy=1,dimy
            do ix=1,dimx
                sum_curJ = sum_curJ + J(ix,iy,izmidplane,3) 
            enddo
        enddo    
        
        current = sum_curJ *( delta **2) * valence
        
    end function current_I_ion

    ! computes grad of chemical potential

    subroutine grad_chem_pot(grad_mu,xsol,xion,psi,iontype)

        use system, only : delta, dimx, dimy, dimz
        use moleculeslist, only : vol, zval, mumin, mumax, xvolmin, xvolmax
        use moleculeslist, only : get_value_moleclist

        ! input arguments 

        real*8, intent(inout) :: grad_mu(:,:,:,:)
        real*8, intent(in) :: xion(:,:,:), xsol(:,:,:)
        real*8, intent(in) :: psi(:,:,:)
        character(len=*), intent(in)  :: iontype

        ! local variables

        integer :: ix, iy, iz, k
        integer ::  id, idxpls, idxmin, idypls, idymin, idzpls, idzmin
        character(len=5):: key
        real*8 :: volum, valence
        real*8 :: mu_zmin, mu_zpls, xvol_zmin, xvol_zpls
        real*8 :: mu(0:dimx+1,0:dimy+1,0:dimz+1)
        real*8 :: xvol(0:dimx+1,0:dimy+1,0:dimz+1)

        key = trim(iontype)

        volum   = get_value_moleclist(vol,key)
        valence = get_value_moleclist(zval,key)
        mu_zmin = get_value_moleclist(mumin,key)
        mu_zpls = get_value_moleclist(mumax,key)
        xvol_zmin = get_value_moleclist(xvolmin,key)
        xvol_zpls = get_value_moleclist(xvolmax,key)
        
        if(DEBUG_ST) grad_mu = 1234567.89d0

        !call chem_potential(mu,xsol,xvol,psi,iontype) 

        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx
                    xvol(ix,iy,iz) = xion(ix,iy,iz) ! xvol and xion not same range !!!
                    mu(ix,iy,iz) = log(xvol(ix,iy,iz)/volum) -log(xsol(ix,iy,iz))*volum + valence * psi(ix,iy,iz)
                enddo
            enddo
        enddo    
    
        ! apply boundary condtions 

        call bc_flux(mu,xvol,iontype)

        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx

                    grad_mu(ix,iy,iz,1)  = mu(ix+1,iy,iz) - mu(ix-1,iy,iz)
                    grad_mu(ix,iy,iz,2)  = mu(ix,iy+1,iz) - mu(ix,iy-1,iz)
                    grad_mu(ix,iy,iz,3)  = mu(ix,iy,iz+1) - mu(ix,iy,iz-1)
                   
                enddo
            enddo
        enddo    

        ! divided by 2 delta 
        
        grad_mu = grad_mu /(2.0d0*delta)
        

    end subroutine  grad_chem_pot

    ! computes flux J

    subroutine calculate_fluxJ() 

        use fields_fkfun, only  : xsol=> xh, psi
        use results, only  : xpos, xneg, xHplus, xOHmin
    
        ! local arguments
        integer :: t

        Jvec = 0.0d0

        do t=1,niontypes
           ! if(isionselfconsistent(t)) then
                if(iontype(t)=="pos")   call fluxJ(Jvec(:,:,:,:,t),xsol,xpos,  psi,iontype(t))
                if(iontype(t)=="neg")   call fluxJ(Jvec(:,:,:,:,t),xsol,xneg,  psi,iontype(t))
                if(iontype(t)=="Hplus") call fluxJ(Jvec(:,:,:,:,t),xsol,xHplus,psi,iontype(t))
                if(iontype(t)=="OHmin") call fluxJ(Jvec(:,:,:,:,t),xsol,xOHmin,psi,iontype(t))
           ! endif 
        enddo

    end subroutine calculate_fluxJ


    subroutine calculate_mu_ion() 

        use fields_fkfun, only  : xsol=>xh, psi
        use results, only  : xpos, xneg, xHplus, xOHmin
    
        ! local arguments
        integer :: t

        mu_ion = 0.0d0

        do t=1,niontypes
           ! if(isionselfconsistent(t)) then
                                
                if(iontype(t)=="pos")   call chem_potential(mu_ion(:,:,:,t),xsol,xpos,  psi,iontype(t))
                if(iontype(t)=="neg")   call chem_potential(mu_ion(:,:,:,t),xsol,xneg,  psi,iontype(t))
                if(iontype(t)=="Hplus") call chem_potential(mu_ion(:,:,:,t),xsol,xHplus,psi,iontype(t))
                if(iontype(t)=="OHmin") call chem_potential(mu_ion(:,:,:,t),xsol,xOHmin,psi,iontype(t))
              
            !endif    
        enddo

    end subroutine calculate_mu_ion

    subroutine unit_test_diff

        use system, only : delta, dimx, dimy, dimz
        use moleculeslist, only : mumin, mumax, xvolmin, xvolmax
        use moleculeslist, only : get_value_moleclist

        integer :: ix,iy,iz

        ! local varailbe not in module just for test
        real*8 :: xion(dimx,dimy,dimz)
        real*8 :: xsol(dimx,dimy,dimz)
        real*8 :: psi(0:dimx+1,0:dimy+1,0:dimz+1)
        real*8 :: divJ(dimx, dimy, dimz)

        real*8 ::  mu_zmin, mu_zpls, xvol_zmin, xvol_zpls
        character(len=5):: key
        real*8 :: slope, intercept, psizmin, psizmax, sumdivJ

        print*,"Warning need to get surface potential imported."
        
        psizmax=0.0d0 
        psizmin=1.0d0

        key = "pos"

        mu_zmin = get_value_moleclist(mumin,key)
        mu_zpls = get_value_moleclist(mumax,key)
        xvol_zmin = get_value_moleclist(xvolmin,key)
        xvol_zpls = get_value_moleclist(xvolmax,key)

        slope = (psizmax-psizmin)/((dimz+1.0d0)*delta)
        intercept = psizmin + slope * delta/2.0d0

        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx
                    xsol(ix,iy,iz) =1.0d0 - 2.0d0 * xvol_zmin 
                    xion(ix,iy,iz) = xvol_zmin
                    psi(ix,iy,iz) =  slope * (iz - 0.5d0) * delta  + intercept
                enddo
            enddo
        enddo        

        call div_flux(divJ,xsol,xion,psi,key)

        sumdivJ = 0.0d0

        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx
                    sumdivJ=sumdivJ+divJ(ix,iy,iz)**2
                enddo
            enddo
        enddo      

        print*,"unit test: sumdivJ=",sumdivJ

    end subroutine unit_test_diff

end module flux
