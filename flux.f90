! module computate divergence of diffusive (flux div_flux)
! and related quantities

module flux

    implicit none

    real*8, allocatable, dimension(:,:,:,:) :: divJ    ! divergence of flux range divJ(nsize,niontypes)
    real*8, allocatable, dimension(:,:,:)   :: mu      ! chemical potential range mu(nsize+bc)
    real*8, allocatable, dimension(:,:,:)   :: mu_ex   ! excess or non ideal chemical potential part range mu(nsize+bc)
    real*8, allocatable, dimension(:,:,:,:,:) :: Jvec  ! vector flux in range J(nsize,3,niontypes)
    real*8, allocatable, dimension(:,:,:,:) :: mu_ion  ! chemical potential range mu_ion(nsize,niontypes) 
    

    logical, parameter  :: DEBUG_ST=.true.
   
    character(len=5), parameter :: iontype(4)=(/"pos  ","neg  ","Hplus","OHmin"/) 
    integer, parameter :: niontypes=4 
    integer, parameter :: xpls=1, xmin=2, ypls=3, ymin=4, zpls=5, zmin=6 ! number edges of cell

    private 
    public :: divJ, mu, Jvec, mu_ion, iontype, niontypes, bc_psi
    public :: div_flux, allocate_flux_var, init_flux_var ,unit_test_divJ, linear_interpolation
    public :: xpls, xmin, ypls, ymin, zpls, zmin

contains

    subroutine allocate_flux_var()
        
        use system, only :ST_bctype

        call allocate_divJ()
        call allocate_mu()
        call allocate_Jvec()
        call allocate_mu_ion()

        if(ST_bctype.eq.3) call allocate_mu_ex()

    end subroutine allocate_flux_var

    subroutine allocate_divJ()

        use system, only : dimx, dimy, dimz
        
        allocate(divJ(dimx, dimy, dimz, niontypes))

    end subroutine allocate_divJ

    subroutine allocate_mu()

        use system, only : dimx, dimy, dimz

        allocate(mu(0:dimx+1, 0:dimy+1, 0:dimz+1))

    end subroutine allocate_mu

    subroutine allocate_mu_ex()

        use system, only : dimx, dimy, dimz

        allocate(mu_ex(0:dimx+1, 0:dimy+1, 0:dimz+1))

    end subroutine allocate_mu_ex


    subroutine allocate_Jvec()

        use system, only : dimx, dimy, dimz

        allocate(Jvec(dimx, dimy, dimz, 3, niontypes))

    end subroutine allocate_Jvec

    subroutine allocate_mu_ion()

        use system, only : dimx, dimy, dimz

        allocate(mu_ion(0:dimx+1, 0:dimy+1, 0:dimz+1, niontypes))

    end subroutine allocate_mu_ion

    subroutine init_flux_var

        use moleculeslist, only : vol, zval, mumin, mumax, xvolmin, xvolmax, Diffcoeff
        use moleculeslist, only : muexmin, muexmax
        use moleculeslist, only : rho_tilde_min, rho_tilde_max, Diffcoeff_tilde_min, Diffcoeff_tilde_max
        use moleculeslist, only : init_diffusion_coeff, init_vol, init_zval, init_xvolmin
        use inputtemp, only : psizmin, psizmax 
        use molecules, only : vsol
        use system, only : ST_bctype

        real*8 :: pibulk 

        call init_diffusion_coeff()
        call init_vol()
        call init_zval()
        call init_xvolmin()

        xvolmax = xvolmin

        pibulk = -log(xvolmin%sol) 

        mumin%pos = log(xvolmin%pos / vol%pos) +pibulk*vol%pos + psizmin*zval%pos
        mumin%neg = log(xvolmin%neg / vol%neg) +pibulk*vol%neg + psizmin*zval%neg
        mumin%Hplus = log(xvolmin%Hplus / vol%Hplus) +pibulk*vol%Hplus + psizmin*zval%Hplus
        mumin%OHmin = log(xvolmin%OHmin / vol%OHmin) +pibulk*vol%OHmin + psizmin*zval%OHmin

        mumax%pos = log(xvolmax%pos / vol%pos) +pibulk*vol%pos + psizmax*zval%pos
        mumax%neg = log(xvolmax%neg / vol%neg) +pibulk*vol%neg + psizmax*zval%neg
        mumax%Hplus = log(xvolmax%Hplus / vol%Hplus) +pibulk*vol%Hplus + psizmax*zval%Hplus
        mumax%OHmin = log(xvolmax%OHmin / vol%OHmin) +pibulk*vol%OHmin + psizmin*zval%OHmin

       ! if(ST_bctype.eq.2) then ! flux using bc of imposedelectric field 
       ! to be done 
       ! endif
        
        if(ST_bctype.eq.3) then ! flux via Slotboom transformation 
        
            ! zmin 
            
            muexmin%pos = pibulk*vol%pos + psizmin*zval%pos
            muexmin%neg = pibulk*vol%neg + psizmin*zval%neg
            muexmin%Hplus = pibulk*vol%Hplus + psizmin*zval%Hplus
            muexmin%OHmin = pibulk*vol%OHmin + psizmin*zval%OHmin

            rho_tilde_min%pos = xvolmin%pos * exp ( muexmin%pos) /(vol%pos *vsol)
            rho_tilde_min%neg = xvolmin%neg * exp ( muexmin%neg) /(vol%neg *vsol)
            rho_tilde_min%Hplus = xvolmin%Hplus * exp ( muexmin%Hplus) /(vol%Hplus *vsol)
            rho_tilde_min%OHmin = xvolmin%OHmin * exp ( muexmin%OHmin) /(vol%OHmin *vsol)

            Diffcoeff_tilde_min%pos = exp ( -muexmin%pos)
            Diffcoeff_tilde_min%neg = exp ( -muexmin%neg)
            Diffcoeff_tilde_min%Hplus = exp ( -muexmin%Hplus)
            Diffcoeff_tilde_min%OHmin = exp ( -muexmin%OHmin)

            ! zmax 

            muexmax%pos = pibulk*vol%pos + psizmax*zval%pos
            muexmax%neg = pibulk*vol%neg + psizmax*zval%neg
            muexmax%Hplus = pibulk*vol%Hplus + psizmax*zval%Hplus
            muexmax%OHmin = pibulk*vol%OHmin + psizmax*zval%OHmin

            rho_tilde_max%pos = xvolmin%pos * exp ( muexmax%pos) /(vol%pos *vsol)
            rho_tilde_max%neg = xvolmin%neg * exp ( muexmax%neg) /(vol%neg *vsol)
            rho_tilde_max%Hplus = xvolmin%Hplus * exp ( muexmax%Hplus) /(vol%Hplus *vsol)
            rho_tilde_max%OHmin = xvolmin%OHmin * exp ( muexmax%OHmin) /(vol%OHmin *vsol)

            Diffcoeff_tilde_max%pos = exp ( -muexmax%pos)
            Diffcoeff_tilde_max%neg = exp ( -muexmax%neg)
            Diffcoeff_tilde_max%Hplus = exp ( -muexmax%Hplus)
            Diffcoeff_tilde_max%OHmin = exp ( -muexmax%OHmin)
            
        endif

    end subroutine init_flux_var

   ! position dependent chem pot  
   ! \beta \mu_i(r) = \ln(\rho_i(r) v_w) + \beta \pi(r) * v_i + \beta * q_i \psi(r)
   ! \beta \mu_i(r) = \ln(x_i(r) v_w/v_i) - \ln(x_w(r)) * v_i/v_w +  z_i *  e  * \beta \psi(r)
   ! input : real*8 xsol,xvol,psi,iontype)    

    subroutine chem_potential(mu,xsol,xvol,psi,iontype)

        use system, only : dimx, dimy, dimz
        use moleculeslist, only :  vol, zval, get_value_moleclist
        use const,only : stdout 

        ! input arguments 
        real*8, intent(inout) :: mu(0:dimx+1,0:dimy+1,0:dimz+1)
        real*8, intent(in) :: xvol(0:dimx+1,0:dimy+1,0:dimz+1), xsol(:,:,:)
        real*8, intent(in) :: psi(:,:,:)
        character(len=*), intent(in)  :: iontype

       
        ! local variables

        integer :: ix, iy,iz
        character(len=5):: key
        real*8 :: volum, valence
        
        key = trim(iontype)
        volum   = get_value_moleclist(vol,key)
        valence = get_value_moleclist(zval,key)
        
        if(DEBUG_ST) then 
            write(stdout,*)"chem_potential: key=",key  
            write(stdout,*)"chem_potential: volum=",volum
            write(stdout,*)"chem_potential: valence=",valence
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
        real*8, intent(inout) :: xvol(0:dimx+1,0:dimy+1,0:dimz+1)
        real*8, intent(in) :: mu(0:dimx+1,0:dimy+1,0:dimz+1),xsol(:,:,:)
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
        real*8, intent(inout) :: mu(0:dimx+1,0:dimy+1,0:dimz+1)
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


    end subroutine bc_flux

    ! == boundary conditions excess chemical potential and xvol 

    subroutine bc_flux_ex(mu_ex,rho_tilde,Diffcoeff_tilde,iontype) 
        
        use system, only : dimx, dimy, dimz  
        use moleculeslist, only : muexmin, muexmax, rho_tilde_min, rho_tilde_max
        use moleculeslist, only : Diffcoeff_tilde_min,Diffcoeff_tilde_max
        use moleculeslist, only : get_value_moleclist

        ! input arguments 
        real*8, intent(inout) :: mu_ex(0:dimx+1,0:dimy+1,0:dimz+1)
        real*8, intent(inout) :: rho_tilde(0:dimx+1,0:dimy+1,0:dimz+1)
        real*8, intent(inout) :: Diffcoeff_tilde(0:dimx+1,0:dimy+1,0:dimz+1)
        character(len=*), intent(in) :: iontype
          
        ! local variables

        character(len=5):: key
        real*8 :: mu_ex_zmin, mu_ex_zmax, rho_tilde_zmin, rho_tilde_zmax, Diffcoeff_tilde_zmin, Diffcoeff_tilde_zmax

        key = trim(iontype)
        mu_ex_zmin = get_value_moleclist(muexmin,key)
        mu_ex_zmax = get_value_moleclist(muexmax,key)
        rho_tilde_zmin = get_value_moleclist(rho_tilde_min,key)
        rho_tilde_zmax = get_value_moleclist(rho_tilde_max,key)
        Diffcoeff_tilde_zmin = get_value_moleclist(Diffcoeff_tilde_min,key)
        Diffcoeff_tilde_zmax = get_value_moleclist(Diffcoeff_tilde_max,key)

        
        ! boundary conditions excess chemical potential and xvol 

        ! x = 0
        mu_ex(0,:,:) = mu_ex(1,:,:) 
        rho_tilde(0,:,:) = rho_tilde(1,:,:) 

        ! x = dimx
        mu_ex(dimx+1,:,:) = mu_ex(dimx,:,:)  
        rho_tilde(dimx+1,:,:) = rho_tilde(dimx,:,:)  

        ! y = 0
        mu_ex(:,0,:) = mu_ex(:,1,:) 
        rho_tilde(:,0,:) = rho_tilde(:,1,:) 

        ! y = dimy
        mu_ex(:,dimy+1,:) = mu_ex(:,dimy,:) 
        rho_tilde(:,dimy+1,:) = rho_tilde(:,dimy,:) 

        ! z = 0
        mu_ex(:,:,0) = mu_ex_zmin
        rho_tilde(:,:,0) = rho_tilde_zmin

        ! z = dimz
        mu_ex(:,:,dimz+1) = mu_ex_zmax
        rho_tilde(:,:,dimz+1) = rho_tilde_zmax


    end subroutine bc_flux_ex

    ! boundary conditions for potential  

    subroutine bc_psi(psi)
        
        use system, only : dimx, dimy, dimz  
        use inputtemp, only : psizmax, psizmin 
    

        ! input arguments 
        real*8, intent(inout) :: psi(0:dimx+1,0:dimy+1,0:dimz+1)
        
        ! local variables

        ! boundary conditions 

        ! x = 0
        psi(0,:,:) = psi(1,:,:) 
        
        ! x = dimx
        psi(dimx+1,:,:) = psi(dimx,:,:)  

        ! y = 0
        psi(:,0,:) = psi(:,1,:) 

        ! y = dimy
        psi(:,dimy+1,:) = psi(:,dimy,:) 
        
        ! z = 0
        psi(:,:,0) = psizmin
        
        ! z = dimz
        psi(:,:,dimz+1) = psizmax

    end subroutine bc_psi

    subroutine div_flux(divJ,xsol,xion,psi,iontype)

        use system, only : dimx, dimy, dimz
        ! input arguments 
        real*8, intent(inout) :: divJ(:,:,:)
        real*8, intent(in) :: xion(:,:,:), xsol(:,:,:)
        real*8, intent(in) :: psi(0:dimx+1,0:dimy+1,0:dimz+1)
        character(len=*), intent(in)  :: iontype

        call div_flux_channel(divJ,xsol,xion,psi,iontype)
                
    end subroutine div_flux
        
     ! Computes div.J with numerical scheme using Gauss's theorem 

    subroutine div_flux_channel(divJ,xsol,xion,psi,iontype)

        use system, only : dimx, dimy, dimz
        use moleculeslist, only : vol, zval, mumin, mumax, xvolmin, xvolmax
        use moleculeslist, only : get_value_moleclist
        use const, only : stdout
        use MPI, only : rank
        use ematrix, only : fvstdint

        ! input arguments 

        real*8, intent(inout) :: divJ(:,:,:)
        real*8, intent(in) :: xion(:,:,:), xsol(:,:,:)
        real*8, intent(in) ::  psi(0:dimx+1,0:dimy+1,0:dimz+1)
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
            if(rank.eq.0) then 
                write(stdout,*)"div_flux_channel: key= ",key  
                write(stdout,*)"div_flux_channel: volum= ",volum
                write(stdout,*)"div_flux_channel: valence= ",valence
                write(stdout,*)"div_flux_channel: mu_zmin= ",mu_zmin
                write(stdout,*)"div_flux_channel: mu_zpls= ",mu_zpls
                write(stdout,*)"div_flux_channel: xvol_zpls= ",xvol_zpls
                write(stdout,*)"div_flux_channel: xvol_zmin= ",xvol_zmin
                write(stdout,*)"div_flux_channel: size(xvol)= ",size(xvol),&
                    "expected size= ",(dimx+2)*(dimy+2)*(dimz+2)
            endif
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
        if(DEBUG_ST) divJ= 123435600.0000d0 ! == used to detect unassinged values of divJ 
        
    
        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx

                    if(fvstdint(ix,iy,iz).eq.0) then

                        divJ(ix,iy,iz) =  xvol(ix,iy,iz) - 0.0d0 ! == inside channel divJ and xvol zero !!
                    
                    else 
                        
                        if(fvstdint(ix+1,iy,iz).eq.0) then 
                            Jdotxpls = 0.0d0
                        else     
                            Jdotxpls = (xvol(ix+1,iy,iz) + xvol(ix  ,iy,iz))*(mu(ix+1,iy,iz) - mu(ix,iy,iz)  ) ! * Ar (ix,iy,iz,xpls)
                        endif

                        if(fvstdint(ix-1,iy,iz).eq.0) then 
                            Jdotxmin = 0.0d0
                        else
                            Jdotxmin = (xvol(ix  ,iy,iz) + xvol(ix-1,iy,iz))*(mu(ix ,iy,iz)  - mu(ix-1,iy,iz)) ! * Ar (ix,iy,iz,xmin)
                        endif
                         
                        if(fvstdint(ix,iy+1,iz).eq.0) then 
                            Jdotypls = 0.0d0
                        else
                            Jdotypls = (xvol(ix,iy+1,iz) + xvol(ix,iy  ,iz))*(mu(ix,iy+1,iz) - mu(ix,iy,iz)  ) ! * Ar (ix,iy,iz,ypls)
                        endif
                        
                        if(fvstdint(ix,iy-1,iz).eq.0) then
                            Jdotymin = 0.0d0
                        else        
                            Jdotymin = (xvol(ix,iy  ,iz) + xvol(ix,iy-1,iz))*(mu(ix,iy  ,iz) - mu(ix,iy-1,iz)) ! * Ar (ix,iy,iz,ymin)
                        endif

                        if(fvstdint(ix,iy,iz+1).eq.0) then 
                            Jdotzpls = 0.0d0
                        else
                            Jdotzpls = (xvol(ix,iy,iz+1) + xvol(ix,iy,iz  ))*(mu(ix,iy,iz+1) - mu(ix,iy,iz)  ) ! * Ar (ix,iy,iz,zpls)
                        endif
                        if(fvstdint(ix,iy,iz-1).eq.0) then
                            Jdotzmin = 0.0d0
                        else      
                            Jdotzmin = (xvol(ix,iy,iz  ) + xvol(ix,iy,iz-1))*(mu(ix,iy,iz  ) - mu(ix,iy,iz-1)) ! * Ar (ix,iy,iz,zmin)
                        endif    
                
                        divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
                        divJ(ix,iy,iz) = - coeff_scaled * divJtmp
                    
                    endif 

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

    end subroutine div_flux_channel

    ! Computes div.J with numerical scheme using Gauss's theorem 

    subroutine div_flux_cubic(divJ,xsol,xion,psi,iontype)

        use system, only : dimx, dimy, dimz
        use moleculeslist, only : vol, zval, mumin, mumax, xvolmin, xvolmax
        use moleculeslist, only : get_value_moleclist
        use const, only : stdout
        use MPI, only : rank

        ! input arguments 

        real*8, intent(inout) :: divJ(:,:,:)
        real*8, intent(in) :: xion(:,:,:), xsol(:,:,:)
        real*8, intent(inout) ::  psi(0:dimx+1,0:dimy+1,0:dimz+1)
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
            if(rank.eq.0) then 
                write(stdout,*)"div_flux_cubic: key= ",key  
                write(stdout,*)"div_flux_cubic: volum= ",volum
                write(stdout,*)"div_flux_cubic: valence= ",valence
                write(stdout,*)"div_flux_cubic: mu_zmin= ",mu_zmin
                write(stdout,*)"div_flux_cubic: mu_zpls= ",mu_zpls
                write(stdout,*)"div_flux_cubic: xvol_zpls= ",xvol_zpls
                write(stdout,*)"div_flux_cubic: xvol_zmin= ",xvol_zmin
                write(stdout,*)"div_flux_cubic: size(xvol)= ",size(xvol),&
                    "expected size= ",(dimx+2)*(dimy+2)*(dimz+2)
            endif
        endif    

        
        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx
                    xvol(ix,iy,iz) = xion(ix,iy,iz) ! xvol and xion not same range !!!
                    mu(ix,iy,iz) = log(xvol(ix,iy,iz)/volum) -log(xsol(ix,iy,iz))*volum + valence * psi(ix,iy,iz)
                enddo
            enddo
        enddo    
      
        ! == apply boundary conditions to mu and xvol 

        call bc_flux(mu,xvol,iontype)
        ! call bc_psi(psi)  ! need to be applied outside check with bc in fkfun !!!
        
    
        divJ=0.0d0
        if(DEBUG_ST) divJ= 123435600.0000d0 ! == used to detect unassinged values of divJ 
        
        
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
                        if(divJ(ix,iy,iz)== 123435600.0000d0) then 
                            print*,"divJ unassiged in ix=",ix,"iy=",iy,'iy=',iy
                        endif    
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
        real*8, intent(in) ::  psi(0:dimx+1,0:dimy+1,0:dimz+1), xsol(:,:,:)
        character(len=*), intent(in)  :: iontype

        ! local arguments

        real*8 :: grad_mu(dimx,dimy,dimz,3)
        real*8 ::  volum, Diffconst, Jvec0
        integer :: i, ix, iy, iz, idx 
        character(len=5) :: key
 
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
        use results, only : xpos, xneg, xHplus, xOHmin
        use fields_fkfun, only : xsol=>xh, psi
       
        ! return arguments 

        real*8 :: currI

        ! local arguments

        real*8 :: J(dimx, dimy, dimz,3)
        integer :: t
        real*8 :: current(niontypes)
        
        current = 0.0d0

        do t=1,niontypes
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
        enddo

        currI = sum(current)

    end function current_I 


    ! Computes  current through plane z=nz/2 delta  for given flux density J 
    ! for given charge velence zval 

    function current_I_ion(J,iontype) result(current)

        use system, only : dimx, dimy, dimz, delta
        use moleculeslist, only : zval, get_value_moleclist
        

        ! input arguments 

        real*8, intent(in) :: J(:,:,:,:)
        character(len=*), intent(in)  :: iontype

        ! return arguments

        real*8 :: current

        ! local variables

        character(len=5) :: key
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
        real*8, intent(in) ::  psi(0:dimx+1,0:dimy+1,0:dimz+1)
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
        !call bc_psi(psi)  ! need to be done outside bc in fkfun !!!

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

    subroutine unit_test_divJ(info)

        use system, only : delta, dimx, dimy, dimz
        use moleculeslist, only : mumin, mumax, xvolmin, xvolmax
        use moleculeslist, only : get_value_moleclist
        use bulk, only : xsolbulk
        use const, only : stdout
        use inputtemp, only : psizmax, psizmin 
        use MPI, only : rank 

        integer, intent(inout) :: info

        ! local arguments

        integer :: ix,iy,iz
        integer :: un_flux, ios
        character(len=8) :: outfilename

        ! local variable  not in module just for testing
        real*8 :: xion(dimx,dimy,dimz)
        real*8 :: xsol(dimx,dimy,dimz)
        real*8 :: psi(0:dimx+1,0:dimy+1,0:dimz+1)
        real*8 :: divJ(dimx, dimy, dimz)

        real*8 ::  mu_zmin, mu_zpls, xvol_zmin, xvol_zpls
        character(len=5):: key
        real*8 :: slope, intercept , sumdivJ        

        key = "pos"

        mu_zmin = get_value_moleclist(mumin,key)
        mu_zpls = get_value_moleclist(mumax,key)
        xvol_zmin = get_value_moleclist(xvolmin,key)
        xvol_zpls = get_value_moleclist(xvolmax,key)

        if(DEBUG_ST) then
            if(rank.eq.0) then 
                write(stdout,*)"unit_test_divJ : key = ",key  
                write(stdout,*)"unit_test_divJ : mu_zmin = ",mu_zmin
                write(stdout,*)"unit_test_divJ : mu_zpls = ",mu_zpls
                write(stdout,*)"unit_test_divJ : xvol_zpls = ",xvol_zpls
                write(stdout,*)"unit_test_divJ : xvol_zmin = ",xvol_zmin
            endif
        endif    


        slope = (psizmax-psizmin)/((dimz+1.0d0)*delta)
        intercept = psizmin + slope * delta/2.0d0

        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx
                    xsol(ix,iy,iz) = xsolbulk
                    xion(ix,iy,iz) = xvol_zmin
                    psi(ix,iy,iz) =  slope * (iz - 0.5d0) * delta  + intercept
                enddo
            enddo
        enddo        

        ! need to apply boundary. conditions of psi
        call bc_psi(psi) 

        call div_flux_cubic(divJ,xsol,xion,psi,key)
   
        sumdivJ = 0.0d0

        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx
                    sumdivJ=sumdivJ+divJ(ix,iy,iz)**2
                enddo
            enddo
        enddo      

        if(DEBUG_ST) then 
            if(rank.eq.0) then 
                outfilename = "flux.out"
                open(newunit=un_flux,file=outfilename, iostat=ios, action="write")
                do iz=1,dimz
                    do iy=1,dimy
                        do ix=1,dimx
                            write(un_flux,*)ix,iy,iz,divJ(ix,iy,iz)
                        enddo
                    enddo
                enddo          
                close(un_flux)

                outfilename = "psi.out"
                open(newunit=un_flux,file=outfilename, iostat=ios, action="write")
                do iz=0,dimz+1
                    do iy=0,dimy+1
                        do ix=0,dimx+1
                            write(un_flux,*)ix,iy,iz,psi(ix,iy,iz)
                        enddo
                    enddo
                enddo          
                close(un_flux)


            endif
        endif        

        info=0 ! ok 
        if((sumdivJ/(xvol_zmin**2))>0.0001d0) info=1

        write(stdout,*)"unit_test_divJ: unit test: sumdivJ=",sumdivJ

    end subroutine unit_test_divJ


    subroutine linear_interpolation(fcn_interp,fcn_begin,fcn_end)

        use system, only : delta, dimz, ST_bctype

        real*8, intent(inout) :: fcn_interp(:)
        real*8, intent(in) :: fcn_begin, fcn_end

        ! local variable 
        real*8 :: slope, intercept
        integer :: k
        logical :: interpolshift


        if(ST_bctype==2) interpolshift=.true. ! not solution test convergence
        if(ST_bctype==1) interpolshift=.false.

        if(interpolshift) then     
            slope = (fcn_end-fcn_begin)/((dimz+1.0d0)*delta) 
            intercept = fcn_begin+slope * delta/2.0d0
        else 
            slope = ((fcn_end-fcn_begin)/(dimz*delta))
            intercept = fcn_begin 
        endif

        do k=1,dimz 
            fcn_interp(k) = slope * (k - 0.5d0) * delta  +  intercept  ! middle of lattice cell in z-direction
        enddo  

    end subroutine linear_interpolation
    

    ! Computes  a laplace equaton that  indirect representaton div J via a Slotboomn transformations 
    ! div . J  = div. ( D rho grad beta mu ) )<=> div.( Dtilde grad  rhotilde )=0 
    ! with Dtilde = D exp( \beta -mu^ex) and rhotilde = rho exp( + beta mu^ex) 
    ! \beta \mu^ex(r) =  \beta pi(r) v_i +\beta psi(r) q_i :  non-dial part of chemincal potential =
    ! \beta \mu(r) = \mu_0 + \ln(rho_i(r) v_w) + \beta pi(r) v_i +\beta psi(r) q_i

    subroutine div_flux_channel_slotboom(divJ,xsol,xion,psi,iontype)

        use system, only : dimx, dimy, dimz, delta
        use molecules, only : vsol
        use moleculeslist, only : vol, zval, mumin, mumax, xvolmin, xvolmax, Diffcoeff
        use moleculeslist, only : get_value_moleclist
        use const, only : stdout
        use MPI, only : rank
        use ematrix, only : fvstdint
        use inputtemp, only : psizmax, psizmin 

        ! input arguments 

        real*8, intent(inout) :: divJ(:,:,:)
        real*8, intent(in) :: xion(:,:,:), xsol(:,:,:)
        real*8, intent(in) ::  psi(0:dimx+1,0:dimy+1,0:dimz+1)
        character(len=*), intent(in) :: iontype
        

        ! local variables
        integer :: ix, iy, iz
        integer ::  id, idxpls, idxmin, idypls, idymin, idzpls, idzmin
        character(len=5):: key
        real*8 :: volum, valence, coeff_scaled, Diffconst, J0
        real*8 :: mu_zmin, mu_zpls, xvol_zmin, xvol_zpls
        real*8 :: mu_ex_zmin, mu_ex_zpls
        real*8 :: rho_tilde_zmin, rho_tilde_zpls, Diffcoeff_tilde_zmin,Diffcoeff_tilde_zpls
        real*8 :: Jdotxpls, Jdotxmin, Jdotypls, Jdotymin, Jdotzpls, Jdotzmin, divJtmp 


        real*8 :: mu_ex(0:dimx+1, 0:dimy+1, 0:dimz+1)
        real*8 :: xvol(0:dimx+1, 0:dimy+1, 0:dimz+1)
        real*8 :: rho_tilde(0:dimx+1, 0:dimy+1, 0:dimz+1)
        real*8 :: Diffcoeff_tilde(0:dimx+1, 0:dimy+1, 0:dimz+1)

        key = trim(iontype)

        volum   = get_value_moleclist(vol,key)
        valence = get_value_moleclist(zval,key)
        mu_zmin = get_value_moleclist(mumin,key)
        mu_zpls = get_value_moleclist(mumax,key)
        xvol_zmin = get_value_moleclist(xvolmin,key)
        xvol_zpls = get_value_moleclist(xvolmax,key) 
        Diffconst = get_value_moleclist(Diffcoeff,key)

        ! pre factor in flux factor (1.0e-9)^2 arise from vsol in unit of nm^3 and grad_mu in 1/nm
        ! unit of Diffccoeff m^2/sec
        ! J in 1/( nm^2) s)
        
        J0 = Diffconst / ( volum * vsol * (1.0d-9)**2) 
        
        coeff_scaled =  1.0d0/(volum*delta*2.0d0)
        coeff_scaled =  1.0d0
        !coeff_scaled = J0/(delta*2.0_dp)

        !  volum from conversion of volumefraction to density  
        !  in cylinder coordinates : 2 delta^2 from division by 2 \pi delta^2 {\bar r}_i 
        !  in cubic coordiantes :  2 delta  :  from interpolation of density

        if(DEBUG_ST) then
            print*,"DEBUG_ST=",DEBUG_ST
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
                    mu_ex(ix,iy,iz) =-log(xsol(ix,iy,iz))*volum + valence * psi(ix,iy,iz)
                    rho_tilde(ix,iy,iz) = xvol(ix,iy,iz) * exp ( mu_ex(ix,iy,iz)) /(volum *vsol)
                    Diffcoeff_tilde(ix,iy,iz)= exp ( -mu_ex(ix,iy,iz)) 
                enddo
            enddo 
        enddo        

        mu_ex_zmin = -log(xvolmin%sol)*volum + valence * psizmin
        mu_ex_zpls = -log(xvolmax%sol)*volum + valence * psizmax
        rho_tilde_zmin = xvol_zmin * exp ( mu_ex_zmin) /(volum *vsol)
        rho_tilde_zpls = xvol_zpls * exp ( mu_ex_zpls) /(volum *vsol)
        Diffcoeff_tilde_zmin = exp ( -mu_ex_zmin) 
        Diffcoeff_tilde_zpls = exp ( -mu_ex_zpls) 

        ! need to apply PBC to mu_ex


        ! call bc_flux_laplace(mu_ex,rho_tilde) 

        divJ=0.0d0
        
        if(DEBUG_ST) divJ= 123435600.0d0 ! used to detect unassinged values of divJ 
        
        
        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx

                    if(fvstdint(ix,iy,iz).eq.0) then

                        divJ(ix,iy,iz) =  rho_tilde(ix,iy,iz) - 0.0d0 ! == inside channel divJ and xvol zero !!
                    
                    else 
                        
                        if(fvstdint(ix+1,iy,iz).eq.0) then 
                            Jdotxpls = 0.0d0
                        else     
                            Jdotxpls = (Diffcoeff_tilde(ix+1,iy,iz) + Diffcoeff_tilde(ix  ,iy,iz))*(rho_tilde(ix+1,iy,iz) - rho_tilde(ix,iy,iz)  ) ! * Ar (ix,iy,iz,xpls)
                        endif

                        if(fvstdint(ix-1,iy,iz).eq.0) then 
                            Jdotxmin = 0.0d0
                        else
                            Jdotxmin = (Diffcoeff_tilde(ix  ,iy,iz) + Diffcoeff_tilde(ix-1,iy,iz))*(rho_tilde(ix ,iy,iz)  - rho_tilde(ix-1,iy,iz)) ! * Ar (ix,iy,iz,xmin)
                        endif
                         
                        if(fvstdint(ix,iy+1,iz).eq.0) then 
                            Jdotypls = 0.0d0
                        else
                            Jdotypls = (Diffcoeff_tilde(ix,iy+1,iz) + Diffcoeff_tilde(ix,iy  ,iz))*(rho_tilde(ix,iy+1,iz) - rho_tilde(ix,iy,iz)  ) ! * Ar (ix,iy,iz,ypls)
                        endif
                        
                        if(fvstdint(ix,iy-1,iz).eq.0) then
                            Jdotymin = 0.0d0
                        else        
                            Jdotymin = (Diffcoeff_tilde(ix,iy  ,iz) + Diffcoeff_tilde(ix,iy-1,iz))*(rho_tilde(ix,iy  ,iz) - rho_tilde(ix,iy-1,iz)) ! * Ar (ix,iy,iz,ymin)
                        endif

                        if(fvstdint(ix,iy,iz+1).eq.0) then 
                            Jdotzpls = 0.0d0
                        else
                            Jdotzpls = (Diffcoeff_tilde(ix,iy,iz+1) + Diffcoeff_tilde(ix,iy,iz  ))*(rho_tilde(ix,iy,iz+1) - rho_tilde(ix,iy,iz)  ) ! * Ar (ix,iy,iz,zpls)
                        endif
                        if(fvstdint(ix,iy,iz-1).eq.0) then
                            Jdotzmin = 0.0d0
                        else      
                            Jdotzmin = (Diffcoeff_tilde(ix,iy,iz  ) + Diffcoeff_tilde(ix,iy,iz-1))*(rho_tilde(ix,iy,iz  ) - rho_tilde(ix,iy,iz-1)) ! * Ar (ix,iy,iz,zmin)
                        endif    
                
                        divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
                        divJ(ix,iy,iz) = - coeff_scaled * divJtmp
                    
                    endif 

                enddo
            enddo
        enddo    

    end subroutine div_flux_channel_slotboom

end module flux
