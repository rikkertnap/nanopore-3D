module fcnpointer

    implicit none

    abstract interface

        subroutine fcn(x,f,ierfcn)
            implicit none
            
            integer*4, intent(inout) :: ierfcn
            real*8, intent(in) :: x(*)
            real*8, intent(inout) :: f(*)

        end subroutine fcn
    end interface

    procedure(fcn), pointer :: fcnptr => null()

end module fcnpointer

module fcnmod
   
    implicit none

contains

subroutine fcn_fkfun(x,f,ier2)

    use system
    use chainsdat
    use molecules
    use const
    use results
    use bulk
    use kai
    use MPI
    use fields_fkfun
    use kinsol
    use conformations
    use ematrix
    use ellipsoid
    use transform
    use kaist
    use mparameters_monomer
    use flux, only : divJ, div_flux, niontypes, iontype

    implicit none

    ! input arguments 
    integer*4, intent(inout) :: ier2
    real*8, intent(in) :: x(*)
    real*8, intent(inout) :: f(*)

    ! local arguments

    !  integer :: ncells
    real*8 :: protemp
    integer :: i,j, ix, iy, iz, ii, ax, ay, az
    integer :: im, ip
    integer :: jx, jy, jz, jj
    real*8 :: xpot(dimx, dimy, dimz, N_monomer)
    integer :: id, noffset ! == indices used for flux contrubution to f
    real*8 :: normvol, normel 

    ! Charge
    real*8 :: psitemp
    real*8 :: MV(3),MU(3),MW(3)
    real*8 :: MVV,MUU,MWW,MVU,MVW,MUW
    real*8 :: psivv,psiuu,psiww, psivu,psivw,psiuw
    real*8 :: psiv(3), epsv(3)
    real*8 :: xtotalsum(dimx,dimy,dimz)

    integer, external :: PBCSYMI, PBCREFI

    ! poor solvent 
    real*8 :: sttemp
    ! MPI
    integer :: tag
    parameter(tag = 0)
    integer :: err
    real*8 :: avpol_temp(dimx,dimy,dimz,N_monomer)
    real*8 :: q_tosend
    real*8 :: gradpsi2
    real*8 :: fv

    ! hamiltonian inception
    real*8 :: hfactor, hd
    real*8 :: hds(100)


    hds = -1

    !-----------------------------------------------------
    ! Common variables

    shift = 1.0

    ! ncells = dimx*dimy*dimz ! numero de celdas == number of cells

    ! Jefe

    if(rank.eq.0) then ! llama a subordinados y pasa vector x  == calls subordinates and passes vector x
        flagsolver = 1
        CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
        CALL MPI_BCAST(x, neqs , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
    endif

    !------------------------------------------------------
    ! DEBUG
    !      if(iter.gt.2000) then
    !      do i = 1, n
    !      write(stdout,*)i, x(i)
    !      enddo
    !      endif


    ! Recupera xh y psi desde x() 
    ! == Retrieve xh and psi from x()

    !  psi = 0.0 
    psi = 0.0d0 ! == without d0 significant number loss can ocur but value overwritten below !!!

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                xh(ix,iy,iz)=x(ix+dimx*(iy-1)+dimx*dimy*(iz-1))  !fraccion solvente == solvent 

                do ip = 1, N_poorsol
                    xtotal(ix,iy,iz,ip) = x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ ip*ncells) !fraccion polimero de tipo ip
                enddo
                if(electroflag.eq.1) psi(ix,iy,iz)=x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsol+1)*ncells)   !potencial electrostatico
               
                if(fluxflag.eq.1) then 

                    xpos(ix, iy, iz)  =  x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsol+1)*ncells +     ncells)    
                    xneg(ix, iy, iz)  =  x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsol+1)*ncells + 2 * ncells) 
                    xHplus(ix, iy, iz)=  x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsol+1)*ncells + 3 * ncells)    
                    xOHmin(ix, iy, iz)=  x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsol+1)*ncells + 4 * ncells)    
                endif
            enddo
        enddo
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    ! Boundary conditions electrostatic potential
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Reflection or PBC, (PBC = 1 or 3)
    
    do jx = 0, dimx+1
        do jy = 0, dimy+1
            do jz = 0, dimz+1

                ix=jx
                iy=jy
                iz=jz ! these lines are necessary for PBC = 0 or 2

                if (PBC(1).eq.1)ix = PBCSYMI(jx,dimx)
                if (PBC(3).eq.1)iy = PBCSYMI(jy,dimy)
                if (PBC(5).eq.1)iz = PBCSYMI(jz,dimz)

                if (PBC(1).eq.3)ix = PBCREFI(jx,dimx)
                if (PBC(3).eq.3)iy = PBCREFI(jy,dimy)
                if (PBC(5).eq.3)iz = PBCREFI(jz,dimz)

                psi(jx, jy, jz) = psi(ix, iy, iz)
            enddo
        enddo
    enddo

    ! Bulk or Wall, PBC = 0 or 2

    select case (PBC(1)) ! x = 0
    case(0) ! set bulk 
        psi(0,:,:) = 0.0d0        ! == added  d0 
    case(2)
        psi(0,:,:) = psi(1,:,:) ! zero charge
    end select

    select case (PBC(2)) ! x = dimx
    case(0) ! set bulk 
        psi(dimx+1,:,:) = 0.0d0  
    case(2)
        psi(dimx+1,:,:) = psi(dimx,:,:) ! zero charge
    end select

    select case (PBC(3)) ! y = 0
    case(0) ! set bulk 
        psi(:,0,:) = 0.0d0  
    case(2)
        psi(:,0,:) = psi(:,1,:) ! zero charge
    end select

    select case (PBC(4)) ! y = dimy
    case(0) ! set bulk 
        psi(:,dimy+1,:) = 0.0d0
    case(2)
        psi(:,dimy+1,:) = psi(:,dimy,:) ! zero charge
    end select

    select case (PBC(5)) ! z = 0
    case(0) ! set bulk 
        psi(:,:,0) = 0.0d0  
    case(2)
        psi(:,:,0) = psi(:,:,1) ! zero charge
    end select

    select case (PBC(6)) ! z = dimz
    case(0) ! set bulk 
        psi(:,:,dimz+1) = 0.0d0
    case(2)
        psi(:,:,dimz+1) = psi(:,:,dimz) ! zero charge
    end select

    ! volume fraction and frdir

    fdis = 0.0d0    ! == added  d0 
    avpol = 0.0d0

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                if(fluxflag.eq.0) then ! Equilibrium  

                    xpos(ix, iy, iz)   = expmupos*(xh(ix, iy, iz)**vsalt)*dexp(-psi(ix, iy, iz)*zpos) ! ion plus volume fraction vsalt=vsal/vsv
                    xneg(ix, iy, iz)   = expmuneg*(xh(ix, iy, iz)**vsalt)*dexp(-psi(ix, iy, iz)*zneg) ! ion neg volume fraction
                    xHplus(ix, iy, iz) = expmuHplus*(xh(ix, iy, iz))*dexp(-psi(ix, iy, iz))           ! H+ volume fraction
                    xOHmin(ix, iy, iz) = expmuOHmin*(xh(ix,iy,iz))*dexp(+psi(ix,iy,iz))               ! OH-  volume fraction
                
                endif 

                do im =1,N_monomer

                    if (zpol(im).eq.1) then !BASE
                        fdis(ix,iy,iz,im) = 1.0d0 /(1.0d0 + xOHmin(ix,iy,iz)/(K0(im)*xh(ix,iy,iz))) !k0 k en fraccion de volumen
                    else if (zpol(im).eq.-1) then !ACID
                        fdis(ix,iy,iz,im) = 1.0d0 /(1.0d0 + xHplus(ix,iy,iz)/(K0(im)*xh(ix,iy,iz)))
                    endif

                enddo

            enddo
        enddo  
    enddo



    ! Compute dielectric permitivity

    xtotalsum = 0.0d0 ! sum of all polymers
    do ip = 1, N_poorsol
        xtotalsum(:,:,:) = xtotalsum(:,:,:) + xtotal(:,:,:,ip)
    enddo
    
    call dielectfcn(xtotalsum,volprot,epsfcn,Depsfcn)

    !------------------------------------------------------------------------
    ! PDFs polimero // == polymer 
    !------------------------------------------------------------------------
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! PARALELO: Cada procesador trabaja sobre una cadena... /== Each processor works on a chain...
    !         
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Calcula xpot

    sttemp = st/(vpol*vsol)

    do im = 1, N_monomer ! loop over different monomer types

        do ix=1,dimx
            do iy=1,dimy
                do iz=1,dimz

                    if(hguess .eq. 0) then

                        hd = sqrt(float((2*ix-dimx)**2+(2*iy-dimy)**2))/2.0*delta
                        hd = hd**2+(oval*float(2*iz-dimz)/2.0*delta)**2
                        hfactor = dexp(-(kp**2)*hd)

                    elseif(hguess .eq. 1) then

                        hd = sqrt(float((2*ix-dimx)**2+(2*iy-dimy)**2))/2.0*delta-hring
                        hd = hd**2+(oval*float(2*iz-dimz)/2.0*delta)**2
                        hfactor = dexp(-(kp**2)*hd)

                    else

                        do i=1,hguess
                            hds(i) = (float(2*ix-dimx)-2*cos(i*2*pi/hguess)*hring/delta)**2+(float(2*iy-dimy)-&
                                2*sin(i*2*pi/hguess)*hring/delta)**2
                            hds(i) = hds(i)/4.0*(delta**2)+(oval*float(2*iz-dimz)/2.0*delta)**2
                        end do
                        hd = minval(hds, mask = hds .gt.0)
                        hfactor = dexp(-(kp**2)*hd)

                    end if

                    !   xpot=exp(-Uj(rj)) para P(alpha)
                    !   == xpot is the exponent of the effective Energy function in P(alpha)

                    fv = (1.0 - volprot(ix,iy,iz)) 
                    !   fraccion de volumen de la celda que es sc volprot->fraccion pared
                    !   == volume fraction of the cell that is sc volprot->wall fraction             

                    xpot(ix, iy, iz, im) = xh(ix,iy,iz)**vpol 
                    
                    ! im:tipo de segmento, término de presion osmotica
                    ! == im: segment type, osmotic pressure term

                    xpot(ix, iy, iz, im) = xpot(ix,iy,iz, im)*dexp(voleps(ix,iy,iz))  
                    
                    ! termino de interaccion con sup de la particula
                    ! == interaction term with the particle's sup =surface ???

                    ! Electrostatics

                    if(zpol(im).ne.0.0) then
                        xpot(ix,iy,iz,im) =  xpot(ix,iy,iz,im)/fdis(ix,iy,iz,im)*dexp(-psi(ix,iy,iz)*zpol(im))  
                        ! fdis: por eq ac. base...  
                    endif
        
                    ! Dielectrics

                    gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy,iz))**2+&
                        (psi(ix,iy,iz+1)-psi(ix,iy,iz))**2 

                    !     gradpsi2 = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy-1,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz-1))**2 
                    !     xpot(ix, iy, iz) = xpot(ix,iy,iz)*exp(-Depsfcn(ix,iy,iz)*(gradpsi2)*constqE)

                    xpot(ix,iy,iz,im) = xpot(ix,iy,iz,im)*exp(Depsfcn(ix,iy,iz)*(gradpsi2)/constq/2.0d0*vpol/fv)

                    ! Poor solvent depende de la grilla donde esta y de sus vecinos
                    ! == Poor solvent depends on the grid where it is and its neighbors



                    if(hydroph(im).ne.0) then

                    protemp=0.0

                    do ax = -Xulimit,Xulimit 
                        do ay = -Xulimit,Xulimit
                            do az = -Xulimit,Xulimit

                                jx = ix+ax
                                jy = iy+ay
                                jz = iz+az

                                if(jx.lt.1) then
                                    if(PBC(1).eq.1)jx = PBCSYMI(jx,dimx)
                                    if(PBC(1).eq.3)jx = PBCREFI(jx,dimx)
                                endif

                                if(jx.gt.dimx) then
                                    if(PBC(2).eq.1)jx = PBCSYMI(jx,dimx)
                                    if(PBC(2).eq.3)jx = PBCREFI(jx,dimx)
                                endif

                                if(jy.lt.1) then
                                    if(PBC(3).eq.1)jy = PBCSYMI(jy,dimy)
                                    if(PBC(3).eq.3)jy = PBCREFI(jy,dimy)
                                endif

                                if(jy.gt.dimy) then
                                    if(PBC(4).eq.1)jy = PBCSYMI(jy,dimy)
                                    if(PBC(4).eq.3)jy = PBCREFI(jy,dimy)
                                endif


                                if(jz.lt.1) then
                                    if(PBC(5).eq.1)jz = PBCSYMI(jz,dimz)
                                    if(PBC(5).eq.3)jz = PBCREFI(jz,dimz)
                                endif

                                if(jz.gt.dimz) then
                                    if(PBC(6).eq.1)jz = PBCSYMI(jz,dimz)
                                    if(PBC(6).eq.3)jz = PBCREFI(jz,dimz)
                                endif


                                if((jx.ge.1).and.(jx.le.dimx)) then
                                    if((jy.ge.1).and.(jy.le.dimy)) then
                                        if((jz.ge.1).and.(jz.le.dimz)) then
                                            fv = (1.0-volprot(jx,jy,jz))

                                            do ip = 1, N_poorsol
                                                protemp = protemp + hfactor*Xu(ax,ay,az)*&
                                                st_matrix(hydroph(im),ip)*sttemp*xtotal(jx,jy,jz,ip)*fv
                                            enddo ! ip

                                        endif
                                    endif
                                endif

                            enddo
                        enddo
                    enddo

                    xpot(ix,iy,iz,im) = xpot(ix,iy,iz,im)*dexp(protemp)

                    endif ! hydrph

!                    write(567,*)xpot(ix,iy,iz,1)

                enddo ! ix
            enddo ! iy
        enddo !iz

    enddo ! N_monomer


    !!!!!!!!!!!!!!!!!!!!!! Calculate pro from xpot !!!!!!!!!!!!!
    call calcavpol(xpot)
    
    !!!!!!!!!!! IMPORTANTE, LOS SUBORDINADOS TERMINAN ACA... 
    if(rank.ne.0)goto 3333
    !!!!!!!!!!!!!!!!!!!!!!! FIN MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !----------------------------------------------------------------------------------------------
    !   Construye Ecuaciones a resolver 
    !    == Build Equations to solve
    !----------------------------------------------------------------------------------------------

    ! Qtot


    qtot = 0.0d0 ! == added  d0 

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz
            
            fv = (1.0-volprot(ix,iy,iz))     ! fv == free volume 

            qtot(ix, iy, iz) =  (zpos*xpos(ix, iy, iz)+zneg*xneg(ix, iy, iz))/vsalt + &
                xHplus(ix, iy, iz) - xOHmin(ix, iy, iz)

            do im = 1, N_monomer
                qtot(ix, iy, iz) =  qtot(ix,iy,iz) + avpol(ix,iy,iz,im)*zpol(im)/vpol*fdis(ix,iy,iz,im)
            enddo

            qtot(ix, iy,iz) = qtot(ix,iy,iz)*fv + volq(ix,iy,iz)*vsol    ! OJO

            enddo
        enddo
    enddo

    ! Volume fraction
    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                f(ix+dimx*(iy-1)+dimx*dimy*(iz-1))= xh(ix,iy,iz) + &
                    xneg(ix, iy, iz) + xpos(ix, iy, iz) + xHplus(ix, iy, iz) + &
                    xOHmin(ix, iy, iz) -1.000000d0  !packing iones+sv

                do im = 1, N_monomer
                    f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) = f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) +&
                         avpol(ix,iy,iz,im) !packing ...+polimero
                enddo

                ! write(123,*)ix,iy,iz,avpol(ix,iy,iz,1),xh(ix,iy,iz)
            enddo
        enddo
    enddo

    ! Poor solvent

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                do ip = 1, N_poorsol
                    f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells) = xtotal(ix,iy,iz,ip)

                    do im = 1, N_monomer
                        if(hydroph(im).eq.ip) then 
                            f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells) = f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells) - &
                                avpol(ix,iy,iz,im)
                        endif
                    enddo ! im
                enddo ! ip

            enddo ! ix
        enddo ! iy
    enddo ! iz


    if(electroflag.eq.1) then

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Poisson equation
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        ! Some auxialiary variables, see Notes Poisson eq. non-cubic grid
        !

        MV(1) = MAT(1,1)
        MV(2) = MAT(1,2)  
        MV(3) = MAT(1,3)

        MU(1) = MAT(2,1)
        MU(2) = MAT(2,2)  
        MU(3) = MAT(2,3)

        MW(1) = MAT(3,1)
        MW(2) = MAT(3,2)  
        MW(3) = MAT(3,3)

        MVV = DOT_PRODUCT(MV,MV)
        MUU = DOT_PRODUCT(MU,MU)
        MWW = DOT_PRODUCT(MW,MW)

        MVU = DOT_PRODUCT(MV,MU)
        MVW = DOT_PRODUCT(MV,MW)
        MUW = DOT_PRODUCT(MU,MW)

        do ix=1,dimx
            do iy=1,dimy
                do iz=1,dimz

                    psivv = psi(ix+1,iy,iz)-2*psi(ix,iy,iz)+psi(ix-1,iy,iz)
                    psiuu = psi(ix,iy+1,iz)-2*psi(ix,iy,iz)+psi(ix,iy-1,iz)
                    psiww = psi(ix,iy,iz+1)-2*psi(ix,iy,iz)+psi(ix,iy,iz-1)

                    psivu = (psi(ix+1,iy+1,iz)+psi(ix-1,iy-1,iz)-psi(ix+1,iy-1,iz)-psi(ix-1,iy+1,iz))/4.0d0
                    psivw = (psi(ix+1,iy,iz+1)+psi(ix-1,iy,iz-1)-psi(ix+1,iy,iz-1)-psi(ix-1,iy,iz+1))/4.0d0
                    psiuw = (psi(ix,iy+1,iz+1)+psi(ix,iy-1,iz-1)-psi(ix,iy+1,iz-1)-psi(ix,iy-1,iz+1))/4.0d0

                    psiv(1) = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))/2.0d0
                    psiv(2) = (psi(ix,iy+1,iz)-psi(ix,iy-1,iz))/2.0d0
                    psiv(3) = (psi(ix,iy,iz+1)-psi(ix,iy,iz-1))/2.0d0

                    epsv(1) = (epsfcn(ix+1,iy,iz)-epsfcn(ix-1,iy,iz))/2.0d0
                    epsv(2) = (epsfcn(ix,iy+1,iz)-epsfcn(ix,iy-1,iz))/2.0d0
                    epsv(3) = (epsfcn(ix,iy,iz+1)-epsfcn(ix,iy,iz-1))/2.0d0

                    psitemp = epsfcn(ix,iy,iz)*&
                        (MVV*psivv+MUU*psiuu+MWW*psiww+2.0d0*MVU*psivu+2.0d0*MVW*psivw+2.0d0*MUW*psiuw)
                    psitemp = psitemp + DOT_PRODUCT(MATMUL(TMAT,epsv),MATMUL(TMAT,psiv))

                    ! OJO CHECK!!!!

                    f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsol+1)*ncells)=&
                        (psitemp + qtot(ix, iy, iz)*constq)/(-2.0d0)

                enddo
            enddo
        enddo

    endif ! electroflag

    if(fluxflag.eq.1) then
        do i=1,niontypes 
            select case (iontype(i))
            case ("Hplus")
                call div_flux(divJ(:,:,:,i),xh,xHplus,psi,"Hplus")
            case ( "OHmin") 
                call div_flux(divJ(:,:,:,i),xh,xOHmin,psi,"OHmin") 
            case("pos")
                call div_flux(divJ(:,:,:,i),xh,xpos,psi,"pos")
            case("neg")
                call div_flux(divJ(:,:,:,i),xh,xneg,psi,"neg")
            case default
                write(stdout,*)"fkfun: error reached unspecified iontype"
                stop 
            end select
        enddo    

        noffset=(N_poorsol+2)*ncells 
        

        do i=1,niontypes
            do iz=1,dimz
                do iy=1,dimy
                    do ix=1,dimx
                        id=ix+dimx*(iy-1)+dimx*dimy*(iz-1)+noffset +(i-1)*ncells 
                        f(id)  =  divJ(ix,iy,iz,i)
                    enddo
                enddo    
            enddo    
        enddo
        
    endif !fluxflag 
 
    ! == norma = 0.0 
    norma = 0.0d0      ! == added  d0  without d0 significant number loss can occur
    normvol = 0.0d0
    normel = 0.0d0

    do i = 1, eqs*ncells
        norma = norma + (f(i))**2
    enddo

    do i = 1,ncells
        normvol= normvol +f(i)**2
    enddo
    if(electroflag.eq.1) then 
        noffset=(N_poorsol+1)*ncells 
        do i= 1, ncells    
            normel = normel +f(i+noffset)**2
        enddo       
    endif
    
    iter = iter + 1
    if(verbose.ge.3) then
        if(rank.eq.0) write(stdout,*)'fkfun:', iter, sqrt(norma), sqrt(normvol), sqrt(normel), q(1)
    endif

    3333 continue
    ier2 = 0.0 

    return

end subroutine fcn_fkfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! == fcn for flux : interation in ion densities 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine fcn_flux(x,f,ier2)

    use system
    use chainsdat
    use molecules
    use const
    use results
    use bulk
    use kai
    use MPI
    use fields_fkfun
    use kinsol
    use conformations
    use ematrix
    use ellipsoid
    use transform
    use kaist
    use mparameters_monomer
    use flux, only : divJ, div_flux, niontypes, iontype, bc_psi 
    use moleculeslist, only : xvolmin, get_value_moleclist
    use inputtemp, only : psizmin, psizmax

    implicit none

    ! input arguments 
    integer*4, intent(inout) :: ier2
    real*8, intent(in) :: x(*)
    real*8, intent(inout) :: f(*)

    ! local arguments

    ! integer :: ncells
    real*8 :: protemp
    integer :: i,j, ix, iy, iz, ii, ax, ay, az
    integer :: idx, im, ip
    integer :: jx, jy, jz, jj
    real*8 :: xpot(dimx, dimy, dimz, N_monomer)
    integer :: id, noffset ! == indices used for flux contribution to f
    real*8 :: normvol, normel,  normflux , normfluxpos

    ! Charge
    real*8 :: psitemp
    real*8 :: MV(3),MU(3),MW(3)
    real*8 :: MVV,MUU,MWW,MVU,MVW,MUW
    real*8 :: psivv,psiuu,psiww, psivu,psivw,psiuw
    real*8 :: psiv(3), epsv(3)
    real*8 :: xtotalsum(dimx,dimy,dimz)

    integer, external :: PBCSYMI, PBCREFI

    ! poor solvent 
    real*8 :: sttemp
    ! MPI
    integer :: tag
    parameter(tag = 0)
    integer :: err

    real*8 :: avpol_temp(dimx,dimy,dimz,N_monomer)
    real*8 :: q_tosend
    real*8 :: gradpsi2
    real*8 :: fv

    ! hamiltonian inception
    real*8 :: hfactor, hd
    real*8 :: hds(100)
    ! flux 
    character(len=5):: key
    real*8 :: xvolbulk

    hds = -1

    ! Common variables

    shift = 1.0d0

    ! ncells = dimx*dimy*dimz !  == number of cells

    ! == head nore 

    if(rank.eq.0) then ! == calls subordinates and passes vector x
        flagsolver = 1
        CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
        CALL MPI_BCAST(x, eqs*ncells , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
    endif

    
    ! == Retrieve xh and psi from x()

    !  psi = 0.0 
    psi = 0.0d0 ! == without d0 significant number loss can ocur but value overwritten below !!!

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz
                idx = ix+dimx*(iy-1)+dimx*dimy*(iz-1)

                xh(ix,iy,iz)=x(idx)  !  == solvent 

                do ip = 1, N_poorsol
                    xtotal(ix,iy,iz,ip) = x(idx+ ip*ncells) ! == fraction  polymer of type ip
                enddo
                
                psi(ix,iy,iz)=x(idx+(N_poorsol+1)*ncells)   ! == potential 
               
                xpos(ix, iy, iz)  =  x(idx+(N_poorsol+2)*ncells )    
                xneg(ix, iy, iz)  =  x(idx+(N_poorsol+3)*ncells ) 
                xHplus(ix, iy, iz)=  x(idx+(N_poorsol+4)*ncells )    
                xOHmin(ix, iy, iz)=  x(idx+(N_poorsol+5)*ncells )    
                
            enddo
        enddo
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    ! Boundary conditions electrostatic potential
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Reflection or PBC, (PBC = 1 or 3)
    
    if(fluxflag.eq.1) then 

        call bc_psi(psi)

    else 

        do jx = 0, dimx+1
            do jy = 0, dimy+1
                do jz = 0, dimz+1

                    ix=jx
                    iy=jy
                    iz=jz ! these lines are necessary for PBC = 0 or 2

                    if (PBC(1).eq.1)ix = PBCSYMI(jx,dimx)
                    if (PBC(3).eq.1)iy = PBCSYMI(jy,dimy)
                    if (PBC(5).eq.1)iz = PBCSYMI(jz,dimz)

                    if (PBC(1).eq.3)ix = PBCREFI(jx,dimx)
                    if (PBC(3).eq.3)iy = PBCREFI(jy,dimy)
                    if (PBC(5).eq.3)iz = PBCREFI(jz,dimz)

                    psi(jx, jy, jz) = psi(ix, iy, iz)
                enddo
            enddo
        enddo
    
        ! Bulk or Wall, PBC = 0 or 2

        select case (PBC(1)) ! x = 0
        case(0) ! set bulk 
            psi(0,:,:) = 0.0d0        ! == added  d0 
        case(2)
            psi(0,:,:) = psi(1,:,:) ! zero charge
        end select

        select case (PBC(2)) ! x = dimx
        case(0) ! set bulk 
            psi(dimx+1,:,:) = 0.0d0  
        case(2)
            psi(dimx+1,:,:) = psi(dimx,:,:) ! zero charge
        end select

        select case (PBC(3)) ! y = 0
        case(0) ! set bulk 
            psi(:,0,:) = 0.0d0  
        case(2)
            psi(:,0,:) = psi(:,1,:) ! zero charge
        end select

        select case (PBC(4)) ! y = dimy
        case(0) ! set bulk 
            psi(:,dimy+1,:) = 0.0d0
        case(2)
            psi(:,dimy+1,:) = psi(:,dimy,:) ! zero charge
        end select

        select case (PBC(5)) ! z = 0
        case(0) ! set bulk 
            psi(:,:,0) = 0.0d0  
        case(2)
            psi(:,:,0) = psi(:,:,1) ! zero charge
        end select

        select case (PBC(6)) ! z = dimz
        case(0) ! set bulk 
            psi(:,:,dimz+1) = 0.0d0
        case(2)
            psi(:,:,dimz+1) = psi(:,:,dimz) ! zero charge
        end select

    endif    

    ! volume fraction and fdis

    fdis = 0.0d0    ! == added  d0 
    avpol = 0.0d0

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                do im =1,N_monomer

                    if (zpol(im).eq.1) then !BASE
                        fdis(ix,iy,iz,im) = 1.0d0 /(1.0d0 + xOHmin(ix,iy,iz)/(K0(im)*xh(ix,iy,iz))) !k0 k en fraccion de volumen
                    else if (zpol(im).eq.-1) then !ACID
                        fdis(ix,iy,iz,im) = 1.0d0 /(1.0d0 + xHplus(ix,iy,iz)/(K0(im)*xh(ix,iy,iz)))
                    endif

                enddo

            enddo
        enddo  
    enddo



    ! Compute dielectric permitivity

    xtotalsum = 0.0d0 ! sum of all polymers
    do ip = 1, N_poorsol
        xtotalsum(:,:,:) = xtotalsum(:,:,:) + xtotal(:,:,:,ip)
    enddo
    
    call dielectfcn(xtotalsum,volprot,epsfcn,Depsfcn)

    !------------------------------------------------------------------------
    ! PDF polymer 
    !------------------------------------------------------------------------
    !
    ! parallel Each processor works on a chain
    !
    !-----------------------------------------------------------------------         
   

    ! Calculate xpot

    sttemp = st/(vpol*vsol)

    do im = 1, N_monomer ! loop over different monomer types

        do ix=1,dimx
            do iy=1,dimy
                do iz=1,dimz

                    if(hguess .eq. 0) then

                        hd = sqrt(float((2*ix-dimx)**2+(2*iy-dimy)**2))/2.0*delta
                        hd = hd**2+(oval*float(2*iz-dimz)/2.0*delta)**2
                        hfactor = dexp(-(kp**2)*hd)

                    elseif(hguess .eq. 1) then

                        hd = sqrt(float((2*ix-dimx)**2+(2*iy-dimy)**2))/2.0*delta-hring
                        hd = hd**2+(oval*float(2*iz-dimz)/2.0*delta)**2
                        hfactor = dexp(-(kp**2)*hd)

                    else

                        do i=1,hguess
                            hds(i) = (float(2*ix-dimx)-2*cos(i*2*pi/hguess)*hring/delta)**2+(float(2*iy-dimy)-&
                                2*sin(i*2*pi/hguess)*hring/delta)**2
                            hds(i) = hds(i)/4.0*(delta**2)+(oval*float(2*iz-dimz)/2.0*delta)**2
                        end do
                        hd = minval(hds, mask = hds .gt.0)
                        hfactor = dexp(-(kp**2)*hd)

                    end if

                    !   == xpot = exp(-Uj(rj)) para P(alpha)
                    !   == xpot is the exponent of the effective Energy function in P(alpha)

                    fv = (1.0 - volprot(ix,iy,iz)) 
    
                    !   == volume fraction of the cell that is sc volprot->wall fraction             

                    xpot(ix, iy, iz, im) = xh(ix,iy,iz)**vpol 
                    
                    ! == im: segment type, osmotic pressure term

                    xpot(ix, iy, iz, im) = xpot(ix,iy,iz, im)*dexp(voleps(ix,iy,iz))  
                    
                    ! == interaction term with the particle's sup =surface 

                    ! == Electrostatics

                    if(zpol(im).ne.0.0) then
                        xpot(ix,iy,iz,im) =  xpot(ix,iy,iz,im)/fdis(ix,iy,iz,im)*dexp(-psi(ix,iy,iz)*zpol(im))  
                        ! == fdis 
                    endif
        
                    ! == Dielectrics

                    gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy,iz))**2+&
                        (psi(ix,iy,iz+1)-psi(ix,iy,iz))**2 

                    ! ==  gradpsi2 = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy-1,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz-1))**2 
                    ! ==  xpot(ix, iy, iz) = xpot(ix,iy,iz)*exp(-Depsfcn(ix,iy,iz)*(gradpsi2)*constqE)

                    xpot(ix,iy,iz,im) = xpot(ix,iy,iz,im)*exp(Depsfcn(ix,iy,iz)*(gradpsi2)/constq/2.0d0*vpol/fv)

                    ! == Poor solvent depends on the grid where it is and its neighbors

                    if(hydroph(im).ne.0) then

                    protemp=0.0

                    do ax = -Xulimit,Xulimit 
                        do ay = -Xulimit,Xulimit
                            do az = -Xulimit,Xulimit

                                jx = ix+ax
                                jy = iy+ay
                                jz = iz+az

                                if(jx.lt.1) then
                                    if(PBC(1).eq.1)jx = PBCSYMI(jx,dimx)
                                    if(PBC(1).eq.3)jx = PBCREFI(jx,dimx)
                                endif

                                if(jx.gt.dimx) then
                                    if(PBC(2).eq.1)jx = PBCSYMI(jx,dimx)
                                    if(PBC(2).eq.3)jx = PBCREFI(jx,dimx)
                                endif

                                if(jy.lt.1) then
                                    if(PBC(3).eq.1)jy = PBCSYMI(jy,dimy)
                                    if(PBC(3).eq.3)jy = PBCREFI(jy,dimy)
                                endif

                                if(jy.gt.dimy) then
                                    if(PBC(4).eq.1)jy = PBCSYMI(jy,dimy)
                                    if(PBC(4).eq.3)jy = PBCREFI(jy,dimy)
                                endif


                                if(jz.lt.1) then
                                    if(PBC(5).eq.1)jz = PBCSYMI(jz,dimz)
                                    if(PBC(5).eq.3)jz = PBCREFI(jz,dimz)
                                endif

                                if(jz.gt.dimz) then
                                    if(PBC(6).eq.1)jz = PBCSYMI(jz,dimz)
                                    if(PBC(6).eq.3)jz = PBCREFI(jz,dimz)
                                endif


                                if((jx.ge.1).and.(jx.le.dimx)) then
                                    if((jy.ge.1).and.(jy.le.dimy)) then
                                        if((jz.ge.1).and.(jz.le.dimz)) then
                                            fv = (1.0-volprot(jx,jy,jz))

                                            do ip = 1, N_poorsol
                                                protemp = protemp + hfactor*Xu(ax,ay,az)*&
                                                st_matrix(hydroph(im),ip)*sttemp*xtotal(jx,jy,jz,ip)*fv
                                            enddo ! ip

                                        endif
                                    endif
                                endif

                            enddo
                        enddo
                    enddo

                    xpot(ix,iy,iz,im) = xpot(ix,iy,iz,im)*dexp(protemp)

                    endif ! hydrph

!                    write(567,*)xpot(ix,iy,iz,1)

                enddo ! ix
            enddo ! iy
        enddo !iz

    enddo ! N_monomer


    !!!!!!!!!!!!!!!!!!!!!! Calculate pro from xpot !!!!!!!!!!!!!
    call calcavpol(xpot)
    
    !!!!!!!!!!! MPORTANT, THE SUBORDINATES' TASKS END HERE ... 
    if(rank.ne.0) goto 2222
    !!!!!!!!!!!!!!!!!!!!!!! END  MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !----------------------------------------------------------------------------------------------
    !    == Build Equations to solve
    !----------------------------------------------------------------------------------------------

    ! Qtot


    qtot = 0.0d0 ! == added  d0 

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz
            
                fv = (1.0-volprot(ix,iy,iz))     ! fv == free volume 

                qtot(ix, iy, iz) =  (zpos*xpos(ix, iy, iz)+zneg*xneg(ix, iy, iz))/vsalt + &
                    xHplus(ix, iy, iz) - xOHmin(ix, iy, iz)

                do im = 1, N_monomer
                    qtot(ix, iy, iz) = qtot(ix,iy,iz) + avpol(ix,iy,iz,im)*zpol(im)/vpol*fdis(ix,iy,iz,im)
                enddo

                qtot(ix, iy,iz) = qtot(ix,iy,iz)*fv + volq(ix,iy,iz)*vsol   

            enddo
        enddo
    enddo

    ! Volume fraction
    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                f(ix+dimx*(iy-1)+dimx*dimy*(iz-1))= xh(ix,iy,iz) + &
                    xneg(ix, iy, iz) + xpos(ix, iy, iz) + xHplus(ix, iy, iz) + &
                    xOHmin(ix, iy, iz) -1.0d0   ! packing iones+sv

                do im = 1, N_monomer
                    f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) = f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) +&
                         avpol(ix,iy,iz,im)          ! packing ...+polymer
                enddo

                ! write(123,*)ix,iy,iz,avpol(ix,iy,iz,1),xh(ix,iy,iz)
            enddo
        enddo
    enddo

    ! Poor solvent

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                do ip = 1, N_poorsol
                    idx = ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells

                    f(idx) = xtotal(ix,iy,iz,ip)

                    do im = 1, N_monomer
                        if(hydroph(im).eq.ip) then 
                            f(idx) = f(idx) - avpol(ix,iy,iz,im)
                        endif
                    enddo ! im
                enddo ! ip

            enddo ! ix
        enddo ! iy
    enddo ! iz


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Poisson equation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !
    ! Some auxialiary variables, see Notes Poisson eq. non-cubic grid
    !

    MV(1) = MAT(1,1)
    MV(2) = MAT(1,2)  
    MV(3) = MAT(1,3)

    MU(1) = MAT(2,1)
    MU(2) = MAT(2,2)  
    MU(3) = MAT(2,3)

    MW(1) = MAT(3,1)
    MW(2) = MAT(3,2)  
    MW(3) = MAT(3,3)

    MVV = DOT_PRODUCT(MV,MV)
    MUU = DOT_PRODUCT(MU,MU)
    MWW = DOT_PRODUCT(MW,MW)

    MVU = DOT_PRODUCT(MV,MU)
    MVW = DOT_PRODUCT(MV,MW)
    MUW = DOT_PRODUCT(MU,MW)

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                psivv = psi(ix+1,iy,iz)-2*psi(ix,iy,iz)+psi(ix-1,iy,iz)
                psiuu = psi(ix,iy+1,iz)-2*psi(ix,iy,iz)+psi(ix,iy-1,iz)
                psiww = psi(ix,iy,iz+1)-2*psi(ix,iy,iz)+psi(ix,iy,iz-1)

                psivu = (psi(ix+1,iy+1,iz)+psi(ix-1,iy-1,iz)-psi(ix+1,iy-1,iz)-psi(ix-1,iy+1,iz))/4.0d0
                psivw = (psi(ix+1,iy,iz+1)+psi(ix-1,iy,iz-1)-psi(ix+1,iy,iz-1)-psi(ix-1,iy,iz+1))/4.0d0
                psiuw = (psi(ix,iy+1,iz+1)+psi(ix,iy-1,iz-1)-psi(ix,iy+1,iz-1)-psi(ix,iy-1,iz+1))/4.0d0

                psiv(1) = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))/2.0d0
                psiv(2) = (psi(ix,iy+1,iz)-psi(ix,iy-1,iz))/2.0d0
                psiv(3) = (psi(ix,iy,iz+1)-psi(ix,iy,iz-1))/2.0d0

                epsv(1) = (epsfcn(ix+1,iy,iz)-epsfcn(ix-1,iy,iz))/2.0d0
                epsv(2) = (epsfcn(ix,iy+1,iz)-epsfcn(ix,iy-1,iz))/2.0d0
                epsv(3) = (epsfcn(ix,iy,iz+1)-epsfcn(ix,iy,iz-1))/2.0d0

                psitemp = epsfcn(ix,iy,iz)*&
                    (MVV*psivv+MUU*psiuu+MWW*psiww+2.0d0*MVU*psivu+2.0d0*MVW*psivw+2.0d0*MUW*psiuw)
                psitemp = psitemp + DOT_PRODUCT(MATMUL(TMAT,epsv),MATMUL(TMAT,psiv))

                ! OJO CHECK!!!!

                f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsol+1)*ncells)=&
                    (psitemp + qtot(ix, iy, iz)*constq)/(-2.0d0)

            enddo
        enddo
    enddo

    ! flux equations 

    do i=1,niontypes 
        select case (iontype(i))
        case ("Hplus")
            call div_flux(divJ(:,:,:,i),xh,xHplus,psi,"Hplus")           
        case ("OHmin") 
            call div_flux(divJ(:,:,:,i),xh,xOHmin,psi,"OHmin") 
        case("pos")
            call div_flux(divJ(:,:,:,i),xh,xpos,psi,"pos")
        case("neg") 
            call div_flux(divJ(:,:,:,i),xh,xneg,psi,"neg")
        case default
            write(stdout,*)"fkfun: error reached unspecified iontype"
            stop 
        end select
    enddo    

    noffset=(N_poorsol+2)*ncells 
        

    do i=1,niontypes
        
        key = trim(iontype(i))
        xvolbulk = get_value_moleclist(xvolmin,key)

        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx
                    id=ix+dimx*(iy-1)+dimx*dimy*(iz-1)+noffset +(i-1)*ncells 
                    f(id)  =  divJ(ix,iy,iz,i)
                enddo
            enddo    
        enddo    
    enddo
    
    ! end flux
 
    ! == norma = 0.0 
    norma = 0.0d0      ! == added  d0  without d0 significant number loss can occur
    normvol = 0.0d0
    normel = 0.0d0
    normflux = 0.0d0
    normfluxpos = 0.0d0

    do i = 1, eqs*ncells
        norma = norma + f(i)**2
    enddo

    do i = 1,ncells
        normvol= normvol + f(i)**2
    enddo
   
    noffset=(N_poorsol+1)*ncells 
    do i= 1, ncells    
        normel = normel + f(i+noffset)**2
    enddo

    noffset=(N_poorsol+2)*ncells 
    do i= 1, 4*ncells    
        normflux = normflux + f(i+noffset)**2
    enddo 

    noffset=(N_poorsol+2)*ncells 
    do i= 1, ncells    
        normfluxpos = normfluxpos + f(i+noffset)**2
    enddo 

    print*,"hello"
    iter = iter + 1
    if(verbose.ge.3) then
        if(rank.eq.0) write(stdout,*)'fkfun:', iter, sqrt(norma), sqrt(normvol), sqrt(normel), sqrt(normflux), q(1)
        write(stdout,*)'rank=', rank, 'norm fluxpos=',normfluxpos
    endif

    2222 continue
    ier2 = 0.0 

    return

end subroutine fcn_flux




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
! ! == fcn for flux : interation in reduced ion densities via the Slotboom transformation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fcn_flux_ST(x,f,ier2)

    use system  ! use system, only : neqs
    use chainsdat
    use molecules
    use const
    use results
    use bulk
    use kai
    use MPI
    use fields_fkfun
    use kinsol
    use conformations
    use ematrix
    use ellipsoid
    use transform
    use kaist
    use mparameters_monomer
    use flux, only : divJ, div_flux, niontypes, iontype,  bc_psi 
    use flux, only : xvol_slotboom, div_flux_channel_slotboom

    implicit none

    ! input arguments 
    integer*4, intent(inout) :: ier2
    real*8, intent(in) :: x(*)
    real*8, intent(inout) :: f(*)

    ! local arguments

!    integer :: ncells
    real*8 :: protemp
    integer :: i, j, ix, iy, iz, ii, ax, ay, az
    integer :: im, ip
    integer :: jx, jy, jz, jj
    real*8 :: xpot(dimx, dimy, dimz, N_monomer)
    integer :: idx, id, noffset         ! == indices used for flux contribution to f
    real*8 :: normvol, normel, normflux, sumdivJ(niontypes), qres

    ! Charge
    real*8 :: psitemp
    real*8 :: MV(3),MU(3),MW(3)
    real*8 :: MVV,MUU,MWW,MVU,MVW,MUW
    real*8 :: psivv,psiuu,psiww, psivu,psivw,psiuw
    real*8 :: psiv(3), epsv(3)
    real*8 :: xtotalsum(dimx,dimy,dimz)

    integer, external :: PBCSYMI, PBCREFI

    ! poor solvent 
    real*8 :: sttemp
    ! MPI
    integer :: tag
    parameter(tag = 0)
    integer :: err
    real*8 :: avpol_temp(dimx,dimy,dimz,N_monomer)
    real*8 :: rho_tilde(0:dimx+1, 0:dimy+1, 0:dimz+1,niontypes) 
    real*8 :: q_tosend
    real*8 :: gradpsi2
    real*8 :: fv

    ! hamiltonian inception
    real*8 :: hfactor, hd
    real*8 :: hds(100)

    hds = -1

    !-----------------------------------------------------
    ! Common variables

    shift = 1.0d0

!    ncells = dimx*dimy*dimz ! == number of cells

    ! == HEAD node

    if(rank.eq.0) then   ! == calls subordinates and passes vector x
        flagsolver = 1
        CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
        !CALL MPI_BCAST(x, eqs*ncells , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
        CALL MPI_BCAST(x, neqs , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
    endif

    ! == Retrieve xh and psi from x()

    !  psi = 0.0 
    psi = 0.0d0 ! == without d0 significant number loss can ocur but value overwritten below !!!

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz
                idx = ix+dimx*(iy-1)+dimx*dimy*(iz-1)

                xh(ix,iy,iz)=x(idx)  !  == solvent 

                do ip = 1, N_poorsol
                    xtotal(ix,iy,iz,ip) = x(idx + ip*ncells) ! == fraction  polymer of type ip
                enddo
                
                psi(ix,iy,iz)=x(idx+(N_poorsol+1)*ncells)   ! == potential 
                
                do ip = 1, niontypes 
                    rho_tilde(ix,iy,iz,ip) = x(idx +(N_poorsol+1)*ncells + ip*ncells)  ! == check stride !!
                enddo      
                
            enddo
        enddo
    enddo

    ! construct ion volume fraction 

     do i=1,niontypes 
        select case (iontype(i))
        case ("Hplus")
            print*,"Hplus",iontype(i),i
            call xvol_slotboom(rho_tilde(:,:,:,i),xh,xHplus,psi,"Hplus")           
        case ("OHmin") 
            print*,"OHmin",iontype(i),i
            call xvol_slotboom(rho_tilde(:,:,:,i),xh,xOHmin,psi,"OHmin") 
        case("pos")
            print*,"pos=",iontype(i),i
            call xvol_slotboom(rho_tilde(:,:,:,i),xh,xpos,psi,"pos")
        case("neg") 
            print*,"neg=",iontype(i),i
            call xvol_slotboom(rho_tilde(:,:,:,i),xh,xneg,psi,"neg")
        case default
            write(stdout,*)"fkfun: error reached unspecified iontype"
            stop 
        end select
    enddo    
    
    if(fluxflag.eq.1) then 

        call bc_psi(psi)

    else 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        ! Boundary conditions electrostatic potential
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Reflection or PBC, (PBC = 1 or 3)
        
        do jx = 0, dimx+1
            do jy = 0, dimy+1
                do jz = 0, dimz+1

                    ix=jx
                    iy=jy
                    iz=jz ! these lines are necessary for PBC = 0 or 2

                    if (PBC(1).eq.1)ix = PBCSYMI(jx,dimx)
                    if (PBC(3).eq.1)iy = PBCSYMI(jy,dimy)
                    if (PBC(5).eq.1)iz = PBCSYMI(jz,dimz)

                    if (PBC(1).eq.3)ix = PBCREFI(jx,dimx)
                    if (PBC(3).eq.3)iy = PBCREFI(jy,dimy)
                    if (PBC(5).eq.3)iz = PBCREFI(jz,dimz)

                    psi(jx, jy, jz) = psi(ix, iy, iz)
                enddo
            enddo
        enddo

        ! Bulk or Wall, PBC = 0 or 2

        select case (PBC(1)) ! x = 0
        case(0) ! set bulk 
            psi(0,:,:) = 0.0d0        ! == added  d0 
        case(2)
            psi(0,:,:) = psi(1,:,:) ! zero charge
        end select

        select case (PBC(2)) ! x = dimx
        case(0) ! set bulk 
            psi(dimx+1,:,:) = 0.0d0  
        case(2)
            psi(dimx+1,:,:) = psi(dimx,:,:) ! zero charge
        end select

        select case (PBC(3)) ! y = 0
        case(0) ! set bulk 
            psi(:,0,:) = 0.0d0  
        case(2)
            psi(:,0,:) = psi(:,1,:) ! zero charge
        end select

        select case (PBC(4)) ! y = dimy
        case(0) ! set bulk 
            psi(:,dimy+1,:) = 0.0d0
        case(2)
            psi(:,dimy+1,:) = psi(:,dimy,:) ! zero charge
        end select

        select case (PBC(5)) ! z = 0
        case(0) ! set bulk 
            psi(:,:,0) = 0.0d0  
        case(2)
            psi(:,:,0) = psi(:,:,1) ! zero charge
        end select

        select case (PBC(6)) ! z = dimz
        case(0) ! set bulk 
            psi(:,:,dimz+1) = 0.0d0
        case(2)
            psi(:,:,dimz+1) = psi(:,:,dimz) ! zero charge
        end select

    endif

    ! volume fraction and fdis

    fdis = 0.0d0    ! == added  d0 
    avpol = 0.0d0

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                do im =1,N_monomer

                    if (zpol(im).eq.1) then !BASE
                        fdis(ix,iy,iz,im) = 1.0d0 /(1.0d0 + xOHmin(ix,iy,iz)/(K0(im)*xh(ix,iy,iz))) !k0 k en fraccion de volumen
                    else if (zpol(im).eq.-1) then !ACID
                        fdis(ix,iy,iz,im) = 1.0d0 /(1.0d0 + xHplus(ix,iy,iz)/(K0(im)*xh(ix,iy,iz)))
                    endif

                enddo

            enddo
        enddo  
    enddo


    ! Compute dielectric permitivity

    xtotalsum = 0.0d0 ! sum of all polymers
    do ip = 1, N_poorsol
        xtotalsum(:,:,:) = xtotalsum(:,:,:) + xtotal(:,:,:,ip)
    enddo
    
    call dielectfcn(xtotalsum,volprot,epsfcn,Depsfcn)

    !------------------------------------------------------------------------
    ! PDF polymer 
    !------------------------------------------------------------------------
    !
    ! parallel Each processor works on a chain
    !
    !-----------------------------------------------------------------------         
   

    ! Calculate xpot

    sttemp = st/(vpol*vsol)

    do im = 1, N_monomer ! loop over different monomer types

        do ix=1,dimx
            do iy=1,dimy
                do iz=1,dimz

                    if(hguess .eq. 0) then

                        hd = sqrt(float((2*ix-dimx)**2+(2*iy-dimy)**2))/2.0*delta
                        hd = hd**2+(oval*float(2*iz-dimz)/2.0*delta)**2
                        hfactor = dexp(-(kp**2)*hd)

                    elseif(hguess .eq. 1) then

                        hd = sqrt(float((2*ix-dimx)**2+(2*iy-dimy)**2))/2.0*delta-hring
                        hd = hd**2+(oval*float(2*iz-dimz)/2.0*delta)**2
                        hfactor = dexp(-(kp**2)*hd)

                    else

                        do i=1,hguess
                            hds(i) = (float(2*ix-dimx)-2*cos(i*2*pi/hguess)*hring/delta)**2+(float(2*iy-dimy)-&
                                2*sin(i*2*pi/hguess)*hring/delta)**2
                            hds(i) = hds(i)/4.0*(delta**2)+(oval*float(2*iz-dimz)/2.0*delta)**2
                        end do
                        hd = minval(hds, mask = hds .gt.0)
                        hfactor = dexp(-(kp**2)*hd)

                    end if

                    ! == xpot = exp(-Uj(rj)) para P(alpha)
                    ! == xpot is the exponent of the effective Energy function in P(alpha)

                    fv = (1.0 - volprot(ix,iy,iz)) 
    
                    ! == volume fraction of the cell that is sc volprot->wall fraction             

                    xpot(ix, iy, iz, im) = xh(ix,iy,iz)**vpol 
                    
                    ! == im: segment type, osmotic pressure term

                    xpot(ix, iy, iz, im) = xpot(ix,iy,iz, im)*dexp(voleps(ix,iy,iz))  
                    
                    ! == interaction term with the particle's sup =surface 

                    ! == Electrostatics

                    if(zpol(im).ne.0.0) then
                        xpot(ix,iy,iz,im) =  xpot(ix,iy,iz,im)/fdis(ix,iy,iz,im)*dexp(-psi(ix,iy,iz)*zpol(im))  
                        ! == fdis 
                    endif
        
                    ! == Dielectrics

                    gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy,iz))**2+&
                        (psi(ix,iy,iz+1)-psi(ix,iy,iz))**2 

                    ! ==  gradpsi2 = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy-1,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz-1))**2 
                    ! ==  xpot(ix, iy, iz) = xpot(ix,iy,iz)*exp(-Depsfcn(ix,iy,iz)*(gradpsi2)*constqE)

                    xpot(ix,iy,iz,im) = xpot(ix,iy,iz,im)*exp(Depsfcn(ix,iy,iz)*(gradpsi2)/constq/2.0d0*vpol/fv)

                    ! == Poor solvent depends on the grid where it is and its neighbors

                    if(hydroph(im).ne.0) then

                    protemp=0.0

                    do ax = -Xulimit,Xulimit 
                        do ay = -Xulimit,Xulimit
                            do az = -Xulimit,Xulimit

                                jx = ix+ax
                                jy = iy+ay
                                jz = iz+az

                                if(jx.lt.1) then
                                    if(PBC(1).eq.1)jx = PBCSYMI(jx,dimx)
                                    if(PBC(1).eq.3)jx = PBCREFI(jx,dimx)
                                endif

                                if(jx.gt.dimx) then
                                    if(PBC(2).eq.1)jx = PBCSYMI(jx,dimx)
                                    if(PBC(2).eq.3)jx = PBCREFI(jx,dimx)
                                endif

                                if(jy.lt.1) then
                                    if(PBC(3).eq.1)jy = PBCSYMI(jy,dimy)
                                    if(PBC(3).eq.3)jy = PBCREFI(jy,dimy)
                                endif

                                if(jy.gt.dimy) then
                                    if(PBC(4).eq.1)jy = PBCSYMI(jy,dimy)
                                    if(PBC(4).eq.3)jy = PBCREFI(jy,dimy)
                                endif


                                if(jz.lt.1) then
                                    if(PBC(5).eq.1)jz = PBCSYMI(jz,dimz)
                                    if(PBC(5).eq.3)jz = PBCREFI(jz,dimz)
                                endif

                                if(jz.gt.dimz) then
                                    if(PBC(6).eq.1)jz = PBCSYMI(jz,dimz)
                                    if(PBC(6).eq.3)jz = PBCREFI(jz,dimz)
                                endif


                                if((jx.ge.1).and.(jx.le.dimx)) then
                                    if((jy.ge.1).and.(jy.le.dimy)) then
                                        if((jz.ge.1).and.(jz.le.dimz)) then
                                            fv = (1.0-volprot(jx,jy,jz))

                                            do ip = 1, N_poorsol
                                                protemp = protemp + hfactor*Xu(ax,ay,az)*&
                                                st_matrix(hydroph(im),ip)*sttemp*xtotal(jx,jy,jz,ip)*fv
                                            enddo ! ip

                                        endif
                                    endif
                                endif

                            enddo
                        enddo
                    enddo

                    xpot(ix,iy,iz,im) = xpot(ix,iy,iz,im)*dexp(protemp)

                    endif ! hydrph

!                    write(567,*)xpot(ix,iy,iz,1)

                enddo ! ix
            enddo ! iy
        enddo !iz

    enddo ! N_monomer


    !!!!!!!!!!!!!!!!!!!!!! Calculate pro from xpot !!!!!!!!!!!!!
    call calcavpol(xpot)
    
    !!!!!!!!!!! MPORTANT, THE SUBORDINATES' TASKS END HERE ... 
    if(rank.ne.0)goto 4444
    !!!!!!!!!!!!!!!!!!!!!!! END  MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !----------------------------------------------------------------------------------------------
    !  == Build Equations to solve
    !----------------------------------------------------------------------------------------------

    ! Qtot


    qtot = 0.0d0 ! == added  d0 
    qres =0.0d0

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz
            
                fv = (1.0-volprot(ix,iy,iz))     ! fv == free volume 

                qtot(ix, iy, iz) =  (zpos*xpos(ix, iy, iz)+zneg*xneg(ix, iy, iz))/vsalt + &
                    xHplus(ix, iy, iz) - xOHmin(ix, iy, iz)

                do im = 1, N_monomer
                    qtot(ix, iy, iz) =  qtot(ix,iy,iz) + &
                        avpol(ix,iy,iz,im)*zpol(im)/vpol*fdis(ix,iy,iz,im)
                enddo

                qtot(ix, iy,iz) = qtot(ix,iy,iz)*fv + volq(ix,iy,iz)*vsol  
               
                qres = qres + qtot(ix, iy,iz) 

            enddo
        enddo
    enddo

    ! Volume fraction
    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                idx = ix+dimx*(iy-1)+dimx*dimy*(iz-1) ! index 

                f(idx) = xh(ix,iy,iz) + xneg(ix, iy, iz) + xpos(ix, iy, iz) + xHplus(ix, iy, iz) + &
                    xOHmin(ix, iy, iz) -1.0d0   ! packing iones+sv

                do im = 1, N_monomer
                    f(idx) = f(idx) + avpol(ix,iy,iz,im)          ! packing ...+polymer
                enddo

                ! write(123,*)ix,iy,iz,avpol(ix,iy,iz,1),xh(ix,iy,iz)
            enddo
        enddo
    enddo

    ! Poor solvent

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz
                !idx= imap(ix,iy,iz)
                do ip = 1, N_poorsol

                    idx = ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells ! index
                
                    f(idx) = xtotal(ix,iy,iz,ip)

                    do im = 1, N_monomer
                        if(hydroph(im).eq.ip) then 
                            f(idx) = f(idx) - avpol(ix,iy,iz,im)
                        endif
                    enddo ! im
                enddo ! ip

            enddo ! ix
        enddo ! iy
    enddo ! iz


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Poisson equation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !
    ! Some auxialiary variables, see Notes Poisson eq. non-cubic grid
    !

    MV(1) = MAT(1,1)
    MV(2) = MAT(1,2)  
    MV(3) = MAT(1,3)

    MU(1) = MAT(2,1)
    MU(2) = MAT(2,2)  
    MU(3) = MAT(2,3)

    MW(1) = MAT(3,1)
    MW(2) = MAT(3,2)  
    MW(3) = MAT(3,3)

    MVV = DOT_PRODUCT(MV,MV)
    MUU = DOT_PRODUCT(MU,MU)
    MWW = DOT_PRODUCT(MW,MW)

    MVU = DOT_PRODUCT(MV,MU)
    MVW = DOT_PRODUCT(MV,MW)
    MUW = DOT_PRODUCT(MU,MW)

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                psivv = psi(ix+1,iy,iz)-2*psi(ix,iy,iz)+psi(ix-1,iy,iz)
                psiuu = psi(ix,iy+1,iz)-2*psi(ix,iy,iz)+psi(ix,iy-1,iz)
                psiww = psi(ix,iy,iz+1)-2*psi(ix,iy,iz)+psi(ix,iy,iz-1)

                psivu = (psi(ix+1,iy+1,iz)+psi(ix-1,iy-1,iz)-psi(ix+1,iy-1,iz)-psi(ix-1,iy+1,iz))/4.0d0
                psivw = (psi(ix+1,iy,iz+1)+psi(ix-1,iy,iz-1)-psi(ix+1,iy,iz-1)-psi(ix-1,iy,iz+1))/4.0d0
                psiuw = (psi(ix,iy+1,iz+1)+psi(ix,iy-1,iz-1)-psi(ix,iy+1,iz-1)-psi(ix,iy-1,iz+1))/4.0d0

                psiv(1) = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))/2.0d0
                psiv(2) = (psi(ix,iy+1,iz)-psi(ix,iy-1,iz))/2.0d0
                psiv(3) = (psi(ix,iy,iz+1)-psi(ix,iy,iz-1))/2.0d0

                epsv(1) = (epsfcn(ix+1,iy,iz)-epsfcn(ix-1,iy,iz))/2.0d0
                epsv(2) = (epsfcn(ix,iy+1,iz)-epsfcn(ix,iy-1,iz))/2.0d0
                epsv(3) = (epsfcn(ix,iy,iz+1)-epsfcn(ix,iy,iz-1))/2.0d0

                psitemp = epsfcn(ix,iy,iz)*&
                    (MVV*psivv+MUU*psiuu+MWW*psiww+2.0d0*MVU*psivu+2.0d0*MVW*psivw+2.0d0*MUW*psiuw)
                psitemp = psitemp + DOT_PRODUCT(MATMUL(TMAT,epsv),MATMUL(TMAT,psiv))

                ! OJO CHECK!!!!

                f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsol+1)*ncells)=&
                    (psitemp + qtot(ix, iy, iz)*constq)/(-2.0d0)

            enddo
        enddo
    enddo

    ! flux equations 

    do i=1,niontypes 
        select case (iontype(i))
        case ("Hplus")
            call div_flux_channel_slotboom(divJ(:,:,:,i),rho_tilde(:,:,:,i), xh, psi,"Hplus")           
        case ( "OHmin") 
            call div_flux_channel_slotboom(divJ(:,:,:,i),rho_tilde(:,:,:,i), xh, psi,"OHmin") 
        case("pos")
            call div_flux_channel_slotboom(divJ(:,:,:,i),rho_tilde(:,:,:,i), xh, psi,"pos")
        case("neg") 
            call div_flux_channel_slotboom(divJ(:,:,:,i),rho_tilde(:,:,:,i), xh, psi,"neg")
        case default
            write(stdout,*)"fkfun: error reached unspecified iontype"
            stop 
        end select
    enddo    

    noffset=(N_poorsol+2)*ncells 
        
    sumdivJ =0.0d0

    do i=1,niontypes
        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx
                    ! id = imap(ix,iy,iz) +(i-1)*ncells 
                    id = ix+dimx*(iy-1)+dimx*dimy*(iz-1)+noffset +(i-1)*ncells 

                    f(id)  =  divJ(ix,iy,iz,i)

                    sumdivJ(i) = sumdivJ(i)+divJ(ix,iy,iz,i)**2

                    ! test out 
                    if(iter==0) then 
                        if(rank.eq.0) then 
                            select case (iontype(i))
                            case ("Hplus")
                                write(100+i,*)ix,iy,iz, divJ(ix,iy,iz,i),rho_tilde(ix,iy,iz,i),xHplus(ix, iy, iz)
                            case ( "OHmin") 
                                write(100+i,*)ix,iy,iz, divJ(ix,iy,iz,i),rho_tilde(ix,iy,iz,i),xOHmin(ix, iy, iz)
                            case("pos")
                                write(100+i,*)ix,iy,iz, divJ(ix,iy,iz,i),rho_tilde(ix,iy,iz,i),xpos(ix, iy, iz)
                            case("neg") 
                                write(100+i,*)ix,iy,iz, divJ(ix,iy,iz,i),rho_tilde(ix,iy,iz,i),xneg(ix, iy, iz)
                            end select
                        endif
                            
                    endif
                enddo
            enddo    
        enddo    
    enddo
    
    ! end flux
 
    ! == norma = 0.0 
    norma = 0.0d0      ! == added  d0  without d0 significant number loss can occur
    normvol = 0.0d0
    normel = 0.0d0
    normflux = 0.0d0

    do i = 1, neqs
        norma = norma + f(i)**2
    enddo

    do i = 1, ncells
        normvol= normvol + f(i)**2
    enddo
   
    noffset = (N_poorsol+1)*ncells 
    do i= 1, ncells    
        normel = normel + f(i+noffset)**2
    enddo

    noffset = (N_poorsol+2)*ncells 
    do i= 1, niontypes * ncells    
        normflux = normflux + f(i+noffset)**2
    enddo 
    
    iter = iter + 1
    if(verbose.ge.3) then
        if(rank.eq.0) write(stdout,*)'fkfun:', iter, sqrt(norma), sqrt(normvol), sqrt(normel), sqrt(normflux), q(1)
        if(rank.eq.0) write(stdout,*)'fkfun:', "sumdivJ=",sumdivJ
        if(rank.eq.0) write(stdout,*)'fkfun:', "sqrt(sumdivJ)=",sqrt(sum(sumdivJ))
        if(rank.eq.0) write(stdout,*)'fkfun:', "niontypes =",niontypes," ", iontype
        if(rank.eq.0) write(stdout,*)'fkfun:', "qres = ", qres
    endif

    4444 continue
    ier2 = 0.0 

    return

end subroutine fcn_flux_ST



! =====================================================================================   
! == fcn for flux : interation in reduced ion densities via the Slotboom transformation
! == sed extra boundary condition to fix electostatic potential of one reservoir and 
! == determined self-consistently the electostatic potential of tte other reservoir 
! =====================================================================================  

subroutine fcn_flux_ST_bc(x,f,ier2)

    use system  ! use system, only : neqs
    use chainsdat
    use molecules
    use const
    use results
    use bulk
    use kai
    use MPI
    use fields_fkfun
    use kinsol
    use conformations
    use ematrix
    use ellipsoid
    use transform
    use kaist
    use mparameters_monomer
    use flux, only : divJ, div_flux, niontypes, iontype,  bc_psi 
    use flux, only : xvol_slotboom, div_flux_channel_slotboom

    implicit none

    ! input arguments 
    integer*4, intent(inout) :: ier2
    real*8, intent(in) :: x(*)
    real*8, intent(inout) :: f(*)

    ! local arguments

!    integer :: ncells
    real*8 :: protemp
    integer :: i, j, ix, iy, iz, ii, ax, ay, az
    integer :: im, ip
    integer :: jx, jy, jz, jj
    real*8 :: xpot(dimx, dimy, dimz, N_monomer)
    integer :: idx, id, noffset         ! == indices used for flux contribution to f
    real*8 :: normvol, normel, normflux, sumdivJ(niontypes), qres

    ! Charge
    real*8 :: psitemp
    real*8 :: MV(3),MU(3),MW(3)
    real*8 :: MVV,MUU,MWW,MVU,MVW,MUW
    real*8 :: psivv,psiuu,psiww, psivu,psivw,psiuw
    real*8 :: psiv(3), epsv(3)
    real*8 :: xtotalsum(dimx,dimy,dimz)

    integer, external :: PBCSYMI, PBCREFI

    ! poor solvent 
    real*8 :: sttemp
    ! MPI
    integer :: tag
    parameter(tag = 0)
    integer :: err
    real*8 :: avpol_temp(dimx,dimy,dimz,N_monomer)
    real*8 :: rho_tilde(0:dimx+1, 0:dimy+1, 0:dimz+1,niontypes) 
    real*8 :: q_tosend
    real*8 :: gradpsi2
    real*8 :: fv

    ! hamiltonian inception
    real*8 :: hfactor, hd
    real*8 :: hds(100)

    hds = -1

    !-----------------------------------------------------
    ! Common variables

    shift = 1.0d0

!    ncells = dimx*dimy*dimz ! == number of cells

    ! == HEAD node

    if(rank.eq.0) then   ! == calls subordinates and passes vector x
        flagsolver = 1
        CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
        !CALL MPI_BCAST(x, eqs*ncells , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
        CALL MPI_BCAST(x, neqs , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
    endif

    ! == Retrieve xh and psi from x()

    !  psi = 0.0 
    psi = 0.0d0 ! == without d0 significant number loss can ocur but value overwritten below !!!

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz
                idx = ix+dimx*(iy-1)+dimx*dimy*(iz-1)

                xh(ix,iy,iz)=x(idx)  !  == solvent 

                do ip = 1, N_poorsol
                    xtotal(ix,iy,iz,ip) = x(idx + ip*ncells) ! == fraction  polymer of type ip
                enddo
                
                psi(ix,iy,iz)=x(idx+(N_poorsol+1)*ncells)   ! == potential 
                
                do ip = 1, niontypes 
                    rho_tilde(ix,iy,iz,ip) = x(idx +(N_poorsol+1)*ncells + ip*ncells)  ! == check stride !!
                enddo      
                
            enddo
        enddo
    enddo

    ! construct ion volume fraction 

     do i=1,niontypes 
        select case (iontype(i))
        case ("Hplus")
            print*,"Hplus",iontype(i),i
            call xvol_slotboom(rho_tilde(:,:,:,i),xh,xHplus,psi,"Hplus")           
        case ("OHmin") 
            print*,"OHmin",iontype(i),i
            call xvol_slotboom(rho_tilde(:,:,:,i),xh,xOHmin,psi,"OHmin") 
        case("pos")
            print*,"pos=",iontype(i),i
            call xvol_slotboom(rho_tilde(:,:,:,i),xh,xpos,psi,"pos")
        case("neg") 
            print*,"neg=",iontype(i),i
            call xvol_slotboom(rho_tilde(:,:,:,i),xh,xneg,psi,"neg")
        case default
            write(stdout,*)"fkfun: error reached unspecified iontype"
            stop 
        end select
    enddo    
    
    if(fluxflag.eq.1) then 

        call bc_psi(psi)

    else 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        ! Boundary conditions electrostatic potential
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Reflection or PBC, (PBC = 1 or 3)
        
        do jx = 0, dimx+1
            do jy = 0, dimy+1
                do jz = 0, dimz+1

                    ix=jx
                    iy=jy
                    iz=jz ! these lines are necessary for PBC = 0 or 2

                    if (PBC(1).eq.1)ix = PBCSYMI(jx,dimx)
                    if (PBC(3).eq.1)iy = PBCSYMI(jy,dimy)
                    if (PBC(5).eq.1)iz = PBCSYMI(jz,dimz)

                    if (PBC(1).eq.3)ix = PBCREFI(jx,dimx)
                    if (PBC(3).eq.3)iy = PBCREFI(jy,dimy)
                    if (PBC(5).eq.3)iz = PBCREFI(jz,dimz)

                    psi(jx, jy, jz) = psi(ix, iy, iz)
                enddo
            enddo
        enddo

        ! Bulk or Wall, PBC = 0 or 2

        select case (PBC(1)) ! x = 0
        case(0) ! set bulk 
            psi(0,:,:) = 0.0d0        ! == added  d0 
        case(2)
            psi(0,:,:) = psi(1,:,:) ! zero charge
        end select

        select case (PBC(2)) ! x = dimx
        case(0) ! set bulk 
            psi(dimx+1,:,:) = 0.0d0  
        case(2)
            psi(dimx+1,:,:) = psi(dimx,:,:) ! zero charge
        end select

        select case (PBC(3)) ! y = 0
        case(0) ! set bulk 
            psi(:,0,:) = 0.0d0  
        case(2)
            psi(:,0,:) = psi(:,1,:) ! zero charge
        end select

        select case (PBC(4)) ! y = dimy
        case(0) ! set bulk 
            psi(:,dimy+1,:) = 0.0d0
        case(2)
            psi(:,dimy+1,:) = psi(:,dimy,:) ! zero charge
        end select

        select case (PBC(5)) ! z = 0
        case(0) ! set bulk 
            psi(:,:,0) = 0.0d0  
        case(2)
            psi(:,:,0) = psi(:,:,1) ! zero charge
        end select

        select case (PBC(6)) ! z = dimz
        case(0) ! set bulk 
            psi(:,:,dimz+1) = 0.0d0
        case(2)
            psi(:,:,dimz+1) = psi(:,:,dimz) ! zero charge
        end select

    endif

    ! volume fraction and fdis

    fdis = 0.0d0    ! == added  d0 
    avpol = 0.0d0

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                do im =1,N_monomer

                    if (zpol(im).eq.1) then !BASE
                        fdis(ix,iy,iz,im) = 1.0d0 /(1.0d0 + xOHmin(ix,iy,iz)/(K0(im)*xh(ix,iy,iz))) !k0 k en fraccion de volumen
                    else if (zpol(im).eq.-1) then !ACID
                        fdis(ix,iy,iz,im) = 1.0d0 /(1.0d0 + xHplus(ix,iy,iz)/(K0(im)*xh(ix,iy,iz)))
                    endif

                enddo

            enddo
        enddo  
    enddo


    ! Compute dielectric permitivity

    xtotalsum = 0.0d0 ! sum of all polymers
    do ip = 1, N_poorsol
        xtotalsum(:,:,:) = xtotalsum(:,:,:) + xtotal(:,:,:,ip)
    enddo
    
    call dielectfcn(xtotalsum,volprot,epsfcn,Depsfcn)

    !------------------------------------------------------------------------
    ! PDF polymer 
    !------------------------------------------------------------------------
    !
    ! parallel Each processor works on a chain
    !
    !-----------------------------------------------------------------------         
   

    ! Calculate xpot

    sttemp = st/(vpol*vsol)

    do im = 1, N_monomer ! loop over different monomer types

        do ix=1,dimx
            do iy=1,dimy
                do iz=1,dimz

                    if(hguess .eq. 0) then

                        hd = sqrt(float((2*ix-dimx)**2+(2*iy-dimy)**2))/2.0*delta
                        hd = hd**2+(oval*float(2*iz-dimz)/2.0*delta)**2
                        hfactor = dexp(-(kp**2)*hd)

                    elseif(hguess .eq. 1) then

                        hd = sqrt(float((2*ix-dimx)**2+(2*iy-dimy)**2))/2.0*delta-hring
                        hd = hd**2+(oval*float(2*iz-dimz)/2.0*delta)**2
                        hfactor = dexp(-(kp**2)*hd)

                    else

                        do i=1,hguess
                            hds(i) = (float(2*ix-dimx)-2*cos(i*2*pi/hguess)*hring/delta)**2+(float(2*iy-dimy)-&
                                2*sin(i*2*pi/hguess)*hring/delta)**2
                            hds(i) = hds(i)/4.0*(delta**2)+(oval*float(2*iz-dimz)/2.0*delta)**2
                        end do
                        hd = minval(hds, mask = hds .gt.0)
                        hfactor = dexp(-(kp**2)*hd)

                    end if

                    ! == xpot = exp(-Uj(rj)) para P(alpha)
                    ! == xpot is the exponent of the effective Energy function in P(alpha)

                    fv = (1.0 - volprot(ix,iy,iz)) 
    
                    ! == volume fraction of the cell that is sc volprot->wall fraction             

                    xpot(ix, iy, iz, im) = xh(ix,iy,iz)**vpol 
                    
                    ! == im: segment type, osmotic pressure term

                    xpot(ix, iy, iz, im) = xpot(ix,iy,iz, im)*dexp(voleps(ix,iy,iz))  
                    
                    ! == interaction term with the particle's sup =surface 

                    ! == Electrostatics

                    if(zpol(im).ne.0.0) then
                        xpot(ix,iy,iz,im) =  xpot(ix,iy,iz,im)/fdis(ix,iy,iz,im)*dexp(-psi(ix,iy,iz)*zpol(im))  
                        ! == fdis 
                    endif
        
                    ! == Dielectrics

                    gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy,iz))**2+&
                        (psi(ix,iy,iz+1)-psi(ix,iy,iz))**2 

                    ! ==  gradpsi2 = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy-1,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz-1))**2 
                    ! ==  xpot(ix, iy, iz) = xpot(ix,iy,iz)*exp(-Depsfcn(ix,iy,iz)*(gradpsi2)*constqE)

                    xpot(ix,iy,iz,im) = xpot(ix,iy,iz,im)*exp(Depsfcn(ix,iy,iz)*(gradpsi2)/constq/2.0d0*vpol/fv)

                    ! == Poor solvent depends on the grid where it is and its neighbors

                    if(hydroph(im).ne.0) then

                    protemp=0.0

                    do ax = -Xulimit,Xulimit 
                        do ay = -Xulimit,Xulimit
                            do az = -Xulimit,Xulimit

                                jx = ix+ax
                                jy = iy+ay
                                jz = iz+az

                                if(jx.lt.1) then
                                    if(PBC(1).eq.1)jx = PBCSYMI(jx,dimx)
                                    if(PBC(1).eq.3)jx = PBCREFI(jx,dimx)
                                endif

                                if(jx.gt.dimx) then
                                    if(PBC(2).eq.1)jx = PBCSYMI(jx,dimx)
                                    if(PBC(2).eq.3)jx = PBCREFI(jx,dimx)
                                endif

                                if(jy.lt.1) then
                                    if(PBC(3).eq.1)jy = PBCSYMI(jy,dimy)
                                    if(PBC(3).eq.3)jy = PBCREFI(jy,dimy)
                                endif

                                if(jy.gt.dimy) then
                                    if(PBC(4).eq.1)jy = PBCSYMI(jy,dimy)
                                    if(PBC(4).eq.3)jy = PBCREFI(jy,dimy)
                                endif


                                if(jz.lt.1) then
                                    if(PBC(5).eq.1)jz = PBCSYMI(jz,dimz)
                                    if(PBC(5).eq.3)jz = PBCREFI(jz,dimz)
                                endif

                                if(jz.gt.dimz) then
                                    if(PBC(6).eq.1)jz = PBCSYMI(jz,dimz)
                                    if(PBC(6).eq.3)jz = PBCREFI(jz,dimz)
                                endif


                                if((jx.ge.1).and.(jx.le.dimx)) then
                                    if((jy.ge.1).and.(jy.le.dimy)) then
                                        if((jz.ge.1).and.(jz.le.dimz)) then
                                            fv = (1.0-volprot(jx,jy,jz))

                                            do ip = 1, N_poorsol
                                                protemp = protemp + hfactor*Xu(ax,ay,az)*&
                                                st_matrix(hydroph(im),ip)*sttemp*xtotal(jx,jy,jz,ip)*fv
                                            enddo ! ip

                                        endif
                                    endif
                                endif

                            enddo
                        enddo
                    enddo

                    xpot(ix,iy,iz,im) = xpot(ix,iy,iz,im)*dexp(protemp)

                    endif ! hydrph

!                    write(567,*)xpot(ix,iy,iz,1)

                enddo ! ix
            enddo ! iy
        enddo !iz

    enddo ! N_monomer


    !!!!!!!!!!!!!!!!!!!!!! Calculate pro from xpot !!!!!!!!!!!!!
    call calcavpol(xpot)
    
    !!!!!!!!!!! MPORTANT, THE SUBORDINATES' TASKS END HERE ... 
    if(rank.ne.0)goto 4444
    !!!!!!!!!!!!!!!!!!!!!!! END  MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !----------------------------------------------------------------------------------------------
    !  == Build Equations to solve
    !----------------------------------------------------------------------------------------------

    ! Qtot


    qtot = 0.0d0 ! == added  d0 
    qres =0.0d0

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz
            
                fv = (1.0-volprot(ix,iy,iz))     ! fv == free volume 

                qtot(ix, iy, iz) =  (zpos*xpos(ix, iy, iz)+zneg*xneg(ix, iy, iz))/vsalt + &
                    xHplus(ix, iy, iz) - xOHmin(ix, iy, iz)

                do im = 1, N_monomer
                    qtot(ix, iy, iz) =  qtot(ix,iy,iz) + &
                        avpol(ix,iy,iz,im)*zpol(im)/vpol*fdis(ix,iy,iz,im)
                enddo

                qtot(ix, iy,iz) = qtot(ix,iy,iz)*fv + volq(ix,iy,iz)*vsol  
               
                qres = qres + qtot(ix, iy,iz) 

            enddo
        enddo
    enddo

    ! Volume fraction
    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                idx = ix+dimx*(iy-1)+dimx*dimy*(iz-1) ! index 

                f(idx) = xh(ix,iy,iz) + xneg(ix, iy, iz) + xpos(ix, iy, iz) + xHplus(ix, iy, iz) + &
                    xOHmin(ix, iy, iz) -1.0d0   ! packing iones+sv

                do im = 1, N_monomer
                    f(idx) = f(idx) + avpol(ix,iy,iz,im)          ! packing ...+polymer
                enddo

                ! write(123,*)ix,iy,iz,avpol(ix,iy,iz,1),xh(ix,iy,iz)
            enddo
        enddo
    enddo

    ! Poor solvent

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz
                !idx= imap(ix,iy,iz)
                do ip = 1, N_poorsol

                    idx = ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells ! index
                
                    f(idx) = xtotal(ix,iy,iz,ip)

                    do im = 1, N_monomer
                        if(hydroph(im).eq.ip) then 
                            f(idx) = f(idx) - avpol(ix,iy,iz,im)
                        endif
                    enddo ! im
                enddo ! ip

            enddo ! ix
        enddo ! iy
    enddo ! iz


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Poisson equation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !
    ! Some auxialiary variables, see Notes Poisson eq. non-cubic grid
    !

    MV(1) = MAT(1,1)
    MV(2) = MAT(1,2)  
    MV(3) = MAT(1,3)

    MU(1) = MAT(2,1)
    MU(2) = MAT(2,2)  
    MU(3) = MAT(2,3)

    MW(1) = MAT(3,1)
    MW(2) = MAT(3,2)  
    MW(3) = MAT(3,3)

    MVV = DOT_PRODUCT(MV,MV)
    MUU = DOT_PRODUCT(MU,MU)
    MWW = DOT_PRODUCT(MW,MW)

    MVU = DOT_PRODUCT(MV,MU)
    MVW = DOT_PRODUCT(MV,MW)
    MUW = DOT_PRODUCT(MU,MW)

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                psivv = psi(ix+1,iy,iz)-2*psi(ix,iy,iz)+psi(ix-1,iy,iz)
                psiuu = psi(ix,iy+1,iz)-2*psi(ix,iy,iz)+psi(ix,iy-1,iz)
                psiww = psi(ix,iy,iz+1)-2*psi(ix,iy,iz)+psi(ix,iy,iz-1)

                psivu = (psi(ix+1,iy+1,iz)+psi(ix-1,iy-1,iz)-psi(ix+1,iy-1,iz)-psi(ix-1,iy+1,iz))/4.0d0
                psivw = (psi(ix+1,iy,iz+1)+psi(ix-1,iy,iz-1)-psi(ix+1,iy,iz-1)-psi(ix-1,iy,iz+1))/4.0d0
                psiuw = (psi(ix,iy+1,iz+1)+psi(ix,iy-1,iz-1)-psi(ix,iy+1,iz-1)-psi(ix,iy-1,iz+1))/4.0d0

                psiv(1) = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))/2.0d0
                psiv(2) = (psi(ix,iy+1,iz)-psi(ix,iy-1,iz))/2.0d0
                psiv(3) = (psi(ix,iy,iz+1)-psi(ix,iy,iz-1))/2.0d0

                epsv(1) = (epsfcn(ix+1,iy,iz)-epsfcn(ix-1,iy,iz))/2.0d0
                epsv(2) = (epsfcn(ix,iy+1,iz)-epsfcn(ix,iy-1,iz))/2.0d0
                epsv(3) = (epsfcn(ix,iy,iz+1)-epsfcn(ix,iy,iz-1))/2.0d0

                psitemp = epsfcn(ix,iy,iz)*&
                    (MVV*psivv+MUU*psiuu+MWW*psiww+2.0d0*MVU*psivu+2.0d0*MVW*psivw+2.0d0*MUW*psiuw)
                psitemp = psitemp + DOT_PRODUCT(MATMUL(TMAT,epsv),MATMUL(TMAT,psiv))

                ! OJO CHECK!!!!

                f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsol+1)*ncells)=&
                    (psitemp + qtot(ix, iy, iz)*constq)/(-2.0d0)

            enddo
        enddo
    enddo

    ! flux equations 

    do i=1,niontypes 
        select case (iontype(i))
        case ("Hplus")
            call div_flux_channel_slotboom(divJ(:,:,:,i),rho_tilde(:,:,:,i), xh, psi,"Hplus")           
        case ( "OHmin") 
            call div_flux_channel_slotboom(divJ(:,:,:,i),rho_tilde(:,:,:,i), xh, psi,"OHmin") 
        case("pos")
            call div_flux_channel_slotboom(divJ(:,:,:,i),rho_tilde(:,:,:,i), xh, psi,"pos")
        case("neg") 
            call div_flux_channel_slotboom(divJ(:,:,:,i),rho_tilde(:,:,:,i), xh, psi,"neg")
        case default
            write(stdout,*)"fkfun: error reached unspecified iontype"
            stop 
        end select
    enddo    

    noffset=(N_poorsol+2)*ncells 
        
    sumdivJ =0.0d0

    do i=1,niontypes
        do iz=1,dimz
            do iy=1,dimy
                do ix=1,dimx
                    ! id = imap(ix,iy,iz) +(i-1)*ncells 
                    id = ix+dimx*(iy-1)+dimx*dimy*(iz-1)+noffset +(i-1)*ncells 

                    f(id)  =  divJ(ix,iy,iz,i)

                    sumdivJ(i) = sumdivJ(i)+divJ(ix,iy,iz,i)**2

                    ! test out 
                    if(iter==0) then 
                        if(rank.eq.0) then 
                            select case (iontype(i))
                            case ("Hplus")
                                write(100+i,*)ix,iy,iz, divJ(ix,iy,iz,i),rho_tilde(ix,iy,iz,i),xHplus(ix, iy, iz)
                            case ( "OHmin") 
                                write(100+i,*)ix,iy,iz, divJ(ix,iy,iz,i),rho_tilde(ix,iy,iz,i),xOHmin(ix, iy, iz)
                            case("pos")
                                write(100+i,*)ix,iy,iz, divJ(ix,iy,iz,i),rho_tilde(ix,iy,iz,i),xpos(ix, iy, iz)
                            case("neg") 
                                write(100+i,*)ix,iy,iz, divJ(ix,iy,iz,i),rho_tilde(ix,iy,iz,i),xneg(ix, iy, iz)
                            end select
                        endif
                            
                    endif
                enddo
            enddo    
        enddo    
    enddo
    
    ! end flux
 
    ! == norma = 0.0 
    norma = 0.0d0      ! == added  d0  without d0 significant number loss can occur
    normvol = 0.0d0
    normel = 0.0d0
    normflux = 0.0d0

    do i = 1, neqs
        norma = norma + f(i)**2
    enddo

    do i = 1, ncells
        normvol= normvol + f(i)**2
    enddo
   
    noffset = (N_poorsol+1)*ncells 
    do i= 1, ncells    
        normel = normel + f(i+noffset)**2
    enddo

    noffset = (N_poorsol+2)*ncells 
    do i= 1, niontypes * ncells    
        normflux = normflux + f(i+noffset)**2
    enddo 
    
    iter = iter + 1
    if(verbose.ge.3) then
        if(rank.eq.0) write(stdout,*)'fkfun:', iter, sqrt(norma), sqrt(normvol), sqrt(normel), sqrt(normflux), q(1)
        if(rank.eq.0) write(stdout,*)'fkfun:', "sumdivJ=",sumdivJ
        if(rank.eq.0) write(stdout,*)'fkfun:', "sqrt(sumdivJ)=",sqrt(sum(sumdivJ))
        if(rank.eq.0) write(stdout,*)'fkfun:', "niontypes =",niontypes," ", iontype
        if(rank.eq.0) write(stdout,*)'fkfun:', "qres = ", qres
    endif

    4444 continue
    ier2 = 0.0 

    return

end subroutine fcn_flux_ST_bc


subroutine calc_std(xpot)

    use MPI
    use fields_fkfun
    use chainsdat
    use conformations
    use molecules
    use ematrix
    use kaist
    use mparameters_monomer
    use results

    implicit none

    real*8, intent(in) :: xpot(dimx, dimy, dimz, N_monomer)

    ! local variables 

    real*8 :: avpol_tosend(dimx,dimy,dimz, N_monomer)
    real*8 :: fv
    real*8 :: q_tosend
    real*8 :: avpol_temp(dimx,dimy,dimz,N_monomer)
    integer :: im,jj,i,j, ix, iy, iz, ii, ax, ay, az
    ! MPI
    integer :: tag
    parameter(tag = 0)
    integer :: err

    shift = 1.0d0 ! == added  d0  ! uniform shift in P(alpha)
    avpol_tosend = 0.0d0
    q = 0.0d0

    do jj = 1, cpp(rank+1)          ! == loop graft  on processor??
        ii = cppini(rank+1)+jj      ! == graft point ii ??

        q_tosend=0.0d0
        avpol_temp = 0.0d0

        do i=1,newcuantas(ii)       ! loop of chains
       
            pro(i, jj)= shift
            
            do j=1,long
                ax = px(i, j, jj) ! cada uno para su cadena... == each to his own chain...
                ay = py(i, j, jj)
                az = pz(i, j, jj)
                pro(i, jj) = pro(i, jj) * xpot(ax, ay, az, segtype(j))
            enddo
            
            pro(i, jj) = pro(i, jj) * dexp(-benergy*ngauche(i,ii)) ! energy of gauche bonds
            pro(i, jj) = pro(i, jj) * dexp(-fz*zfinal(i,jj))  ! termino Fz

            do j=1,long
                fv = fvstd(px(i,j, jj),py(i,j, jj),pz(i,j, jj))
                im = segtype(j)
                avpol_temp(px(i,j, jj),py(i,j, jj),pz(i,j, jj),im)= &
                avpol_temp(px(i,j, jj),py(i,j, jj),pz(i,j, jj),im)+&
                pro(i, jj)*vpol*vsol/(delta**3)/fv*ngpol(ii)*sc ! ngpol(ii) has the number of chains grafted to the point ii
            enddo

            q_tosend=q_tosend+pro(i, jj)

        enddo ! i
        ! norma 
        do im = 1, N_monomer
            do ix=1,dimx
                do iy=1,dimy
                    do iz=1,dimz
                        avpol_tosend(ix,iy,iz,im)=avpol_tosend(ix, iy, iz,im) + &
                            avpol_temp(ix,iy,iz,im)/q_tosend
                    enddo
                enddo
            enddo
        enddo

        q(ii) = q_tosend ! no la envia ahora

    enddo ! jj
    !------------------ MPI ----------------------------------------------
    !1. Todos al jefe


    call MPI_Barrier(MPI_COMM_WORLD, err)

    ! Junta avpol       
    call MPI_REDUCE(avpol_tosend, avpol, dimx*dimy*dimz*N_monomer, MPI_DOUBLE_PRECISION, &
        MPI_SUM,0, MPI_COMM_WORLD, err)
    
end subroutine calc_std


subroutine calcavpol(xpot)
    use mparameters_monomer
    use mkl
    use system
    implicit none

    real*8 xpot(dimx, dimy, dimz, N_monomer)

    if(flagmkl.eq.0)call calc_std(xpot)
#ifdef _MKL
    if(flagmkl.eq.1)call calc_mkl(xpot)
    if(flagmkl.eq.2)call calc_mkl_map(xpot)
#endif

end subroutine calcavpol


! selects appropiate  fcn function 
! ST_bctype == 1 : div flux scf determined of density 
! ST_bctype == 2 : 1 + different boundary conditions  
! ST_bctype == 3 : div flux scft determined by Slotboom transform of densities 
! ST_bctype == 4 : 3 + different boundary conditions 

subroutine set_fcn

    use fcnpointer
    use system, only : fluxflag, ST_bctype

    if(fluxflag.eq.0) then

        fcnptr => fcn_fkfun
    
    else if(fluxflag.eq.1) then    
    
        if(ST_bctype.eq.1) fcnptr => fcn_flux      
        !if(ST_bctype.eq.2) fcnptr => fcn_flux_bc
        
        if(ST_bctype.eq.3) fcnptr => fcn_flux_ST
        
        if(ST_bctype.eq.4) fcnptr => fcn_flux_ST_bc

    endif    
         
end subroutine set_fcn


end module fcnmod

