subroutine fkfun(x,f,ier2)

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

    integer :: ncells
    real*8 :: proAtemp, proBtemp
    integer :: i,j, ix, iy, iz, ii, ax, ay, az
    integer :: im, ip
    integer :: jx, jy, jz, jj
    real*8 :: xpotA(dimx, dimy, dimz, N_monomerA)
    real*8 :: xpotB(dimx, dimy, dimz, N_monomerB)

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
    real*8 :: avpolA_temp(dimx,dimy,dimz,N_monomerA) 
    real*8 :: avpolB_temp(dimx,dimy,dimz,N_monomerB)
    real*8 :: qA_tosend, qB_tosend
    real*8 :: gradpsi2
    real*8 :: fv

    ! hamiltonian inception
    real*8 :: hfactor, hd
    real*8 :: hds(100)


    hds = -1

    !-----------------------------------------------------
    ! Common variables

    shiftA = 1.0d0

    ncells = dimx*dimy*dimz ! numero de celdas == number of cells

    ! Jefe

    if(rank.eq.0) then ! llama a subordinados y pasa vector x  == calls subordinates and passes vector x
        flagsolver = 1
        CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
        CALL MPI_BCAST(x, eqs*ncells , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
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

                do ip = 1, N_poorsolA
                    xtotalA(ix,iy,iz,ip) = x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ ip*ncells) !fraccion polimero de tipo ip
                enddo

                do ip = 1, N_poorsolB
                    xtotalB(ix,iy,iz,ip) = x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ (ip+N_poorsolA)*ncells) !fraccion polimero de tipo ip
                enddo

                if(electroflag.eq.1) psi(ix,iy,iz)=&
                    x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsolA+N_poorsolB+1)*ncells)   !potencial electrostatico
               
                if(fluxflag.eq.1) then 

                    xpos(ix, iy, iz)  =  x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsolA+N_poorsolB+2)*ncells )    
                    xneg(ix, iy, iz)  =  x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsolA+N_poorsolB+3)*ncells ) 
                    xHplus(ix, iy, iz)=  x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsolA+N_poorsolB+4)*ncells )    
                    xOHmin(ix, iy, iz)=  x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsolA+N_poorsolB+5)*ncells )   

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

    fdisA = 0.0d0    ! == added  d0 
    avpolA = 0.0d0
    fdisB = 0.0d0    
    avpolB = 0.0d0

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                if(fluxflag.eq.0) then ! Equilibrium  

                    xpos(ix, iy, iz)   = expmupos*(xh(ix, iy, iz)**vsalt)*dexp(-psi(ix, iy, iz)*zpos) ! ion plus volume fraction vsalt=vsal/vsv
                    xneg(ix, iy, iz)   = expmuneg*(xh(ix, iy, iz)**vsalt)*dexp(-psi(ix, iy, iz)*zneg) ! ion neg volume fraction
                    xHplus(ix, iy, iz) = expmuHplus*(xh(ix, iy, iz))*dexp(-psi(ix, iy, iz))           ! H+ volume fraction
                    xOHmin(ix, iy, iz) = expmuOHmin*(xh(ix,iy,iz))*dexp(+psi(ix,iy,iz))               ! OH-  volume fraction
                
                endif 

                do im =1,N_monomerA

                    if (zpolA(im).eq.1) then !BASE
                        fdisA(ix,iy,iz,im) = 1.0d0 /(1.0d0 + xOHmin(ix,iy,iz)/(K0A(im)*xh(ix,iy,iz))) !k0 k en fraccion de volumen
                    else if (zpolA(im).eq.-1) then !ACID
                        fdisA(ix,iy,iz,im) = 1.0d0 /(1.0d0 + xHplus(ix,iy,iz)/(K0A(im)*xh(ix,iy,iz)))
                    endif

                enddo

                do im =1,N_monomerB

                    if (zpolB(im).eq.1) then !BASE
                        fdisB(ix,iy,iz,im) = 1.0d0 /(1.0d0 + xOHmin(ix,iy,iz)/(K0B(im)*xh(ix,iy,iz))) !k0 k en fraccion de volumen
                    else if (zpolB(im).eq.-1) then !ACID
                        fdisB(ix,iy,iz,im) = 1.0d0 /(1.0d0 + xHplus(ix,iy,iz)/(K0B(im)*xh(ix,iy,iz)))
                    endif

                enddo

            enddo
        enddo  
    enddo



    ! Compute dielectric permitivity

    xtotalsum = 0.0d0 
    ! sum of all A polymers
    do ip = 1, N_poorsolA
        xtotalsum(:,:,:) = xtotalsum(:,:,:) + xtotalA(:,:,:,ip)
    enddo
    ! sum of all B polymers
    do ip = 1, N_poorsolB
        xtotalsum(:,:,:) = xtotalsum(:,:,:) + xtotalB(:,:,:,ip)
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

    ! Calcula xpotA

    sttemp = st/(vpolA*vsol)

    do im = 1, N_monomerA ! loop over different monomer types

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

                    fv = (1.0d0 - volprot(ix,iy,iz)) 
                    !   fraccion de volumen de la celda que es sc volprot->fraccion pared
                    !   == volume fraction of the cell that is sc volprot->wall fraction             

                    xpotA(ix, iy, iz, im) = xh(ix,iy,iz)**vpolA 
                    
                    ! im:tipo de segmento, término de presion osmotica
                    ! == im: segment type, osmotic pressure term

                    xpotA(ix, iy, iz, im) = xpotA(ix,iy,iz, im)*dexp(volepsA(ix,iy,iz))  
                    
                    ! termino de interaccion con sup de la particula
                    ! == interaction term with the particle's sup =surface ???

                    ! Electrostatics

                    if(zpolA(im).ne.0.0) then
                        xpotA(ix,iy,iz,im) =  xpotA(ix,iy,iz,im)/fdisA(ix,iy,iz,im)*dexp(-psi(ix,iy,iz)*zpolA(im))  
                        ! fdis: por eq ac. base...  
                    endif
        
                    ! Dielectrics

                    gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy,iz))**2+&
                        (psi(ix,iy,iz+1)-psi(ix,iy,iz))**2 

                    !     gradpsi2 = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy-1,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz-1))**2 
                    !     xpot(ix, iy, iz) = xpot(ix,iy,iz)*exp(-Depsfcn(ix,iy,iz)*(gradpsi2)*constqE)

                    xpotA(ix,iy,iz,im) = xpotA(ix,iy,iz,im)*exp(Depsfcn(ix,iy,iz)*(gradpsi2)/constq/2.0d0*vpolA/fv)

                    ! Poor solvent depende de la grilla donde esta y de sus vecinos
                    ! == Poor solvent depends on the grid where it is and its neighbors



                    if(hydrophA(im).ne.0) then

                    proAtemp=0.0

                    do ax = -XulimitA,XulimitA 
                        do ay = -XulimitA,XulimitA
                            do az = -XulimitA,XulimitA

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
                                            fv = (1.0d0-volprot(jx,jy,jz))

                                            do ip = 1, N_poorsolA
                                                proAtemp = proAtemp + hfactor*XuA(ax,ay,az)*&
                                                st_matrixA(hydrophA(im),ip)*sttemp*xtotalA(jx,jy,jz,ip)*fv
                                            enddo ! ip

                                        endif
                                    endif
                                endif

                            enddo
                        enddo
                    enddo

                    xpotA(ix,iy,iz,im) = xpotA(ix,iy,iz,im)*dexp(proAtemp)

                    endif ! hydrph

!                    write(567,*)xpot(ix,iy,iz,1)

                enddo ! ix
            enddo ! iy
        enddo !iz

    enddo ! N_monomer

 ! Calcula xpotA

    sttemp = st/(vpolB*vsol)

    do im = 1, N_monomerB ! loop over different monomer types

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

                    fv = (1.0d0 - volprot(ix,iy,iz)) 
                    !   fraccion de volumen de la celda que es sc volprot->fraccion pared
                    !   == volume fraction of the cell that is sc volprot->wall fraction             

                    xpotB(ix, iy, iz, im) = xh(ix,iy,iz)**vpolB 
                    
                    ! im:tipo de segmento, término de presion osmotica
                    ! == im: segment type, osmotic pressure term

                    xpotB(ix, iy, iz, im) = xpotB(ix,iy,iz, im)*dexp(volepsB(ix,iy,iz))  
                    
                    ! termino de interaccion con sup de la particula
                    ! == interaction term with the particle's sup =surface ???

                    ! Electrostatics

                    if(zpolA(im).ne.0.0) then
                        xpotB(ix,iy,iz,im) =  xpotB(ix,iy,iz,im)/fdisB(ix,iy,iz,im)*dexp(-psi(ix,iy,iz)*zpolB(im))  
                        ! fdis: por eq ac. base...  
                    endif
        
                    ! Dielectrics

                    gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy,iz))**2+&
                        (psi(ix,iy,iz+1)-psi(ix,iy,iz))**2 

                    !     gradpsi2 = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy-1,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz-1))**2 
                    !     xpot(ix, iy, iz) = xpot(ix,iy,iz)*exp(-Depsfcn(ix,iy,iz)*(gradpsi2)*constqE)

                    xpotB(ix,iy,iz,im) = xpotB(ix,iy,iz,im)*exp(Depsfcn(ix,iy,iz)*(gradpsi2)/constq/2.0d0*vpolA/fv)

                    ! Poor solvent depende de la grilla donde esta y de sus vecinos
                    ! == Poor solvent depends on the grid where it is and its neighbors



                    if(hydrophB(im).ne.0) then

                    proBtemp=0.0

                    do ax = -XulimitB,XulimitB 
                        do ay = -XulimitB,XulimitB
                            do az = -XulimitB,XulimitB

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
                                            fv = (1.0d0-volprot(jx,jy,jz))

                                            do ip = 1, N_poorsolB
                                                proBtemp = proBtemp + hfactor*XuB(ax,ay,az)*&
                                                st_matrixB(hydrophB(im),ip)*sttemp*xtotalB(jx,jy,jz,ip)*fv
                                            enddo ! ip

                                        endif
                                    endif
                                endif

                            enddo
                        enddo
                    enddo

                    xpotB(ix,iy,iz,im) = xpotB(ix,iy,iz,im)*dexp(proBtemp)

                    endif ! hydrph

!                    write(567,*)xpot(ix,iy,iz,1)

                enddo ! ix
            enddo ! iy
        enddo !iz

    enddo ! N_monomer

 



    !!!!!!!!!!!!!!!!!!!!!! Calculate pro from xpot !!!!!!!!!!!!!
    call calcavpolA(xpotA)

    call calcavpolB(xpotB)


    
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

            do im = 1, N_monomerA
                qtot(ix, iy, iz) =  qtot(ix,iy,iz) + avpolA(ix,iy,iz,im)*zpolA(im)/vpolA*fdisA(ix,iy,iz,im)
            enddo

            do im = 1, N_monomerB
                qtot(ix, iy, iz) =  qtot(ix,iy,iz) + avpolB(ix,iy,iz,im)*zpolB(im)/vpolB*fdisB(ix,iy,iz,im)
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

                do im = 1, N_monomerA
                    f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) = f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) +&
                         avpolA(ix,iy,iz,im) !packing ...+polimero
                enddo

                do im = 1, N_monomerB
                    f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) = f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) +&
                         avpolB(ix,iy,iz,im) !packing ...+polimero
                enddo

                ! write(123,*)ix,iy,iz,avpol(ix,iy,iz,1),xh(ix,iy,iz)
            enddo
        enddo
    enddo

    ! Poor solvent chain A 

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                do ip = 1, N_poorsolA
                    f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells) = xtotalA(ix,iy,iz,ip)

                    do im = 1, N_monomerA
                        if(hydrophA(im).eq.ip) then 
                            f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells) = f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells) - &
                                avpolA(ix,iy,iz,im)
                        endif
                    enddo ! im
                enddo ! ip

            enddo ! ix
        enddo ! iy
    enddo ! iz

    !  Poor solvent chain B 

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                do ip = 1, N_poorsolB
                    f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells) = xtotalB(ix,iy,iz,ip)

                    do im = 1, N_monomerB
                        if(hydrophB(im).eq.ip) then 
                            f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(ip+N_poorsolA)*ncells) = &
                            f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(ip+N_poorsolB)*ncells) -avpolB(ix,iy,iz,im)
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

                    f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsolA+N_poorsolB+1)*ncells)=&
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

        noffset=(N_poorsolA+N_poorsolB+2)*ncells 
        

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
        noffset=(N_poorsolA+N_poorsolB+1)*ncells 
        do i= 1, ncells    
            normel = normel +f(i+noffset)**2
        enddo       
    endif
    
    iter = iter + 1
    if(verbose.ge.3) then
        if(rank.eq.0) write(stdout,*)'fkfun:', iter, sqrt(norma), sqrt(normvol), sqrt(normel), qA(1)
    endif

    3333 continue
    ier2 = 0.0 

    return

end subroutine fkfun


subroutine calc_stdA(xpotA)

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

    real*8, intent(in) :: xpotA(dimx, dimy, dimz, N_monomerA)

    ! local variables 

    real*8 :: avpolA_tosend(dimx,dimy,dimz, N_monomerA)
    real*8 :: fv
    real*8 :: qA_tosend
    real*8 :: avpolA_temp(dimx,dimy,dimz,N_monomerA)
    integer :: im,jj,i,j, ix, iy, iz, ii, ax, ay, az
    ! MPI
    integer :: tag
    parameter(tag = 0)
    integer :: err

    shiftA = 1.0d0                   ! == added  d0  ! uniform shift in P(alpha)
    avpolA_tosend = 0.0d0
    qA = 0.0d0

    do jj = 1, cpp(rank+1)          ! == graft point on each node
        ii = cppini(rank+1)+jj      ! == actual graft point ii 

        qA_tosend = 0.0d0
        avpolA_temp = 0.0d0

        if(hasGraftA(ii)) then 

            do i=1,newcuantasA(ii)       ! == loop of chains conformation for graft point ii 
        
                proA(i, jj)= shiftA       ! == pro of conf i belong to graftpoint jj 
                

                do j=1,nsegA           
                    ax = pxA(i, j, jj)   ! == each to his own chain
                    ay = pyA(i, j, jj)
                    az = pzA(i, j, jj)
                    proA(i, jj) = proA(i, jj) * xpotA(ax, ay, az, segtypeA(j))
                enddo
                
                proA(i, jj) = proA(i, jj) * dexp(-benergyA*ngaucheA(i,ii)) ! == energy of gauche bonds
                proA(i, jj) = proA(i, jj) * dexp(-fz*zfinalA(i,jj))       ! == terminal end energy Fz

                do j=1,nsegA

                    fv = fvstd(pxA(i,j, jj),pyA(i,j, jj),pzA(i,j, jj))
                    im = segtypeA(j)
                    avpolA_temp(pxA(i,j, jj),pyA(i,j, jj),pzA(i,j, jj),im)= &
                    avpolA_temp(pxA(i,j, jj),pyA(i,j, jj),pzA(i,j, jj),im)+&
                        proA(i, jj)*vpolA*vsol/(delta**3)/fv*ngpol(ii)*sc ! ngpol(ii) has the number of chains grafted to the point ii
                
                enddo

                qA_tosend=qA_tosend+proA(i, jj)
            
            enddo ! == end loop conf
      
            ! norma 
            do im = 1, N_monomerA
                do iz=1,dimz
                    do iy=1,dimy
                        do ix=1,dimx
                            avpolA_tosend(ix,iy,iz,im) = avpolA_tosend(ix, iy, iz,im) + &
                                avpolA_temp(ix,iy,iz,im)/qA_tosend
                        enddo
                    enddo
                enddo
            enddo

        endif   

        qA(ii) = qA_tosend ! no la envia ahora

    enddo ! jj
    !------------------ MPI ----------------------------------------------
    !1. Todos al jefe


    call MPI_Barrier(MPI_COMM_WORLD, err)

    ! Junta avpol       
    call MPI_REDUCE(avpolA_tosend, avpolA, dimx*dimy*dimz*N_monomerA, MPI_DOUBLE_PRECISION, &
        MPI_SUM,0, MPI_COMM_WORLD, err)
        
end subroutine calc_stdA



subroutine calc_stdB(xpotB)

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

    real*8, intent(in) :: xpotB(dimx, dimy, dimz, N_monomerB)

    ! local variables 

    real*8 :: avpolB_tosend(dimx,dimy,dimz, N_monomerB)
    real*8 :: fv
    real*8 :: qB_tosend
    real*8 :: avpolB_temp(dimx,dimy,dimz,N_monomerB)
    integer :: im,jj,i,j, ix, iy, iz, ii, ax, ay, az
    ! MPI
    integer :: tag
    parameter(tag = 0)
    integer :: err

    shiftB = 1.0d0                   ! == added  d0  ! uniform shift in P(alpha)
    avpolB_tosend = 0.0d0
    qB = 0.0d0

    do jj = 1, cpp(rank+1)          ! == graft point on each node
        ii = cppini(rank+1)+jj      ! == actual graft point ii 

        qB_tosend = 0.0d0
        avpolB_temp = 0.0d0

        if(hasGraftB(ii)) then 

            do i=1,newcuantasB(ii)       ! == loop of chains conformation for graft point ii 
        
                proB(i, jj)= shiftB       ! == pro of conf i belong to graftpoint jj 
                
                do j=1,nsegB           
                    ax = pxB(i, j, jj)   ! == each to his own chain
                    ay = pyB(i, j, jj)
                    az = pzB(i, j, jj)
                    proB(i, jj) = proB(i, jj) * xpotB(ax, ay, az, segtypeB(j))
                enddo
                
                proB(i, jj) = proB(i, jj) * dexp(-benergyB*ngaucheB(i,ii)) ! == energy of gauche bonds
                proB(i, jj) = proB(i, jj) * dexp(-fz*zfinalB(i,jj))       ! == terminal end energy Fz

                do j=1,nsegB

                    fv = fvstd(pxB(i,j, jj),pyB(i,j, jj),pzB(i,j, jj))
                    im = segtypeB(j)
                    avpolB_temp(pxB(i,j, jj),pyB(i,j, jj),pzB(i,j, jj),im)= &
                    avpolB_temp(pxB(i,j, jj),pyB(i,j, jj),pzB(i,j, jj),im)+&
                        proB(i, jj)*vpolB*vsol/(delta**3)/fv*ngpol(ii)*sc ! ngpol(ii) has the number of chains grafted to the point ii
                
                enddo

                qB_tosend=qB_tosend+proB(i, jj)

            enddo ! == end loop conf
            ! norma 
            do im = 1, N_monomerB
                do iz=1,dimz
                    do iy=1,dimy
                        do ix=1,dimx
                            avpolB_tosend(ix,iy,iz,im) = avpolB_tosend(ix, iy, iz,im) + &
                                avpolB_temp(ix,iy,iz,im)/qB_tosend
                        enddo
                    enddo
                enddo
            enddo

        endif 

        qB(ii) = qB_tosend ! no la envia ahora

    enddo ! jj
    !------------------ MPI ----------------------------------------------
    !1. Todos al jefe


    call MPI_Barrier(MPI_COMM_WORLD, err)

    ! Junta avpol       
    call MPI_REDUCE(avpolB_tosend, avpolB, dimx*dimy*dimz*N_monomerB, MPI_DOUBLE_PRECISION, &
        MPI_SUM,0, MPI_COMM_WORLD, err)
        
end subroutine calc_stdB


subroutine calcavpolA(xpotA)

    use mparameters_monomer
    use mkl
    use system
    implicit none

    real*8 :: xpotA(dimx, dimy, dimz, N_monomerA)

    if(flagmkl.eq.0)call calc_stdA(xpotA)
#ifdef _MKL
    if(flagmkl.eq.1)call calc_mklA(xpotA)
    if(flagmkl.eq.2)call calc_mklA_map(xpotA)
#endif

end subroutine calcavpolA



subroutine calcavpolB(xpotB)

    use mparameters_monomer
    use mkl
    use system
    implicit none

    real*8 :: xpotB(dimx, dimy, dimz, N_monomerB)

    if(flagmkl.eq.0)call calc_stdB(xpotB)
#ifdef _MKL
    if(flagmkl.eq.1)call calc_mklB(xpotB)
    if(flagmkl.eq.2)call calc_mklB_map(xpotB)
#endif

end subroutine calcavpolB


