!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!    Free Energy Calculation...
!
!
!
subroutine Free_Energy_Calc(looped)

    use system
    use const
    use fields_fkfun
    use MPI
    use molecules
    use kai
    use bulk
    use results
    use ematrix
    use montecarlo
    use ellipsoid
    use transform
    use kaist
    use conformations
    use mparameters_monomer
    use mkl
    implicit none

    ! == input argement 
    integer, intent(inout) ::  looped  !! input : is case number ??

    ! == local variables

    real*8 :: qA_tosend(ncha), sumgaucheA_tosend(ncha)
    real*8 :: qA0(ncha), sumgaucheA0(ncha)
    integer :: newcuantasA0(ncha)
    real*8 :: F_Mix_s, F_Mix_pos
    real*8 :: F_Mix_neg, F_Mix_Hplus
    real*8 :: Free_energy2, sumpi, sumrho, sumel, sumdiel, suma, mupolA
    real*8 :: temp
    real*8 :: F_Mix_OHmin, F_gaucheA, F_ConfA, F_EqA, F_vdWA, F_epsA, F_electro
    real*8 :: proA0(cuantasA, maxcpp)         ! pro(cuantas, maxcpp))
    real*8 :: entropyA(dimx,dimy,dimz)
    character*5 :: title
    real*8 :: xtotalAsum(dimx,dimy,dimz)
    
    ! MPI
    !integer :: stat(MPI_STATUS_SIZE) 
    type(MPI_Status) :: stat

    integer :: source
    integer :: dest
    integer :: tag
    parameter(tag = 0)
    integer :: err

    ! Dummies
    integer :: ix, iy, iz, i, ii, ax, ay, az, jj
    integer :: jx, jy, jz,iii
    integer :: im, ip, ipp
    real*8 :: gradpsi2
    real*8 :: fv, fv2

    integer :: counter
    real*8 :: psiv(3)

    integer, external :: PBCSYMI, PBCREFI

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
    !
    !  Recupera pro(i) de todos los procesos para calculo de F
    !

    ! Subordinados

    entropyA = 0.0d0
    qA0 = 0.0d0
    qA_tosend = 0.0d0
    sumgaucheA_tosend = 0.0d0

    if(flagmkl.ne.0) then
        proA = 0.0
        do jj = 1, cpp(rank+1)
            iii = cppini(rank+1)+jj
            do i = 1, newcuantasA(iii)
                proA(i,jj) = promkl(iii)%pro(i)
            enddo ! i
        enddo ! jj
    endif

    if(rank.ne.0) then
        dest = 0
        ! Envia q

        do jj = 1, cpp(rank+1)
            iii = cppini(rank+1)+jj
            qA_tosend(iii) = q(iii)
        enddo

        call MPI_REDUCE(qA_tosend, qA0, ncha, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

        ! newcuantas
        
        call MPI_REDUCE(newcuantasA, newcuantasA0, ncha, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, err)

        ! Envia pro

        ! pro(cuantas, maxcpp)) 
        CALL MPI_SEND(proA, cuantasA*cpp(rank+1) , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD,err)
        
        ! sum gauche

        do jj = 1, cpp(rank+1)
            iii = cppini(rank+1)+jj
            sumgaucheA_tosend(iii) = 0.0
            do i = 1, newcuantasA(iii)
                sumgaucheA_tosend(iii) = sumgaucheA_tosend(iii)+ ngaucheA(i,iii)*proA(i,jj)/qA(iii)
            enddo ! i
        enddo ! jj

        call MPI_REDUCE(sumgaucheA_tosend, sumgaucheA0, ncha, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)


        goto 888

    endif

    Free_Energy = 0.0
    Free_Energy2 = 0.0

    ! 1. Mezcla solvente

    F_Mix_s = 0.0 

    do ix = 1, dimx
        do iy = 1, dimy
            do iz = 1, dimz
                fv=(1.0-volprot(ix,iy,iz))
                F_Mix_s = F_Mix_s + xh(ix, iy,iz)*(dlog(xh(ix, iy, iz))-1.0)*fv
                F_Mix_s = F_Mix_s - xsolbulk*(dlog(xsolbulk)-1.0)*fv
            enddo      
        enddo      
    enddo      
    F_Mix_s = F_Mix_s * delta**3/vsol
    Free_Energy = Free_Energy + F_Mix_s

    ! 2. Mezcla ion positivo

    F_Mix_pos = 0.0 

    do ix = 1, dimx
        do iy = 1, dimy
            do iz = 1, dimz
      
                fv=(1.0-volprot(ix,iy,iz))

                F_Mix_pos = F_Mix_pos + xpos(ix, iy,iz) &
                    *(dlog(xpos(ix, iy, iz)/vsalt)-1.0-dlog(expmupos) + dlog(vsalt))*fv

                F_Mix_pos = F_Mix_pos - xposbulk &
                    *(dlog(xposbulk/vsalt)-1.0-dlog(expmupos) + dlog(vsalt))*fv

            enddo
        enddo
    enddo
    F_Mix_pos = F_Mix_pos * delta**3/vsol/vsalt
    Free_Energy = Free_Energy + F_Mix_pos

    ! 3. Mezcla ion negativo

    F_Mix_neg = 0.0

    do ix = 1, dimx
        do iy = 1, dimy
            do iz = 1, dimz

                fv=(1.0-volprot(ix,iy,iz))

                F_Mix_neg = F_Mix_neg + xneg(ix, iy,iz) &
                    *(dlog(xneg(ix, iy, iz)/vsalt)-1.0- dlog(expmuneg) + dlog(vsalt))*fv

                F_Mix_neg = F_Mix_neg - xnegbulk &
                    *(dlog(xnegbulk/vsalt)-1.0- dlog(expmuneg) + dlog(vsalt))*fv

            enddo 
        enddo 
    enddo 
    F_Mix_neg = F_Mix_neg * delta**3/vsol/vsalt
    Free_Energy = Free_Energy + F_Mix_neg

    ! 4. Mezcla protones

    F_Mix_Hplus = 0.0

    do ix = 1, dimx
        do iy = 1, dimy
            do iz = 1, dimz

                fv=(1.0-volprot(ix,iy,iz))

                F_Mix_Hplus = F_Mix_Hplus &
                +xHplus(ix, iy, iz)*(dlog(xHplus(ix,iy,iz))-1.0 -dlog(expmuHplus))*fv

                F_Mix_Hplus = F_Mix_Hplus &
                    -xHplusbulk*(dlog(xHplusbulk)-1.0 -dlog(expmuHplus))*fv

            enddo
        enddo
    enddo
    F_Mix_Hplus = F_Mix_Hplus * delta**3/vsol
    Free_Energy = Free_Energy + F_Mix_Hplus

    ! 5. Mezcla hidroxilos

    F_Mix_OHmin = 0.0

    do ix = 1, dimx
        do iy = 1, dimy
            do iz = 1, dimz

                fv=(1.0-volprot(ix,iy,iz))

                F_Mix_OHmin = F_Mix_OHmin + &
                    xOHmin(ix, iy,iz)*(dlog(xOHmin(ix, iy, iz))-1.0-dlog(expmuOHmin))*fv

                F_Mix_OHmin = F_Mix_OHmin - xOHminbulk*(dlog(xOHminbulk)-1.0-dlog(expmuOHmin))*fv

            enddo
        enddo
    enddo
    F_Mix_OHmin = F_Mix_OHmin * delta**3/vsol
    Free_Energy = Free_Energy + F_Mix_OHmin

    ! 6. Entropia interna polimero

    F_ConfA = 0.0

    ! Jefe

    if (rank.eq.0) then ! Igual tiene que serlo, ver arriba

        do jj = 1, cpp(rank+1)
            iii = jj
            qA_tosend(iii) = qA(iii)
        enddo

        call MPI_REDUCE(qA_tosend, qA0, ncha, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

        call MPI_REDUCE(newcuantasA, newcuantasA0, ncha, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, err)

        do jj = 1, cpp(rank+1)
            do i = 1, newcuantasA0(jj)
                iii = jj
      
                F_ConfA = F_ConfA + (proA(i, jj)/qA0(iii)) &
                    *dlog((proA(i, jj))/qA0(iii))*ngpol(iii)

                entropyA(p0(iii,1),p0(iii,2),p0(iii,3)) =  - dlog(qA0(iii)/shiftA) 
            enddo
        enddo 

        do ii = 2, size ! loop sobre los procesadores restantes 
                        ! == loop over the remaining processors

            source = ii-1
           ! pro(cuantas, maxcpp))
            call MPI_RECV(proA0, cuantasA*cpp(ii), MPI_DOUBLE_PRECISION, source, tag, MPI_COMM_WORLD, stat, err)

            do jj = 1, cpp(ii)
                !       write(stdout,*) ii, jj, pro0(10,jj)
                iii = cppini(ii)+jj
                do i = 1, newcuantasA0(iii)

                    F_ConfA = F_ConfA + (proA0(i, jj)/qA0(iii))*dlog((proA0(i, jj))/qA0(iii))*ngpol(iii)

                    entropyA(p0(iii,1),p0(iii,2),p0(iii,3)) =  - dlog(qA0(iii)/shiftA) 

                enddo
            enddo

        enddo ! ii

    endif ! rank  == end rank ==0 

    Free_Energy = Free_Energy + F_ConfA

    if(rank.eq.0) then

    !      title = 'entpy'
    !      call savetodisk(entropyA, title, looped)
 
        open (unit=8, file='entropyA.out', form='unformatted')
        write(8)dimx,dimy,dimz
        write(8)entropyA
        close(8)

    endif

    ! 6.5 Energy of gauche bonds

    F_gaucheA = 0.0

    ! Jefe

    if (rank.eq.0) then ! Igual tiene que serlo, ver arriba == Same process see above. Has to do the same thing 

        do jj = 1, cpp(rank+1)
            iii = cppini(rank+1)+jj
            sumgaucheA_tosend(iii) = 0.0
            do i = 1, newcuantasA(iii)
                sumgaucheA_tosend(iii) = sumgaucheA_tosend(iii)+ ngaucheA(i,iii)*proA(i,jj)/qA(iii)
            enddo ! i
        enddo ! jj

        call MPI_REDUCE(sumgaucheA_tosend, sumgaucheA0, ncha, &
            MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

        do ii = 1, ncha
            F_gaucheA = F_gaucheA + sumgaucheA0(ii)*ngpol(ii)*benergyA
        enddo  

    endif ! rank

    Free_Energy = Free_Energy + F_gaucheA

    if(rank.eq.0) then
        title = 'entpy'
        call savetodisk(entropyA, title, looped)
 
        open (unit=8, file='entropyA.out', form='unformatted')
        write(8)dimx,dimy,dimz
        write(8)entropyA
        close(8)
    endif

    
    ! 7. Chemical Equilibrium
    F_EqA = 0.0 
            

    do ix  = 1, dimx
        do iy  = 1, dimy
            do iz  = 1, dimz

                do im = 1, N_monomerA
      
                    fv=(1.0-volprot(ix,iy,iz))

                    if(zpolA(im).ne.0) then

                        F_EqA = F_EqA + fdisA(ix,iy,iz,im)*dlog(fdisA(ix,iy,iz,im)) &
                            *avpolA(ix,iy,iz,im)/vpolA*fv

                        F_EqA = F_EqA + (1.0d0-fdisA(ix,iy,iz,im)) &
                            *dlog(1.0d0-fdisA(ix,iy,iz,im))*avpolA(ix,iy,iz,im)/vpolA*fv

                        F_EqA = F_EqA + (1.0-fdisA(ix,iy,iz,im))*dlog(K0A(im))*avpolA(ix,iy,iz,im)/vpolA*fv

                        select case (zpolA(im))
                        case (-1) ! acid
                            F_EqA= F_EqA + &
                                (1.0d0-fdisA(ix,iy,iz,im))*(-dlog(expmuHplus))*avpolA(ix,iy,iz,im)/vpolA*fv
                        case (1) ! base
                            F_EqA = F_EqA + &
                                (1.0-fdisA(ix,iy,iz,im))*(-dlog(expmuOHmin))*avpolA(ix,iy,iz,im)/vpolA*fv
                        end select

                    endif ! zpol

                enddo ! im
   
            enddo
        enddo
    enddo

    F_eqA = F_eqA * delta**3/vsol

    Free_Energy = Free_Energy + F_EqA

    ! 8.vdW ! Ojo, los kai son negativos => atraccion

    F_vdWA = 0.0

    do ix = 1, dimx
        do iy = 1, dimy
            do iz = 1, dimz

                fv=(1.0-volprot(ix,iy,iz))

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
                                fv2 = (1.0-volprot(jx,jy,jz)) 

                                do ip = 1, N_poorsolA
                                do ipp = 1, N_poorsolA
                    
                                    F_vdWA = F_vdWA - 0.5d0*delta**3*xtotalA(ix,iy,iz,ip) &
                        *xtotalA(jx,jy,jz,ipp)*XuA(ax, ay, az)*st*st_matrixA(ip,ipp)*fv*fv2/(vpolA*vpolA*vsol*vsol)
                    
                                enddo ! ip
                                enddo ! ipp

                            endif
                            endif
                            endif

                        enddo
                    enddo
                enddo

            enddo
        enddo
    enddo

    Free_Energy = Free_Energy + F_vdWA

    ! 9. Electrostatic ! OJO

    F_electro = 0.0    

    do ix  = 1, dimx
        do iy  = 1, dimy
            do iz  = 1, dimz
                F_electro = F_electro & 
                    + delta**3*psi(ix, iy, iz)*qtot(ix, iy, iz)/2.0d0/vsol

            enddo
        enddo
    enddo
  
    print*, F_electro

    Free_Energy = Free_Energy + F_electro

    ! 10. Pol-prot

    F_epsA = 0.0 

    do ix = 1, dimx
        do iy = 1, dimy
            do iz = 1, dimz
                fv=(1.0-volprot(ix,iy,iz))
                do im = 1, N_monomerA
                    F_epsA = F_epsA - avpolA(ix,iy,iz,im)*volepsA(ix,iy,iz)*(delta**3)/vpolA/vsol*fv
                enddo
            enddo
        enddo
    enddo

    Free_Energy = Free_Energy + F_epsA

    if (verbose.ge.1) then
        write(stdout,*) 'Free_Energy_Calc: Free energy(1) = ', Free_energy
    endif

    ! minimal F

    Free_Energy2 = 0.0d0

    xtotalAsum = 0.0d0
    do ip = 1, N_poorsolA
        xtotalAsum(:,:,:)= xtotalAsum(:,:,:)+xtotalA(:,:,:,ip)
    enddo


    sumpi = 0.0d0
    sumrho=0.0d0
    sumel=0.0d0
    sumdiel = 0.0d0

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz

                fv=(1.0d0-volprot(ix,iy,iz))

                sumpi = sumpi+dlog(xh(ix, iy, iz))*fv     
                sumpi = sumpi-dlog(xsolbulk)*fv
     
                sumrho = sumrho + ( - xh(ix, iy, iz) -xHplus(ix, iy, iz) &
                - xOHmin(ix, iy, iz) - (xpos(ix, iy, iz)+xneg(ix, iy, iz))/vsalt)*fv! sum over  rho_i i=+,-,s


                sumrho = sumrho - ( - xsolbulk -xHplusbulk &
                    -xOHminbulk - (xposbulk+xnegbulk)/vsalt)*fv ! sum over  rho_i i=+,-,s

                sumel = sumel - qtot(ix, iy, iz)*psi(ix, iy, iz)/2.0d0 
      
                sumel = sumel + volq(ix,iy,iz)*psi(ix,iy,iz)*vsol                   


                psiv(1) = psi(ix+1,iy,iz)-psi(ix,iy,iz)
                psiv(2) = psi(ix,iy+1,iz)-psi(ix,iy,iz)
                psiv(3) = psi(ix,iy,iz+1)-psi(ix,iy,iz)

                gradpsi2 = DOT_PRODUCT(MATMUL(TMAT, psiv), MATMUL(TMAT, psiv))
         
                sumdiel = sumdiel + 0.5d0/constq*xtotalAsum(ix,iy,iz)*gradpsi2*Depsfcn(ix,iy,iz)

            enddo
        enddo
    enddo
         
    sumpi = (delta**3/vsol)*sumpi
    sumrho = (delta**3/vsol)*sumrho
    sumel = (delta**3/vsol)*sumel
    sumdiel = (delta**3/vsol)*sumdiel


    suma = sumpi + sumrho + sumel + sumdiel


    do ii = 1, ncha
        Free_Energy2 = Free_Energy2-dlog(qA0(ii)/shiftA)*ngpol(ii) 
    enddo

    Free_Energy2 = Free_Energy2 + suma - F_vdWA

    if (verbose.ge.1) then
        write(stdout,*) 'Free_Energy_Calc: Free energy(2) = ', Free_energy2, sumdiel
    endif

    ! Guarda energia libre


    mupolA = 0.0
    do ii = 1, ncha
        mupolA = mupolA - dlog(qA0(ii)/shiftA)*ngpol(ii)
    enddo

    temp = sum(ngpol)
    mupolA = mupolA/temp

    if(rank.eq.0) then

        write(301,*)looped, Free_energy
        flush(301)
        write(302,*)looped, F_Mix_s 
        write(303,*)looped, F_Mix_pos
        write(304,*)looped, F_Mix_neg
        write(305,*)looped, F_Mix_Hplus
        write(306,*)looped, F_Mix_OHmin
	    write(3071,*)looped, F_gaucheA
        write(307,*)looped, F_ConfA
        write(308,*)looped, F_EqA
        write(309,*)looped, F_vdWA
        write(410,*)looped, F_epsA
        write(311,*)looped, F_electro

        write(312,*)looped, Free_energy2

        write(313,*)looped, mupolA

    endif
    
    
888     call MPI_BCAST(free_energy, 1, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, err)

return

end




