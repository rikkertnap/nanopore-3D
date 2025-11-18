
subroutine initmpi

    use MPI
!    use chainsdat

    implicit none

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

end subroutine initmpi

subroutine initconst

    use const
    use molecules
    use ellipsoid
    use mparameters_monomer
    use MPI, only : rank

    implicit none

    real*8 :: pi_my=3.141592653589793d0

    !pi = acos(-1.0)  
    pi = acos(-1.0d0)

    if(abs(pi-pi_my)>0.0d0) then ! check accuracy 
        if(rank.eq.0)write(stdout,*) 'init const: pi',pi,' pi_my', pi_my
    endif    

    lb = 0.714d0 ! bjerrum lenght in nm
    zpos = 1.0d0
    zneg = -1.0d0
    vsol = vsol0
    vsalt = ((4.0d0/3.0d0)*pi*(0.2d0)**3)/vsol  ! volume salt in units of vsol 0.2=radius salt  
    constq = delta*delta*4.0d0*pi*lb/vsol   ! multiplicative factor in poisson eq  
    pKw = 14.0d0
    Kw = 10.0d0**(-pKw)
    error = 1e-4 ! para comparar con la norma... ! == to compare to norm
    errel = 1d-6
    itmax = 200

    ! == eqs number of equation in unit of lattice size

    if(electroflag.eq.0) eqs = (1+N_poorsolA+N_poorsolB) 
    if(electroflag.eq.1) then 
        if (fluxflag.eq.0) then 
            eqs = (2+N_poorsolA+N_poorsolB)    !== Equilbrium 
        else if(fluxflag.eq.1) then
            eqs = (2+4+N_poorsolA)       !== Steady state : 4 iontypes
        endif
    endif     
         
end subroutine

! == init Input-dependent variables
 
subroutine initall

    use molecules
    use const
    use bulk
    use MPI
    use ellipsoid
    use chainsdat
    use inputtemp
    use mparameters_monomer
    use flux, only : init_flux_var, allocate_flux_var
    
    implicit none
    integer :: im

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Open common files
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    if(rank.eq.0) then
        open(unit=301, file='F_tot.dat', access='APPEND')
        open(unit=302, file='F_mixs.dat',  access='APPEND')
        open(unit=303, file='F_mixpos.dat',  access='APPEND')
        open(unit=304, file='F_mixneg.dat',  access='APPEND')
        open(unit=305, file='F_mixH.dat',  access='APPEND')
        open(unit=306, file='F_mixOH.dat',  access='APPEND')
        open(unit=307, file='F_conf.dat',  access='APPEND')
        open(unit=3071, file='F_gauche.dat',  access='APPEND')
        open(unit=308, file='F_eq.dat',  access='APPEND')
        open(unit=309, file='F_vdW.dat',  access='APPEND')
        open(unit=410, file='F_eps.dat',  access='APPEND')
        open(unit=311, file='F_electro.dat',  access='APPEND')
        open(unit=312, file='F_tot2.dat',  access='APPEND')
        open(unit=314, file='F_mixpos2.dat',  access='APPEND')
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Input-dependent variables
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    vpolA = vpolA/vsol ! vpol in units of vsol
    vpolB = vpolB/vsol
    constqE = vpolA/(2.0d0*constq)
    print*,"Warning constqE involvs vpol"

    dielW = 78.54d0
    dielPr = dielP/dielW
    dielSr = dielS/dielW

    call initbulk

    do im = 1, N_monomerA
        KaA(im)=10.0d0**(-pKaA(im))
        select case (zpolA(im))
        case (-1) ! acid
            K0A(im) = (KaA(im)*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant, Ka
        case (1) ! base
            K0A(im) = ((Kw/KaA(im))*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant, Kb 
        end select
    enddo

     do im = 1, N_monomerB
        KaB(im)=10.0d0**(-pKaB(im))
        select case (zpolB(im))
        case (-1) ! acid
            K0B(im) = (KaB(im)*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant, Ka
        case (1) ! base
            K0B(im) = ((Kw/KaB(im))*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant, Kb 
        end select
    enddo

    if(fluxflag.eq.1) then
        ! allocate flux variables
        call allocate_flux_var() 
        ! init  flux/steady state related variables
        call init_flux_var()
    endif    

end subroutine initall


subroutine initbulk 

    use molecules, only : vsol, vsalt, zpos,zneg  
    use inputtemp, only : xsalt, pHbulk, csalt    
    use const, only : pi, Na, pKw
    use bulk, only : expmupos, expmuneg, expmuHplus, expmuOHmin
    use bulk, only : xsolbulk, xposbulk, xnegbulk, xHplusbulk, xOHminbulk

    cHplus = 10.0d0**(-pHbulk)                ! concentration H+ in bulk
    xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
    pOHbulk = pKw -pHbulk
    cOHmin = 10.0d0**(-pOHbulk)               ! concentration OH- in bulk
    xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol  
    xsalt =(csalt*Na/(1.0d24))*(vsalt*vsol)    ! volume fraction salt,csalt in mol/l 

    if(pHbulk.le.7) then  ! pH<= 7
        xposbulk = xsalt/zpos
        xnegbulk = -xsalt/zneg+(xHplusbulk -xOHminbulk) *vsalt ! NaCl+ HCl  
    else                  ! pH >7 
        xposbulk = xsalt/zpos +(xOHminbulk -xHplusbulk) *vsalt ! NaCl+ NaOH   
        xnegbulk = -xsalt/zneg 
    endif

    xsolbulk = 1.0d0 -xHplusbulk -xOHminbulk -xnegbulk -xposbulk 

    expmupos = xposbulk /xsolbulk**vsalt
    expmuneg = xnegbulk /xsolbulk**vsalt
    expmuHplus = xHplusbulk /xsolbulk   ! vsol = vHplus 
    expmuOHmin = xOHminbulk /xsolbulk   ! vsol = vOHmin 

end subroutine initbulk

! == close all open 301.. 313 file and ends mpi and stop program
subroutine endall
    use MPI
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!
    ! Close common files
    !!!!!!!!!!!!!!!!!!!!!!

    close(301)
    close(302)
    close(303)
    close(304)
    close(305)
    close(306)
    close(307)
    close(3071)
    close(308)
    close(309)
    close(310)
    close(311)
    close(312)
    close(313)

    call MPI_FINALIZE(ierr) ! finaliza MPI    
    stop

end subroutine



subroutine savedata(cccc)
    use system
    use results
    use const
    use molecules
    use chainsdat
    use kai
    use ematrix
    use fields_fkfun
    use MPI
    use kinsol
    use kaist
    use mparameters_monomer
    use channelcurved, only :radiusC, radiusL, Lengthchannel, total_surface_area_curv, lengthsection, nsections
    use inputtemp, only : csalt, pHbulk

    implicit none

    integer, intent(in) :: cccc

    character*20 :: filename
    character*5  :: title
    real*8 :: temp(dimx,dimy,dimz)
    real*8 :: sumpolA, sumpolB
    real*8 :: avfdisA(N_monomerA), avfdisB(N_monomerB)
    real*8 :: sumavpolA, sumavpolB 
    integer :: ix,iy,iz, im
    real*8 :: fv
    real*8 :: area

    !----------------------------------------------------------
    !  OUTPUT
    !----------------------------------------------------------

    if(rank.eq.0) then 

     ! == save files 
        ! == polymer 

        temp = 0.0d0
        do im = 1, N_monomerA
            temp(:,:,:) =  temp(:,:,:) + avpolA(:,:,:, im)*(1.0d0 - volprot(:,:,:))
        enddo

        title = 'xpolA' ! 'avpol' 
        call savetodisk(temp, title, cccc)

        temp = 0.0d0
        do im = 1, N_monomerB
            temp(:,:,:) =  temp(:,:,:) + avpolB(:,:,:, im)*(1.0d0 - volprot(:,:,:))
        enddo

        title = 'xpolB' ! 'avpol' 
        call savetodisk(temp, title, cccc)

        ! == polymer, by type 
        
        do im = 1, N_monomerA
            temp(:,:,:) = avpolA(:,:,:,im)*(1.0d0 - volprot(:,:,:))
            write(title,'(A3, I2.2)')'xAp',im      ! avp 
            call savetodisk(temp, title, cccc)
        enddo

        do im = 1, N_monomerB
            temp(:,:,:) = avpolB(:,:,:,im)*(1.0d0 - volprot(:,:,:))
            write(title,'(A3, I2.2)')'xBp',im      ! avp 
            call savetodisk(temp, title, cccc)
        enddo

        ! == solvent
        temp(:,:,:) = xh(:,:,:)*(1.0 - volprot(:,:,:))

        title = 'avsol'
        call savetodisk(temp, title, cccc)
        
        
        ! Cationes
        ! title = 'avpos'
        ! call savetodisk(xpos, title, cccc)
        ! Aniones
        ! title = 'avneg'
        ! call savetodisk(xneg, title, cccc)
        ! H+
        !  title = 'avHpl'
        !  call savetodisk(xHplus, title, cccc)
        ! OH-
        !  title = 'avOHm'
        !  call savetodisk(xOHmin, title, cccc)
        ! fdis
        
        title = 'fdisA' 
        temp(1:dimx,1:dimy, 1:dimz) = fdisA(1:dimx,1:dimy, 1:dimz,1)
        call savetodisk(temp, title, cccc)

        ! polymer charge

        !temp = 0.0d0

        !do ix=1,dimx
        !    do iy=1,dimy
        !        do iz=1,dimz
        !           fv = (1.0d0-volprot(ix,iy,iz))
        !            do im = 1, N_monomer
        !                temp(ix, iy, iz) = temp(ix,iy,iz) + &
        !                    avpol(ix,iy,iz,im)*zpol(im)/vpol/vsol*fdis(ix,iy,iz,im)! units of |e|/nm^3 
        !            enddo
        !        enddo
        !    enddo
        ! enddo

        !  title = 'qpol_'
        !  call savetodisk(temp, title, cccc)

        ! electostatic potential 

        temp(1:dimx,1:dimy, 1:dimz) = psi(1:dimx,1:dimy, 1:dimz)

        title = 'poten'
        call savetodisk(temp, title, cccc)


        !  Particle
        !  title = 'avpar'
        !  call savetodisk(volprot, title, cccc)

        ! save volprot for supercell

        if(rank.eq.0) then
            open (unit=8, file='out.par', form='unformatted')
            do ix=1,dimx
                do iy=1,dimy
                    do iz=1,dimz
                        xpar(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) = volprot(ix,iy,iz)
                    enddo
                enddo
            enddo
            write(8)xpar
            close(8)
        endif

        ! == system

        ! == total number of segments of chain A
        sumpolA = 0.0d0  
        do im = 1, N_monomerA
            do iz = 1, dimz
                do iy = 1, dimy
                    do ix= 1, dimx
                        sumpolA = sumpolA + avpolA(ix,iy,iz,im)*(delta**3)*(1.0d0-volprot(ix,iy,iz))/vpolA/vsol
                    enddo
                enddo
            enddo
        enddo

        ! == total number of segments of chain B
        sumpolB = 0.0d0  
        do im = 1, N_monomerB
            do iz = 1, dimz
                do iy = 1, dimy
                    do ix= 1, dimx
                        sumpolB = sumpolB + avpolB(ix,iy,iz,im)*(delta**3)*(1.0d0-volprot(ix,iy,iz))/vpolB/vsol
                    enddo
                enddo
            enddo
        enddo


        ! == average fraction of charged monomers of type im of chain A
      
        do im = 1, N_monomerA
            if (zpolA(im).ne.0) then 
                avfdisA(im) = 0.0d0
                sumavpolA = 0.0d0
                do iz=1,dimz
                    do iy=1,dimy
                        do ix=1,dimz
                            fv = (1.0d0-volprot(ix,iy,iz))
                            avfdisA(im)= avfdisA(im)+ avpolA(ix,iy,iz,im)*fv*zpolA(im)/vpolA/vsol*fdisA(ix,iy,iz,im) ! units of |e|/nm^3 
                            sumavpolA= sumavpolA+avpolA(ix,iy,iz,im)*fv/vpolA/vsol      
                        enddo
                    enddo
                enddo
                avfdisA(im)= avfdisA(im)/sumavpolA
            else  
                avfdisA(im) = 0.0d0
            endif      
        enddo

        ! == average fraction of charged monomers of type im of chain A
      
        do im = 1, N_monomerB
            if (zpolB(im).ne.0) then 
                avfdisB(im) = 0.0d0
                sumavpolB = 0.0d0
                do iz=1,dimz
                    do iy=1,dimy
                        do ix=1,dimz
                            fv = (1.0d0-volprot(ix,iy,iz))
                            avfdisB(im)= avfdisB(im)+ avpolB(ix,iy,iz,im)*fv*zpolB(im)/vpolB/vsol*fdisB(ix,iy,iz,im) ! units of |e|/nm^3 
                            sumavpolB = sumavpolB+avpolB(ix,iy,iz,im)*fv/vpolB/vsol      
                        enddo
                    enddo
                enddo
                avfdisB(im)= avfdisB(im)/sumavpolB
            else  
                avfdisB(im) = 0.0d0
            endif      
        enddo


        if(curvedflag==0) area=dimx*dimy*delta*delta      ! == straight nanopore 
        if(curvedflag==1) then 
            area = total_surface_area_curv(radiusL,radiusC,lengthsection,nsections) !== curved nanopore
        endif

    

        ! system
        if(curvedflag==0) then 
            area=dimx*dimy*delta*delta      ! == straight nanopore 
            print*,"Warning area update needed for curvflag=0"
        endif    
        if(curvedflag==1) then 
            area = total_surface_area_curv(radiusL,radiusC,lengthsection,nsections) !== curved nanopore
        endif

        write(filename,'(A7, I3.3, A4)')'system.', cccc, '.dat'
        
        open (unit=310, file=filename)
        write(310,*)'GIT version = ', _VERSION
        write(310,*)'fnorm       = ',norma      ! residual size of iteration vector
        write(310,*)'lsegA       = ',lsegA       
        write(310,*)'delta       = ',delta
        write(310,*)'dimx        = ',dimx
        write(310,*)'dimy        = ',dimy
        write(310,*)'dimz        = ',dimz
        write(310,*)'vsol        = ',vsol
        write(310,*)'vsalt       = ',vsalt*vsol
        write(310,*)'vpolA       = ',vpolA*vsol
        write(310,*)'vpolB       = ',vpolB*vsol
        write(310,*)'pKw         = ',pKw
        write(310,*)'zpos        = ',zpos
        write(310,*)'zneg        = ',zneg
        write(310,*)'nsegA       = ',nsegA
        write(310,*)'nsegB       = ',nsegB
        write(310,*)'csalt       = ',csalt
        write(310,*)'pH          = ',pHbulk
        write(310,*)'iterations  = ',iter
        write(310,*)'sigma cad/nm2 = ',ncha/area
        write(310,*)'kaiA        = ', XuA
        write(310,*)'kaiB        = ', XuB
        write(310,*)'st          = ',st  
        write(310,*)'Number of A segments =          ', sumpolA
        write(310,*)'Number of B segments =          ', sumpolB
        do im = 1, N_monomerA
            write(310,*)'avfdisA(',im,') = ',avfdisA(im) 
        enddo 
        do im = 1, N_monomerB
            write(310,*)'avfdisB(',im,') = ',avfdisB(im) 
        enddo
        
        
        close(310)

    endif ! == if(rank==)  

end subroutine

subroutine store2disk(counter) ! saves state to disk

    use ellipsoid
    use kinsol
    use montecarlo
    use ematrix
    use results
    use MPI
    use const
    
    implicit none
    integer :: counter
    character*20 :: filename

    if(rank.eq.0) then
        open (unit=8, file='out.out', form='unformatted')
        write(8)counter
        write(8)free_energy
        write(8)xflag
        close(8)
    endif

    if(rank.eq.0) then
        write(filename,'(A4, I3.3, A4)')'out.', counter, '.dat'
        open(unit=8, file=filename, form='unformatted')
        write(8)counter
        write(8)free_energy
        write(8)xflag
        close(8)
    endif

end subroutine store2disk

subroutine retrivefromdisk(counter) ! saves state to disk

    use ellipsoid
    use kinsol
    use montecarlo
    use ematrix
    use results
    use const

    implicit none
    integer :: counter

    open (unit=8, file='in.in', form='unformatted')
    read(8)counter
    read(8)free_energy
    read(8)xflag
    close(8)

end subroutine retrivefromdisk 

subroutine mirror
    use const
    use kinsol
    use mparameters_monomer
    implicit none

    real*8 :: xh(dimx,dimy,dimz), psi(dimx,dimy,dimz)
    real*8 :: xtotalA(dimx,dimy,dimz,N_poorsolA)
    real*8 :: xtotalB(dimx,dimy,dimz,N_poorsolB)
    integer :: ip
    integer :: ix,iy,iz
    real*8 :: temp
    integer :: ncells


    ncells = dimx*dimy*dimz

    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz
                xh(ix,iy,iz)=xflag(ix+dimx*(iy-1)+dimx*dimy*(iz-1))

                do ip = 1, N_poorsolA
                    xtotalA(ix,iy,iz,ip) = xflag(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells)
                enddo

                do ip = 1, N_poorsolB
                    xtotalB(ix,iy,iz,ip) = xflag(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(ip+N_poorsolA)*ncells)
                enddo


                if(electroflag.eq.1) psi(ix,iy,iz)=xflag(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+&
                    (N_poorsolA+N_poorsolB+1)*ncells)
            enddo
        enddo
    enddo 
   
    do ix=1,int(dimx/2)
        do iy=1,dimy
            do iz=1,dimz
                temp = xh(ix,iy,iz)
                xh(ix,iy,iz) = xh(dimx-ix,iy,iz)
                xh(dimx-ix,iy,iz) = temp

                do ip = 1, N_poorsolA
                    temp = xtotalA(ix,iy,iz,ip)
                    xtotalA(ix,iy,iz,ip) = xtotalA(dimx-ix,iy,iz,ip)
                    xtotalA(dimx-ix,iy,iz,ip) = temp
                enddo

                do ip = 1, N_poorsolB
                    temp = xtotalB(ix,iy,iz,ip)
                    xtotalB(ix,iy,iz,ip) = xtotalB(dimx-ix,iy,iz,ip)
                    xtotalB(dimx-ix,iy,iz,ip) = temp
                enddo

                if(electroflag.eq.1) then
                    temp = psi(ix,iy,iz)
                    psi(ix,iy,iz) = psi(dimx-ix,iy,iz)
                    psi(dimx-ix,iy,iz) = temp
                endif
  
            enddo
        enddo
    enddo
    
  
    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz
                xflag(ix+dimx*(iy-1)+dimx*dimy*(iz-1))= xh(ix,iy,iz)

                do ip = 1, N_poorsolA
                    xflag(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells) =  xtotalA(ix,iy,iz,ip)
                enddo

                do ip = 1, N_poorsolB
                    xflag(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(ip+N_poorsolA)*ncells) =  xtotalB(ix,iy,iz,ip)
                enddo

                if(electroflag.eq.1) xflag(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsolA+N_poorsolB+1)*ncells)= &
                    psi(ix,iy,iz)
            enddo
        enddo
    enddo

end subroutine mirror





