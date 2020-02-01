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

integer looped
real*8  q_tosend(ncha), sumgauche_tosend(ncha)
real*8  q0(ncha), sumgauche0(ncha)
integer newcuantas0(ncha)
real*8 F_Mix_s, F_Mix_pos
real*8 F_Mix_neg, F_Mix_Hplus
real*8 Free_energy2, sumpi, sumrho, sumel, sumdiel, suma, mupol
real*8 temp
real*8 F_Mix_OHmin, F_gauche, F_Conf, F_Eq, F_vdW, F_eps, F_electro
real*8 pro0(cuantas, maxcpp)
real*8 entropy(dimx,dimy,dimz)
character*5  title
real*8 xtotalsum(dimx,dimy,dimz)
 
! MPI
integer stat(MPI_STATUS_SIZE) 
integer source
integer dest
integer tag
parameter(tag = 0)
integer err

! Dummies
integer ix, iy, iz, i, ii, ax, ay, az, jj
integer jx, jy, jz,iii
integer im, ip, ipp
real*8 gradpsi2
real*8 fv, fv2

integer counter
real*8 psiv(3)

integer, external :: PBCSYMI, PBCREFI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!
!  Recupera pro(i) de todos los procesos para calculo de F
!

! Subordinados

entropy = 0.0
q0 = 0.0
q_tosend = 0.0
sumgauche_tosend = 0.0

if(flagmkl.eq.1) then
       pro = 0.0
       do jj = 1, cpp(rank+1)
       iii = cppini(rank+1)+jj
       do i = 1, newcuantas(iii)
        pro(i,jj) = promkl(iii)%pro(i)
       enddo ! i
       enddo ! jj
endif

if(rank.ne.0) then
       dest = 0
! Envia q

       do jj = 1, cpp(rank+1)
       iii = cppini(rank+1)+jj
       q_tosend(iii) = q(iii)
       enddo

        call MPI_REDUCE(q_tosend, q0, ncha, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

! newcuantas
        call MPI_REDUCE(newcuantas, newcuantas0, ncha, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, err)

! Envia pro

        CALL MPI_SEND(pro, cuantas*cpp(rank+1) , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD,err)

! sum gauche

       do jj = 1, cpp(rank+1)
       iii = cppini(rank+1)+jj
       sumgauche_tosend(iii) = 0.0
       do i = 1, newcuantas(iii)
       sumgauche_tosend(iii) = sumgauche_tosend(iii)+ ngauche(i,iii)*pro(i,jj)/q(iii)
       enddo ! i
       enddo ! jj

        call MPI_REDUCE(sumgauche_tosend, sumgauche0, ncha, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)


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

      F_Mix_OHmin = F_Mix_OHmin + xOHmin(ix, iy,iz)*(dlog(xOHmin(ix, iy, iz))-1.0-dlog(expmuOHmin))*fv

      F_Mix_OHmin = F_Mix_OHmin - xOHminbulk*(dlog(xOHminbulk)-1.0-dlog(expmuOHmin))*fv

      enddo
      enddo
      enddo
      F_Mix_OHmin = F_Mix_OHmin * delta**3/vsol
      Free_Energy = Free_Energy + F_Mix_OHmin

! 6. Entropia interna polimero

      F_Conf = 0.0

! Jefe

       if (rank.eq.0) then ! Igual tiene que serlo, ver arriba

       do jj = 1, cpp(rank+1)
       iii = jj
       q_tosend(iii) = q(iii)
       enddo

        call MPI_REDUCE(q_tosend, q0, ncha, &
        MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

        call MPI_REDUCE(newcuantas, newcuantas0, ncha, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, err)

       do jj = 1, cpp(rank+1)
       do i = 1, newcuantas0(jj)
       iii = jj
      
         F_Conf = F_Conf + (pro(i, jj)/q0(iii)) &
      *dlog((pro(i, jj))/q0(iii))*ngpol(iii)

       entropy(p0(iii,1),p0(iii,2),p0(iii,3)) =  - dlog(q0(iii)/shift) 
       enddo
       enddo 

         do ii = 2, size ! loop sobre los procesadores restantes

        source = ii-1

        call MPI_RECV(pro0, cuantas*cpp(ii), &
        MPI_DOUBLE_PRECISION, source, tag, MPI_COMM_WORLD,stat, err)


       do jj = 1, cpp(ii)
!       write(stdout,*) ii, jj, pro0(10,jj)
       iii = cppini(ii)+jj
       do i = 1, newcuantas0(iii)

         F_Conf = F_Conf + (pro0(i, jj)/q0(iii))*dlog((pro0(i, jj))/q0(iii))*ngpol(iii)

       entropy(p0(iii,1),p0(iii,2),p0(iii,3)) =  - dlog(q0(iii)/shift) 

       enddo
       enddo

         enddo ! ii

       endif ! rank

      Free_Energy = Free_Energy + F_Conf

if(rank.eq.0) then

!      title = 'entpy'
!      call savetodisk(entropy, title, looped)
 
      open (unit=8, file='entropy.out', form='unformatted')
      write(8)dimx,dimy,dimz
      write(8)entropy
      close(8)
endif

! 6.5 Energy of gauche bonds

      F_gauche = 0.0

! Jefe

       if (rank.eq.0) then ! Igual tiene que serlo, ver arriba

       do jj = 1, cpp(rank+1)
       iii = cppini(rank+1)+jj
       sumgauche_tosend(iii) = 0.0
       do i = 1, newcuantas(iii)
       sumgauche_tosend(iii) = sumgauche_tosend(iii)+ ngauche(i,iii)*pro(i,jj)/q(iii)
       enddo ! i
       enddo ! jj

        call MPI_REDUCE(sumgauche_tosend, sumgauche0, ncha, &
        MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

       do ii = 1, ncha
       F_gauche = F_gauche + sumgauche0(ii)*ngpol(ii)*benergy
       enddo  

       endif ! rank

      Free_Energy = Free_Energy + F_gauche

if(rank.eq.0) then
      title = 'entpy'
      call savetodisk(entropy, title, looped)
 
      open (unit=8, file='entropy.out', form='unformatted')
      write(8)dimx,dimy,dimz
      write(8)entropy
      close(8)
endif

     

! 7. Chemical Equilibrium
      F_Eq = 0.0 
            

      do ix  = 1, dimx
      do iy  = 1, dimy
      do iz  = 1, dimz

      do im = 1, N_monomer
      
      fv=(1.0-volprot(ix,iy,iz))

      if(zpol(im).ne.0) then

      F_Eq = F_Eq + fdis(ix,iy,iz,im)*dlog(fdis(ix,iy,iz,im)) &
      *avpol(ix,iy,iz,im)/vpol*fv

      F_Eq = F_Eq + (1.0-fdis(ix,iy,iz,im)) &
      *dlog(1.0-fdis(ix,iy,iz,im))*avpol(ix,iy,iz,im)/vpol*fv

      F_Eq = F_Eq + (1.0-fdis(ix,iy,iz,im))*dlog(K0(im))*avpol(ix,iy,iz,im)/vpol*fv

      select case (zpol(im))
      case (-1) ! acid
       F_Eq = F_Eq + (1.0-fdis(ix,iy,iz,im))*(-dlog(expmuHplus))*avpol(ix,iy,iz,im)/vpol*fv
      case (1) ! base
       F_Eq = F_Eq + (1.0-fdis(ix,iy,iz,im))*(-dlog(expmuOHmin))*avpol(ix,iy,iz,im)/vpol*fv
      end select

      endif ! zpol

      enddo ! im
   
      enddo
      enddo
      enddo

      F_eq = F_eq *delta**3/vsol

      Free_Energy = Free_Energy + F_Eq
! 8.vdW ! Ojo, los kai son negativos => atraccion

       F_vdW = 0.0

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz

      fv=(1.0-volprot(ix,iy,iz))

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
                fv2 = (1.0-volprot(jx,jy,jz)) 

            do ip = 1, N_poorsol
            do ipp = 1, N_poorsol
 
                F_vdW = F_vdW - 0.5000*delta**3*xtotal(ix,iy,iz,ip) &
        *xtotal(jx,jy,jz,ipp)*Xu(ax, ay, az)*st*st_matrix(ip,ipp)*fv*fv2/(vpol*vpol*vsol*vsol)
 
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

      Free_Energy = Free_Energy + F_vdW

! 9. Electrostatic ! OJO

      F_electro = 0.0    

      do ix  = 1, dimx
      do iy  = 1, dimy
      do iz  = 1, dimz

      F_electro = F_electro &
       + delta**3*psi(ix, iy, iz)*qtot(ix, iy, iz)/2.0/vsol

      enddo
      enddo
      enddo
  
      print*, F_electro

      Free_Energy = Free_Energy + F_electro

! 10. Pol-prot

      F_eps = 0.0 

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz
      fv=(1.0-volprot(ix,iy,iz))
      do im = 1, N_monomer
      F_eps = F_eps - avpol(ix,iy,iz,im)*voleps(ix,iy,iz)*(delta**3)/vpol/vsol*fv
      enddo
      enddo
      enddo
      enddo

      Free_Energy = Free_Energy + F_eps

      if (verbose.ge.1) then
      write(stdout,*) 'Free_Energy_Calc: Free energy(1) = ', Free_energy
      endif

! minimal F

      Free_Energy2 = 0.0

      xtotalsum = 0.0
      do ip = 1, N_poorsol
      xtotalsum(:,:,:)= xtotalsum(:,:,:)+xtotal(:,:,:,ip)
      enddo


        sumpi = 0.0
        sumrho=0.0
        sumel=0.0
        sumdiel = 0.0

        do ix=1,dimx
        do iy=1,dimy
        do iz=1,dimz

      fv=(1.0-volprot(ix,iy,iz))

           sumpi = sumpi+dlog(xh(ix, iy, iz))*fv     
           sumpi = sumpi-dlog(xsolbulk)*fv
     
           sumrho = sumrho + ( - xh(ix, iy, iz) -xHplus(ix, iy, iz) &
        - xOHmin(ix, iy, iz) - (xpos(ix, iy, iz)+xneg(ix, iy, iz))/vsalt)*fv! sum over  rho_i i=+,-,s


           sumrho = sumrho - ( - xsolbulk -xHplusbulk &
       -xOHminbulk - (xposbulk+xnegbulk)/vsalt)*fv ! sum over  rho_i i=+,-,s

         sumel = sumel - qtot(ix, iy, iz)*psi(ix, iy, iz)/2.0 
      
         sumel = sumel + volq(ix,iy,iz)*psi(ix,iy,iz)*vsol                   


         psiv(1) = psi(ix+1,iy,iz)-psi(ix,iy,iz)
         psiv(2) = psi(ix,iy+1,iz)-psi(ix,iy,iz)
         psiv(3) = psi(ix,iy,iz+1)-psi(ix,iy,iz)

         gradpsi2 = DOT_PRODUCT(MATMUL(TMAT, psiv), MATMUL(TMAT, psiv))
         sumdiel = sumdiel + 0.5/constq*xtotalsum(ix,iy,iz)*gradpsi2*Depsfcn(ix,iy,iz)

         enddo
         enddo
         enddo
         
         sumpi = (delta**3/vsol)*sumpi
         sumrho = (delta**3/vsol)*sumrho
         sumel = (delta**3/vsol)*sumel
         sumdiel = (delta**3/vsol)*sumdiel


         suma = sumpi + sumrho + sumel + sumdiel


         do ii = 1, ncha
         Free_Energy2 = Free_Energy2-dlog(q0(ii)/shift)*ngpol(ii) 
         enddo

         Free_Energy2 = Free_Energy2 + suma - F_vdW

      if (verbose.ge.1) then
      write(stdout,*) 'Free_Energy_Calc: Free energy(2) = ', Free_energy2, sumdiel
      endif

! Guarda energia libre


        mupol = 0.0
        do ii = 1, ncha
        mupol = mupol - dlog(q0(ii)/shift)*ngpol(ii)
        enddo

        temp = sum(ngpol)
        mupol = mupol/temp

          if(rank.eq.0) then

         write(301,*)looped, Free_energy
          flush(301)
         write(302,*)looped, F_Mix_s 
         write(303,*)looped, F_Mix_pos
         write(304,*)looped, F_Mix_neg
         write(305,*)looped, F_Mix_Hplus
         write(306,*)looped, F_Mix_OHmin
	 write(3071,*)looped, F_gauche
         write(307,*)looped, F_Conf
         write(308,*)looped, F_Eq
         write(309,*)looped, F_vdW
         write(410,*)looped, F_eps
         write(311,*)looped, F_electro

         write(312,*)looped, Free_energy2

         write(313,*)looped, mupol

         endif
 
 888     call MPI_BCAST(free_energy, 1, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, err)

         return

         end




