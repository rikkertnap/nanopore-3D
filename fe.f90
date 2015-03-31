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
implicit none

integer looped
real*8  q_tosend(ncha)
real*8  q0(ncha)
integer newcuantas0(ncha)
real*8 F_Mix_s, F_Mix_pos
real*8 F_Mix_neg, F_Mix_Hplus
real*8 Free_energy2, sumpi, sumrho, sumel, sumdiel, suma, mupol
real*8 temp
real*8 F_Mix_OHmin, F_Conf, F_Eq, F_vdW, F_eps, F_electro
real*8 pro0(cuantas, maxcpp)
 
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

real*8 gradpsi2
real*8 fv, fv2

integer counter
real*8 psiv(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!
!  Recupera pro(i) de todos los procesos para calculo de F
!

! Subordinados

q0 = 0
q_tosend = 0

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

       enddo
       enddo 

         do ii = 2, size ! loop sobre los procesadores restantes

        source = ii-1

        call MPI_RECV(pro0, cuantas*cpp(ii), &
        MPI_DOUBLE_PRECISION, source, tag, MPI_COMM_WORLD,stat, err)

       do jj = 1, cpp(rank+1)
!       print*, ii, jj, pro0(10,jj)
       iii = cppini(ii)+jj
       do i = 1, newcuantas0(iii)

         F_Conf = F_Conf + (pro0(i, jj)/q0(iii))*dlog((pro0(i, jj))/q0(iii))*ngpol(iii)

       enddo
       enddo

         enddo ! ii

       endif ! rank

      Free_Energy = Free_Energy + F_Conf

! 7. Chemical Equilibrium
      F_Eq = 0.0 
            

      do ix  = 1, dimx
      do iy  = 1, dimy
      do iz  = 1, dimz
      
      fv=(1.0-volprot(ix,iy,iz))

      F_Eq = F_Eq + fdis(ix, iy, iz)*dlog(fdis(ix, iy, iz)) &
      *avpol(ix,iy,iz)/vpol*fv

      F_Eq = F_Eq + (1.0-fdis(ix,iy,iz)) &
      *dlog(1.0-fdis(ix, iy,iz))*avpol(ix, iy,iz)/vpol*fv

      F_Eq = F_Eq + (1.0-fdis(ix, iy,iz))*dlog(K0)*avpol(ix, iy,iz)/vpol*fv

      F_Eq = F_Eq + (1.0-fdis(ix, iy,iz)) &
      *(-dlog(expmuHplus))*avpol(ix, iy, iz)/vpol*fv

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

            jx = mod(jx-1+5*dimx, dimx) + 1
            jy = mod(jy-1+5*dimy, dimy) + 1
 
     select case(PBC)
     case(0)  
        if(((iz+az).ge.1).and.((iz+az).le.dimz)) then 
        fv2 = (1.0-volprot(jx,jy,iz+az)) 
        F_vdW = F_vdW - 0.5000*delta**3*xtotal(ix, iy, iz) &
        *xtotal(jx, jy,iz+az)*Xu(ax, ay, az)*st*fv*fv2
       endif 
    case(1)
        jz = iz+az
        jz = mod(jz-1+5*dimz, dimz) + 1

        fv2 = (1.0-volprot(jx,jy,jz))
        F_vdW = F_vdW - 0.5000*delta**3*xtotal(ix, iy, iz) &
        *xtotal(jx, jy,jz)*Xu(ax, ay, az)*st*fv*fv2
     endselect

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


      Free_Energy = Free_Energy + F_electro

! 10. Pol-prot

      F_eps = 0.0 

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz
      fv=(1.0-volprot(ix,iy,iz))
      F_eps = F_eps - avpol(ix,iy,iz)*voleps(ix,iy,iz)*(delta**3)/vpol/vsol*fv
      enddo
      enddo
      enddo

      Free_Energy = Free_Energy + F_eps

      if (verbose.ge.1) then
      print*, 'Free_Energy_Calc: Free energy(1) = ', Free_energy
      endif

! minimal F

      Free_Energy2 = 0.0

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
         sumdiel = sumdiel + 0.5/constq*xtotal(ix,iy,iz)*gradpsi2*Depsfcn(ix,iy,iz)

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
      print*, 'Free_Energy_Calc: Free energy(2) = ', Free_energy2, sumdiel
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




