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
implicit none

integer*4 ier2
integer ncells
real*8 x(*),f(*)
real*8 protemp
integer i,j, ix, iy, iz, ii, ax, ay, az
integer im, ip
integer jx, jy, jz, jj
real*8 xpot(dimx, dimy, dimz, N_monomer)
! Charge
real*8 psitemp
real*8 MV(3),MU(3),MW(3)
real*8 MVV,MUU,MWW,MVU,MVW,MUW
real*8 psivv,psiuu,psiww, psivu,psivw,psiuw
real*8 psiv(3), epsv(3)
real*8 xtotalsum(dimx,dimy,dimz)

integer, external :: PBCSYMI, PBCREFI

! poor solvent 
real*8 sttemp
! MPI
integer tag
parameter(tag = 0)
integer err
real*8 avpol_tosend(dimx,dimy,dimz,N_monomer)
real*8 avpol_temp(dimx,dimy,dimz,N_monomer)
real*8 q_tosend, sumgauche_tosend
real*8 gradpsi2
real*8 fv

! hamiltonian inception
real*8 hfactor, hd
real*8, allocatable :: hds(:)
ALLOCATE(hds(100))
hds = -1

!-----------------------------------------------------
! Common variables

shift = 1.0d-100

ncells = dimx*dimy*dimz ! numero de celdas

! Jefe

if(rank.eq.0) then ! llama a subordinados y pasa vector x
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

psi = 0.0
do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
     xh(ix,iy,iz)=x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) !fraccion solvente

     do ip = 1, N_poorsol
      xtotal(ix,iy,iz,ip) = x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ ip*ncells) !fraccion polimero de tipo ip
     enddo
     if(electroflag.eq.1)psi(ix,iy,iz)=x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsol+1)*ncells)   !potencial electrostatico

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
   psi(0,:,:) = 0.0 
case(2)
   psi(0,:,:) = psi(1,:,:) ! zero charge
endselect

select case (PBC(2)) ! x = dimx
case(0) ! set bulk 
   psi(dimx+1,:,:) = 0.0  
case(2)
   psi(dimx+1,:,:) = psi(dimx,:,:) ! zero charge
endselect

select case (PBC(3)) ! y = 0
case(0) ! set bulk 
   psi(:,0,:) = 0.0  
case(2)
   psi(:,0,:) = psi(:,1,:) ! zero charge
endselect

select case (PBC(4)) ! y = dimy
case(0) ! set bulk 
   psi(:,dimy+1,:) = 0.0
case(2)
   psi(:,dimy+1,:) = psi(:,dimy,:) ! zero charge
endselect

select case (PBC(5)) ! z = 0
case(0) ! set bulk 
   psi(:,:,0) = 0.0  
case(2)
   psi(:,:,0) = psi(:,:,1) ! zero charge
endselect

select case (PBC(6)) ! z = dimz
case(0) ! set bulk 
   psi(:,:,dimz+1) = 0.0
case(2)
   psi(:,:,dimz+1) = psi(:,:,dimz) ! zero charge
endselect

! volume fraction and frdir

fdis = 0.0
avpol = 0.0

do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
    xpos(ix, iy, iz) = expmupos*(xh(ix, iy, iz)**vsalt)*dexp(-psi(ix, iy, iz)*zpos) ! ion plus volume fraction vsalt=vsal/vsv
    xneg(ix, iy, iz) = expmuneg*(xh(ix, iy, iz)**vsalt)*dexp(-psi(ix, iy, iz)*zneg) ! ion neg volume fraction
    xHplus(ix, iy, iz) = expmuHplus*(xh(ix, iy, iz))*dexp(-psi(ix, iy, iz))           ! H+ volume fraction
    xOHmin(ix, iy,iz) = expmuOHmin*(xh(ix,iy,iz))*dexp(+psi(ix,iy,iz))           ! OH-  volume fraction

     do im =1,N_monomer
        if (zpol(im).eq.1) then !BASE
          fdis(ix,iy,iz,im) = 1.0 /(1.0 + xOHmin(ix,iy,iz)/(K0(im)*xh(ix,iy,iz))) !k0 k en fraccion de volumen
        else if (zpol(im).eq.-1) then !ACID
          fdis(ix,iy,iz,im) = 1.0 /(1.0 + xHplus(ix,iy,iz)/(K0(im)*xh(ix,iy,iz)))
        endif
     enddo

   enddo
 enddo  
enddo



! Compute dielectric permitivity

xtotalsum = 0.0 ! sum of all polymers
do ip = 1, N_poorsol
xtotalsum(:,:,:) = xtotalsum(:,:,:) + xtotal(:,:,:,ip)
enddo
 
call dielectfcn(xtotalsum,volprot,epsfcn,Depsfcn)

!------------------------------------------------------------------------
! PDFs polimero
!------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! PARALELO: Cada procesador trabaja sobre una cadena...
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
       hds(i) = (float(2*ix-dimx)-2*cos(i*2*pi/hguess)*hring/delta)**2+(float(2*iy-dimy)-2*sin(i*2*pi/hguess)*hring/delta)**2
       hds(i) = hds(i)/4.0*(delta**2)+(oval*float(2*iz-dimz)/2.0*delta)**2
      end do
      hd = minval(hds, mask = hds .gt.0)
      hfactor = dexp(-(kp**2)*hd)

     end if
!xpot=exp(-Uj(rj)) para P(alpha)

     fv = (1.0 - volprot(ix,iy,iz)) !fraccion de volumen de la celda que es sc volprot->fraccion pared
     xpot(ix, iy, iz, im) = xh(ix,iy,iz)**vpol ! im:tipo de segmento, término de presion osmotica
     xpot(ix, iy, iz, im) = xpot(ix,iy,iz, im)*dexp(voleps(ix,iy,iz))  !termino de interaccion con sup de la particula

 ! Electrostatics

     if(zpol(im).ne.0.0) then
         xpot(ix,iy,iz,im) =  xpot(ix,iy,iz,im)/fdis(ix,iy,iz,im)*dexp(-psi(ix,iy,iz)*zpol(im))  !fdis: por eq ac. base...  
     endif
  
 ! Dielectrics

     gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz))**2 

!     gradpsi2 = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy-1,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz-1))**2 
!     xpot(ix, iy, iz) = xpot(ix,iy,iz)*exp(-Depsfcn(ix,iy,iz)*(gradpsi2)*constqE)

     xpot(ix,iy,iz,im) = xpot(ix,iy,iz,im)*exp(Depsfcn(ix,iy,iz)*(gradpsi2)/constq/2.0*vpol/fv)

 ! Poor solvent depende de la grilla donde esta y de sus vecinos

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
               protemp=protemp + hfactor*Xu(ax,ay,az)*st_matrix(hydroph(im),ip)*sttemp*xtotal(jx,jy,jz,ip)*fv
               enddo ! ip

            endif
            endif
            endif

       enddo
      enddo
     enddo

     xpot(ix,iy,iz,im) = xpot(ix,iy,iz,im)*dexp(protemp)

     endif ! hydrph

   enddo ! ix
  enddo ! iy
enddo !iz

enddo ! N_monomer

avpol_tosend = 0.0
q = 0.0
sumgauche = 0.0

do jj = 1, cpp(rank+1)    !punto de anclaje
   ii = cppini(rank+1)+jj

   q_tosend=0.0
   sumgauche_tosend = 0.0
   avpol_temp = 0.0

 do i=1,newcuantas(ii)  !conformaciones 
   pro(i, jj)=shift
   do j=1,long   ! segmentos
    ax = px(i, j, jj) ! cada uno para su cadena...
    ay = py(i, j, jj)
    az = pz(i, j, jj)         
    pro(i, jj) = pro(i, jj) * xpot(ax, ay, az, segtype(j))
    pro(i, jj) = pro(i, jj) * dexp(-fz*zfinal(i,j,jj)*zpol(segtype(j)))  ! termino Fz 
   enddo
    pro(i,jj) = pro(i,jj)*exp(-benergy*ngauche(i,ii)) ! energy of gauche bonds P(alpha).q-> sin normalizar

   do j=1,long
   fv = (1.0-volprot(px(i,j, jj),py(i,j, jj),pz(i,j, jj)))
   im = segtype(j)
    avpol_temp(px(i,j, jj),py(i,j, jj),pz(i,j, jj), im)= &
    avpol_temp(px(i,j, jj),py(i,j, jj),pz(i,j, jj), im)+pro(i, jj)*vpol*vsol/(delta**3)/fv* &
    ngpol(ii)*sc !ngpol(ii) has the number of chains grafted to the point ii. fraccion de vol de pol sin norm
   enddo

   q_tosend=q_tosend+pro(i, jj)
   sumgauche_tosend = sumgauche_tosend+ngauche(i, ii)*pro(i,jj)

 enddo ! i
! norma 
    
avpol_tosend=avpol_tosend + avpol_temp/q_tosend  !normalizar por q

q(ii) = q_tosend ! no la envia ahora
sumgauche(ii) = sumgauche_tosend/q_tosend

!write(stdout,*) rank+1,jj,ii,q(ii)
enddo ! jj

!------------------ MPI ----------------------------------------------
!1. Todos al jefe


call MPI_Barrier(MPI_COMM_WORLD, err)

! Jefe
if (rank.eq.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol, ncells*N_monomer, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
endif
! Subordinados
if(rank.ne.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol, ncells*N_monomer, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err) 
!!!!!!!!!!! IMPORTANTE, LOS SUBORDINADOS TERMINAN ACA... SINO VER !MPI_allreduce!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  goto 3333
endif

!!!!!!!!!!!!!!!!!!!!!!! FIN MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------------------------------
!   Construye Ecuaciones a resolver 
!----------------------------------------------------------------------------------------------

! Qtot


qtot = 0.0

do ix=1,dimx
do iy=1,dimy
do iz=1,dimz
  
 fv = (1.0-volprot(ix,iy,iz))

 qtot(ix, iy, iz) =  (zpos*xpos(ix, iy, iz)+zneg*xneg(ix, iy, iz))/vsalt + xHplus(ix, iy, iz) - xOHmin(ix, iy, iz)

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
  f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) = f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)) + avpol(ix,iy,iz,im) !packing ...+polimero
 enddo

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
    f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells) = f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ip*ncells) - avpol(ix,iy,iz,im)
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

psivu = (psi(ix+1,iy+1,iz)+psi(ix-1,iy-1,iz)-psi(ix+1,iy-1,iz)-psi(ix-1,iy+1,iz))/4.0
psivw = (psi(ix+1,iy,iz+1)+psi(ix-1,iy,iz-1)-psi(ix+1,iy,iz-1)-psi(ix-1,iy,iz+1))/4.0
psiuw = (psi(ix,iy+1,iz+1)+psi(ix,iy-1,iz-1)-psi(ix,iy+1,iz-1)-psi(ix,iy-1,iz+1))/4.0

psiv(1) = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))/2.0
psiv(2) = (psi(ix,iy+1,iz)-psi(ix,iy-1,iz))/2.0
psiv(3) = (psi(ix,iy,iz+1)-psi(ix,iy,iz-1))/2.0

epsv(1) = (epsfcn(ix+1,iy,iz)-epsfcn(ix-1,iy,iz))/2.0
epsv(2) = (epsfcn(ix,iy+1,iz)-epsfcn(ix,iy-1,iz))/2.0
epsv(3) = (epsfcn(ix,iy,iz+1)-epsfcn(ix,iy,iz-1))/2.0

psitemp = epsfcn(ix,iy,iz)*(MVV*psivv+MUU*psiuu+MWW*psiww+2.0*MVU*psivu+2.0*MVW*psivw+2.0*MUW*psiuw)
psitemp = psitemp + DOT_PRODUCT(MATMUL(TMAT,epsv),MATMUL(TMAT,psiv))

! OJO CHECK!!!!

f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+(N_poorsol+1)*ncells)=(psitemp + qtot(ix, iy, iz)*constq)/(-2.0)

enddo
enddo
enddo

endif ! electroflag
 
norma = 0.0

do i = 1, eqs*ncells
  norma = norma + (f(i))**2
enddo

iter = iter + 1
if(verbose.ge.3) then
if(rank.eq.0)write(stdout,*)'fkfun:', iter, norma, q(1)
endif

3333 continue
ier2 = 0.0 

return
end
