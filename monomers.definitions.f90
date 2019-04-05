subroutine monomer_definitions

use MPI
use mparameters_monomer


implicit none
integer i,j

!read epsilon.in: N_poorsol and st_matrix
open(file='epsilon.in', unit=333)
read(333,*) N_poorsol
ALLOCATE (st_matrix(N_poorsol, N_poorsol)) ! interaction between monomer types in fraction of st, scaled by st-scale during running....
do i = 1, N_poorsol
 read(333,*)(st_matrix(i,j), j = 1, i)
 do j = 1, i
  st_matrix(j,i) = st_matrix(i,j)
 enddo
enddo
close(333)

!read monomer.in : N_monomer i zpol, hydroph pka
open(file='monomer.in', unit=333)
read(333,*) N_monomer
ALLOCATE (zpol(N_monomer))    ! charge of monomer segment: 1: base, -1: acid, 0:neutral
ALLOCATE (hydroph(N_monomer)) ! 0: hydrophilic, 1 < x < N_poorsol, type of poor solvent
ALLOCATE (pKa(N_monomer), Ka(N_monomer), K0(N_monomer))
do j =1, N_monomer
 read(333,*) i
 read(333,*) zpol(i)
 read(333,*) hydroph(i)
 read(333,*) pka(i)
enddo
close(333)

!implicit none
!integer i

!N_poorsol = 1 ! number of different kais
!N_monomer = 1

!ALLOCATE (st_matrix(N_poorsol, N_poorsol)) ! interaction between monomer types in fraction of st, scaled by st-scale during running....
!ALLOCATE (zpol(N_monomer))    ! charge of monomer segment: 1: base, -1: acid, 0:neutral
!ALLOCATE (hydroph(N_monomer)) ! 0: hydrophilic, 1 < x < N_poorsol, type of poor solvent
!ALLOCATE (pKa(N_monomer), Ka(N_monomer), K0(N_monomer))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!st_matrix(1,1)=1.0

!Segmento carga+1
!i=1
!zpol(i) = 1 !carga positiva
!hydroph(i)= 1 ! sv pobre
!pka(i)=8

! Segment type 1 for NPC, positive base, hydrophilic
!i = 1
!zpol(i) = 0
!hydroph(i) = 1
!pKa(i) = 10.0

!i = 2
!zpol(i) = 0
!hydroph(i) = 1
!pKa(i) = 11.0

!! Segment type 2 for NPC, negative , hydrophilic
!
!zpol(2) = -1
!hydroph(2) = 0
!pKa(2) = 5.0
!
!! Segment type 3 for NPC, neutral , hydrophilic
!
!zpol(3) = 0
!hydroph(3) = 0
!pKa(3) = 1 ! set any number if zpol = 0...
!
!! Segment type 4 for NPC , neutral, hydrophobic, 1
!
!zpol(4) = 0
!hydroph(4) = 1
!pKa(4) = 1 ! set any number if zpol = 0....
!
end

