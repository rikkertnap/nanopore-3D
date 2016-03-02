subroutine monomer_definitions

use MPI
use mparameters_monomer

implicit none


N_poorsol = 1 ! number of different kais
N_monomer = 4 

ALLOCATE (st_matrix(N_poorsol, N_poorsol)) ! interaction between monomer types in fraction of st, scaled by st-scale during running....
ALLOCATE (zpol(N_monomer))    ! charge of monomer segment: 1: base, -1: acid, 0:neutral
ALLOCATE (hydroph(N_monomer)) ! 0: hydrophilic, 1 < x < N_poorsol, type of poor solvent
ALLOCATE (pKa(N_monomer), Ka(N_monomer), K0(N_monomer))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

st_matrix(1,1)=1.0

! Segment type 1 for NPC, positive base, hydrophilic

zpol(1) = 1
hydroph(1) = 0
pKa(1) = 11.0

! Segment type 2 for NPC, negative , hydrophilic

zpol(2) = -1
hydroph(2) = 0
pKa(2) = 5.0

! Segment type 3 for NPC, neutral , hydrophilic

zpol(3) = 0
hydroph(3) = 0
pKa(3) = 1 ! set any number if zpol = 0...

! Segment type 4 for NPC , neutral, hydrophobic, 1

zpol(4) = 0
hydroph(4) = 1
pKa(4) = 1 ! set any number if zpol = 0....

end

