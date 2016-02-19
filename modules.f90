!EDITAR:

!OK - chequear limites de psi y xtotal en fkfun
!OK - Poisson eq en fkfun
!OK - mapear particula en el nuevo sistema
!OK - check free energy

! 
! Crystal v1.0
! 
! BASED on MCPARv4.0
! + transformation in non-cubic lattices
!


module channel
real*8 rchannel
real*8 originc(2)
endmodule

module system 
real*8 delta
real*8 dx,dy,dz
real*8 cdiva
integer  dimx 
integer  dimy 
integer  dimz 
integer PBC(6)
integer vtkflag
integer electroflag
endmodule

module s2d
integer scx,scy,scz
endmodule

module mkinsol
double precision, allocatable :: pp(:)
endmodule

module montecarlo
real*8 free_energy
endmodule

module chainsdat
integer cuantas 
integer, allocatable :: newcuantas(:)
integer long 
integer ncha 
real*8, ALLOCATABLE :: in1(:,:)  ! segment positions 
integer ing ! number of gauches in current chain
real*8, ALLOCATABLE :: posicion(:,:) ! posicion graft de la cadena ncha
real*8, ALLOCATABLE :: ngpol(:) ! posicion graft de la cadena ncha
integer, ALLOCATABLE :: cpp(:)
integer, ALLOCATABLE :: cppini(:)
integer maxcpp
real*8 lseg
integer readchains
endmodule

module molecules
use system
real*8 vsol
real*8 vpol
real*8 vpol0
real*8 vsol0
real*8 vsalt
real*8 zpos,zneg, zpol
real*8 benergy
endmodule

module kaist
integer nst
real*8 st
real*8 sts(100)

integer nsc
real*8 sc
real*8 scs(100)

endmodule

module fields_fkfun
use system
use chainsdat
real*8, allocatable :: xtotal(:, :, :) ! xtotal para poor solvent
real*8, allocatable :: psi(:, :, :) 
real*8, allocatable :: q(:)
real*8, allocatable :: sumgauche(:)
real*8, allocatable :: pro(:,:)
real*8, allocatable :: xh(:, :, :)
real*8 shift
endmodule

module conformations
integer*1, allocatable :: px(:,:,:)
integer*1, allocatable :: py(:,:,:)
integer*1, allocatable :: pz(:,:,:)
integer*1, allocatable :: ngauche(:,:)
endmodule

module MPI
include 'mpif.h' ! librerias MPI
integer rank, size, ierr
integer flagsolver
endmodule

module kinsol
use system
integer iter
integer *4 ier ! Kinsol error flag
integer *8 neq ! Kinsol number of equations
real*8 norma
real*8, ALLOCATABLE :: xflag(:) 
real*8, ALLOCATABLE :: xpar(:) 
endmodule

module const
real*8 dielW, dielP, dielS
real*8 constqE
real*8 dielPr, dielSr
real*8 pKw
real*8 pi 
real*8, parameter :: Na = 6.02d23 
real*8 constq
real*8 lb
integer seed
real*8 error  ! para comparar con la norma...
real*8 errel
integer itmax
integer infile
integer randominput
integer verbose
endmodule

module kai
integer Xulimit
real*8 cutoff
real*8, allocatable :: Xu(:,:,:)
real*8 sumXu
endmodule

module results
use system
real*8, allocatable :: avpol(:,:,:)
real*8, allocatable :: epsfcn(:,:,:)
real*8, allocatable :: Depsfcn(:,:,:)
real*8, allocatable :: xpos(:,:,:) ! pos ion
real*8, allocatable :: xneg(:,:,:) ! neg ioni
real*8, allocatable :: qtot(:,:,:) ! Carga total
real*8, allocatable :: xHplus(:,:,:) ! H+
real*8, allocatable :: xOHmin(:,:,:) ! OH-
real*8, allocatable :: fdis(:,:,:)
endmodule

module bulk
real*8 expmupos,expmuneg,expmuHplus,expmuOHmin
real*8 K0
real*8 xsolbulk, xposbulk, xnegbulk, xHplusbulk,xOHminbulk
endmodule


module ellipsoid
integer systemtype
integer NNN
real*8, allocatable :: rotmatrix(:,:,:)
real*8, allocatable :: Aell(:,:)
real*8, allocatable :: AellS(:,:)
real*8, allocatable :: AellL(:,:)
real*8, allocatable :: AellX(:,:)
real*8, allocatable :: AAA(:,:,:)
real*8, allocatable :: AAAS(:,:,:)
real*8, allocatable :: AAAL(:,:,:)
real*8, allocatable :: AAAX(:,:,:)
real*8, allocatable :: Rell(:,:)
real*8, allocatable :: Rellf(:,:)
real*8, allocatable :: orient(:,:)
real*8, allocatable :: echarge(:)
real*8, allocatable :: sigma(:)
real*8, allocatable :: eeps(:)
end module

module ematrix
use system
real*8, allocatable :: volprot(:,:,:)
real*8, allocatable :: volprot1(:,:,:)
real*8, allocatable :: voleps(:,:,:)
real*8, allocatable :: voleps1(:,:,:)
real*8, allocatable :: volq(:,:,:)
real*8, allocatable :: volq1(:,:,:)

integer, parameter :: maxvolx = 10000
real*8 volx(maxvolx)
real*8 com(maxvolx,3)
integer p0(maxvolx,3)
end module

module inputtemp
real*8 xsalt
real*8 pHbulk
real*8 pOHbulk
real*8 pKa, Ka
real*8 csalt
real*8 cHplus, cOHmin
end module

module transform
real*8 gama0
real*8 MAT(3,3)
real*8 TMAT(3,3)
real*8 IMAT(3,3)
endmodule

