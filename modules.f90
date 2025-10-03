

module mkl

    type compressed_matrix

        real*8, allocatable :: inc_values(:)
        integer, allocatable :: inc_columns(:)
        integer, allocatable :: pntrb(:)
        integer, allocatable :: pntre(:)
        integer :: nonzeros
        integer :: gidx
        integer, allocatable :: localmap(:)
        integer, allocatable :: localmapr(:)
        integer :: mapsize
        real*8, allocatable :: lnxpot(:)
        real*8, allocatable :: avpol_tmp(:)

    endtype

    type compressed_prob
        real*8, allocatable :: pro(:)
        real*8, allocatable :: lnpro(:)
    endtype

    type(compressed_prob), allocatable :: promkl(:)

    type(compressed_matrix), allocatable :: csrm(:,:) ! csr matrix, first index is grafing point, second index is monomer type
    integer :: flagmkl
    integer :: sumnewcuantas

end module

module maps                                
    integer, allocatable :: imap(:,:,:)               ! == hash table form 3D to 1D 
    integer, allocatable :: mapx(:), mapy(:), mapz(:) ! == inverse 1D to 3D
endmodule
 

module mparameters_monomer

    integer :: N_poorsol                    ! number of different kais
    integer :: N_monomer                    ! number of different monomer types
    real*8, allocatable  :: st_matrix(:,:)  ! interaction between monomer types in fraction of st, scaled by st-scale during running....
    integer, allocatable :: zpol(:)         ! charge of monomer segment: 1: base, -1: acid, 0:neutral
    integer, allocatable :: hydroph(:)      ! 0: hydrophilic, 1 < x < N_poorsol, type of poor solvent
    real*8, allocatable  ::  pKa(:), Ka(:), K0(:)

endmodule mparameters_monomer


module branches
    integer longb(3), longbb
    integer branched
    integer indexncha
endmodule

module system 
    integer :: systemtype   ! == descriptor of system, see parser.f90
    integer :: vscan        ! == select type of loop of VdW variable ?? 
    real*8 :: delta         ! == unit of length volume cell  in nm ?? 
    real*8 :: dx,dy,dz      ! == 
    real*8 :: cdiva         ! == cdiva ???
    integer :: dimx         ! == number of volume cell  in x-directions
    integer :: dimy 
    integer :: dimz       
    integer :: PBC(6)       ! == Periodic bondary conditions 
    integer :: vtkflag      ! == if flag ==1 make vtk output formatted file ?? 
    integer :: electroflag  ! == if flag ==1 use electrostatics i.e. solve Poisson Eq
    integer :: eqs          ! == number of set of equations, total number of non-linear equatios eqs * (nsize = dimx *dimy *dimz)  
endmodule



module ematrix
    use system

    real*8, allocatable :: volprot(:,:,:)
    real*8, allocatable :: volprot1(:,:,:)
    real*8, allocatable :: voleps(:,:,:)
    real*8, allocatable :: voleps1(:,:,:)
    real*8, allocatable :: volq(:,:,:)
    real*8, allocatable :: volq1(:,:,:)
    integer, parameter :: maxvolx = 50000
    real*8 volx(maxvolx)
    real*8 com(maxvolx,3)
    integer p0(maxvolx,3)
    real*8, allocatable :: fvstd(:,:,:)
    real*8, allocatable :: fvmkl(:)

end module

module rotchain
    use ematrix, only : maxvolx
    real*8 :: rotangle(maxvolx)
endmodule

module channel

    real*8 :: rchannel                   ! == radius nanochannel
    real*8 :: originc(2)                  
    real*8 :: echargec, sigmac, eepsc, sigmar
    integer :: NBRUSH
    integer :: RdimZ                     ! size of reservoirs in delta units
    integer :: Nrings                    ! number of rings for systemtype = 42
    real*8, allocatable :: ringpos(:)    ! position along the pore
    integer :: Npolx, Npoly

endmodule

module s2d        
    integer scx,scy,scz                  ! == ranges in vtk file 
endmodule

module mkinsol
    double precision, allocatable :: pp(:)    ! == pre-condition variable in kinsol
endmodule

module montecarlo
    real*8 :: free_energy
endmodule

module chainsdat
    integer :: cuantas                    ! == number of conformations
    integer, allocatable :: newcuantas(:) ! == ???
    integer :: long                       ! == length of polymer chain /number of segments
    integer, allocatable :: segtype(:)    ! sequence of the chain 
    integer :: ncha 
    real*8, ALLOCATABLE :: in1(:,:)       ! segment positions 
    integer :: ing                        ! number of gauches in current chain
    real*8, ALLOCATABLE :: posicion(:,:)  ! posicion graft de la cadena ncha
    real*8, ALLOCATABLE :: ngpol(:)       ! posicion graft de la cadena ncha
    integer, ALLOCATABLE :: cpp(:)
    integer, ALLOCATABLE :: cppini(:)
    integer :: maxcpp                     ! == ??
    real*8 :: lseg                        ! == length segment 
    integer ::  readchains                ! == variable that select reading stored conformation 
endmodule

module molecules
    use system
    real*8 :: vsol                        ! == volume solvent 
    real*8 :: vpol                        ! == volume polymer segment 
    real*8 :: vpol0
    real*8 :: vsol0
    real*8 :: vsalt                       ! == volume salt ions
    real*8 :: zpos,zneg                   ! == valence ions  
    real*8 :: benergy                     ! == energy gauche bond ??
    real*8 :: fz                          ! == ?  
endmodule

module kaist
    integer hguess
    real*8 hring
    real*8 oval
    integer nkp
    real*8 kp
    real*8 kps(100)

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
    real*8, allocatable :: xtotal(:, :, :, :) ! xtotal para poor solvent
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
    real*4 , allocatable :: zfinal(:,:)
endmodule

module MPI
    !include 'mpif.h' ! librerias MPI
    use mpi_f08  
    integer :: rank       ! local rank node
    integer :: size       ! number of nodes, size of mpi size override intrinic function size 
    integer :: ierr       ! output flag  
    integer :: flagsolver ! continuation stop flag of  
endmodule

module kinsol
    use system
    integer :: iter                 ! == number fkfun evals /iterations           
    integer*4 :: ier                ! Kinsol error flag
    integer*8 :: neq                ! Kinsol number of equations
    real*8 :: norma                 ! l2norm of fkfun     
    real*8, ALLOCATABLE :: xflag(:) 
    real*8, ALLOCATABLE :: xpar(:) 
endmodule

module const
    real*8 :: dielW, dielP, dielS    ! == dielectric contant of water, polymer and surfce    
    real*8 :: constqE                ! == pre factor in Poisson Equation 
    real*8 :: dielPr, dielSr         ! == relative dielectric constant ?? 
    real*8 :: pKw, Kw
    real*8 :: pi 
    real*8, parameter :: Na = 6.02d23                    ! == Avogadro's number  
   ! real*8, parameter :: Na = 6.022140857e23_dp          ! == Avogadro's number unit 
    real*8 :: constq                 ! == pre factor in Poisson Eq 
    real*8 :: lb                     ! == Bjerrum length 
    integer :: seed                                 
    integer :: seed2
    real*8 :: error                  ! para comparar con la norma...
    real*8 :: errel
    integer :: itmax
    integer :: infile
    integer :: randominput
    integer :: epstype
    integer :: verbose     
    integer :: stdout             ! == unit number of write of stdout
endmodule

module kai
    integer :: Xulimit
    real*8 :: cutoff
    real*8, allocatable :: Xu(:,:,:)
    real*8 :: sumXu
endmodule

module results
    use system
    real*8, allocatable :: avpol(:,:,:,:)   ! == volume fraction polymer : indices ix iy iz im : im type monomer
    real*8, allocatable :: epsfcn(:,:,:)    ! == dielectric constant at ix iy iz   
    real*8, allocatable :: Depsfcn(:,:,:)   ! == derivative of dielectric constant ?? 
    real*8, allocatable :: xpos(:,:,:)      ! == volume fraction pos ion
    real*8, allocatable :: xneg(:,:,:)      ! == volumer fraction neg ion
    real*8, allocatable :: qtot(:,:,:)      ! == total charge density 
    real*8, allocatable :: xHplus(:,:,:)    ! == volume fraction of H+
    real*8, allocatable :: xOHmin(:,:,:)    ! == volume fraction of OH-
    real*8, allocatable :: fdis(:,:,:,:)    ! == degree of dissociation of momomer type im : indices ix, iy ,iz im 
endmodule

module bulk                                 ! = bulk chemical potentails and volume fractions
    real*8 :: expmupos,expmuneg,expmuHplus,expmuOHmin
    real*8 :: xsolbulk, xposbulk, xnegbulk, xHplusbulk,xOHminbulk
endmodule


module ellipsoid
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

module inputtemp     
    real*8 :: xsalt            ! == volume fraction salt in reservoir
    real*8 :: pHbulk           ! == pH  resevoir 
    real*8 :: pOHbulk          ! == pOH 
    real*8 :: csalt            ! == concentration of salt in Mol/L 
    real*8 :: cHplus, cOHmin   ! == concentration of H+ and OH- 
end module

module transform               ! == coordinate transformation 
    real*8 :: gama0            ! == angle between basis vector in oblique or prism coordinates 
    real*8 :: MAT(3,3)
    real*8 :: TMAT(3,3)
    real*8 :: IMAT(3,3)
endmodule

