subroutine allocation

    use system
    use fields_fkfun
    use conformations
    use chainsdat
    use kinsol
    use results
    use ematrix
    use mkinsol
    use ellipsoid
    use MPI
    use kai
    use mparameters_monomer

    implicit none

    ! fields_fkfun
    !ALLOCATE(xtotal(1-Xulimit:dimx+Xulimit, 1-Xulimit:dimy+Xulimit, 1-Xulimit:dimz+Xulimit)) ! xtotal para poor solvent
    ALLOCATE(xtotalA(dimx, dimy, dimz, N_poorsolA))
    ALLOCATE(xtotalB(dimx, dimy, dimz, N_poorsolB))
    ALLOCATE(psi(0:dimx+1, 0:dimy+1, 0:dimz+1))
    ALLOCATE(xh(dimx, dimy, dimz))

    ! kinsol
    ALLOCATE (xflag(eqs*dimx*dimy*dimz))
    ALLOCATE (xpar(dimx*dimy*dimz))

    ! results
    ALLOCATE (avpolA(dimx, dimy, dimz, N_monomerA))   ! == polymer volume fraction
    ALLOCATE (avpolB(dimx, dimy, dimz, N_monomerB))   ! == polymer volume fraction
    ALLOCATE (xpos(dimx, dimy, dimz))               ! pos ion
    ALLOCATE (xneg(dimx, dimy, dimz))               ! neg ioni
    ALLOCATE (qtot(dimx, dimy, dimz))               ! total charge density 
    ALLOCATE (xHplus(dimx, dimy, dimz))             ! H+
    ALLOCATE (xOHmin(dimx, dimy, dimz))             ! OH-
    ALLOCATE (fdisA(dimx, dimy, dimz, N_monomerA))
    ALLOCATE (fdisB(dimx, dimy, dimz, N_monomerB))
    ALLOCATE (epsfcn(0:dimx+1, 0:dimy+1, 0:dimz+1))
    ALLOCATE (Depsfcn(0:dimx+1, 0:dimy+1, 0:dimz+1))

    ! ematrix
    ALLOCATE (volprot(dimx,dimy,dimz))
    ALLOCATE (volprot1(dimx,dimy,dimz))
    ALLOCATE (voleps(dimx,dimy,dimz))
    ALLOCATE (voleps1(dimx,dimy,dimz))

    ALLOCATE (volepsA(dimx,dimy,dimz))
    ALLOCATE (volepsA1(dimx,dimy,dimz))
    ALLOCATE (volepsB(dimx,dimy,dimz))
    ALLOCATE (volepsB1(dimx,dimy,dimz))


    ALLOCATE (volq(dimx,dimy,dimz))
    ALLOCATE (volq1(dimx,dimy,dimz))
    ALLOCATE (fvstd(dimx,dimy,dimz))
    ALLOCATE (fvmkl(dimx*dimy*dimz))
    ! mkinsol
    ALLOCATE (pp(eqs*dimx*dimy*dimz))

    ! chainsdat
    allocate(inA1(nsegA,3))
    allocate(inB1(nsegB,3))
    allocate(cpp(size))
    allocate(cppini(size))



end subroutine
