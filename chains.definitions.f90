subroutine chains_definitions

    use chainsdat, only : nsegA, nsegB, segtypeA, segtypeB
    use MPI, only : ierr
    ! use branches
    use const, only : stdout

    implicit none
    
    integer :: i, l
 
    ALLOCATE (segtypeA(nsegA)) 
    l=0
    open(file="sequenceA.in",unit=333)
    do i=1,nsegA
        read(333,*) segtypeA(i)
         l=l+1
    enddo

    if(l.ne.nsegA) then
        write(stdout,*) "error in sequenceA.in"
        call MPI_FINALIZE(ierr) ! finaliza MPI
        stop
    endif
    close(333)

    ALLOCATE (segtypeB(nsegB)) 
    l=0
    open(file="sequenceB.in",unit=333)
    do i=1,nsegB
        read(333,*) segtypeB(i)
         l=l+1
    enddo

    if(l.ne.nsegB) then
        write(stdout,*) "error in sequenceB.in"
        call MPI_FINALIZE(ierr) ! finaliza MPI
        stop
    endif
    close(333)




    !if(branched.eq.1) then
    !   segtype(1:longbb+longb(1)) = 2 ! backbone and first branch is hydrophobic
    !endif

    !if(branched.ne.1) segtype = 2

end subroutine chains_definitions




