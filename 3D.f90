subroutine solve(flagcrash)

    use system
    use const
    use kai
    use chainsdat
    use molecules
    use results
    use kinsol
    ! use bulk
    use MPI
    use ellipsoid
    use ematrix
    use mparameters_monomer
    use bulk, only : xsolbulk, xposbulk, xnegbulk, xHplusbulk, xOHminbulk
    use inputtemp, only : psizmax, psizmin 
    use flux, only : linear_interpolation
    use maps, only : mapx, mapy, mapz 

    implicit none

    integer, intent(inout) :: flagcrash
    
    ! ==  local variables 

    integer :: i, ix, iy, iz, ip, k, idx 
    integer :: noffset(4) 
    real*8 :: temp, psitemp

    ! == variable to resolve 

    real*8 :: x1(eqs*dimx*dimy*dimz),xg1(eqs*dimx*dimy*dimz)
    real*8 :: f(eqs*dimx*dimy*dimz) 
    integer :: ncells

    ! == Volume fraction ! local copies
    real*8 :: xh(dimx, dimy, dimz)
    real*8 :: xtotal(dimx,dimy,dimz,N_poorsol)
    real*8 :: psi(dimx, dimy, dimz) ! potential
    real*8 :: psiz( dimz)

    ! == MPI
    integer :: tag, source
    parameter(tag = 0)
    integer :: err
    integer :: ier_tosend
    double  precision :: norma_tosend

    ! == number of equations

    ncells = dimx*dimy*dimz

    ! == stride in x vector storage 

    do k=1,4     
        noffset(k) =(N_poorsol+2)*ncells + (k-1) * ncells
    enddo 

    ! == Initial guess

    if((infile.eq.2).or.(infile.eq.-1).or.(infile.eq.3)) then
        do i = 1, eqs*ncells  
            xg1(i) = xflag(i)     
            x1(i) = xflag(i)
        enddo
    endif

    if(infile.eq.0) then
    
        do i=1,ncells
            xg1(i)=xsolbulk
            x1(i)=xsolbulk
        enddo

        do i = ncells+1,(N_poorsol+1)*ncells
            xg1(i) = 0.0d0
            x1(i) = 0.0d0
        enddo

        if(electroflag.eq.1) then 
        
            if(fluxflag.eq.0) then
         
                do i=(N_poorsol+1)*ncells+1, (N_poorsol+2)*ncells 
                    xg1(i) = 0.0d0   ! potential 
                    x1(i) = 0.0d0
                enddo

            else if(fluxflag.eq.1) then 

                call linear_interpolation(psiz,psizmin,psizmin)
                ! potential 
                do iz=1,dimz
                    psitemp=psiz(iz)
                    do iy=1,dimy
                        do ix=1,dimx
                            idx = ix + dimx*(iy-1) + dimx*dimy*(iz-1) +(N_poorsol+1)*ncells
                            xg1(idx) = psitemp
                            x1(idx)  = psitemp
                        enddo
                    enddo
                enddo

                ! == ion densities 
                do idx=1, ncells

                    ix = mapx(idx)
                    iy = mapy(idx)
                    iz = mapz(idx)

                    if(fvstdint(ix,iy,iz).eq.1) then  
                        xg1(idx+noffset(1)) = xposbulk
                        x1(idx+noffset(1)) = xposbulk
                        xg1(idx+noffset(2)) = xnegbulk
                        x1(idx+noffset(2)) = xnegbulk
                        xg1(idx+noffset(3)) = xHplusbulk
                        x1(idx+noffset(3)) = xHplusbulk
                        xg1(idx+noffset(4)) = xOHminbulk
                        x1(idx+noffset(4)) = xOHminbulk
                    else                            ! not accesible for ions inside membrane 
                        xg1(idx+noffset(1)) = 0.0d0
                        x1(idx+noffset(1)) = 0.0d0 
                        xg1(idx+noffset(2)) = 0.0d0
                        x1(idx+noffset(2)) = 0.0d0
                        xg1(idx+noffset(3)) = 0.0d0
                        x1(idx+noffset(3)) = 0.0d0
                        xg1(idx+noffset(4)) = 0.0d0
                        x1(idx+noffset(4)) = 0.0d0
                    endif    

                enddo   
            endif

        endif ! electroflag     

    endif ! infile

    !--------------------------------------------------------------
    ! Solve               
    !--------------------------------------------------------------

    ! == head node
    if(rank.eq.0) then ! Only the head node calls the solver
        iter = 0
        write(stdout,*) 'solve: Enter solver ', eqs*ncells, ' eqs'

        if(infile.ge.0) then
            call solver(x1, xg1, ier)
        endif

        if(infile.eq.-1) then
            call fkfun(x1, f, ier)
        endif
        flagsolver = 0
        call MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
    endif
  
    ! == compute nodes

    if(rank.ne.0) then
        do
            flagsolver = 0
            source = 0
            CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
            if(flagsolver.eq.1) then
                call call_fkfun(x1) ! == there is still no solution  => fkfun 
            endif ! flagsolver
            if(flagsolver.eq.0) exit ! == stop program for these nodes
        enddo
    endif

    ! == Recover the value of ier and the norm!
    ! == This way, the compute nodes or subordinates find out if the solver converged or 
    ! == if the strategy needs to be changed.

    ! == head node

    if (rank.eq.0) then
        norma_tosend = norma
        ier_tosend = ier ! == different type of integer
        CALL MPI_BCAST(norma_tosend, 1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
        CALL MPI_BCAST(ier_tosend,1, MPI_INTEGER,0,MPI_COMM_WORLD,err)
    endif

    ! == compute nodes

    if (rank.ne.0) then
        CALL MPI_BCAST(norma_tosend, 1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
        CALL MPI_BCAST(ier_tosend, 1, MPI_INTEGER,0,MPI_COMM_WORLD,err)
        norma = norma_tosend
        ier = ier_tosend
    endif

    ! == recover get  xh and psi not common on all nodes
    do ix=1,dimx
        do iy=1,dimy
            do iz=1,dimz
                idx = ix+dimx*(iy-1)+dimx*dimy*(iz-1)
                xh(ix,iy,iz) = x1(idx)

                do ip=1, N_poorsol
                    xtotal(ix,iy,iz,ip)=x1(idx + ip*ncells)
                enddo

                if(electroflag.eq.1) then  
                    if(fluxflag.eq.0) then 
                        psi(ix,iy,iz)=x1(idx + (N_poorsol+1)*ncells)
                    endif
                    if(fluxflag.eq.1) then 
                        psi(ix,iy,iz) = x1(idx +(N_poorsol+1)*ncells)
                      !  xpos(ix,iy,iz) = x1( idx + noffset(1))
                      !  xneg(ix,iy,iz) = x1( idx + noffset(2))
                      !  xHplus(ix,iy,iz) = x1( idx + noffset(3))
                      !  xOHmin(ix,iy,iz) = x1( idx + noffset(4))
                    endif
                endif
            enddo
        enddo  
    enddo

    ! == Check solution  => set value flag 

    if(infile.ne.-1) then
        if((ier.lt.0).or.(.not.((norma.gt.0).or.(norma.lt.0))).or.(norma.gt.error)) then 
            if(rank.eq.0)write(stdout,*) 'solve: Error in solver: ', ier
            if(rank.eq.0)write(stdout,*) 'solve: norma ', norma
            flagcrash = 1
            return
        endif
    endif    

    ! == save xflag
    ! == xflag serves as input for the next iteration

    do i = 1, eqs*ncells
        xflag(i) = x1(i) 
    enddo
    
    infile = 2 ! == It does not read infile again.

    flagcrash = 0

end subroutine solve



