! main program 
! ## Montecarlo - Molecular theory for the adsorption of a particle on a brush


program main

!   == variable and constant declaractions and definition 

    use system
    use MPI
    use ellipsoid
    use kinsol
    use const
    use montecarlo
    use ematrix
    use kaist
    use mkl

    implicit none

    integer :: counter, counterr
    integer :: MCpoints
    integer :: saveevery 
    real*8  :: maxmove
    real*8  :: maxrot
    integer :: moves,rots
    real*8  :: rv(3)
    real*8  :: temp
    real*8  :: theta
    real*8, external :: rands
    logical :: flag              
    character*10 :: filename
    integer :: j, i, ii, iii
    integer :: flagcrash
    real*8  ::  stOK,kpOK
    real*8  :: time0, timeF

    stdout = 6                           ! == unit number defined modules          

    !!!!!!!!!!!! global parameters ...  !!!!!!

    saveevery = 1

    counter = 0
    counterr = 1

    call initmpi

    if(rank.eq.0) write(stdout,*) 'Program Crystal'
    if(rank.eq.0) write(stdout,*) 'GIT Version: ', _VERSION
    if(rank.eq.0) write(stdout,*) 'MPI OK'

    if(rank.eq.0) write(10,*) 'Program Crystal'
    if(rank.eq.0) write(10,*) 'GIT Version: ', _VERSION
    if(rank.eq.0) write(10,*) 'MPI OK'
    
    call readinput        ! == DEFINITIONS.txt file 

    call monomer_definitions
    call chains_definitions
    call makemaps

    call initconst
    call inittransf ! Create transformation matrixes
    call initellpos ! calculate real positions for ellipsoid centers
    call initall
    call allocation

    !!! General files

    if(systemtype.eq.1) then
        do j = 1, NNN
            write(filename,'(A3,I3.3, A4)')'pos',j,'.dat'
            open(file=filename, unit=5000+j)
            write(filename,'(A3,I3.3, A4)')'orn',j,'.dat'
            open(file=filename, unit=6000+j)
            write(filename,'(A3,I3.3, A4)')'rot',j,'.dat'
            open(file=filename, unit=7000+j)
        enddo
    endif

    open(file='free_energy.dat', unit=9000)

    call kais                               ! == Van der Waals chi=kai 
    
    if(rank.eq.0)write(stdout,*) 'Kai OK'

    ! == select system 

    if (systemtype.eq.1) then
        call update_matrix(flag)            ! updates 'the matrix'
    elseif (systemtype.eq.2) then
        call update_matrix_channel(flag)    ! == channel 
    elseif (systemtype.eq.3) then
        call update_matrix_channel_3(flag)  ! updates 'the matrix'
    elseif (systemtype.eq.4) then
        call update_matrix_channel_4(flag)  ! updates 'the matrix'
    elseif (systemtype.eq.41) then
        call update_matrix_channel_4(flag)  ! 
    elseif (systemtype.eq.42) then
        call update_matrix_channel_4(flag)  ! == channel
    elseif (systemtype.eq.52) then
        call update_matrix_channel_4(flag)  ! == rod
    elseif (systemtype.eq.6) then
        call update_matrix_planar(flag)     ! == planar surface
    elseif (systemtype.eq.60) then
        call update_matrix_60(flag)         ! channel + particles
    endif
    
    call calcfv                             ! == calulate variable fv matrix ! IMPORTANT

    if(flag.eqv..true.) then
        write(stdout,*) 'Initial position of particle does not fit in z'
        write(stdout,*) 'or particles collide'
        stop
    else
        if(rank.eq.0) write(stdout,*) 'Particle OK'
    endif

    call graftpoints
    if(rank.eq.0) write(stdout,*) 'Graftpoints OK'

    call creador ! Genera cadenas
    if(rank.eq.0) write(stdout,*) 'Creador OK'


#ifdef _MKL
    if (flagmkl.eq.1) then ! use compressed MKL CSR format to store chains
        call px2csr
    if(rank.eq.0)print*, 'PX2CSR OK'
        elseif (flagmkl.eq.2) then
    call px2csr_map
        if(rank.eq.0)print*, 'PX2CSR MAP OK'
    endif
#endif

    ! == initial guess  what is flag ??

    if(infile.ne.0) then
        call retrivefromdisk(counter)
        if(rank.eq.0) write(stdout,*) 'Load input from file'
        if(rank.eq.0) write(stdout,*) 'Free energy', free_energy
        if(infile.eq.3) call mirror
        if(infile.ne.-1)infile = 2
        if(flag.eqv..true.) then
            write(stdout,*) 'Initial position of particle does not fit in z'
            write(stdout,*) 'or particles collide'
            stop
        endif
    endif

    ii = 1
    sc = scs(ii)  ! == what is scs 

    select case (vscan) ! == variable scan over variabel kp ??
    case (1)
        st = sts(1)
        kp = 1.0d10+kps(1)
        do i = 1, nkp
            do while (kp.ne.kps(i))
                kp = kps(i)
                if(rank.eq.0) write(stdout,*)'Switch to kp = ', kp
                flagcrash = 1
                do while(flagcrash.eq.1)
                    flagcrash = 0
                    call CPU_TIME(time0)

                    call solve(flagcrash)       ! == call to solver system
                    
                    call CPU_TIME(timeF)
                    if(rank.eq.0) print*,'Timer:',timeF-time0
                    if(flagcrash.eq.1) then
                        if(i.eq.1) stop
                        kp = (kp + kpOK)/2.0
                        if(rank.eq.0) write(stdout,*)'Error, switch to kp = ', kp
                    endif
                enddo

                kpOK = kp ! last st solved OK
                if(rank.eq.0) write(stdout,*) 'Solved OK, kp: ', kpOK
    
            enddo

            counterr = counter + i + ii  - 1
            call Free_Energy_Calc(counterr)
            if(rank.eq.0) write(stdout,*) 'Free energy after solving', free_energy

            call savedata(counterr)
            if(rank.eq.0)write(stdout,*) 'Save OK'
            call store2disk(counterr)

        enddo

    case (2)

        kp = 0
        st = 1.0d10+sts(1)
        do i = 1, nst
            do while (st.ne.sts(i))
                st = sts(i)
                if(rank.eq.0) write(stdout,*)'Switch to st = ', st
                flagcrash = 1
                do while(flagcrash.eq.1)
                    flagcrash = 0

                    call solve(flagcrash) ! == call to solver system
                    
                    if(flagcrash.eq.1) then
                        if(i.eq.1) stop
                        st = (st + stOK)/2.0
                        if(rank.eq.0) write(stdout,*)'Error, switch to st = ', st
                    endif
                enddo

                stOK = st ! last st solved OK
                if(rank.eq.0)write(stdout,*) 'Solved OK, st: ', stOK

            enddo

            counterr = counter + i + ii  - 1
            
            call Free_Energy_Calc(counterr)
            if(rank.eq.0) write(stdout,*) 'Free energy after solving', free_energy
            
            call savedata(counterr)
            if(rank.eq.0) write(stdout,*) 'Save OK'
            
            call store2disk(counterr)

        enddo

    end select

    call endall

end program main


