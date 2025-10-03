!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Esta subrutina se encarga de poner a todas los segmentos dentro del slab
! === This subroutine is responsible for placing all the segments within the slab.



subroutine pxs

    use system
    use MPI
    use chainsdat
    use conformations
    use const
    use transform
    use mparameters_monomer
    implicit none
        
    integer j, ii, jj,i
    real*8 pxtemp(3,long)
    real*8 xx(3)
    real*8 x(3)
    real*8 v(3)
    integer testsystem
    integer testsystemr
    integer testsystemc
    real*8 maxx(3)
    integer flag
    integer aa
    real*4 ztemp

    integer, external :: PBCREFI, PBCSYMI

    maxx(1) = float(dimx)*delta       ! == max dimension lattice 
    maxx(2) = float(dimy)*delta
    maxx(3) = float(dimz)*delta

    do jj = 1, cpp(rank+1)            ! == cpp , jj ?? comes from cadenas and subroutine graftpoit
        ii = cppini(rank+1)+jj        ! == cppini ii ??
        flag = 0

        ztemp = 0.0
        do j=1,long
            x(1) = in1(j ,2)         ! == coordianes chain for cadenas/creador
            x(2) = in1(j, 3)
            x(3) = in1(j, 1)

            if((systemtype.eq.2).or.(systemtype.eq.3).or.(systemtype.eq.4).or.(systemtype.eq.41)   &
        .or.(systemtype.eq.42).or.(systemtype.eq.52).or.(systemtype.eq.60)) call rot_chain_cyl(x,ii)
                                                 ! ==  rotation ??   

            x = x + posicion(ii,:)               ! == translate to graft position given by ii

            ztemp = ztemp+ x(3)*zpol(segtype(j)) ! == charge of z position * segtype(j)=segment_number of segment j
                                                 ! == dipole moment ???

            v = MATMUL(MAT,x)                    ! == coordinate transform from xyz to uvz 
            pxtemp(:,j) = v(:)                   ! == assgin v to pxtemp   

            ! == test location 

            select case (systemtype)
            case (1)
 
                if(testsystem(x).eq.-1) then ! if testsystem = -1,  there is a collision with all or particle 
                    flag = -1
                    exit
                endif

                if(testsystem(x).eq.-2) then ! if testsystem = -2, the polymer goes out-of-system
                    write(stdout,*) 'pxs: out-of-system'
                    stop
                endif

            case (6)
       
                if(testsystem(x).eq.-1) then ! if testsystem = -1,  there is a collision with all or particle 
                    flag = -1
                    exit
                endif

                if(testsystem(x).eq.-2) then ! if testsystem = -2, the polymer goes out-of-system
                    write(stdout,*) 'pxs: out-of-system'
                    stop
                endif

            case (60)

                if(testsystemc(x).eq.-1) then ! if testsystem = -1,  there is a collision with channel
                    flag = -1
                    exit
                endif

                if(testsystemc(x).eq.-2) then ! if testsystem = -2, the polymer goes out-of-system
                    write(stdout,*) 'pxs: out-of-system'
                    stop
                endif

                if(testsystem(x).eq.-1) then ! if testsystem = -1,  there is a collision with all or one particle 
                    flag = -1
                    exit
                endif

                if(testsystem(x).eq.-2) then ! if testsystem = -2, the polymer goes out-of-system
                    write(stdout,*) 'pxs: out-of-system'
                    stop
                endif

            case (2, 3, 4, 41, 42)


                if(testsystemc(x).eq.-1) then ! if testsystem = -1,  there is a collision with all or particle 
                    flag = -1
                    exit
                endif

                if(testsystemc(x).eq.-2) then ! if testsystem = -2, the polymer goes out-of-system
                    write(stdout,*) 'pxs: out-of-system'
                    stop
                endif

            case (52)

                if(testsystemr(x).eq.-1) then ! if testsystem = -1,  there is a collision with all or particle 
                    flag = -1
                    exit
                endif

                if(testsystemr(x).eq.-2) then ! if testsystem = -2, the polymer goes out-of-system
                    write(stdout,*) 'pxs: out-of-system'
                    stop
                endif

            end select

        enddo ! j
    
        if(flag.eq.0) then        ! == no conflict :location good 

            newcuantas(ii) = newcuantas(ii)+1. ! == increase newcuats for ii
            zfinal(newcuantas(ii),jj) = ztemp
            ngauche(newcuantas(ii),ii) = ing

            do j = 1, long
                aa = floor(pxtemp(1,j)/delta) + 1
                px(newcuantas(ii),j,jj) = aa
                if(aa.lt.1) then
                    if(PBC(1).eq.1)px(newcuantas(ii),j,jj) = PBCSYMI(aa,dimx)
                    if(PBC(1).eq.3)px(newcuantas(ii),j,jj) = PBCREFI(aa,dimx)
                endif
                if(aa.gt.dimx) then
                    if(PBC(2).eq.1)px(newcuantas(ii),j,jj) = PBCSYMI(aa,dimx)
                    if(PBC(2).eq.3)px(newcuantas(ii),j,jj) = PBCREFI(aa,dimx)
                endif

                aa = floor(pxtemp(2,j)/delta) + 1
                py(newcuantas(ii),j,jj) = aa
                if(aa.lt.1) then
                    if(PBC(3).eq.1)py(newcuantas(ii),j,jj) = PBCSYMI(aa,dimy)
                    if(PBC(3).eq.3)py(newcuantas(ii),j,jj) = PBCREFI(aa,dimy)
                endif
                if(aa.gt.dimy) then
                    if(PBC(4).eq.1)py(newcuantas(ii),j,jj) = PBCSYMI(aa,dimy)
                    if(PBC(4).eq.3)py(newcuantas(ii),j,jj) = PBCREFI(aa,dimy)
                endif

                aa = floor(pxtemp(3,j)/delta) + 1
                pz(newcuantas(ii),j,jj) = aa
                if(aa.lt.1) then
                    if(PBC(5).eq.1)pz(newcuantas(ii),j,jj) = PBCSYMI(aa,dimz)
                    if(PBC(5).eq.3)pz(newcuantas(ii),j,jj) = PBCREFI(aa,dimz)
                endif
                if(aa.gt.dimz) then
                    if(PBC(6).eq.1)pz(newcuantas(ii),j,jj) = PBCSYMI(aa,dimz)
                    if(PBC(6).eq.3)pz(newcuantas(ii),j,jj) = PBCREFI(aa,dimz)
                endif
    
            enddo
        endif ! == flag
    enddo ! jj
    
    return
end subroutine pxs
      


subroutine rot_chain_cyl(x,ii)
    use rotchain
    implicit none
    integer, intent(in) :: ii
    real*8 , intent(inout) :: x(3)

    ! local variable 
    real*8 y(3)
    real*8 t

    t = rotangle(ii)
    y = x
    x(1) = cos(t)*y(1)+sin(t)*y(2)
    x(2) = -sin(t)*y(1)+cos(t)*y(2)
    x(3) = y(3)
end
