subroutine makemaps
    use maps
    use system
    use ematrix

    implicit none
    integer ix,iy,iz,i

    ALLOCATE(imap(dimx,dimy,dimz))
    ALLOCATE(mapx(dimx*dimy*dimz))
    ALLOCATE(mapy(dimx*dimy*dimz))
    ALLOCATE(mapz(dimx*dimy*dimz))

    do ix = 1, dimx
        do iy = 1, dimy
            do iz = 1, dimz
                imap(ix,iy,iz) = ix+dimx*(iy-1)+dimx*dimy*(iz-1)
                mapx(imap(ix,iy,iz))=ix
                mapy(imap(ix,iy,iz))=iy
                mapz(imap(ix,iy,iz))=iz
            enddo
        enddo
    enddo
end subroutine 

subroutine calcfv

    use system
    use ematrix
    use maps

    implicit none
    integer ix,iy,iz,i

    do ix = 1, dimx
        do iy = 1, dimy
            do iz = 1, dimz
                fvstd(ix,iy,iz) = 1.0-volprot(ix,iy,iz)
                fvmkl(imap(ix,iy,iz))=1.0-volprot(ix,iy,iz)
            enddo
        enddo
    enddo

end subroutine calcfv


! == calcutate binarized version of fvstd 
! == if free volume cell = fv < epsfv  then fvstdint = 0  otherwise 1 
! == Indicate whether a cell has accesible volume or not 
! == used  in div_flux_channel to determine which cell have zero flux 
! == pre volprot needs to be computed 

subroutine calcfvint
    use const, only : stdout
    use system, only : dimx, dimy, dimz
    use ematrix, only : volprot, fvstdint
    use MPI, only : rank 
    use channel, only : Rdimz
    
    implicit none

    ! local arguments 
    integer :: ix,iy,iz,i
    real*8 :: fv
    integer :: teller, maxteller
    real*8, parameter :: epsfv = 0.001d0 ! threshold 

    fvstdint = 0 ! init 
    teller = 0

    maxteller = dimx*dimy*(dimz-2*Rdimz)

    do iz = 1, dimz
        do iy = 1, dimy
            do ix = 1, dimx
    
                fv = 1.0d0-volprot(ix,iy,iz) 

                if(fv<epsfv) then         ! carefull float precsion volprot , check value
                    fvstdint(ix,iy,iz) = 0
                    teller = teller +1
                else 
                    fvstdint(ix,iy,iz) = 1
                endif         
            
            enddo
        enddo
    enddo

    ! x = 0
    fvstdint(0,:,:) = fvstdint(1,:,:) 
        
    ! x = dimx
    fvstdint(dimx+1,:,:) = fvstdint(dimx,:,:)  

    ! y = 0
    fvstdint(:,0,:) = fvstdint(:,1,:) 
        
    ! y = dimy
    fvstdint(:,dimy+1,:) = fvstdint(:,dimy,:) 
        
    ! z = 0
    fvstdint(:,:,0) = fvstdint(:,:,1)
        
    ! z = dimz
    fvstdint(:,:,dimz+1) = fvstdint(:,:,dimz)

    if(rank.eq.0) write(stdout,*) 'calcintfv: teller=', teller, " max_teller= ",maxteller

end subroutine calcfvint 


!   T1(x,y,z) = idx with 
!   idx = x + nx*y + nx * ny * z  -nx -nx*ny = x = x +nx*(y-1)+nx * ny(z-1) = 
!        =  x +dimx* (y-1) +dimz*(z-1)
!   T2(ix,iy,iz) = idx with 
!   idx = ix+dimx*(iy-1)+dimx*dimy*(iz-1)
!   T1 identical to T2  

subroutine linearIndexFromCoordinate(x,y,z,idx)
     
    use system, only : dimx, dimy
    implicit none 
      
    integer, intent(in)  :: x,y,z
    integer, intent(out) :: idx
      
    integer :: a,b,c,d

    a = 1
    b = dimx 
    c = dimx * dimy 
    d = 1-a-b-c
    idx = a*x + b*y + c*z + d 

end subroutine linearIndexFromCoordinate


subroutine coordinateFromLinearIndex(idx, x,y, z)

    use system, only : dimx, dimy
    implicit none

    integer, intent(out)  :: x,y,z
    integer, intent(in)   :: idx
    integer :: idxtmp
   
    idxtmp=idx
    x =  mod(idxtmp-1,dimx)+1
    idxtmp =int((idxtmp-1)/dimx)+1
    y = mod(idxtmp-1,dimy)+1
    idxtmp = int((idxtmp-1)/dimy)+1
    z = idxtmp
    
end subroutine coordinateFromLinearIndex

subroutine allocate_hashtable

    use maps, only :  coordtoindex, indextocoord
    use system, only : dimx, dimy, dimz
   
    integer ::  nsize

    nsize = dimx * dimy * dimz

    allocate(coordtoindex(dimx,dimy,dimz))
    allocate(indextocoord(nsize,3))
        
end subroutine allocate_hashtable

subroutine make_hashtable

    use system, only : dimx, dimy, dimz
    use maps, only :  coordtoindex, indextocoord

    integer :: idx, ix, iy, iz
    integer ::  nsize   

    call allocate_hashtable

    nsize = dimx * dimy * dimz
    
    do idx=1,nsize
            
        call coordinateFromLinearIndex(idx, ix, iy, iz)
        coordtoindex(ix,iy,iz) = idx  ! equivalent to imap

        indextocoord(idx,1) = ix       ! equivalent to mapx, mapy, mapz
        indextocoord(idx,2) = iy
        indextocoord(idx,3) = iz

    enddo
        
end subroutine make_hashtable