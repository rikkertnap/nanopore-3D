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
end

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
end


! == calcutate binarized version of fvstd 
! == if fvstd < 0.000001d0 then fvstdint = 0  otherwise 1 
! == indicated whether a cell has accesible volume or not 
! == used  in div_flux_channel to determine which cell have zero flux 
! == pre volprot needs to be computed 

subroutine calcintfv

    use system, only : dimx, dimy, dimz
    use ematrix, only : volprot, fvstdint
    
    implicit none

    ! local arguments 
    integer :: ix,iy,iz,i
    real*8 :: fv

    fvstdint = 0 ! init 

    do iz = 1, dimz
        do iy = 1, dimy
            do ix = 1, dimx
    
                fv = 1.0-volprot(ix,iy,iz) 

                if(fv<0.00001d0) then ! carefull float precsion , check this
                    fvstdint = 0
                else 
                    fvstdint = 1
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
 
end subroutine calcintfv 

