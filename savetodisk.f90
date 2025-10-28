subroutine savetodisk(array, title, counter)

    use system, only : dimx,dimy,dimz, vtkflag
    use const, only : stdout

    implicit none

    real*8, intent(inout) :: array(dimx, dimy, dimz) 
    character*5, intent(in) :: title
    integer, intent(in) :: counter

    select case (vtkflag)
    case(1)
        call savetodisk_binary(array, title, counter)
    case(2) 
        call savetodisk_ascii(array, title, counter)
    case(3)  
        call savetodisk_xyz(array, title, counter)
    case default
        write(stdout,*)"savetodisk failed vtkflag :",vtkflag        
    end select    

end subroutine savetodisk


subroutine savetodisk_binary(array, title, counter)

    use system
    use transform
    use s2d
    implicit none

    real*8, intent(inout) :: array(dimx, dimy, dimz) 
    character*5, intent(in) ::  title
    integer, intent(in) :: counter

    ! local arguments

    integer :: npoints
    integer :: ix, iy, iz, jx, jy, jz
    real*8 :: arraytemp(dimx, dimy, dimz)
    character*21 :: filename, tempc
    real*4 :: singlepres
    real*8 :: v(3), x(3)
    real*4 :: xx(3)
    integer, external :: PBCSYMI, PBCREFI
    integer :: i
    character*1 :: lf
    character*8 :: str1, str2, str3

    npoints = (dimx*scx+1)*(dimy*scy+1)*(dimz*scz+1)

    lf = char(10) ! line feed character

    !-----  coordinates -------------------------
    ! Variables

    arraytemp = 0
    do iz = 1, dimz
        arraytemp(:,:,iz) = array(:,:,iz)
    enddo
    
    array = -1000
    do iz = 1, dimz
        array(:,:,iz) = arraytemp(:,:,iz)
    enddo

    ! VTK for Paraview


    if(vtkflag.eq.1) then

        write(filename,'(A5, A1,I3.3, A4)')title,'.', counter, '.vtk'
        open (unit=45, file=filename, access='stream',form='unformatted', convert='BIG_ENDIAN')

        write(45)'# vtk DataFile Version 2.0', lf
        write(45)title, lf
        write(45)'BINARY', lf
        write(45)'DATASET STRUCTURED_GRID ', lf
        write(str1(1:8),'(i8)') scz*dimz+1
        write(str2(1:8),'(i8)') scy*dimy+1
        write(str3(1:8),'(i8)') scx*dimx+1

        write(45)'DIMENSIONS', str1, ' ', str2, ' ',str3, lf

        write(str1(1:8),'(i8)') npoints
        write(45)'POINTS ',str1,' float', lf


        i = 1
        do ix = 0, scx*dimx
            do iy = 0, scy*dimy
                do iz = 0, scz*dimz

                    v(1) = (ix-float(dimx)/2.)*delta ! transformed coordinates
                    v(2) = (iy-float(dimy)/2.)*delta
                    v(3) = (iz-float(dimz)/2.)*delta

                    x = MATMUL(IMAT,v)
                    xx = x
                    write(45) xx(:)
                enddo
            enddo
        enddo

        write(str1(1:8),'(i8)') scx*scy*scz*dimx*dimy*dimz
        write(45)lf, 'CELL_DATA ', str1, lf
        tempc = 'SCALARS ' // title // ' float 1'
        write(45)tempc, lf
        write(45)'LOOKUP_TABLE default', lf

        do ix = 1, scx*dimx
            do iy = 1, scy*dimy
                do iz = 1, scz*dimz

                    if(PBC(2).ne.3)jx = PBCSYMI(ix,dimx)
                    if(PBC(2).eq.3)jx = PBCREFI(ix,dimx)

                    if(PBC(4).ne.3)jy = PBCSYMI(iy,dimy)
                    if(PBC(4).eq.3)jy = PBCREFI(iy,dimy)

                    if(PBC(6).ne.3)jz = PBCSYMI(iz,dimz)
                    if(PBC(6).eq.3)jz = PBCREFI(iz,dimz)

                    singlepres = array(jx, jy, jz) ! Lo necesito en single presicion
        
                    write(45)singlepres
                enddo
            enddo
        enddo
        write(45)lf
        close (45)

    endif

end subroutine savetodisk_binary

subroutine savetodisk_xyz(array, title, counter)


    use system, only : dimx,dimy,dimz

    !use transform
    !use s2d
    implicit none

    real*8, intent(inout) :: array(dimx, dimy, dimz) 
    character*5, intent(in) :: title
    integer, intent(in) :: counter

    ! local arguments
    integer :: ix, iy, iz
    real*8 :: arraytemp(dimx, dimy, dimz)
    character*21 :: filename
   

    arraytemp = 0
    do iz = 1, dimz
        arraytemp(:,:,iz) = array(:,:,iz)
    enddo
    
    array = -1000
    do iz = 1, dimz
        array(:,:,iz) = arraytemp(:,:,iz)
    enddo


    write(filename,'(A5, A1, I3.3, A4)') title,'.', counter, '.raw' 
    open(unit=45, file=filename)
 
    do iz = 1, dimz
        do iy = 1, dimy
            do ix = 1, dimx
                write(45,*)array(ix,iy,iz)
            enddo
        enddo
    enddo
end subroutine savetodisk_xyz

! outine to save vtk files in ASCII format

subroutine savetodisk_ascii(array, title, counter)

    use system
    use transform
    use s2d
    implicit none
 
    real*8, intent(inout) :: array(dimx, dimy, dimz)
    character*5, intent(in) :: title
    integer, intent(in)  :: counter
    
    ! local variables
    integer :: ix, iy, iz, jx, jy, jz
    real*8 :: arraytemp(dimx, dimy, dimz)
    real*8 :: arrayz(dimz)
    character*21 :: filename, tempc
    real*4 :: singlepres
    real*8 :: x(3), v(3)
    integer, external :: PBCSYMI, PBCREFI

    !-----  coordenadas -------------------------
    ! Variables

    arraytemp = 0
    do iz = 1, dimz
        arraytemp(:,:,iz) = array(:,:,iz)
    enddo
    
    array = -1000
    do iz = 1, dimz
        array(:,:,iz) = arraytemp(:,:,iz)
    enddo

   
    if(vtkflag.ne.1) then

        write(filename,'(A5, A1,I3.3, A4)')title,'.', counter, '.vtk'
        open (unit=45, file=filename)
        write(45,'(A)')'# vtk DataFile Version 2.0'
        write(45,'(A)')title
        write(45,'(A)')'ASCII'
        write(45,'(A)')'DATASET STRUCTURED_GRID '
        write(45,'(A, I5, A1, I5, A1, I5)')'DIMENSIONS', scz*dimz+1, ' ', scy*dimy+1, ' ',scx*dimx+1
        write(45,'(A, I8, A)')'POINTS ',(scx*dimx+1)*(scy*dimy+1)*(scz*dimz+1),' float'
        do ix = 0, scx*dimx
            do iy = 0, scy*dimy
                do iz = 0, scz*dimz

                v(1) = (ix-float(dimx)/2.)*delta ! transformed coordinates
                v(2) = (iy-float(dimy)/2.)*delta
                v(3) = (iz-float(dimz)/2.)*delta

                x = MATMUL(IMAT,v)
        
                    write(45, *) x(1),'   ',  x(2), '   ', x(3) ! Cartesian coordinates 
                enddo
            enddo
        enddo

        write(45,'(A, I8)')'CELL_DATA ', scx*scy*scz*dimx*dimy*dimz
        tempc = 'SCALARS ' // title // ' float 1'
        write(45,'(A)')tempc
        write(45,'(A)')'LOOKUP_TABLE default'

        do ix = 1, scx*dimx
            do iy = 1, scy*dimy
                do iz = 1, scz*dimz

                    if(PBC(2).ne.3)jx = PBCSYMI(ix,dimx)
                    if(PBC(2).eq.3)jx = PBCREFI(ix,dimx)

                    if(PBC(4).ne.3)jy = PBCSYMI(iy,dimy)
                    if(PBC(4).eq.3)jy = PBCREFI(iy,dimy)

                    if(PBC(6).ne.3)jz = PBCSYMI(iz,dimz)
                    if(PBC(6).eq.3)jz = PBCREFI(iz,dimz)

                    singlepres = array(jx, jy, jz) ! Lo necesito en single presicion
        
                    write(45,*)singlepres
                enddo
            enddo
        enddo
        close (45)
    endif

end subroutine savetodisk_ascii
 

