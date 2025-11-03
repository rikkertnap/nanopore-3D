module graftpoint 

    implicit none 

contains 

    subroutine read_graftpoint(rS,rL,lenChannel,ngraft,positiongraft,info)
        
        use const, only : stdout
        use MPI, only : rank, ierr

        implicit none

        ! ==input arguments

        real*8, intent(in) :: rS                ! == smallest radius channel  
        real*8, intent(in) :: rL                ! == largest  radius channel   
        real*8, intent(in) :: lenChannel        ! == length  channel   
        integer, intent(in) :: ngraft           ! == number of graft points
        real*8, intent(inout) :: positiongraft(:,:)     ! == locatation graft point  => integer goes to p1
        integer, intent(inout) :: info

        ! == local argument 

        character(len=14) :: fname
        integer :: un, ngraftin, ios, i, unnew, infotmp
        real*8 :: rSin, rLin, lenChannelin
        character(len=80) :: text
        character(len=100) :: io_msg  
        logical :: isReadGood, isWrite
        character :: comment

        info=0
        infotmp=0
        isWrite=.false.

        ! == reading in of graft point from file
        write(fname,'(A14)')'graftpoints.in'
        open(newunit=un,file=fname,iostat=ios,status='old',iomsg=io_msg)
        if(ios >0 ) then
            write(stdout,*)'Error opening graftpoints.in file : iostat =', ios
            infotmp = 1
        endif
    
        ! == read preamble

        read(un,*,iostat=ios)rSin
        if(ios/=0) isReadGood=.false.
        read(un,*,iostat=ios)rLin 
        if(ios/=0) isReadGood=.false.   
        read(un,*,iostat=ios)lenChannelin 
        if(ios/=0) isReadGood=.false.   
        read(un,*,iostat=ios)ngraftin
            
        if(abs(rSin-rS)>0.0000001) then 
            write(stdout,*)"RadiusS graft file: ",Rsin," not equal inputted internal Rs: ",rS
            infotmp=2
        endif
        if(abs(rLin-rL)>0.0000001) then 
            write(stdout,*)"RadiusL graft file: ",RLin," not equal inputted internal RL: ",rL
            infotmp=3
        endif
        if(abs(lenChannel-lenChannelin)>0.0000001) then 
            write(stdout,*)" Lenchannel graft file not equal inputted internal value"
            infotmp=4
        endif
        if(abs(ngraft-ngraftin)>0.0000001) then 
            write(stdout,*)"ngraft graft file: ",ngraftin," not equal inputted internal ngraft: ",ngraft
            infotmp=5
        endif

        read(un,*,iostat=ios)comment

        do i=1,ngraft
            read(un,*,iostat=ios)positiongraft(i,1),positiongraft(i,2),positiongraft(i,3)
            if(ios/=0) then 
                isReadGood=.false.
                infotmp=6
            endif    
        enddo 

        if(isReadGood.eqv..false.) info=infotmp
        
        close(un)

        if(info.ne.0) then
            write(stdout,*) 'failure to read graftpoint info:',info
            call MPI_FINALIZE(ierr) ! finaliza MPI
            stop
        endif    

        if(isWrite) call write_graftpoint(rS,rL,lenChannel,ngraft,positiongraft,info)

    end subroutine read_graftpoint

    ! write grafing points 
 

    subroutine write_graftpoint(rS,rL,lenChannel,ngraft,positiongraft,info)

        use const, only : stdout
        use MPI, only : rank
        implicit none

        ! ==input arguments

        real*8, intent(in) :: rS                ! == smallest radius channel  
        real*8, intent(in) :: rL                ! == largest  radius channel   
        real*8, intent(in) :: lenChannel        ! == length  channel   
        integer, intent(in) :: ngraft           ! == number of graft points
        real*8, intent(in) :: positiongraft(:,:)    ! == locatation graft point  => integer goes to p1
        integer, intent(inout) :: info          !== output information 
       
        ! == local argument 

        character(len=15) :: fname
        integer :: un, ios, i, unnew
        character(len=80) :: text
        character(len=100) :: io_msg  
        logical :: isWriteGood
      
        info=0

        if(rank==0) then 

            ! == writing in of graft point from file
            write(fname,'(A15)')'graftpoints.out'
            open(newunit=un,file=fname,iostat=ios,status='replace',iomsg=io_msg)
            if(ios >0 ) then
                write(stdout,*)'Error opening graftpoint-run.out file : iostat =', ios
                info = ios
            endif
        
            ! == write preamble
            write(un,*,iostat=ios)rS
            if(ios/=0) isWriteGood=.false.
            write(un,*,iostat=ios)rL
            if(ios/=0) isWriteGood=.false.     
            write(un,*,iostat=ios)Lenchannel
            if(ios/=0) isWriteGood=.false. 
            write(un,*,iostat=ios)ngraft
            if(ios/=0) isWriteGood=.false. 
            write(un,*,iostat=ios)"#"

            do i=1,ngraft
                write(un,*,iostat=ios)positiongraft(i,1),positiongraft(i,2),positiongraft(i,3)
                if(ios/=0) isWriteGood=.false.   
            enddo 

            if(isWriteGood.eqv..false.) info=1
          
            close(un)
     
        endif    
    
    end subroutine write_graftpoint

end module graftpoint