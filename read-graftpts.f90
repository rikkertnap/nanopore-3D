module graftpoint 

    implicit none 

contains 

    ! Reads in of ngraft graftpoints from file graftpoints.in that are grafted on a surface 
    ! defined by rS, rl , LenChannel 
    ! retrun real*8 positiongraft(:,:) range (ngraft,3) 
    !        integer info : =0 is succesfull nonzero stops program
 
    subroutine read_graftpoint(rS,rL,lenChannel,ngraft,positiongraft,info)
        
        use const, only : stdout
        use MPI, only : rank

        implicit none

        ! ==input arguments

        real*8, intent(in) :: rS                ! == smallest radius channel  
        real*8, intent(in) :: rL                ! == largest  radius channel   
        real*8, intent(in) :: lenChannel        ! == length channel   
        integer, intent(in) :: ngraft           ! == number of graft points
        real*8, intent(inout) :: positiongraft(:,:)     ! == locatation graft point  => integer goes to p1
        integer, intent(inout) :: info           ! return status

        ! == local argument 

        character(len=14) :: fname
        integer :: un, ngraftin, ios, i, unnew 
        real*8 :: rSin, rLin, lenChannelin
        character(len=80) :: text
        character(len=100) :: io_msg  
        logical :: isReadGood, exist, isWrite
        character :: comment

        info=0
        isWrite=.False.
    
        ! == reading in of graft point from file
        write(fname,'(A14)')'graftpoints.in'
        inquire(file=fname,exist=exist)
        if(exist) then 
            open(newunit=un,file=fname,iostat=ios,status='old',iomsg=io_msg)
            if(ios >0 ) then
                if(rank.eq.0) write(stdout,*)'Error opening file graftpoints.in iostat =', ios
                info = 1
                return
            endif
        else
            if(rank.eq.0) write(stdout,*)'Error file graftpoints.in does not exist.'
            info= 7
            return
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
            if(rank.eq.1)write(stdout,*)"RadiusS graft file: ",Rsin," not equal inputted internal Rs: ",rS
            info=2
            return
        endif
        if(abs(rLin-rL)>0.0000001) then 
            if(rank.eq.0)write(stdout,*)"RadiusL graft file: ",RLin," not equal inputted internal RL: ",rL
            info=3
            return
        endif
        if(abs(lenChannel-lenChannelin)>0.0000001) then 
            if(rank.eq.0)write(stdout,*)" Lenchannel graft file not equal inputted internal value"
            info=4
            return
        endif
        if(abs(ngraft-ngraftin)>0.0000001) then 
            if(rank.eq.0)write(stdout,*)"ngraft graft file: ",ngraftin," not equal inputted internal ngraft: ",ngraft
            info=5
            return
        endif

        read(un,*,iostat=ios)comment

        do i=1,ngraft
            read(un,*,iostat=ios)positiongraft(i,1),positiongraft(i,2),positiongraft(i,3)
            if(ios/=0) then 
                isReadGood=.false.
                if(rank.eq.0)write(stdout,*)"Reading error in reading graftpoint.in: ios:",ios 
                info=6
                return
            endif    
        enddo 

        close(un)

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
                if(rank.eq.0) write(stdout,*)'Error opening graftpoint-run.out file : iostat =', ios
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

    ! read in graft point sequence 

    subroutine read_pattern_grafts(ngraft,hasGraftA,hasGraftB,info)

        use const, only : stdout
        use MPI, only : rank

        implicit none

        ! ==input arguments

        integer, intent(in) :: ngraft                       ! == number of graft points
        logical, intent(inout) :: hasGraftA(:),hasGraftB(:) ! == logical variable indicating that graft point i has a A and or B  chain grafted     integer, intent(inout) :: info                      ! == return status
        integer, intent(inout) :: info 

        ! == local argumets

        logical :: exist
        integer :: un, i, line, maxline, ios , ix, iy
        character(len=80) :: str 
        character(len=16) :: fname
        character(len=100) :: io_msg  

        info = 0 

        ! == reading sequence of graft point from file
        write(fname,'(A16)')'graftsequence.in'
        open(newunit=un,file=fname,iostat=ios,status='old',iomsg=io_msg)
        inquire(file=fname,exist=exist)
        if(exist) then 
            open(newunit=un,file=fname,iostat=ios,status='old',iomsg=io_msg)
            if(ios >0 ) then
                if(rank.eq.0)write(stdout,*)'Error opening graftsequence.in file : iostat =', ios
                info=1
                return
            endif
        else
            if(rank.eq.0)write(stdout,*)'Error graftsequence.indoes not exist.'
            info = 7
            return
        endif
      

        ! init
        hasGraftA = .False.
        hasGraftB = .False.

        line = 0
        ios = 0
        maxline = ngraft
        
        do while (line<maxline.and.ios==0)
            line=line+1
            read(un,*,iostat=ios)ix,iy
            if(ix.eq.1) hasgraftA(line)=.True.
            if(iy.eq.1) hasgraftB(line)=.True.
        enddo
        
        close(un)

         if(line/=maxline.or.ios/=0) then 
            str="reached end of file before all elements read or ios error"
            if(rank.eq.0)write(stdout,*)str
            str="read file "//trim(adjustl(fname))//" failed"
            if(rank.eq.0) write(stdout,*)str
            info = 1
            return
        endif
        
    end subroutine read_pattern_grafts


end module graftpoint