module moleculeslist

    implicit none

    type moleclist
        real*8 :: sol
        real*8 :: pos
        real*8 :: neg
        real*8 :: Hplus
        real*8 :: OHmin
    end type moleclist

    type(moleclist) :: vol, zval, mumin, mumax, xvolmin, xvolmax, Diffcoeff
    
    type(moleclist) :: muexmin, muexmax, rho_tilde_min, rho_tilde_max, Diffcoeff_tilde_min, Diffcoeff_tilde_max

contains

    subroutine init_zero_moleclist(struct)
    
        type(moleclist), intent(inout) :: struct

        struct%sol = 0.0d0
        struct%pos = 0.0d0
        struct%neg = 0.0d0
        struct%Hplus = 0.0d0
        struct%OHmin = 0.0d0
 
    end subroutine init_zero_moleclist


    function get_value_moleclist(struct,member) result(val)

        type(moleclist) , intent(in) :: struct 
        character(len=5),  intent(in) :: member
        real*8 :: val

        select case (member)
            case ("sol") 
                val=struct%sol
            case ("pos") 
                val=struct%pos
            case ("neg") 
                val=struct%neg
            case ("Hplus") 
                val=struct%Hplus
            case ("OHmin") 
                val=struct%OHmin
            case default
                print*,"Wrong value member molecule list :  ",member
                stop
        end select      
            
    end function get_value_moleclist


    function sum_value_moleclist(struct) result(val)

        type(moleclist) , intent(in) :: struct
        real*8 :: val

        val = struct%sol + struct%pos +struct%neg + struct%Hplus + struct%OHmin
        
    end function sum_value_moleclist   


    subroutine init_diffusion_coeff()

        ! values from "Diffusion- Mass transfer in Fuild System" by E.L.Cussler page 143
        ! at T= 25 C
        ! all values are in 10-5 cm^2/sec 10-9 m^2/sec

        Diffcoeff%pos   = 1.33d-9 ! Na
        Diffcoeff%neg   = 2.03d-9 ! Cl 
        Diffcoeff%Hplus = 9.31d-9
        Diffcoeff%OHmin = 5.28d-9  
        
        ! assign value not needed for calcualting conductivity  not charged 
        Diffcoeff%sol   = 0.0d-9   
        
    end subroutine  init_diffusion_coeff


    subroutine init_vol()

        use molecules, only : vsol, vsalt

        vol%sol   = 1.0d0
        vol%pos   = vsalt ! == deivide by vsol
        vol%neg   = vsalt 
        vol%Hplus = 1.0d0
        vol%OHmin = 1.0d0

    end subroutine init_vol    


    subroutine init_zval()

        use molecules, only : zpos, zneg 

        zval%sol   = 0.0d0
        zval%pos   = 1.0d0*zpos
        zval%neg   = 1.0d0*zneg 
        zval%Hplus = 1.0d0
        zval%OHmin = -1.0d0

    end subroutine init_zval

    subroutine init_xvolmin()
        
        use bulk, only : xsolbulk, xposbulk, xnegbulk, xHplusbulk,xOHminbulk

        xvolmin%sol   = xsolbulk
        xvolmin%pos   = xposbulk
        xvolmin%neg   = xnegbulk
        xvolmin%Hplus = xHplusbulk
        xvolmin%OHmin = xOHminbulk
 
    end subroutine init_xvolmin



end module moleculeslist

