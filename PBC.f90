integer function PBCSYMI(i,dimi) ! returns the PBC cell coordinate 
integer i, dimi
PBCSYMI = mod(i-1+5*dimi, dimi) + 1
end function

integer function PBCREFI(i,dimi) ! returns the reflection cell coordinate 
integer i, dimi, iaux, iaux2
iaux = mod(i-1+5*dimi,dimi)+1
iaux2 = mod(int((i-1)/dimi),2)
PBCREFI = iaux+(dimi-2*iaux+1)*iaux2
end function


real function PBCSYMR(i,dimi) ! returns the PBC cell coordinate 
real*8 i, dimi
PBCSYMR = mod(i-1.0+5.0*dimi, dimi) + 1.0
end function

real function PBCREFR(i,dimi) ! returns the reflection cell coordinate 
real*8 i, dimi, iaux
integer iaux2
iaux = mod(i-1.0+5.0*dimi,dimi)+1.0
iaux2 = mod(int((i-1.0)/dimi),2)
PBCREFR = iaux+(dimi-2.0*iaux+1.0)*float(iaux2)
end function
