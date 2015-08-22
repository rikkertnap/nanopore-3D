integer function PBCSYMI(i,dimi) ! returns the PBC cell coordinate 
integer i, dimi
PBCSYMI = mod(i-1+5*dimi, dimi) + 1
end function

integer function PBCREFI(i,dimi) ! returns the reflection cell coordinate 
integer i, dimi, iaux, iaux2
iaux = mod(i-1+5*dimi,dimi)+1
iaux2 = mod(int(float(i-1)/float(dimi)),2)
PBCREFI = iaux+(dimi-2*iaux+1)*iaux2
end function


real function PBCSYMR(i,dimi) ! returns the PBC cell coordinate 
real*8 i, dimi
PBCSYMR = mod(i+5.0*dimi,dimi)
end function

real function PBCREFR(i,dimi) ! returns the reflection cell coordinate 
real*8 i, dimi, iaux
integer iaux2
iaux = mod(i+5.0*dimi,dimi)
iaux2 = mod(floor(i/dimi),2)
PBCREFR = iaux+(dimi-2.0*iaux)*float(iaux2)
end function
