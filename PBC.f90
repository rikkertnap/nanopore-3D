integer function PBCSYMI(i,dimi) ! returns the PBC cell coordinate 
integer i, dimi, p
p = abs(i) + 10
PBCSYMI = mod(i-1+p*dimi, dimi) + 1
end function

integer function PBCREFI(i,dimi) ! returns the reflection cell coordinate 
integer i, dimi, iaux, iaux2, p
p = abs(i) + 10
iaux = mod(i-1+p*dimi,dimi)+1
iaux2 = abs(mod(floor(float(i-1)/float(dimi)),2))
PBCREFI = iaux+(dimi-2*iaux+1)*iaux2
end function


real function PBCSYMR(i,dimi) ! returns the PBC cell coordinate 
real*8 i, dimi, p
p = abs(i) +  10.0
PBCSYMR = mod(i+p*dimi,dimi)
end function

real function PBCREFR(i,dimi) ! returns the reflection cell coordinate 
real*8 i, dimi, iaux, p
integer iaux2
p = abs(i) + 10.0
iaux = mod(i+p*dimi,dimi)
iaux2 = abs(mod(floor(i/dimi),2))
PBCREFR = iaux+(dimi-2.0*iaux)*float(iaux2)
end function
