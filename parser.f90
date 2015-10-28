
subroutine readinput
use molecules
use const
use bulk
use MPI
use ellipsoid
use chainsdat
use inputtemp
use transform
use system
use kai
use kaist
use s2d
implicit none

! Input related variables
character (len=100)  buffer,label
integer pos
integer, parameter :: fh = 15
integer ios
integer line, linemax
integer i, j
character(len=50) :: filename = 'DEFINITIONS.txt'
character basura
integer ndi
real*8 ndr

! not defined variables, change if any variable can take the value

PBC = 1

ndi = -1e5
ndr = -1.0d10

verbose = 5

electroflag = 1 ! system with electrostatics?

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Check validity of input
!

scx = ndi
scy = ndi
scz = ndi
vtkflag = ndi
kaptype = ndi
dimx = ndi
dimy = ndi
dimz = ndi
long = ndi
cuantas = ndi
readchains = ndi
infile = ndi
randominput = ndi
cutoff = ndr
lseg = ndr
nst = ndi
dielS = ndr
pHbulk = ndr
dielP = ndr
delta = ndr
dx = ndr
dy = ndr
dz = ndr
cdiva = ndr
csalt = ndr
zpol = ndr
pKa = ndr
vpol0 = ndr
vsol0 = ndr
gama0 = ndr

nsc = 1
scs(1) = 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Control file variables

line = 0
ios = 0

open(fh, file=filename)

if(rank.eq.0)print*, 'parser:', 'Reading parameters from ', filename

! ios is negative  if an end of record condition is encountered or if
! an endfile condition was detected.  It is positive  if an error was
! detected.  ios is zero otherwise.

do while (ios == 0)

 read(fh, '(A)', iostat=ios) buffer
 if (ios == 0) then
 line = line + 1

! Find the first instance of whitespace.  Split label and data.

 pos = scan(buffer, ' ')

 label = buffer(1:pos)
 buffer = buffer(pos+1:)


 select case (label)

 case ('PBC')
   read(buffer, *, iostat=ios) PBC(1),PBC(2),PBC(3),PBC(4),PBC(5),PBC(6)
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

   do j = 1,5,2
    if((PBC(j).eq.1).and.(PBC(j+1).ne.1)) then 
      print*, 'parser:', 'Error in PBC'
      stop
    endif
    if((PBC(j+1).eq.1).and.(PBC(j).ne.1)) then
      print*, 'parser:', 'Error in PBC'
      stop
    endif
   enddo

 case ('verbose')
   read(buffer, *, iostat=ios) verbose
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('vtkflag')
   read(buffer, *, iostat=ios) vtkflag
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('electroflag')
   read(buffer, *, iostat=ios) electroflag
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('randominput')
   read(buffer, *, iostat=ios) randominput
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('readchains')
   read(buffer, *, iostat=ios) readchains
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dimx')
   read(buffer, *, iostat=ios) dimx
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('scx')
   read(buffer, *, iostat=ios) scx
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('scy')
   read(buffer, *, iostat=ios) scy
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('scz')
   read(buffer, *, iostat=ios) scz
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)


 case ('delta')
   read(buffer, *, iostat=ios) delta
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dx')
   read(buffer, *, iostat=ios) dx
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dy')
   read(buffer, *, iostat=ios) dy
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dz')
   read(buffer, *, iostat=ios) dz
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('cdiva')
   read(buffer, *, iostat=ios) cdiva
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)


 case ('dimy')
   read(buffer, *, iostat=ios) dimy
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dimz')
   read(buffer, *, iostat=ios) dimz
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('long')
   read(buffer, *, iostat=ios) long
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('cuantas')
   read(buffer, *, iostat=ios) cuantas
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('lseg')
   read(buffer, *, iostat=ios) lseg
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)


 case ('dielP')
   read(buffer, *, iostat=ios) dielP
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dielS')
   read(buffer, *, iostat=ios) dielS
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('csalt')
   read(buffer, *, iostat=ios) csalt
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('zpol')
   read(buffer, *, iostat=ios) zpol
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('vsol')
   read(buffer, *, iostat=ios) vsol0
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('vpol')
   read(buffer, *, iostat=ios) vpol0
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('pKa')
   read(buffer, *, iostat=ios) pKa
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('gama')
   read(buffer, *, iostat=ios) gama0
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('pHbulk')
   read(buffer, *, iostat=ios) pHbulk
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('infile')
   read(buffer, *, iostat=ios) infile
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('nst')
   read(buffer, *, iostat=ios) nst
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)
  
   do i = 1, nst
   read(fh,*)sts(i)
   enddo 


 case ('nsc')
   read(buffer, *, iostat=ios) nsc
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)
  
   do i = 1, nsc
   read(fh,*)scs(i)
   enddo 



 case ('Xucutoff')
   read(buffer, *, iostat=ios) cutoff
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('kaptype')
   read(buffer, *, iostat=ios) kaptype
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

   select case (kaptype)
    case(1) 
     read(fh, *) basura
     read(fh, *)NNN

     if(NNN.ne.0) then

     call allocateell
     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), Rellf(1,j), Rellf(2,j), Rellf(3,j)
     if(rank.eq.0)print*, 'parser:','Set particle',j,'pos to',  Rellf(1,j), Rellf(2,j), Rellf(3,j)
     enddo
     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), Aell(1,j), Aell(2,j), Aell(3,j)
     if(rank.eq.0)print*, 'parser:','Set particle',j,'axis to',  Aell(1,j), Aell(2,j), Aell(3,j)
     enddo
     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
     read(fh, *), rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
     read(fh, *), rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
     if(rank.eq.0) then
         print*, 'parser:','Set particle',j,'rotation to:'
         print*, 'parser:', rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
         print*, 'parser:', rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
         print*, 'parser:', rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
     endif
     enddo

     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), sigma(j)
     if(rank.eq.0)print*, 'parser:','Set particle',j,'surface coverage to', sigma(j)
     enddo

     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), echarge(j)
     if(rank.eq.0)print*, 'parser:','Set particle',j,'charge to', echarge(j)
     enddo
     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), eeps(j)
     if(rank.eq.0)print*, 'parser:','Set particle',j,'hydrophobicity to', eeps(j)
     enddo

     endif ! NNN
endselect
endselect

endif

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Check validity of input
! 


if(vtkflag.eq.ndi)call stopundef('vtkflag')
if(dimx.eq.ndi)call stopundef('dimx')
if(scx.eq.ndi)call stopundef('scx')
if(scy.eq.ndi)call stopundef('scy')
if(scz.eq.ndi)call stopundef('scz')
if(dimy.eq.ndi)call stopundef('dimy')
if(dimz.eq.ndi)call stopundef('dimz')
if(ncha.eq.ndi)call stopundef('ncha')
if(long.eq.ndi)call stopundef('long')
if(cuantas.eq.ndi)call stopundef('cuantas')
if(infile.eq.ndi)call stopundef('infile')
if(randominput.eq.ndi)call stopundef('randominput')
if(cutoff.eq.ndr)call stopundef('Xucutoff')
if(readchains.eq.ndi)call stopundef('readchains')
if(kaptype.eq.ndi)call stopundef('kaptype')
if(nst.eq.ndi)call stopundef('nst')

if(delta.eq.ndr)call stopundef('delta')
if(dx.eq.ndr)call stopundef('dx')
if(dy.eq.ndr)call stopundef('dy')
if(dz.eq.ndr)call stopundef('dz')
if(cdiva.eq.ndr)call stopundef('cdiva')
if(dielS.eq.ndr)call stopundef('dielS')
if(dielP.eq.ndr)call stopundef('dielP')
if(lseg.eq.ndr)call stopundef('lseg')
if(csalt.eq.ndr)call stopundef('csalt')
if(pHbulk.eq.ndr)call stopundef('pHbulk')
if(zpol.eq.ndr)call stopundef('zpol')
if(vpol0.eq.ndr)call stopundef('vpol')
if(vsol0.eq.ndr)call stopundef('vsol')
if(gama0.eq.ndr)call stopundef('gama')
if(pKa.eq.ndr)call stopundef('pKa')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine

subroutine stopundef(namevar)
character(len=*) :: namevar
print*, 'parser:', 'Variable ', namevar, ' is undefined '
call MPI_FINALIZE(ierr) ! finaliza MPI
stop
end

