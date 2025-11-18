TARGET = nanopore-mixture

SRC = modules.f90 maps.f90  moleculelist.f90 flux.f90  read-graftpts.f90 channel-curved.f90 SPmain.f90 channel.f90 PBC.f90 parser.f90 init.f90 allocation.f90 allocatencha.f90 allocateell.f90 3D.f90  allocatecpp.f90  cadenas.f90 cadenas_b.f90 cadenas_b2.f90  fe.f90  fkfun.f90  kai.f90 anderson.f90 kinsol.f90 solver.f90  pxs.f90  savetodisk.f90 rands.f90 ellipsoid.f90 dielectric.f90 transform.f90 testsystem.f90 testsystemc.f90 testsystemr.f90 monomers.definitions.f90 chains.definitions.f90 channel-part.f90 

# STANDARD FLAGS
LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib

ifeq ($(MKL),true)
   SRC += csr.f90 csr_map.f90
   MKLCOMMENTED = -D_MKL
endif

HOST=$(shell hostname)
$(info HOST is ${HOST})

ifeq ($(HOST),skay)
LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so -L/usr/local/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial      -Wl,-rpath,/usr/local/lib

# MKL
ifeq ($(MKL),true)
LFLAGS += -I/opt/intel/mkl/include -w -fcray-pointer -cpp -Wl,--no-as-needed -lpthread -ldl  -Wl,--start-group "/opt/intel/mkl/lib/intel64"/libmkl_intel_lp64.a "/opt/intel/mkl/lib/intel64"/libmkl_core.a "/opt/intel/mkl/lib/intel64"/libmkl_gnu_thread.a -Wl,--end-group -L "/opt/intel/mkl/../compiler/lib/intel64" -liomp5 -lm
endif

endif

ifeq ($(HOST),carapa)
LFLAGS = -lsundials_fkinsol -lsundials_fnvecserial -lsundials_kinsol -lsundials_nvecserial -lm
endif

ifeq ($(HOST),juno)
LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib
endif

# some definitions
SHELL = /bin/bash
FFLAGS= -O3 # -fbacktrace -fbounds-check # -O3

ifeq ($(HOST),piluso.rosario-conicet.gov.ar)
LFLAGS = -L/home/mtagliazucchi.inquimae/software/kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s
endif


ifeq ($(HOST),ensalada)
LDFLAGS=-lm /usr/lib/x86_64-linux-gnu/librt.so -L/opt/local/sundials-2.6.1-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial     -Wl,-rpath,/opt/local/sundials-2.6.1-openmpi/lib -llapack -lblas

LFLAGS=-lm /usr/lib/x86_64-linux-gnu/librt.so -L/opt/local/sundials-2.6.1-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial     -Wl,-rpath,/opt/local/sundials-2.6.1-openmpi/lib -llapack -lblas
endif

ifeq ($(HOST),mccbmec02cc21mlvdm)
LDFLAGS=
LFLAGS= -lm  -framework Accelerate
endif


ifeq ($(HOST),pear)
LFLAGS=-L/home/mario/software/KINSOL/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux
endif

ifeq ($(HOST),cnode01)
LFLAGS = -L/home/mtagliazucchi/software/Kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s
endif

ifeq ($(HOST),master) 
LFLAGS = -L/shared/software/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath
LFLAGS += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl


endif

ifeq ($(HOST),mate.bme.northwestern.edu) 
LFLAGS = -L/home/mario/software/kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath
endif


ifeq ($(HOST),quser40)
FFLAGS=  -O3 #-fcheck=all -fbounds-check -Warray-bounds -g -fbacktrace # -Wargument-mismatch -Wpedantic #-Wall

LFLAGS= -lm /usr/lib64/librt.so -L/projects/p31445/sundials/sundials-2.6.1-openmpi-gfortran84/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial     -Wl,-rpath,/projects/p31445/sundials/sundials-2.6.1-openmpi-gfortran84/lib -L/software/lapack/3.10.1/lib64  -llapack
endif

ifeq ($(HOST),quser41)
FFLAGS=  -O3 #-fcheck=all -fbounds-check -Warray-bounds -g -fbacktrace # -Wargument-mismatch -Wpedantic #-Wall

LFLAGS= -lm /usr/lib64/librt.so -L/projects/p31445/sundials/sundials-2.6.1-openmpi-gfortran84/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial     -Wl,-rpath,/projects/p31445/sundials/sundials-2.6.1-openmpi-gfortran84/lib -L/software/lapack/3.10.1/lib64  -llapack
endif

ifeq ($(HOST),quser42)
FFLAGS=  -O3 #-fcheck=all -fbounds-check -Warray-bounds -g -fbacktrace # -Wargument-mismatch -Wpedantic #-Wall

LFLAGS= -lm /usr/lib64/librt.so -L/projects/p31445/sundials/sundials-2.6.1-openmpi-gfortran84/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial     -Wl,-rpath,/projects/p31445/sundials/sundials-2.6.1-openmpi-gfortran84/lib -L/software/lapack/3.10.1/lib64  -llapack
endif

ifeq ($(HOST),quser43)
FFLAGS=  -O3 #-fcheck=all -fbounds-check -Warray-bounds -g -fbacktrace # -Wargument-mismatch -Wpedantic #-Wall

LFLAGS= -lm /usr/lib64/librt.so -L/projects/p31445/sundials/sundials-2.6.1-openmpi-gfortran84/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial     -Wl,-rpath,/projects/p31445/sundials/sundials-2.6.1-openmpi-gfortran84/lib -L/software/lapack/3.10.1/lib64  -llapack
endif


ifeq ($(HOST),quser44) 
FFLAGS=  -O3 #-fcheck=all -fbounds-check -Warray-bounds -g -fbacktrace # -Wargument-mismatch -Wpedantic #-Wall

LFLAGS= -lm /usr/lib64/librt.so -L/projects/p31445/sundials/sundials-2.6.1-openmpi-gfortran84/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial     -Wl,-rpath,/projects/p31445/sundials/sundials-2.6.1-openmpi-gfortran84/lib -L/software/lapack/3.10.1/lib64  -llapack
endif


GIT_VERSION := $(shell git describe --abbrev=6 --dirty --always --tags)
GFLAGS=-cpp -D_VERSION=\"$(GIT_VERSION)\" $(MKLCOMMENTED)

FF = mpif90 #${F90}
VER = ~/bin

all:	$(TARGET)

$(TARGET): $(SRC:.f90=.o)
	$(FF) $(GFLAGS) -o $(TARGET) $(SRC:.f90=.o) $(LFLAGS)
	cp $(TARGET) $(VER)

$(SRC:.f90=.o): $(SRC)
	${FF} $(GFLAGS) -c ${FFLAGS}  $(SRC) $(LFLAGS)

install: all
	cp $(TARGET) $(VER)

clean:	
	@rm -f $(SRC:.f90=.o) $(SRC:.f90=.d) $(TARGET) *~ *.mod

realclean: clean
	@rm -f .depend

depend dep:
	@$(FF)  $(CFLAGS) -MM $(SRC) > .depend 

ifeq (.depend, $(wildcard .depend))
include .depend
endif

















































