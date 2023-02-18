CC=icc
F90=ifort
F77=ifort
#

D0=-Dtype1h -DFULLINEXP
DEFS?= -DC12PT1 -DC22PT1 -DADDC3 -DC23PT1 -DC33PT1 -DADC2 -DC11PT2 -DC12PT2 -DC22PT2 #ADC(2,2)
TARG?=main.exe
#
SINT?=-i4
#
ifeq ($(SINT),-i8)
	MKLINT=-lmkl_intel_ilp64
else
	MKLINT=-lmkl_intel_lp64
endif
#
OPTS=-xHost -O2 -qopenmp -qopt-prefetch -assume buffered_io $(DEFS) $(SINT) -fpp
COPTS=-xHost -O2 -qopenmp -qopt-prefetch $(DEFS)
#COPTS= -g -O0 -qopenmp $(DEFS)
#
LIBGBTO=/home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/
#
INC= \
    -I${LIBGBTO}/mod
#
LIB= \
    $(MKLINT) -lmkl_intel_thread -lmkl_core -qopenmp -lpthread \
    -L${LIBGBTO}/lib -lGBTO #\
#
OBJECTS= globalmod.o phismod.o main.o
#
all: $(TARG) 

globalmod.o: globalmod.f90 Makefile
phismod.o: $(LIBGBTO)/lib/libGBTO.a globalmod.o phismod.f90 Makefile
main.o: globalmod.o phismod.o $(LIBGBTO)/lib/libGBTO.a main.f90 Makefile

%.o: %.f90
	$(F90) $(OPTS) $(INC) -c $<

%.o: %.f
	$(F77) $(OPTS) $(INC) -c $<

%.o: %.c
	$(CC) $(COPTS) $(INC) -c $<

$(TARG): $(OBJECTS) Makefile
	$(F90) $(OPTS) $(D0) $(DEFS) $(OBJECTS) $(LIB) -o $(TARG)

clean:
	rm -f *.o *.mod *__genmod.f90

distclean:
	rm -f *.o *.mod *__genmod.f90 *.exe *~
