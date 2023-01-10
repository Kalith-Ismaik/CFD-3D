# The makefile compiles (type make) and makes an executable called mse-bs

EXEPATH = .

EXE = $(EXEPATH)/dem.x

F90 = gfortran

OBJS = Main.o CFD-init.o DiscreteCFD.o Navier-stokes.o Poisson.o CBound.o Swrite.o

# compile and load
default:
	@echo " "
	@echo "Compiling Code MD"
	@echo "Version 2.0"
	@echo "FORTRAN 90"
	$(MAKE) $(EXE)

$(EXE):	$(OBJS)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $(EXE)  $(OBJS)

.SUFFIXES: .f90 .o
.f90.o:
	$(F90) $(F90FLAGS) -c $*.f90


clean:
	rm -f *.o
	rm -f dem.x