FC=gfortran
FOPT=-Og -g -Wall -pedantic -fbounds-check -fbacktrace -fimplicit-none -static -static-libgfortran -static-libgcc -static-libstdc++
RM=del

SRC_utils   = types.f90 utils.f90 lapack.f90 linalg.f90 
SRC_modules = modinit.f90 moduniinv.f90 modfunc.f90 modcore.f90
SRC_main    = subroutines.f90 singlecluster.f90
DLL         = lapack.dll
SRC = $(SRC_utils) $(SRC_modules) $(SRC_main)
OBJ = $(SRC:.f90=.o)
EXE = singlecluster.exe

%.o: %.f90
	$(FC) $(FOPT) -o $@ -c $<

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FOPT) $(OBJ) -o $@ $(DLL)

clean:
	-$(RM) -f *.o *.mod *~	