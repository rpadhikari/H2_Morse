FC=gfortran -c
LD=gfortran
SRC=fact.f90 chgm.f90 lnkx.f90 psi.f90 main.f90
OBJ=fact.o chgm.o lnkx.o psi.o main.o
energy:
	$(FC) $(SRC)
	$(LD) $(OBJ) -o function.x -llapack -lblas
	rm -rf *.o
clean:
	rm -rf *.o *.x *.pdf

