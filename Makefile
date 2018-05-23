FC=gfortran -c
LD=gfortran
SRC=fact.f90 chgm.f90 main.f90
OBJ=fact.o chgm.o main.o
energy:
	$(FC) $(SRC)
	$(LD) $(OBJ) -o function.x
	rm -rf *.o
clean:
	rm -rf *.o *.x *.pdf

