main: main.o modd.o
	gfortran main.o modd.o -o main
main.o: main.f95 modd.mod
	gfortran -c main.f95
modd.mod modd.o: modd.f95
	gfortran -c modd.f95
clear:
	rm -f *.o *.mod main
result: main
	./main
