OBJ=main.o pythia-6.4.25.o \
	hadronization.o
GF=gfortran
CLEAR=clear
task1:$(OBJ)
	$(GF) -o main -O $(OBJ)
main.o: main.f
	$(GF) -c -O  main.f
hadronization.o:hadronization.f
	$(GF) -c hadronization.f
#if you don't need to compile pythia, comment below two lines...
pythia-6.4.25.o: pythia-6.4.25.f
	$(GF) -c pythia-6.4.25.f
clean: 
	$(RM) ana*.dat main.o hadronization.o main *~ |$(CLEAR)

