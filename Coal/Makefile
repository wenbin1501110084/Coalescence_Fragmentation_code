OBJ=main.o \
	hadronization.o
GF=gfortran
CLEAR=clear
task1:$(OBJ)
	$(GF) -o main -O $(OBJ)
main.o: main.f
	$(GF) -c -O  main.f
hadronization.o:hadronization.f
	$(GF) -c hadronization.f
clean: 
	$(RM) ana*.dat main.o hadronization.o main *~ |$(CLEAR)

