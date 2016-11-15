IC=icc
TARGET= Poisson_Equation
OBJECT= Poisson_Equation.o Poisson_Solver.o Jacobi.o SOR.o Conjugate_Gradient.o

FCFLAGS= -O2 -L$(MPI_HOME)/lib -I$(MPI_HOME)/include
LDFLAGS= -lmpi

all : $(TARGET)
$(TARGET) : $(OBJECT)
	$(IC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

.SUFFIXES. : .o .c

%.o : %.c
	$(IC) $(FCFLAGS) -c $< $(LDFLAGS)

clean :
	rm -f *.o
	rm -f ./RESULT/*.plt
