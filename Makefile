IC=icc
TARGET= Poisson_Equation
OBJECT= Poisson_Equation.o Poisson_Solver.o Jacobi.o SOR.o Conjugate_Gradient.o
FCFLAGS= -O2 -openmp
LDFLAGS= -lrt

all : $(TARGET)
$(TARGET) : $(OBJECT)
	$(IC) $(FCFLAGS) -o $@ $^ -lrt

.SUFFIXES. : .o .c

%.o : %.c
	$(IC) $(FCFLAGS) -c $< -lrt

clean :
	rm -f *.o
	rm -f ./RESULT/*.plt
