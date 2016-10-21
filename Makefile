IC=icc
TARGET= Poisson_Equation
OBJECT= Poisson_Equation.o Poisson_Solver.o Jacobi.o SOR.o Conjugate_Gradient.o
FCFLAGS= -O2 -openmp

all : $(TARGET)
$(TARGET) : $(OBJECT)
	$(IC) $(FCFLAGS) -o $@ $^

.SUFFIXES. : .o .c

%.o : %.c
	$(IC) $(FCFLAGS) -c $<

clean :
	rm -f *.o
	rm -f ./RESULT/*.plt
