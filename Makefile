IC=icc
TARGET= Poisson_Equation
OBJECT= Poisson_Equation.o Poisson_Solver.o Jacobi.o SOR.o Conjugate_Gradient.o

FCFLAGS= -O2 -I${MKL_HOME}/include -L${MKL_HOME}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -openmp
LDFLAGS= -lrt

all : $(TARGET)
$(TARGET) : $(OBJECT)
	$(IC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

.SUFFIXES. : .o .c

%.o : %.c
	$(IC) $(FCFLAGS) -c $< $(LDFLAGS)

clean :
	rm -f *.o
	rm -f ./RESULT/*.plt
