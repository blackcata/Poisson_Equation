IC=icc
TARGET= Poisson_Equation
OBJECT= Poisson_Equation.o Poisson_Solver.o Jacobi.o SOR.o Conjugate_Gradient.o
FCFLAGS= -O2 -qopenmp
#LDFLAGS= -lrt
LDFLAGS= 
# LDLFAGS have to be determined to case
# At linux it must has -lrt option
# But in mac it must remove -lrt option


all : $(TARGET)
$(TARGET) : $(OBJECT)
	$(IC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

.SUFFIXES. : .o .c

%.o : %.c
	$(IC) $(FCFLAGS) -c $< $(LDFLAGS)

clean :
	rm -f *.o
	rm -f ./RESULT/*.plt
