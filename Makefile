IC=icc
TARGET= Poisson_Equation
OBJECT= Poisson_Equation.o Jacobi.c

all : $(TARGET)
$(TARGET) : $(OBJECT)
	$(IC) -o $@ $^

.SUFFIXES. : .o .c

%.o : %.c
	$(IC) -c $<

clean :
	rm -f *.o
