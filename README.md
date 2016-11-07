This program is used for solving Poisson Equation.

###1.
Poisson Equation is first solved by 3 methods.
  - Jacobi method
  - SOR method
  - Conjugate Gradient method
  
###2.
CSR(Compressed Sparse Row) matrix is used to CG method
  - At new 'CSR' branch
  
###3.
OpenMP technique is applied to each method except SOR
  - At new branch 'OpenMP' branch
    - Jacobi method is parallelized
    - Conjugate gradient method is parallelized
    
  - At 'CSR' branch
    - MKL libraries are used 
      - to calculate matrix vector product, vector vector product, norm2
      - and also tested with changing OPEN_NUM_THREADS
      
###4.
MPI technique is applied to each method
  - At new branch 'MPI' branch