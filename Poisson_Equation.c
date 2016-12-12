//
//
//  Programming for 2D Poisson Equation
//                    By. Kyung Min Noh 2016
//
//  Solving 2D Poisson Equation
//
//  Laplce(u(x,y)) = f(x,y) for (x,y) in domain
//  on the boundary boundary_D and boundary_N
//
//  u(x,y) = g(x,y) on boundary_D
//  du/dn  = h(x,y) on boundary N
//
//
//  Domian
//      [x,y] is in [0,1] X [0,1]
//      f(x,y) = sin(pi*x) * cos(pi *y)
//      Analytic solution is
//          u^a(x,y) = -1/(2*pi^2) * sin(pi*x) * cos(pi*y)
//
//  Boundary Condition
//      Case 1
//          Diriclet bondary condition u(x,y) = 0 in x direction q
//          Neumann boundary condition du/dn = 0 in y direction
//
//      Case 2
//          Dirichlet boundary conditioins using analytical solution
//          both x and y directions.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "def.h"

//----------------------------------------------------------------------------//
//                             Productivity tools                             //
//----------------------------------------------------------------------------//
void poisson_solver(double **u, double **u_anal, double tol, double omega,
                    int BC, int method, int mpi_xsize, int mpi_ysize,
                    int write_type,char *dir_name);

//----------------------------------------------------------------------------//
//                                                                            //
//                          Main Poisson Equation                             //
//                                                                            //
//----------------------------------------------------------------------------//
int main(int argc, char* argv[])
{
    double **u;
    double **u_anal;

    char *dir_name ;

    int i, method, BC, myrank, nprocs, write_type, mpi_xsize, mpi_ysize;
    double tol, omega;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    //------------------------------------------------------------------------//
    //                           Memory allocation                            //
    //------------------------------------------------------------------------//
    u      = (double **) malloc(ROW *sizeof(double));
    u_anal = (double **) malloc(ROW *sizeof(double));

    for (i=0;i<ROW;i++)
    {
      u[i]      = (double *) malloc(COL * sizeof(double));
      u_anal[i] = (double *) malloc(COL * sizeof(double));
    }

    //------------------------------------------------------------------------//
    //                             Initial setting                            //
    //------------------------------------------------------------------------//
    tol = 1e-6;
    omega = 1.9;
    dir_name = "./RESULT/";

    //------------------------------------------------------------------------//
    //                        Poisson Solver Type                             //
    //                                                                        //
    // BC         = 1 : Boundary condition Case 1                             //
    //            = 2 : Boundary condition Case 2                             //
    //                                                                        //
    // method     = 1 : Jacobi method                                         //
    //            = 2 : SOR method                                            //
    //            = 3 : Conjugate Gradient method                             //
    //                                                                        //
    // write_type = 1 : Post Reassembly                                       //
    //            = 2 : Single task I/O                                       //
    //            = 3 : MPI I/O                                               //
    //------------------------------------------------------------------------//
    BC         = 1;
    method     = 1;
    write_type = 1;

    if (method != 3){
      if (argc < 3) {
        if (myrank==0){
        printf("%s\n","----------------------------------------------------------");
        printf("%s\n","                 <<< NO-INPUT WARNING >>>                 ");
        printf("%s\n","       You have to determine (mpi_xsize),(mpi_ysize)      ");
        printf("%s\n","       ex) mpirun -np=4 Poisson_Equation 2 2              ");
        printf("%s\n","                                                          ");
        printf("%s\n","         PROGRAM POISSON IS FORCIBLYY TERMINATED          ");
        printf("%s\n","----------------------------------------------------------");
        }
        exit(1);
      }
      mpi_xsize  = atof(argv[1]);
      mpi_ysize  = atof(argv[2]);

      if (nprocs != atoi(argv[1])*atoi(argv[2])) {
        if (myrank ==0) {
          printf("%s\n","----------------------------------------------------------");
          printf("%s\n","                 <<< MISMATCH ERROR >>>                   ");
          printf("%s\n","The number of cores and (mpi_xsize) X (mpi_ysize) mismatch");
          printf("Nprocs : %d, (mpi_xsize) X (mpi_ysize) : %d X %d = %d \n",
                      nprocs,atoi(argv[1]),atoi(argv[2]),atoi(argv[1])*atoi(argv[2]));
          printf("%s\n","                                                          ");
          printf("%s\n","         PROGRAM POISSON IS FORCIBLYY TERMINATED          ");
          printf("%s\n","----------------------------------------------------------");
        }
        exit(1);
      }
    }

    if (myrank == 0 ){
      int make_fold= system("mkdir RESULT");
      printf("\n");
      printf("%s\n","----------------------------------------------------------");
      printf("%s\n","             Poisson Equation Initial Setting             ");
      printf("Nx : %d, Ny : %d\n",ROW,COL);
      printf("Tolerance : %f, Omega : %f \n",tol, omega);
      switch (method) {
        case 1 : printf("Method : Jacobi method \n");
          break;

        case 2 :
          printf("Method : SOR method \n");
          break;

        case 3:
          printf("Method : Conjugate Gradient method \n");
          break;
      }
      printf("%s\n","----------------------------------------------------------");
      printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    poisson_solver(u,u_anal,tol,omega,BC,method,mpi_xsize,mpi_ysize,
                                                           write_type,dir_name);
    MPI_Finalize();

    free(u);
    free(u_anal);

    return 0;
}
