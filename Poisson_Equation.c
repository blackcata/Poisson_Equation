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

//-----------------------------------
//        Productivity tools
//-----------------------------------

void poisson_solver(double **u, double **u_anal, double tol, double omega,
                    int BC, int method, char *dir_name);

int main(int argc, char* argv[])
{
    double **u;
    double **u_anal;

    char *dir_name ;

    int i, method, BC, myrank;
    double tol, omega;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    // --------------------------------------------------------
    //                    Memory allocation
    // --------------------------------------------------------
    u      = (double **) malloc(ROW *sizeof(double));
    u_anal = (double **) malloc(ROW *sizeof(double));

    for (i=0;i<ROW;i++)
    {
      u[i]      = (double *) malloc(COL * sizeof(double));
      u_anal[i] = (double *) malloc(COL * sizeof(double));
    }

    //--------------------
    //   Initial setting
    //--------------------
    tol = 1e-6;
    omega = 1.8;
    dir_name = "./RESULT/";

    if (myrank == 0 ){
      int make_fold= system("mkdir RESULT");
      printf("\n");
      printf("---------------------------------------- \n");
      printf("Nx : %d, Ny : %d\n",ROW,COL);
      printf("Tolerance : %f, Omega : %f \n",tol, omega);
      printf("---------------------------------------- \n");
      printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //----------------------------------------
    //       Poisson Solver Type
    //
    // BC = 1 : Boundary condition Case 1
    //    = 2 : Boundary condition Case 2
    //
    // method = 1 : Jacobi method
    //        = 2 : SOR method
    //        = 3 : Conjugate Gradient method
    //----------------------------------------
    BC = 1;
    method = 3;


    poisson_solver(u,u_anal,tol,omega,BC,method,dir_name);
    MPI_Finalize();

    free(u);
    free(u_anal);


    return 0;
}
