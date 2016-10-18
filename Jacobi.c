#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "def.h"
double func(int i, int j, double dx, double dy);
void initialization(double **p);

void Jacobi(double **p,double dx, double dy, double tol,
                       double *tot_time,int *iter,int BC)
{
    int i,j,k,it;
    int Nx,Ny;
    double beta,rms;
    double SUM1,SUM2;
    double **p_new;
    time_t start_t =0, end_t =0;

    start_t = clock();
    beta = dx/dy;

    p_new = (double **) malloc(ROW *sizeof(double));
    for (i=0;i<ROW;i++)
    {
      p_new[i]      = (double *) malloc(COL * sizeof(double));
    }
    initialization(p_new);

    for (it=1;it<itmax;it++){
        SUM1 = 0;
        SUM2 = 0;

        for (i=1;i<ROW-1;i++){
            for (j=1;j<COL-1;j++){
            p_new[i][j] =  (p[i+1][j]+p[i-1][j]
                            + pow(beta,2) *(p[i][j+1]+p[i][j-1])
                            - dx*dx*func(i,j,dx,dy))/(2*(1+pow(beta,2)));
            }
        }

        //------------------------
        //  Boundary conditions
        //------------------------

        // Boundary - Case 1
        if (BC == 1){
          for (j=0;j<COL;j++){
              p_new[0][j] = 0;
              p_new[ROW-1][j] = 0;
          }

          for (i=0;i<ROW;i++) {
              p_new[i][0] = p_new[i][1];
              p_new[i][COL-1] = p_new[i][COL-2];
          }
        }
        // Boundary - Case 2
        else if (BC ==2){
          for (j=0;j<COL;j++){
              p_new[0][j] = -1/(2*pow(pi,2))*func(0,j,dx,dy);
              p_new[ROW-1][j] = -1/(2*pow(pi,2))*func(ROW-1,j,dx,dy);
          }

          for (i=0;i<ROW;i++) {
              p_new[i][0] = -1/(2*pow(pi,2))*func(i,0,dx,dy);
              p_new[i][COL-1] = -1/(2*pow(pi,2))*func(i,COL-1,dx,dy);
          }
        }

        //------------------------
        //  Convergence Criteria
        //------------------------
        for (i=1;i<ROW-1;i++){
            for (j=1;j<COL-1;j++){
                SUM1 += fabs(p_new[i][j]);
                SUM2 += fabs(p_new[i+1][j] + p_new[i-1][j]
                             + pow(beta,2)*(p_new[i][j+1] + p_new[i][j-1])
                             - (2+2*pow(beta,2))*p_new[i][j]-dx*dx*func(i,j,dx,dy));
            }
        }

        if ( SUM2/SUM1 < tol ){
            free(p_new);
            *iter = it;
            end_t = clock();
            *tot_time = (double)(end_t - start_t)/(CLOCKS_PER_SEC);
            break;
        }

        // printf("Iteration : %d, SUM1 : %f, SUM2 : %f, Ratio : %f \n",it,SUM1,SUM2,SUM2/SUM1);

        //------------------------
        //         Update
        //------------------------
        for (i=0;i<ROW;i++){
            for (j=0;j<COL;j++){
                p[i][j] = p_new[i][j];}}

    }
}
