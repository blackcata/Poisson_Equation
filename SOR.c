#include <stdio.h>
#include <math.h>
#include "def.h"
double func(int i, int j, double dx, double dy);

void SOR(double(*p)[COL],double dx, double dy, double tol, double omega)
{
    int i,j,k,it;
    int Nx,Ny;
    double beta,rms;
    double SUM1,SUM2;
    double p_new[ROW][COL]={0};
    printf(" %d %d \n",COL,ROW);
    
    beta = dx/dy;

    for (it=1;it<itmax;it++){
        SUM1 = 0;
        SUM2 = 0;
        
        
        for (i=1;i<ROW-1;i++){
            for (j=1;j<COL-1;j++){
                p_new[i][j] =  (p[i+1][j]+p_new[i-1][j]
                                + pow(beta,2) *(p[i][j+1]+p_new[i][j-1])
                                - dx*dx*func(i,j,dx,dy))/(2*(1+pow(beta,2)));
                p_new[i][j] = p[i][j] + omega * (p_new[i][j] - p[i][j]);
            }
        }
        
        //------------------------
        //  Boundary conditions
        //------------------------
        for (j=0;j<COL;j++){
            p_new[0][j] = 0;
            p_new[ROW-1][j] = 0;
        }
        
        for (i=0;i<ROW;i++) {
            p_new[i][0] = p_new[i][1];
            p_new[i][COL-1] = p_new[i][COL-2];
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
        
        if ( SUM2/SUM1 < tol )
            break;
        
        //------------------------
        //         Update
        //------------------------
        for (i=0;i<ROW;i++){
            for (j=0;j<COL;j++){
                p[i][j] = p_new[i][j];}}


        
        printf("Iteration : %d, SUM1 : %f, SUM2 : %f, Ratio : %f \n",it,SUM1,SUM2,SUM2/SUM1);
    }
}

