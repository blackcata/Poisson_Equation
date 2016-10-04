#include <stdio.h>
#include <math.h>
#include <string.h>
#include "def.h"

//-----------------------------------
//      Matrix Calculations
//-----------------------------------
double norm_L2(double *a);
double vvdot(double *a, double *b);
void vmdot(double (*A)[ROW*COL],double *x,double *b);

void make_Abx(double (*A)[ROW*COL], double *b, double *x, double(*u)[COL]
              ,double dx, double dy);

//-----------------------------------
//      Mathematical functions
//-----------------------------------
double func(int i, int j, double dx, double dy);

void Conjugate_Gradient(double(*p)[COL],double dx, double dy, double tol)
{
    int i,j,k,it;
    double A[ROW*COL][ROW*COL] = {0};
    
    double x[ROW*COL] = {0};
    double b[ROW*COL] = {0};
    double r[ROW*COL] = {0};
    double z[ROW*COL] = {0};
    double tmp[ROW*COL] = {0};

}

void make_Abx(double (*A)[ROW*COL],double *b,double *x,
              double (*u)[COL],double dx, double dy)
{
    int i,j,k,l;
    //--------------------------------
    //         Make Matrix A
    //--------------------------------
    for (k=1;k<ROW;k++){
        for (l=1;l<COL;l++){
            if (k==l){
                for (i=0;i<ROW;i++){
                    if (i == 0){
                        A[COL*k+i][ROW*l+i]   = -4;
                        A[COL*k+i+1][ROW*l+i] = 1;
                    }
                    else if (i == ROW){
                        A[COL*k+i][ROW*l+i]   = -4;
                        A[COL*k+i-1][ROW*l+i] = 1;
                    }
                    else {
                        A[COL*k+i][ROW*l+i]   = -4;
                        A[COL*k+i-1][ROW*l+i] = 1;
                        A[COL*k+i+1][ROW*l+i] = 1;
                    }
                    
                }
            }
            else if ( abs(k-l) == 1 ){
                for (i=0;i<ROW;i++){
                    A[COL*k+i][ROW*l+i] = 1;
                }
            }
            else{
                for (i=0;i<ROW;i++){
                    for (j=0;j<COL;j++){
                        A[COL*k+i][ROW*l+j] = 0;
                    }
                }
            }
            
        }
    }

    
    //--------------------------------
    //         Make Vector x
    //--------------------------------
    for (i=0;i<ROW;i++){
        for (j=0;j<COL;j++){
            x[ROW*i+j] = u[i][j];
        }
    }
    //--------------------------------
    //        Make Vector b
    //--------------------------------
    for (i=0;i<ROW;i++){
        for (j=0;j<COL;j++){
            b[ROW*i+j] = func(i,j,dx,dy);
        }
    }
}

//------------------------------------------------------------
//              Matrix Calcuation Functions
//------------------------------------------------------------
double norm_L2(double *a)
{
    int i;
    double sum = 0;
    
    for (i=0;i<ROW*COL;i++){
        sum = sum + pow(a[i],2);
    }
    return sqrt(sum);
}

void vmdot(double (*A)[ROW*COL],double *x,double *b)
{
    int i,j;
    
    for (i=0;i<ROW*COL;i++){
        for (j=0;j<ROW*COL;j++){
            b[i] = b[i] + A[i][j]*x[j];
        }
        
    }
}

double vvdot(double *a, double *b)
{
    int i;
    double c = 0;
    
    for (i=1;i<ROW*COL;i++){
        c = c + a[i]*b[i];

    }

    return c;
}

