#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include "def.h"

//-----------------------------------
//      Matrix Calculations
//-----------------------------------
double norm_L2(int num, double *a);
double vvdot(int num, double *a, double *b);
void vmdot(int row,int col,double **A,double *x,double *b);

void make_Abx(double **A, double *b, double *x, double**u
              ,double dx, double dy);

//-----------------------------------
//      Mathematical functions
//-----------------------------------
double func(int i, int j, double dx, double dy);

void Conjugate_Gradient(double **p,double dx, double dy, double tol,
                                   double *tot_time,int *iter, int BC)
{
    int i,j,k,it;
    int nproc,myrank,ista,iend;
    double alpha,beta,ts,te ;

    double **A,**A_tmp;
    double *tmp, *x, *b, *z, *r, *r_new;

    time_t start_t = 0, end_t = 0;

    start_t = clock();
    ts = MPI_Wtime();

    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    printf("ROW X COL : %d, nproc : %d \n",ROW*COL,ROW*COL/nproc);
    
    ista = myrank*(ROW*COL/nproc);
    iend = (myrank+1)*(ROW*COL/nproc) - 1;
    printf("[ista,iend] : [%d,%d]\n \n",ista,iend );

    A = (double **) malloc(ROW*COL * sizeof(double));
    for (i=0;i<ROW*COL;i++) {
      A[i] = (double *) malloc(ROW*COL * sizeof(double));
    }
    tmp    = (double *) malloc(ROW*COL * sizeof(double));
    x      = (double *) malloc(ROW*COL * sizeof(double));
    b      = (double *) malloc(ROW*COL * sizeof(double));
    z      = (double *) malloc(ROW*COL * sizeof(double));
    r      = (double *) malloc(ROW*COL * sizeof(double));
    r_new  = (double *) malloc(ROW*COL * sizeof(double));

    make_Abx(A,b,x,p,dx,dy);
    vmdot(ROW*COL,ROW*COL,A,x,tmp);

   for (i=0;i<ROW;i++){
       for (j=0;j<COL;j++){
           r[COL*i+j] = b[COL*i+j] - tmp[COL*i+j];
           z[COL*i+j] = r[COL*i+j];
       }
   }

   //---------------------------------------
   //   Main Loop of Conjugate_Gradient
   //---------------------------------------
   for (it=0;it<itmax;it++)
   {
       vmdot(ROW*COL,ROW*COL,A,z,tmp);
       alpha = vvdot(ROW*COL,r,r)/vvdot(ROW*COL,z,tmp);


       for (i=0;i<ROW;i++){
           for (j=0;j<COL;j++){
               x[COL*i+j] = x[COL*i+j] + alpha * z[COL*i+j];
               r_new[COL*i+j] = r[COL*i+j] - alpha*tmp[COL*i+j];
           }
       }

       if (norm_L2(ROW*COL,r_new) < tol ){
          // printf("iteration : %d, tol : %f, value : %f\n",it,tol,norm_L2(r_new) );
          //---------------------------------------
          //   Redistribute x vector to array
          //---------------------------------------
          for (i=0;i<ROW;i++)
          {
            for (j=0;j<COL;j++)
            {
              p[i][j] = x[COL*i+j];
            }
          }
          *iter = it;
          free(A);
          free(tmp);
          free(x);
          free(b);
          free(z);
          free(r);
          free(r_new);

          end_t = clock();
          te = MPI_Wtime();
          *tot_time = (double)(end_t - start_t)/(CLOCKS_PER_SEC);
          if(myrank==0) printf("Total time is : %f s \n",te-ts );
          break;
       }

       beta = vvdot(ROW*COL,r_new,r_new)/vvdot(ROW*COL,r,r);
       for (i=0;i<ROW;i++){
           for (j=0;j<COL;j++){
               z[COL*i+j] = r_new[COL*i+j] + beta*z[COL*i+j];
               r[COL*i+j] = r_new[COL*i+j];
           }
       }
   }

}

//------------------------------------------------------------
//             Make Stiffness matrix of CG method
//------------------------------------------------------------
void make_Abx(double **A,double *b,double *x,
              double **u,double dx, double dy)
{
    int i,j,k,l;
    //--------------------------------
    //         Make Matrix A
    //--------------------------------
    for (k=0;k<ROW;k++){
        for (l=0;l<COL;l++){
            if (k==l){
                if (k==0 || k==ROW-1){
                  for (i=0;i<ROW;i++){
                      A[COL*k+i][ROW*l+i]   = 1;
                  }
                }
                else{
                  for (i=0;i<ROW;i++){
                      if (i == 0){
                          A[COL*k+i][ROW*l+i]   = -1;
                          A[COL*k+i+1][ROW*l+i] = 1;
                      }
                      else if (i == ROW-1){
                          A[COL*k+i][ROW*l+i]   = -1;
                          A[COL*k+i-1][ROW*l+i] = 1;
                      }
                      else {
                          A[COL*k+i][ROW*l+i]   = -4;
                          A[COL*k+i-1][ROW*l+i] = 1;
                          A[COL*k+i+1][ROW*l+i] = 1;
                      }
                  }
                }
            }

            else if ( abs(k-l) == 1 && k!=0 && k!=ROW-1){
                for (i=0;i<ROW;i++){
                  if (i==0 || i==ROW-1)
                    A[COL*k+i][ROW*l+i] = 0;
                  else
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
            // printf("i: %d, j :  %d \n",k,l);
            // for (j=0;j<ROW;j++){
            //   printf("%d ",i);
            //   for (i=0;i<COL;i++){
            //       printf("%f ",A[COL*k+i][ROW*l+j]);
            //   }
            //   printf("\n");
            // }
            // printf("\n");
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
          if (i==0 || i==ROW-1 || j==0 || j==COL-1)
            b[ROW*i+j] = 0;//1/(2*pow(pi,2))*func(i,j,dx,dy);
          else
              b[ROW*i+j] = dx*dx*func(i,j,dx,dy);

          // printf(" i:%d, j:%d, b[i][j] : %f\n",i,j,b[ROW*i+j] );
        }
    }
}

//------------------------------------------------------------
//              Matrix Calcuation Functions
//------------------------------------------------------------
double norm_L2(int num, double *a)
{
    int i;
    double sum = 0;

    for (i=0;i<num;i++){
        sum = sum + pow(a[i],2);
    }
    return sqrt(sum);
}

void vmdot(int row,int col,double **A,double *x,double *b)
{
    int i,j;

    for (i=0;i<row;i++){
            b[i] = 0;
    }

    for (i=0;i<row;i++){
        for (j=0;j<col;j++){
            b[i] = b[i] + A[i][j]*x[j];
        }

    }
}

double vvdot(int num, double *a, double *b)
{
    int i;
    double c = 0;

    for (i=0;i<num;i++){
        c = c + a[i]*b[i];
    }

    return c;
}
