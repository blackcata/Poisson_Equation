#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mkl.h>

#include "def.h"

//-----------------------------------
//      Matrix Calculations
//-----------------------------------
double norm_L2(double *a);
double vvdot(double *a, double *b);
void vmdot(double *nnzeros, int *col_ind,int *row_ptr,
           double *x,double *b);

void make_Abx(double *nnzeros, int *col_ind,int *row_ptr,
              double *b, double *x, double **u,
              int nnz,double dx, double dy);

//-----------------------------------
//      Mathematical functions
//-----------------------------------
double func(int i, int j, double dx, double dy);

void Conjugate_Gradient(double **p,double dx, double dy, double tol,
                                   double *tot_time,int *iter, int BC)
{
    int i,j,it,nnz;
    double alpha, beta, alp, bet, diff;
    MKL_INT m,k;

    int *col_ind, *row_ptr, *r_ptr_b, *r_ptr_e;
    double *nnzeros,*tmp,*x,*b,*z,*r,*r_new ;
    char *transa, *matdescra ;
    struct timespec start,end ;

    time_t start_t = 0, end_t = 0;

    start_t = clock();
    clock_gettime(CLOCK_MONOTONIC,&start);
    // printf("%d \n",clock_gettime(CLOCK_REALTIME,&start) );

    nnz = 5*(ROW-2)*(ROW-2) + 2*(ROW-2)*2 + ROW*2;

    col_ind  = (int *) malloc(nnz * sizeof(int));
    row_ptr  = (int *) malloc((ROW*COL+1) * sizeof(int));
    r_ptr_b  = (int *) malloc((ROW*COL) * sizeof(int));
    r_ptr_e  = (int *) malloc((ROW*COL) * sizeof(int));

    nnzeros  = (double *) malloc(nnz * sizeof(double));
    tmp      = (double *) malloc(ROW*COL * sizeof(double));
    x        = (double *) malloc(ROW*COL * sizeof(double));
    b        = (double *) malloc(ROW*COL * sizeof(double));
    z        = (double *) malloc(ROW*COL * sizeof(double));
    r        = (double *) malloc(ROW*COL * sizeof(double));
    r_new    = (double *) malloc(ROW*COL * sizeof(double));

    printf("nnz : %d \n",nnz);

    make_Abx(nnzeros,col_ind,row_ptr,b,x,p,nnz,dx,dy);
    m = ROW*COL;
    k = ROW*COL;
    alp = 1.0;
    bet = 0.0;

    printf("m : %d, k : %d, alp : %f, bet : %f\n",m,k,alp,bet );

    mkl_dcsrmv("N",&m,&k,&alp,"G**C",nnzeros,
                col_ind,row_ptr,row_ptr+1,x,&bet,tmp);

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
       mkl_dcsrmv("N",&m,&k,&alp,"G**C",nnzeros,
                   col_ind,row_ptr,row_ptr+1,z,&bet,tmp);
       alpha = cblas_ddot(ROW*COL,r,1,r,1)/cblas_ddot(ROW*COL,z,1,tmp,1);

       for (i=0;i<ROW;i++){
           for (j=0;j<COL;j++){
               x[COL*i+j] = x[COL*i+j] + alpha * z[COL*i+j];
               r_new[COL*i+j] = r[COL*i+j] - alpha*tmp[COL*i+j];
           }
       }

       if (norm_L2(r_new) < tol ){
          *iter = it;
          //---------------------------------------
          //   Redistribute x vector to array
          //---------------------------------------
          for (i=0;i<ROW;i++){
            for (j=0;j<COL;j++){
              p[i][j] = x[COL*i+j];
            }
          }

          free(col_ind);
          free(row_ptr);
          free(nnzeros);
          free(tmp);
          free(x);
          free(b);
          free(z);
          free(r);
          free(r_new);

          end_t = clock();
          clock_gettime(CLOCK_MONOTONIC,&end);
          diff = 1e+9 * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
          printf("elapsed time = %llu nanoseconds\n", (long long unsigned int) diff);
          
          *tot_time = (double)(end_t - start_t)/(CLOCKS_PER_SEC);
          break;
       }

       beta = cblas_ddot(ROW*COL,r_new,1,r_new,1)/cblas_ddot(ROW*COL,r,1,r,1);

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
void make_Abx(double *nnzeros, int *col_ind,int *row_ptr,
              double *b, double *x, double **u,
              int nnz,double dx, double dy){
  int i,j,row,row_out,row_in,count=0;

  //----------------------------------------
  //         Make Matrix A using CSR
  //----------------------------------------
  for (row_out=0;row_out<ROW;row_out++){
    for (row_in=0;row_in<ROW;row_in++){

        row = row_out*ROW+row_in;
        row_ptr[row] = count  ;

        if(row_out == 0){
          nnzeros[count] = 1; col_ind[count++] = row;
        }
        else if (row_out == ROW-1 ){
          nnzeros[count] = 1; col_ind[count++] = row;
        }
        else{
          if (row_in == 0){
            nnzeros[count] = -1 ; col_ind[count++] = row;
            nnzeros[count] = 1  ; col_ind[count++] = row+1;
          }
          else if (row_in == ROW-1){
            nnzeros[count] = 1  ; col_ind[count++] = row-1;
            nnzeros[count] = -1 ; col_ind[count++] = row;
          }
          else{
            nnzeros[count] = 1;  col_ind[count++] = row-ROW;
            nnzeros[count] = 1;  col_ind[count++] = row-1;
            nnzeros[count] = -4; col_ind[count++] = row;
            nnzeros[count] = 1;  col_ind[count++] = row+1;
            nnzeros[count] = 1;  col_ind[count++] = row+ROW;
          }
        }

    }
  }
  row_ptr[ROW*COL] = count;

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
        sum = sum + a[i]*a[i];
    }
    return sqrt(sum);
}

void vmdot(double *nnzeros, int *col_ind,int *row_ptr,
           double *x,double *b)
{
    int i,j;
    double sum ;

    for (i=0;i<ROW*COL;i++){
            b[i] = 0;
    }

    for (i=0;i<ROW*COL;i++){
        sum = 0;
        for (j=row_ptr[i];j<row_ptr[i+1];j++){
          sum += nnzeros[j] * x[col_ind[j]];
        }
        b[i] = sum;
    }
}

double vvdot(double *a, double *b)
{
    int i;
    double c = 0;

    for (i=0;i<ROW*COL;i++){
        c = c + a[i]*b[i];
    }

    return c;
}
