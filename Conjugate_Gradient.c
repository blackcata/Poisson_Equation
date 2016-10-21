#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
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
    double alpha,beta ;

    int *col_ind, *row_ptr;
    double *nnzeros,*tmp,*x,*b,*z,*r,*r_new;

    time_t start_t =0, end_t =0;

    start_t = clock();
    nnz = 5*(ROW-2)*(ROW-2) + 2*(ROW-2)*2 + ROW*2;

    col_ind  = (int *) malloc(nnz * sizeof(int));
    row_ptr  = (int *) malloc((ROW*COL+1) * sizeof(int));

    nnzeros  = (double *) malloc(nnz * sizeof(double));
    tmp      = (double *) malloc(ROW*COL * sizeof(double));
    x        = (double *) malloc(ROW*COL * sizeof(double));
    b        = (double *) malloc(ROW*COL * sizeof(double));
    z        = (double *) malloc(ROW*COL * sizeof(double));
    r        = (double *) malloc(ROW*COL * sizeof(double));
    r_new    = (double *) malloc(ROW*COL * sizeof(double));
    printf("nnz : %d \n",nnz);

    make_Abx(nnzeros,col_ind,row_ptr,b,x,p,nnz,dx,dy);
    vmdot(nnzeros,col_ind,row_ptr,x,tmp);

   #pragma omp parallel for shared(r,z,b,tmp) private(i,j)
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
       vmdot(nnzeros,col_ind,row_ptr,z,tmp);
       alpha = vvdot(r,r)/vvdot(z,tmp);

       #pragma omp parallel for shared(x,z,r_new,r,tmp) private(i,j)
       for (i=0;i<ROW;i++){
           for (j=0;j<COL;j++){
               x[COL*i+j] = x[COL*i+j] + alpha * z[COL*i+j];
               r_new[COL*i+j] = r[COL*i+j] - alpha*tmp[COL*i+j];
           }
       }

       if (norm_L2(r_new) < tol ){
          // printf("iteration : %d, tol : %f, value : %f\n",it,tol,norm_L2(r_new) );
          *iter = it;
          //---------------------------------------
          //   Redistribute x vector to array
          //---------------------------------------
          #pragma omp parallel for shared(p,x) private(i,j)
          for (i=0;i<ROW;i++)
          {
            for (j=0;j<COL;j++)
            {
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
          *tot_time = (double)(end_t - start_t)/(CLOCKS_PER_SEC);
          break;
       }

       beta = vvdot(r_new,r_new)/vvdot(r,r);
       #pragma omp parallel for shared(z,r_new,r) private(i,j)
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

      // printf("row : %d, count : %d, nnzeros : %f, row in : %d, row_out : %d \n",row,count,nnzeros[count-1],row_in,row_out );
    }
  }
  row_ptr[ROW*COL] = count;

  //--------------------------------
  //         Make Vector x
  //--------------------------------
  #pragma omp parallel for shared(x,u) private(i,j)
  for (i=0;i<ROW;i++){
      for (j=0;j<COL;j++){
          x[ROW*i+j] = u[i][j];
      }
  }
  //--------------------------------
  //        Make Vector b
  //--------------------------------
  #pragma omp parallel for shared(b) private(i,j)
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
