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
void vmdot(int myrank, int nproc,int row,int col,
           double **L, double **A, double **R,double *x,double *b);

void make_Abx(int ista,int iend,double **A,double **L,double **R,
              double *b,double *x,double**u,double dx,double dy);

//-----------------------------------
//      Mathematical functions
//-----------------------------------
double func(int i, int j, double dx, double dy);

void Conjugate_Gradient(double **p,double dx, double dy, double tol,
                                   double *tot_time,int *iter, int BC)
{
    int i,j,k,it,tt;
    int nproc,myrank,ista,iend;
    double alpha,beta,ts,te;
    double rnew_sum,rnew_sum_loc,rr_sum,rr_sum_loc,rn_sum,rn_sum_loc,zAz_sum,zAz_sum_loc ;

    double **A, **L, **R;
    double *tmp, *x, *b, *z, *r, *r_new;
    double *tmp_loc, *r_loc, *r_new_loc, *x_loc, *z_loc;

    time_t start_t = 0, end_t = 0;

    start_t = clock();
    ts = MPI_Wtime();

    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    printf("ROW X COL : %d, nproc : %d \n",ROW*COL,ROW*COL/nproc);

    ista = myrank*(ROW/nproc);
    iend = (myrank+1)*(ROW/nproc)-1;
    printf("[ista,iend] : [%d,%d]\n \n",ista,iend );

    A = (double **) malloc(ROW*COL/nproc * sizeof(double));
    L = (double **) malloc(ROW*COL/nproc * sizeof(double));
    R = (double **) malloc(ROW*COL/nproc * sizeof(double));
    for (i=0;i<ROW*COL/nproc;i++) {
      A[i] = (double *) malloc(ROW*COL/nproc * sizeof(double));
      L[i] = (double *) malloc(ROW*COL/nproc * sizeof(double));
      R[i] = (double *) malloc(ROW*COL/nproc * sizeof(double));
    }

    tmp    = (double *) malloc(ROW*COL * sizeof(double));
    x      = (double *) malloc(ROW*COL * sizeof(double));
    b      = (double *) malloc(ROW*COL * sizeof(double));
    z      = (double *) malloc(ROW*COL * sizeof(double));
    r      = (double *) malloc(ROW*COL * sizeof(double));
    r_new  = (double *) malloc(ROW*COL * sizeof(double));

    x_loc       = (double *) malloc(ROW*COL/nproc * sizeof(double));
    z_loc       = (double *) malloc(ROW*COL/nproc * sizeof(double));
    tmp_loc     = (double *) malloc(ROW*COL/nproc * sizeof(double));
    r_loc       = (double *) malloc(ROW*COL/nproc * sizeof(double));
    r_new_loc   = (double *) malloc(ROW*COL/nproc * sizeof(double));

    MPI_Barrier(MPI_COMM_WORLD);
    make_Abx(ista,iend,A,L,R,b,x,p,dx,dy);

    vmdot(myrank,nproc,ROW*COL/nproc,ROW*COL/nproc,L,A,R,x,tmp_loc);
  //   MPI_Allgather(&tmp_loc[0],ROW*COL/nproc,MPI_DOUBLE,
  //                 tmp,ROW*COL/nproc,MPI_DOUBLE,MPI_COMM_WORLD);
   //
  //  tt = 0;
  //  for (i=ista;i<iend+1;i++){
  //      for (j=0;j<COL;j++){
  //          r_loc[tt] = b[COL*i+j] - tmp[COL*i+j];
  //          z_loc[tt] = r_loc[tt];
  //          tt+=1;
  //      }
  //  }
   //
  //  MPI_Allgather(&r_loc[0],ROW*COL/nproc,MPI_DOUBLE,
  //                r,ROW*COL/nproc,MPI_DOUBLE,MPI_COMM_WORLD);
  //  MPI_Allgather(&z_loc[0],ROW*COL/nproc,MPI_DOUBLE,
  //                z,ROW*COL/nproc,MPI_DOUBLE,MPI_COMM_WORLD);

  //  //---------------------------------------
  //  //   Main Loop of Conjugate_Gradient
  //  //---------------------------------------
  //  for (it=0;it<itmax;it++)
  //  {
  //      vmdot(ROW*COL/nproc,ROW*COL,A,z,tmp_loc);
  //      rr_sum_loc = vvdot(ROW*COL/nproc,r_loc,r_loc);
  //      zAz_sum_loc = vvdot(ROW*COL/nproc,z_loc,tmp_loc);
   //
  //      MPI_Allreduce(&rr_sum_loc,&rr_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  //      MPI_Allreduce(&zAz_sum_loc,&zAz_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  //      alpha = rr_sum/zAz_sum;
   //
  //      tt = 0;
  //      for (i=ista;i<iend+1;i++){
  //          for (j=0;j<COL;j++){
  //              x_loc[tt] = x_loc[tt] + alpha * z_loc[tt];
  //              r_new_loc[tt] = r_loc[tt] - alpha*tmp_loc[tt];
  //              tt+=1;
  //          }
  //      }
   //
  //      rnew_sum_loc = pow(norm_L2(ROW*COL/nproc,r_new_loc),2);
  //      MPI_Allreduce(&rnew_sum_loc,&rnew_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   //
  //      if ( sqrt(rnew_sum) < tol ){
  //         //---------------------------------------
  //         //   Redistribute x vector to array
  //         //---------------------------------------
  //         MPI_Allgather(&x_loc[0],ROW*COL/nproc,MPI_DOUBLE,
  //                     x,ROW*COL/nproc,MPI_DOUBLE,MPI_COMM_WORLD);
  //         if (myrank == 0){
  //           tt = 0;
  //           for (i=0;i<ROW;i++){
  //             for (j=0;j<COL;j++){
  //               p[i][j] = x[tt];
  //               tt += 1;
  //             }
  //           }
  //         }
   //
  //         *iter = it;
  //         free(A);
  //         free(tmp);
  //         free(x);
  //         free(b);
  //         free(z);
  //         free(r);
  //         free(r_new);
   //
  //         free(x_loc);
  //         free(z_loc);
  //         free(tmp_loc);
  //         free(r_loc);
  //         free(r_new_loc);
   //
  //         end_t = clock();
  //         te = MPI_Wtime();
  //         *tot_time = (double)(end_t - start_t)/(CLOCKS_PER_SEC);
  //         if(myrank==0) printf("Total time is : %f s \n",te-ts );
  //         break;
  //      }
   //
  //      rn_sum_loc = vvdot(ROW*COL/nproc,r_new_loc,r_new_loc);
  //      MPI_Allreduce(&rn_sum_loc,&rn_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  //      beta = rn_sum/rr_sum;
   //
  //      tt = 0;
  //      for (i=ista;i<iend+1;i++){
  //          for (j=0;j<COL;j++){
  //              z_loc[tt] = r_new_loc[tt] + beta*z_loc[tt];
  //              r_loc[tt] = r_new_loc[tt];
  //              tt += 1;
  //          }
  //      }
   //
  //      MPI_Allgather(&z_loc[0],ROW*COL/nproc,MPI_DOUBLE,
  //                    z,ROW*COL/nproc,MPI_DOUBLE,MPI_COMM_WORLD);
  //  }

}

//------------------------------------------------------------
//             Make Stiffness matrix of CG method
//------------------------------------------------------------
void make_Abx(int ista,int iend,double **A,double **L,double **R,
              double *b,double *x,double**u,double dx,double dy)
{
    int i,j,k,l,tmp=0;

    //--------------------------------
    //         Make Matrix A
    //--------------------------------
    for (k=ista;k<iend+1;k++){
        for (l=0;l<COL;l++){
            if (k==l){
                if (k==0 || k==ROW-1){
                  for (i=0;i<ROW;i++){
                      A[COL*tmp+i][COL*tmp+i]   = 1;
                  }
                }
                else{
                  for (i=0;i<ROW;i++){
                      if (i == 0){
                          A[COL*tmp+i][COL*tmp+i]   = -1;
                          A[COL*tmp+i+1][COL*tmp+i] = 1;
                      }
                      else if (i == ROW-1){
                          A[COL*tmp+i][COL*tmp+i]   = -1;
                          A[COL*tmp+i-1][COL*tmp+i] = 1;
                      }
                      else {
                          A[COL*tmp+i][COL*tmp+i]   = -4;
                          A[COL*tmp+i-1][COL*tmp+i] = 1;
                          A[COL*tmp+i+1][COL*tmp+i] = 1;
                      }
                  }
                }
            }

            else if ( abs(k-l) == 1 && k!=0 && k!=ROW-1){
                for (i=0;i<ROW;i++){
                  if (i==0 || i==ROW-1){
                    L[COL*tmp+i][COL*tmp+i] = 0;
                    R[COL*tmp+i][COL*tmp+i] = 0;
                  }
                  else{
                    L[COL*tmp+i][COL*tmp+i] = 1;
                    R[COL*tmp+i][COL*tmp+i] = 1;
                  }
                }
            }
        }
        tmp += 1;
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

void vmdot(int myrank, int nproc,int row,int col,
           double **L, double **A, double **R,double *x,double *b)
{
    int i,j;
    double *send_l, *recv_l,*send_r, *recv_r;
    double *b_l, *b_c, *b_r;
    MPI_Request req1, req2, req3, req4;
    MPI_Status stat1, stat2, stat3, stat4;

    send_l = (double *) malloc(row * sizeof(double));
    recv_l = (double *) malloc(row * sizeof(double));
    send_r = (double *) malloc(row * sizeof(double));
    recv_r = (double *) malloc(row * sizeof(double));

    b_l = (double *) malloc(row * sizeof(double));
    b_c = (double *) malloc(row * sizeof(double));
    b_r = (double *) malloc(row * sizeof(double));

    for (i=0;i<row;i++){
        send_r[i] = x[i];
        send_l[i] = x[i];
    }

    // Left exchange
    if (myrank!=0){
      MPI_Irecv(recv_l,row,MPI_DOUBLE,myrank-1,101,MPI_COMM_WORLD,&req1);
    }
    if (myrank!=nproc-1){
      MPI_Isend(send_l,row,MPI_DOUBLE,myrank+1,101,MPI_COMM_WORLD,&req2);
    }

    // Right exchange
    if (myrank!=nproc-1){
      MPI_Irecv(recv_r,row,MPI_DOUBLE,myrank+1,102,MPI_COMM_WORLD,&req3);
    }
    if (myrank!=0){
      MPI_Isend(send_r,row,MPI_DOUBLE,myrank-1,102,MPI_COMM_WORLD,&req4);
    }

    for (i=0;i<row;i++){
        b_l[i] = 0;
        b_c[i] = 0;
        b_r[i] = 0;
    }

    // Center calculation
    for (i=0;i<row;i++){
        for (j=0;j<col;j++){
            b_c[i] = b_c[i] + A[i][j]*x[j];
        }
    }

    //Left calculation
    if (myrank!=0) MPI_Wait(&req1,&stat1);
    if (myrank!=nproc-1) MPI_Wait(&req2,&stat2);

    if (myrank!=nproc-1){
      for (i=0;i<row;i++){
        for (j=0;j<row;j++){
          b_l[i] = b_l[i] + L[i][j]*recv_l[i];
        }
      }
    }

    // Right calculation
    if (myrank!=nproc-1) MPI_Wait(&req3,&stat3);
    if (myrank!=0) MPI_Wait(&req4,&stat4);

    if (myrank!=0){
      for (i=0;i<row;i++){
        for (j=0;j<row;j++){
          b_r[i] = b_r[i] + R[i][j]*recv_r[i];
        }
      }
    }

    for (i=0;i<row;i++){
      b[i] = b_l[i] + b_c[i] + b_r[i];
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
