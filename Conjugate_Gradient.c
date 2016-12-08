#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#include "def.h"

void write_u(char *dir_nm,char *file_nm,int write_type,
             double *p,double dx,double dy);

//----------------------------------------------------------------------------//
//                            Matrix Calculations                             //
//----------------------------------------------------------------------------//
double norm_L2(int num, double *a);
double vvdot(int num, double *a, double *b);
void vmdot(int myrank, int nproc,int row,int col,
           double **L, double **A, double **R,double *x,double *b);

void make_Abx(int ista,int iend,double **A,double **L,double **R,
              double *b,double *x,double**u,double dx,double dy);

//----------------------------------------------------------------------------//
//                          Mathematical functions                            //
//----------------------------------------------------------------------------//
double func(int i, int j, double dx, double dy);

void Conjugate_Gradient(double **p,double dx, double dy, double tol,
                        double *tot_time,int *iter,int BC,
                        char *file_name,char *dir_name,int write_type)
{
    int i,j,k,it;
    int nproc,myrank,ista,iend;
    double alpha,beta,ts,te,time_ws,time_we;
    double rnew_sum,rnew_sum_loc,rr_sum,rr_sum_loc,rn_sum,rn_sum_loc,
           zAz_sum,zAz_sum_loc ;

    double **A, **L, **R;
    double *tmp, *x, *z, *r, *r_new;
    double *tmp_loc, *b_loc, *r_loc, *r_new_loc, *x_loc, *z_loc;

    time_t start_t = 0, end_t = 0;
    char loc_name[50];

    start_t = clock();
    ts = MPI_Wtime();

    //-----------------------------------------------------------------------//
    //                              MPI Setting                              //
    //-----------------------------------------------------------------------//
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    printf("ROW X COL : %d, nproc : %d \n",ROW*COL,ROW*COL/nproc);

    ista = myrank*(ROW/nproc);
    iend = (myrank+1)*(ROW/nproc)-1;
    printf("[ista,iend] : [%d,%d]\n \n",ista,iend );

    //------------------------------------------------------------------------//
    //                           Memory allocation                            //
    //------------------------------------------------------------------------//
    A = (double **) malloc(ROW*COL/nproc * sizeof(double));
    L = (double **) malloc(ROW * sizeof(double));
    R = (double **) malloc(ROW * sizeof(double));
    for (i=0;i<ROW*COL/nproc;i++) {
      A[i] = (double *) malloc(ROW*COL/nproc * sizeof(double));
    }
    for (i=0;i<ROW;i++){
      L[i] = (double *) malloc(COL * sizeof(double));
      R[i] = (double *) malloc(COL * sizeof(double));
    }

    tmp    = (double *) malloc(ROW*COL * sizeof(double));
    x      = (double *) malloc(ROW*COL * sizeof(double));

    z      = (double *) malloc(ROW*COL * sizeof(double));
    r      = (double *) malloc(ROW*COL * sizeof(double));
    r_new  = (double *) malloc(ROW*COL * sizeof(double));

    b_loc       = (double *) malloc(ROW*COL * sizeof(double));
    x_loc       = (double *) malloc(ROW*COL/nproc * sizeof(double));
    z_loc       = (double *) malloc(ROW*COL/nproc * sizeof(double));
    tmp_loc     = (double *) malloc(ROW*COL/nproc * sizeof(double));
    r_loc       = (double *) malloc(ROW*COL/nproc * sizeof(double));
    r_new_loc   = (double *) malloc(ROW*COL/nproc * sizeof(double));

    MPI_Barrier(MPI_COMM_WORLD);
    make_Abx(ista,iend,A,L,R,b_loc,x_loc,p,dx,dy);

    vmdot(myrank,nproc,ROW*COL/nproc,ROW*COL/nproc,L,A,R,x_loc,tmp_loc);

    for (i=0;i<ROW*COL/nproc;i++){
           r_loc[i] = b_loc[i] - tmp_loc[i];
           z_loc[i] = r_loc[i];
    }

   //-------------------------------------------------------------------------//
   //                Main Loop of Conjugate Gradient Method                   //
   //-------------------------------------------------------------------------//
   for (it=0;it<itmax;it++)
   {

       vmdot(myrank,nproc,ROW*COL/nproc,ROW*COL/nproc,L,A,R,z_loc,tmp_loc);
       rr_sum_loc = vvdot(ROW*COL/nproc,r_loc,r_loc);
       zAz_sum_loc = vvdot(ROW*COL/nproc,z_loc,tmp_loc);

       MPI_Allreduce(&rr_sum_loc,&rr_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
       MPI_Allreduce(&zAz_sum_loc,&zAz_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
       alpha = rr_sum/zAz_sum;

       for (i=0;i<ROW*COL/nproc;i++){
               x_loc[i] = x_loc[i] + alpha * z_loc[i];
               r_new_loc[i] = r_loc[i] - alpha*tmp_loc[i];
       }

       //---------------------------------------------------------------------//
       //                        Convergence Criteria                         //
       //---------------------------------------------------------------------//
       rnew_sum_loc = pow(norm_L2(ROW*COL/nproc,r_new_loc),2);
       MPI_Allreduce(&rnew_sum_loc,&rnew_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
       if (myrank==0) printf("cri : %f, tol : %f, \n",sqrt(rnew_sum),tol);

       if ( sqrt(rnew_sum) < tol ){
          //------------------------------------------------------------------//
          //                  Redistribute x vector to array                  //
          //------------------------------------------------------------------//
          *iter = it;

          end_t = clock();
          te = MPI_Wtime();
          *tot_time = (double)(end_t - start_t)/(CLOCKS_PER_SEC);
          if (myrank==0) printf("Total time is : %f s \n",te-ts );

          time_ws = MPI_Wtime();
          switch (write_type) {
            case 1 :
            //-----------------------------------------------------------------//
            //                         Post Reassembly                         //
            //-----------------------------------------------------------------//
              MPI_Allgather(&x_loc[0],ROW*COL/nproc,MPI_DOUBLE,
                          x,ROW*COL/nproc,MPI_DOUBLE,MPI_COMM_WORLD);
              if(myrank==0) write_u(dir_name,file_name,write_type,x,dx,dy);
              break;

            case 2 :
            //-----------------------------------------------------------------//
            //                          Single Task I/O                        //
            //-----------------------------------------------------------------//
              sprintf(loc_name,"%d.%s",myrank,file_name);
              write_u(dir_name,loc_name,write_type,x_loc,dx,dy);
              break;

            case 3:
            //-----------------------------------------------------------------//
            //                            MPI I/O                              //
            //-----------------------------------------------------------------//
              write_u(dir_name,file_name,write_type,x_loc,dx,dy);
              break;
          }
          time_we = MPI_Wtime();
          if (myrank==0) printf("Writing time is : %f s \n",time_we-time_ws);
          free(b_loc);
          free(x_loc);
          free(z_loc);
          free(tmp_loc);
          free(r_loc);
          free(r_new_loc);

          break;
       }

       rn_sum_loc = vvdot(ROW*COL/nproc,r_new_loc,r_new_loc);
       MPI_Allreduce(&rn_sum_loc,&rn_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
       beta = rn_sum/rr_sum;

       //---------------------------------------------------------------------//
       //                               Update                                //
       //---------------------------------------------------------------------//
       for (i=0;i<ROW*COL/nproc;i++){
           z_loc[i] = r_new_loc[i] + beta*z_loc[i];
           r_loc[i] = r_new_loc[i];
       }
   }

}

//----------------------------------------------------------------------------//
//                     Make Stiffness matrix of CG method                     //
//----------------------------------------------------------------------------//
void make_Abx(int ista,int iend,double **A,double **L,double **R,
              double *b,double *x,double**u,double dx,double dy)
{
    int i,j,k,l,tmp=0,tmp2=0;;

    //------------------------------------------------------------------------//
    //                              Make Matrix A                             //
    //------------------------------------------------------------------------//
    for (k=ista;k<iend+1;k++){
        tmp2 = 0;
        for (l=ista;l<iend+1;l++){
            if (k==l){
                if (k==0 || k==ROW-1){
                  for (i=0;i<ROW;i++){
                      A[COL*tmp+i][ROW*tmp2+i]   = 1;
                  }
                }
                else{
                  for (i=0;i<ROW;i++){
                      if (i == 0){
                          A[COL*tmp+i][ROW*tmp2+i]   = -1;
                          A[COL*tmp+i+1][ROW*tmp2+i] = 1;
                      }
                      else if (i == ROW-1){
                          A[COL*tmp+i][ROW*tmp2+i]   = -1;
                          A[COL*tmp+i-1][ROW*tmp2+i] = 1;
                      }
                      else {
                          A[COL*tmp+i][ROW*tmp2+i]   = -4;
                          A[COL*tmp+i-1][ROW*tmp2+i] = 1;
                          A[COL*tmp+i+1][ROW*tmp2+i] = 1;
                      }
                  }
                }
            }

            else if ( abs(k-l) == 1 && k!=0 && k!=ROW-1){
                for (i=0;i<ROW;i++){
                  if (i==0 || i==ROW-1){
                    A[COL*tmp+i][ROW*tmp2+i] = 0;
                  }
                  else{
                    A[COL*tmp+i][ROW*tmp2+i] = 1;
                  }
                }
            }
            tmp2 += 1;
        }
        tmp += 1;
    }
    //------------------------------------------------------------------------//
    //                            Make matrix L,R                             //
    //------------------------------------------------------------------------//
    for (i=0;i<ROW;i++){
      if (i==0 || i==ROW-1){
        L[i][i] = 0;
        R[i][i] = 0;
      }
      else{
        L[i][i] = 1;
        R[i][i] = 1;
      }
    }

    //------------------------------------------------------------------------//
    //                              Make Vector x                             //
    //------------------------------------------------------------------------//
    tmp = 0;
    for (i=ista;i<iend+1;i++){
        for (j=0;j<COL;j++){
            x[tmp] = u[i][j];
            tmp += 1;
        }
    }
    //------------------------------------------------------------------------//
    //                              Make Vector b                             //
    //------------------------------------------------------------------------//
    tmp =0;
    for (i=ista;i<iend+1;i++){
        for (j=0;j<COL;j++){
          if (i==0 || i==ROW-1 || j==0 || j==COL-1)
              b[tmp] = 0;//1/(2*pow(pi,2))*func(i,j,dx,dy);
          else
              b[tmp] = dx*dx*func(i,j,dx,dy);
          tmp +=1;
        }
    }
}

//----------------------------------------------------------------------------//
//                        Matrix Calcuation Functions                         //
//    -L2 norm, vector-matrix production, vector-vector production functions  //
//----------------------------------------------------------------------------//
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
    int i,j,tmp;
    double *send_l, *recv_l,*send_r, *recv_r;

    MPI_Request req1, req2, req3, req4;
    MPI_Status stat1, stat2, stat3, stat4;

    send_l = (double *) malloc(COL * sizeof(double));
    recv_l = (double *) malloc(COL * sizeof(double));
    send_r = (double *) malloc(COL * sizeof(double));
    recv_r = (double *) malloc(COL * sizeof(double));

    for (i=0;i<COL;i++){
        send_l[i] = x[i];
        send_r[i] = x[row-COL+i];
    }

    for (i=0;i<ROW*COL/nproc;i++){
        b[i] = 0;
    }
    //------------------------------------------------------------------------//
    //                            Left exchange                               //
    //------------------------------------------------------------------------//
    if (myrank!=0){
      MPI_Irecv(recv_l,COL,MPI_DOUBLE,myrank-1,101,MPI_COMM_WORLD,&req1);
      MPI_Isend(send_l,COL,MPI_DOUBLE,myrank-1,102,MPI_COMM_WORLD,&req2);
    }
    //------------------------------------------------------------------------//
    //                            Right exchange                              //
    //------------------------------------------------------------------------//
    if (myrank!=nproc-1){
      MPI_Irecv(recv_r,COL,MPI_DOUBLE,myrank+1,102,MPI_COMM_WORLD,&req3);
      MPI_Isend(send_r,COL,MPI_DOUBLE,myrank+1,101,MPI_COMM_WORLD,&req4);
    }

    //------------------------------------------------------------------------//
    //                          Center calculation                            //
    //------------------------------------------------------------------------//
    for (i=0;i<row;i++){
        for (j=0;j<col;j++){
            b[i] = b[i] + A[i][j]*x[j];
        }
    }
    //------------------------------------------------------------------------//
    //                           Left calculation                             //
    //------------------------------------------------------------------------//
    if (myrank!=0) MPI_Wait(&req1,&stat1);
    if (myrank!=nproc-1) MPI_Wait(&req4,&stat4);

    if (myrank!=0){
      tmp = 0;
      for (i=0;i<COL;i++){
        for (j=0;j<COL;j++){
          b[i] = b[i] + L[tmp][j]*recv_l[j];
        }
        tmp += 1;
      }
    }
    //------------------------------------------------------------------------//
    //                          Right calculation                             //
    //------------------------------------------------------------------------//
    if (myrank!=0) MPI_Wait(&req2,&stat2);
    if (myrank!=nproc-1) MPI_Wait(&req3,&stat3);

    if (myrank!=nproc-1){
      tmp = 0;
      for (i=row-COL;i<row;i++){
        for (j=0;j<COL;j++){
          b[i] = b[i] + R[tmp][j]*recv_r[j];
        }
        tmp += 1;
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
