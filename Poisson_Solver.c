#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "def.h"

void initialization(double **p);
void write_u(char *dir_nm,char *file_nm,int write_type,
             double *p,double dx,double dy);

//----------------------------------------------------------------------------//
//                    Poisson Solvers - Jacobi, SOR, CG                       //
//----------------------------------------------------------------------------//
void Jacobi(double **p,double dx, double dy, double tol,
            double *tot_time, int *iter,int BC,
            char* file_name, char* dir_name,int write_type);
void SOR(double **p,double dx, double dy, double tol, double omega,
         double *tot_time, int *iter,int BC,
         char* file_name, char* dir_name,int write_type);
void Conjugate_Gradient(double **p,double dx, double dy, double tol,
                        double *tot_time,int *iter,int BC,
                        char *file_name,char *dir_name,int write_type);

//----------------------------------------------------------------------------//
//                              Mathematical tools                            //
//----------------------------------------------------------------------------//
double func(int i, int j, double dx, double dy);
void func_anal(double **p, int row_num, int col_num, double dx, double dy);
void error_rms(double **p, double **p_anal, double *err);

//----------------------------------------------------------------------------//
//                                                                            //
//                      Main Poisson_Solver Subroutine                        //
//                                                                            //
//----------------------------------------------------------------------------//
void poisson_solver(double **u, double **u_anal, double tol, double omega,
                    int BC, int method, int write_type,char *dir_name){

  char *file_name;

  int iter = 0, myrank;
  double Lx = 1.0, Ly = 1.0;
  double dx, dy, err = 0, tot_time = 0;

  dx = Lx/(ROW-1);
  dy = Ly/(COL-1);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  //--------------------------------------------------------------------------//
  //                            Analytic Solutions                            //
  //--------------------------------------------------------------------------//
  file_name = "Analytic_solution.plt";
  func_anal(u_anal,ROW,COL,dx,dy);
  // write_u(dir_name,file_name,u_anal,dx,dy);

  switch (method) {
    case 1 :
      //----------------------------------------------------------------------//
      //                              Jacobi Method                           //
      //----------------------------------------------------------------------//
      file_name = "Jacobi_result.plt";

      initialization(u);
      Jacobi(u,dx,dy,tol,&tot_time,&iter,BC,file_name,dir_name,write_type);
      error_rms(u,u_anal,&err);

      if(myrank==0) {
        printf("Jacobi Method - Error : %e, Iteration : %d, Time : %f s \n",
                err,iter,tot_time);
      }
      break;

    case 2 :
       //---------------------------------------------------------------------//
       //                             SOR Method                              //
       //---------------------------------------------------------------------//
       file_name = "SOR_result.plt";

       initialization(u);
       SOR(u,dx,dy,tol,omega,&tot_time,&iter,BC,file_name,dir_name,write_type);
       error_rms(u,u_anal,&err);

       if(myrank==0) {
         printf("SOR Method - Error : %e, Iteration : %d, Time : %f s \n",
                err,iter,tot_time);
       }
      break;

    case 3 :
       //---------------------------------------------------------------------//
       //                       Conjugate Gradient Method                     //
       //---------------------------------------------------------------------//
       file_name = "CG_result.plt";

       initialization(u);
       Conjugate_Gradient(u,dx,dy,tol,&tot_time,&iter,BC,
                          file_name,dir_name,write_type);
       error_rms(u,u_anal,&err);

       MPI_Barrier(MPI_COMM_WORLD);
       if(myrank==0){
         printf("CG method - Error : %e, Iteration : %d, Time : %f s \n",
                  err,iter,tot_time);
       }
       break;
  }

}


void initialization(double **p)
{
    int i,j;
    for (i=0;i<ROW;i++){
        for (j=0;j<COL;j++){
            p[i][j] = 0; }}

}

//----------------------------------------------------------------------------//
//                           Mathematical functions                           //
//----------------------------------------------------------------------------//
double func(int i,int j,double dx,double dy)
{
    return sin(pi*i*dx)*cos(pi*j*dy);
}

void error_rms(double **p, double **p_anal, double *err)
{
  int i,j;
  for (i=0;i<ROW;i++){
    for (j=0;j<COL;j++){
      *err = *err + pow(p[i][j] - p_anal[i][j],2);
    }
  }

  *err = sqrt(*err)/(ROW*COL);
}


void func_anal(double **p, int row_num, int col_num, double dx, double dy)
{
    int i,j;
    for (i=0;i<row_num;i++){
        for (j=0;j<col_num;j++){
            p[i][j] = -1/(2*pow(pi,2))*sin(pi*i*dx)*cos(pi*j*dy); }}
}

//----------------------------------------------------------------------------//
//                             Writing functions                              //
//----------------------------------------------------------------------------//
void write_u(char *dir_nm,char *file_nm,int write_type,
             double *p,double dx,double dy)
{
    FILE* stream;
    MPI_File thefile;
    MPI_Offset disp;

    int i,j,myrank,nproc;
    double tmp[3];
    char file_path[50], header[50];

    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    sprintf(file_path,"%s%s",dir_nm,file_nm);

    switch (write_type) {
      case 1:
        //-------------------------------------------------------------------//
        //                           Post Reassembly                         //
        //-------------------------------------------------------------------//
        stream=fopen(file_path,"w");
        fprintf(stream,"ZONE I=%d J=%d \n",ROW,COL);
        for (i=0;i<ROW;i++){
            for(j=0;j<COL;j++){
                fprintf(stream,"%f %f %f \n",i*dx,j*dy,p[i*ROW+j]);
            }
        }
        fclose(stream);
        break;

      case 2:
        //-------------------------------------------------------------------//
        //                           Single Task I/O                         //
        //-------------------------------------------------------------------//
        stream=fopen(file_path,"w");
        fprintf(stream,"ZONE I=%d J=%d \n",COL,ROW/nproc);
        for (i=myrank*(ROW/nproc);i<(myrank+1)*(ROW/nproc);i++){
            for(j=0;j<COL;j++){
                fprintf(stream,"%f %f %f \n",i*dx,j*dy,
                                p[(i-myrank*(ROW/nproc))*ROW+j]);
            }
        }
        fclose(stream);
        break;

      case 3:
        //-------------------------------------------------------------------//
        //                              MPI I/O                              //
        //-------------------------------------------------------------------//
        sprintf(header,"ZONE I=%d J=%d \n",COL,ROW);

        MPI_File_open(MPI_COMM_WORLD,file_path,
                    MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,&thefile);

        MPI_File_write(thefile,header,15,MPI_DOUBLE,MPI_STATUS_IGNORE);
        disp = myrank*(ROW*COL/nproc)*sizeof(double)+15;
        MPI_File_set_view(thefile,disp,MPI_INT,MPI_INT,"native",MPI_INFO_NULL);
        for (i=myrank*(ROW/nproc);i<(myrank+1)*(ROW/nproc);i++){
            for(j=0;j<COL;j++){
            tmp[0] = i*dx;
            tmp[1] = j*dy;
            tmp[2] = p[(i-myrank*(ROW/nproc))*ROW+j];
            MPI_File_write(thefile,tmp,3,MPI_DOUBLE,MPI_STATUS_IGNORE);
         // MPI_File_write(thefile,p,ROW*COL/nproc,MPI_DOUBLE,MPI_STATUS_IGNORE);
          }
        }
        MPI_File_close(&thefile);
        break;

    }

}
