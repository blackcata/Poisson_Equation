#include <stdio.h>
#include <math.h>

#include <string.h>
#include "def.h"

void initialization(double **p);
void write_u(char *dir_nm,char *file_nm, double **p,double dx, double dy);

//-----------------------------------
// Poisson Solvers - Jacobi, SOR, CG
//-----------------------------------
void Jacobi(double **p,double dx, double dy, double tol,
                       double *tot_time, int *iter,int BC);
void SOR(double **p,double dx, double dy, double tol, double omega,
                               double *tot_time,int *iter,int BC);
void Conjugate_Gradient(double **p,double dx, double dy, double tol,
                                   double *tot_time,int *iter,int BC);

//-----------------------------------
//        Mathematical tools
//-----------------------------------
double func(int i, int j, double dx, double dy);
void func_anal(double **p, int row_num, int col_num, double dx, double dy);
void error_rms(double **p, double **p_anal, double *err);
void poisson_solver(double **u, double **u_anal, double tol, double omega,
                    int BC, int method, char *dir_name){

  char *file_name ;

  int iter = 0;
  double Lx = 1.0, Ly = 1.0;
  double dx, dy, err = 0, tot_time = 0;

  dx = Lx/(ROW-1);
  dy = Ly/(COL-1);

  //-----------------------------
  //      Analytic Solutions
  //-----------------------------
  file_name = "Analytic_solution.plt";
  func_anal(u_anal,ROW,COL,dx,dy);
  write_u(dir_name,file_name,u_anal,dx,dy);

  switch (method) {
    case 1 :
      //-----------------------------
      //        Jacobi Method
      //-----------------------------
      initialization(u);
      Jacobi(u,dx,dy,tol,&tot_time,&iter,BC);
      error_rms(u,u_anal,&err);
      printf("Jacobi Method - Error : %f, Iteration : %d, Time : %f s \n",err,iter,tot_time);

      file_name = "Jacobi_result.plt";
      write_u(dir_name,file_name,u,dx,dy);
      break;

    case 2 :
       //-----------------------------
       //         SOR Method
       //-----------------------------
       initialization(u);
       SOR(u,dx,dy,tol,omega,&tot_time,&iter,BC);
       error_rms(u,u_anal,&err);
       printf("SOR Method - Error : %f, Iteration : %d, Time : %f s \n",err,iter,tot_time);

       file_name = "SOR_result.plt";
       write_u(dir_name,file_name,u,dx,dy);
      break;

    case 3 :
       //-----------------------------
       //  Conjugate Gradient Method
       //-----------------------------
       initialization(u);
       Conjugate_Gradient(u,dx,dy,tol,&tot_time,&iter,BC);
       error_rms(u,u_anal,&err);
       printf("CG method - Error : %f, Iteration : %d, Time : %f s \n",err,iter,tot_time);

       file_name = "CG_result.plt";
       write_u(dir_name,file_name,u,dx,dy);
      break;
  }

}

double func(int i,int j,double dx,double dy)
{
    return sin(pi*i*dx)*cos(pi*j*dy);
}

void initialization(double **p)
{
    int i,j;
    for (i=0;i<ROW;i++){
        for (j=0;j<COL;j++){
            p[i][j] = 0; }}

}

void error_rms(double **p, double **p_anal, double *err)
{
  int i,j;
  for (i=0;i<ROW;i++){
    for (j=0;j<COL;j++){
      *err = *err + pow(p[i][j] -p_anal[i][j],2);
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

void write_u(char *dir_nm,char *file_nm, double **p,double dx, double dy)
{
    FILE* stream;
    int i,j;
    char file_path[50];
    sprintf(file_path,"%s%s",dir_nm,file_nm);

    stream=fopen(file_path,"w");
    fprintf(stream,"ZONE I=%d J=%d \n",ROW,COL);
    for (i=0;i<ROW;i++){
        for(j=0;j<COL;j++){
            fprintf(stream,"%f %f %f \n",i*dx,j*dy,p[i][j]); }}
    fclose(stream);
}
