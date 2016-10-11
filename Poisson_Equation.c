//
//
//  Programming for 2D Poisson Equation
//                    By. Kyung Min Noh 2016
//
//  Solving 2D Poisson Equation
//
//  Laplce(u(x,y)) = f(x,y) for (x,y) in domain
//  on the boundary boundary_D and boundary_N
//
//  u(x,y) = g(x,y) on boundary_D
//  du/dn  = h(x,y) on boundary N
//
//
//  Domian
//      [x,y] is in [0,1] X [0,1]
//      f(x,y) = sin(pi*x) * cos(pi *y)
//      Analytic solution is
//          u^a(x,y) = -1/(2*pi^2) * sin(pi*x) * cos(pi*y)
//
//  Boundary Condition
//      Case 1
//          Diriclet bondary condition u(x,y) = 0 in x direction q
//          Neumann boundary condition du/dn = 0 in y direction
//
//      Case 2
//          Dirichlet boundary conditioins using analytical solution
//          both x and y directions.
//


#include <stdio.h>
#include <math.h>
#include <string.h>
#include "def.h"

//-----------------------------------
// Poisson Solvers - Jacobi, SOR, CG
//-----------------------------------
void Jacobi(double **p,double dx, double dy, double tol, int *iter,int BC);
void SOR(double **p,double dx, double dy, double tol, double omega,int *iter,int BC);
void Conjugate_Gradient(double **p,double dx, double dy, double tol, int *iter,int BC);

//-----------------------------------
//        Productivity tools
//-----------------------------------
void print_u_array(double **p);
void initialization(double **p);
void write_u(char *dir_nm,char *file_nm, double **p,double dx, double dy);

//-----------------------------------
//        Mathematical tools
//-----------------------------------
void func_anal(double **p, int row_num, int col_num, double dx, double dy);
double func(int i, int j, double dx, double dy);
void error_rms(double **p, double **p_anal, double *err);

int main(void)
{
    char *dir_name ;
    char *file_name ;
    int i,j,k;
    int Nx, Ny, BC, iter;

    double **u;
    double **u_anal;
    double Lx = 1.0, Ly = 1.0;
    double dx, dy, tol, omega, err;


    u      = (double **) malloc(ROW *sizeof(double));
    u_anal = (double **) malloc(ROW *sizeof(double));
    for (i=0;i<ROW;i++)
    {
      u[i]      = (double *) malloc(COL * sizeof(double));
      u_anal[i] = (double *) malloc(COL * sizeof(double));
    }

    //--------------------
    //   Initial setting
    //--------------------
    iter = 0;
    dir_name = "./RESULT/";

    dx = Lx/(ROW-1);
    dy = Ly/(COL-1);

    tol = 1e-6;
    omega = 1.8;
    BC = 1;

    printf("\n");
    printf("Nx : %d, Ny : %d, dx : %f, dy : %f \n",ROW,COL,dx,dy);
    printf("Tolerance : %f, Omega : %f \n",tol, omega);
    //-----------------------------
    //      Analytic Solutions
    //-----------------------------
    file_name = "Analytic_solution.plt";
    func_anal(u_anal,ROW,COL,dx,dy);
    write_u(dir_name,file_name,u_anal,dx,dy);

   //-----------------------------
   //        Jacobi Method
   //-----------------------------
   initialization(u);
   Jacobi(u,dx,dy,tol,&iter,BC);
   error_rms(u,u_anal,&err);
   printf("Jacobi Method - Error : %f, Iteration : %d \n",err,iter);

   file_name = "Jacobi_result.plt";
   write_u(dir_name,file_name,u,dx,dy);

   //-----------------------------
   //         SOR Method
   //-----------------------------
   initialization(u);
   SOR(u,dx,dy,tol,omega,&iter,BC);
   error_rms(u,u_anal,&err);
   printf("SOR Method - Error : %f, Iteration : %d \n",err,iter);

   file_name = "SOR_result.plt";
   write_u(dir_name,file_name,u,dx,dy);

   //-----------------------------
   //  Conjugate Gradient Method
   //-----------------------------
   initialization(u);
   Conjugate_Gradient(u,dx,dy,tol,&iter,BC);
   error_rms(u,u_anal,&err);
   printf("CG method - Error : %f, Iteration : %d \n",err,iter);

   file_name = "CG_result.plt";
   write_u(dir_name,file_name,u,dx,dy);


    return 0;
}

void print_u_array(double **p)
{
    int i,j;
    for (i=0;i<ROW;i++){
        for (j=0;j<COL;j++){
            printf("%d %d %lf \n",i,j,p[i][j]); }}

}

void initialization(double **p)
{
    int i,j;
    for (i=0;i<ROW;i++){
        for (j=0;j<COL;j++){
            p[i][j] = 0; }}

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


double func(int i,int j,double dx,double dy)
{
    return sin(pi*i*dx)*cos(pi*j*dy);
}

void func_anal(double **p, int row_num, int col_num, double dx, double dy)
{
    int i,j;
    for (i=0;i<row_num;i++){
        for (j=0;j<col_num;j++){
            p[i][j] = -1/(2*pow(pi,2))*sin(pi*i*dx)*cos(pi*j*dy); }}
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
