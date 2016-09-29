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

void Jacobi(double(*p)[COL],double dx, double dy, double tol);
void SOR(double(*p)[COL],double dx, double dy, double tol, double omega);

void print_u_array(double (*p)[COL]);
void initialization(double (*p)[COL]);
void write_u(char *file_nm, double(*p)[COL],double dx, double dy);
void func_anal(double(*p)[COL], int row_num, int col_num, double dx, double dy);

double func(int i, int j, double dx, double dy);


int main(void)
{
    char file_path[50];
    char *dir_name ;
    char *file_name ;
    int i,j,k;
    int Nx, Ny;
    
    double u[ROW][COL] = {0};
    double u_anal[ROW][COL] = {0};
    double Lx = 1.0, Ly = 1.0;
    double dx, dy, tol, omega;
    
    //--------------------
    //   Initial setting 
    //--------------------
    
    dir_name = "./RESULT/";
    Nx = sizeof(u)/sizeof(u[0]) - 1;
    Ny = sizeof(u[0])/sizeof(u[0][0]) - 1;

    dx = Lx/Nx;
    dy = Ly/Ny;

    tol = 1e-6;
    omega = 1.8;
    
    printf("\n");
    printf("Nx : %d, Ny : %d, dx : %f, dy : %f \n",Nx,Ny,dx,dy);
    printf("Tolerance : %f, Omega : %f \n",tol, omega);
    
    //---------------------------
    //     Jacobi Method
    //---------------------------
    initialization(u);
    Jacobi(u,dx,dy,tol);
    
    file_name = "Jacobi_result.plt";
    sprintf(file_path,"%s%s",dir_name,file_name);
    write_u(file_path,u,dx,dy);
    
    //---------------------------
    //     SOR Method
    //---------------------------
    initialization(u);
    SOR(u,dx,dy,tol,omega);
    
    file_name = "SOR_result.plt";
    sprintf(file_path,"%s%s",dir_name,file_name);
    write_u(file_path,u,dx,dy);
    
    //---------------------------
    //     Analytic Solutions
    //---------------------------
    file_name = "Analytic_solution.plt";
    func_anal(u_anal,ROW,COL,dx,dy);
    sprintf(file_path,"%s%s",dir_name,file_name);
    write_u(file_path,u_anal,dx,dy);
    
    return 0;
}

void print_u_array(double (*p)[COL])
{
    int i,j;
    for (i=0;i<ROW;i++){
        for (j=0;j<COL;j++){
            printf("%d %d %lf \n",i,j,p[i][j]); }}

}

void initialization(double (*p)[COL])
{
    int i,j;
    for (i=0;i<ROW;i++){
        for (j=0;j<COL;j++){
            p[i][j] = 0; }}

}

void write_u(char *file_nm, double(*p)[COL],double dx, double dy)
{
    FILE* stream;
    int i,j;
    
    stream=fopen(file_nm,"w");
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

void func_anal(double(*p)[COL], int row_num, int col_num, double dx, double dy)
{
    int i,j;
    for (i=0;i<row_num;i++){
        for (j=0;j<col_num;j++){
            p[i][j] = -1/(2*pow(pi,2))*sin(pi*i*dx)*cos(pi*j*dy); }}
}

