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
#include "def.h"

//#define ROW 11
//#define COL 11
//#define pi 3.141592
//#define itmax 1000000

void initialization(double(*p)[COL]);
void write_u(double(*p)[COL],double dx, double dy);
double func(int i, int j, double dx, double dy);
            
//void test(double(*p)[COL])
//{
//    int Nx,Ny;
//    Nx = sizeof(p)/sizeof(p[0]) - 1;
//    Ny = sizeof(p[0])/sizeof(p[0][0]) - 1;
//    printf("%d %d \n",Nx,Ny);
//}
void Jacobi(double p[][COL],double dx, double dy, double tol)
{
    int i,j,k,it;
    int Nx,Ny;
    double beta,rms;
    double SUM1,SUM2;
    double p_new[ROW][COL]={0};
    printf(" %d %d \n",COL,ROW);
    
    beta = dx/dy;
    
    for (it=1;it<itmax;it++){
        SUM1 = 0;
        SUM2 = 0;
        
        for (i=1;i<ROW-1;i++){
            for (j=1;j<COL-1;j++){
                p_new[i][j] =  (p[i+1][j]+p[i-1][j]
                                + pow(beta,2) *(p[i][j+1]+p[i][j-1])
                                - dx*dx*func(i,j,dx,dy))/(2*(1+pow(beta,2)));
                
//                printf(" %d %d %f \n",i,j,p_new[i][j]);
            }
        }
        
        //------------------------
        //  Boundary conoditions
        //------------------------
        for (j=0;j<COL;j++){
            p_new[0][j] = 0;
            p_new[ROW-1][j] = 0;
        }
        
        for (i=0;i<ROW;i++) {
            p_new[i][0] = p_new[i][1];
            p_new[i][COL-1] = p_new[i][COL-2];
        }
        
        //------------------------
        //  Convergence Criteria
        //------------------------
        for (i=1;i<ROW-1;i++){
            for (j=1;j<COL-1;j++){
                SUM1 += fabs(p_new[i][j]);
                SUM2 += fabs(p_new[i+1][j] + p_new[i-1][j]
                             + pow(beta,2)*(p_new[i][j+1] + p_new[i][j-1])
                             - (2+2*pow(beta,2))*p_new[i][j]-dx*dx*func(i,j,dx,dy));
            }
        }
        
        if ( SUM2/SUM1 < tol )
            break;
        
        //------------------------
        //         Update
        //------------------------
        for (i=0;i<ROW;i++){
            for (j=0;j<COL;j++){
                p[i][j] = p_new[i][j];}}
        
        printf("Iteration : %d, SUM1 : %f, SUM2 : %f, Ratio : %f \n",it,SUM1,SUM2,SUM2/SUM1);
    }
}

int main(void)
{
  
    int i,j,k;
    int Nx, Ny;
    
    double u[ROW][COL] = {0};
    double Lx = 1.0, Ly = 1.0;
    double dx, dy, tol, omega;
    
    //--------------------
    //   Initial setting 
    //--------------------
    Nx = sizeof(u)/sizeof(u[0]) - 1;
    Ny = sizeof(u[0])/sizeof(u[0][0]) - 1;

    dx = Lx/Nx;
    dy = Ly/Ny;

    tol = 1e-6;
    omega = 1.0;

    initialization(u);
    Jacobi(u,dx,dy,tol);

    printf("\n");
    printf("Nx : %d, Ny : %d, dx : %f, dy : %f \n",Nx,Ny,dx,dy);
    printf("Tolerance : %f, Omega : %f \n",tol, omega);
    
    write_u(u,dx,dy);

    return 0;
}


void initialization(double(*p)[COL])
{
    int i,j;
    for (i=0;i<ROW;i++){
        for (j=0;j<COL;j++){
            p[i][j] = 0; }}
}

void write_u(double(*p)[COL],double dx, double dy)
{
    FILE* stream;
    int i,j;
    stream=fopen("data.plt","w");
    fprintf(stream,"ZONE I=%d J=%d \n",ROW,COL);
    for (i=0;i<ROW;i++){
        for(j=0;j<COL;j++){
            fprintf(stream,"%f %f %f \n",i*dx,j*dy,p[i][j]); }}
    fclose(stream);
}


double func(int i, int j, double dx, double dy)
{
    return sin(pi*i*dx)*cos(pi*j*dy);
}
