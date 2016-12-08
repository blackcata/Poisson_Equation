#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "def.h"
double func(int i, int j, double dx, double dy);

//----------------------------------------------------------------------------//
//                     MPI TYPE, MPI setting functions                        //
//----------------------------------------------------------------------------//
typedef struct mympi{
    int rank_sur[4];
    int nprocs;
    int myrank;
    int nx_mpi,ny_mpi;
    int mpisize_x,mpisize_y,mpirank_x,mpirank_y;
} MYMPI;

void mpi_setup(int nx, int ny, int mpi_xsize, int mpi_ysize, MYMPI *mpi_info);
void send_north(double **u, int nx_mpi, int ny_mpi, MYMPI *mpi_info);
void send_south(double **u, int nx_mpi, int ny_mpi, MYMPI *mpi_info);
void send_west(double **u,  int nx_mpi, int ny_mpi, MYMPI *mpi_info);
void send_east(double **u,  int nx_mpi, int ny_mpi, MYMPI *mpi_info);

//----------------------------------------------------------------------------//
//                                                                            //
//                           Main Jacobi Subroutine                           //
//                                                                            //
//----------------------------------------------------------------------------//
void Jacobi(double **p,double dx, double dy, double tol,
            double *tot_time, int *iter,int BC,
            char* file_name, char* dir_name,int write_type)
{
    int i,j,k,it;
    int Nx,Ny,ista,iend,jsta,jend;
    double beta,rms;
    double SUM1 = 0,SUM2 = 0,SUM1_loc,SUM2_loc;
    double *p_tmp;
    double **p_new, **p_loc;
    time_t start_t =0, end_t =0;
    char loc_name[50],file_path[50];
    FILE* stream;
    MYMPI mpi_info;

    start_t = clock();
    beta = dx/dy;

    //----------------------------------------------------------------------//
    //                             MPI Setting                              //
    //----------------------------------------------------------------------//
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_info.nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_info.myrank);

    // This two variables has to be assigned at Poisson_Equation.c
    int mpi_xsize = 2;
    int mpi_ysize = 2;
    mpi_setup(ROW,COL,mpi_xsize,mpi_ysize,&mpi_info);

    ista = mpi_info.mpirank_x*mpi_info.nx_mpi;
    iend = ista + mpi_info.nx_mpi -1;
    jsta = mpi_info.mpirank_y*mpi_info.ny_mpi;
    jend = jsta + mpi_info.ny_mpi -1;

    printf("Myrank : %d\n",mpi_info.myrank);
    printf("mpirank_x : %d, mpirank_y : %d\n",mpi_info.mpirank_x,mpi_info.mpirank_y);
    printf("nx_mpi : %d, ny_mpi : %d\n",mpi_info.nx_mpi,mpi_info.ny_mpi);
    printf("mpisize_x : %d, mpisize_y : %d\n",mpi_info.mpisize_x,mpi_info.mpisize_y);
    printf("(%d,%d) X (%d,%d)\n",ista,iend,jsta,jend);
    printf("e_rank : %d, w_rank : %d, s_rank : %d, n_rank : %d \n",
    mpi_info.rank_sur[0],mpi_info.rank_sur[1],mpi_info.rank_sur[2],mpi_info.rank_sur[3]);
    printf("\n");

    //------------------------------------------------------------------------//
    //                          Memory Allocation                             //
    //------------------------------------------------------------------------//
    p_tmp = (double *) malloc((COL*ROW/mpi_info.nprocs)*sizeof(double));
    p_loc = (double **) malloc((ROW/mpi_xsize+2)*sizeof(double));
    p_new = (double **) malloc((ROW/mpi_xsize+2)*sizeof(double));

    for (i=0;i<ROW/mpi_xsize+2;i++){
      p_loc[i] = (double *) malloc((COL/mpi_xsize+2)*sizeof(double));
      p_new[i] = (double *) malloc((COL/mpi_ysize+2)*sizeof(double));
    }

    for (i=ista;i<iend+1;i++){
      for (j=jsta;j<jend+1;j++){
        p_loc[i-ista+1][j-jsta+1] = p[i][j];
      }
    }

    //------------------------------------------------------------------------//
    //                       Main Loop of Jacobi method                       //
    //------------------------------------------------------------------------//
    for (it=1;it<itmax;it++){
        SUM1_loc = 0;
        SUM2_loc = 0;

        send_north(p_loc, mpi_info.nx_mpi, mpi_info.ny_mpi, &mpi_info);
        send_south(p_loc, mpi_info.nx_mpi, mpi_info.ny_mpi, &mpi_info);
        send_west(p_loc,  mpi_info.nx_mpi, mpi_info.ny_mpi, &mpi_info);
        send_east(p_loc,  mpi_info.nx_mpi, mpi_info.ny_mpi, &mpi_info);

        for (i=1;i<=mpi_info.nx_mpi;i++){
            for (j=1;j<=mpi_info.ny_mpi;j++){
            p_new[i][j] =  (p_loc[i+1][j]+p_loc[i-1][j]
                            + pow(beta,2) *(p_loc[i][j+1]+p_loc[i][j-1])
                            - dx*dx*func(i+ista-1,j+jsta-1,dx,dy))
                                                           /(2*(1+pow(beta,2)));
            }
        }

        //--------------------------------------------------------------------//
        //                        Boundary conditions                         //
        //--------------------------------------------------------------------//

        //--------------------------------------//
        //           Boundary - Case 1          //
        //--------------------------------------//
        if (BC == 1){
           if(ista == 0){
             for (j=1;j<=mpi_info.ny_mpi;j++){
               p_new[1][j] = 0;
             }
           }

           if(iend == ROW-1){
             for (j=1;j<=mpi_info.ny_mpi;j++){
               p_new[mpi_info.nx_mpi][j] = 0;
             }
           }

           if(jsta == 0){
             for (i=1;i<=mpi_info.nx_mpi;i++){
               p_new[i][1] = p_new[i][2];
             }
           }

           if(jend == COL-1){
             for (i=1;i<=mpi_info.nx_mpi;i++){
               p_new[i][mpi_info.ny_mpi] = p_new[i][mpi_info.ny_mpi-1];
             }
           }

        }
        //--------------------------------------//
        //           Boundary - Case 2          //
        //--------------------------------------//
        else if (BC ==2){
          for (j=0;j<COL;j++){
              p_new[0][j] = -1/(2*pow(pi,2))*func(0,j,dx,dy);
              p_new[ROW-1][j] = -1/(2*pow(pi,2))*func(ROW-1,j,dx,dy);
          }

          for (i=0;i<ROW;i++) {
              p_new[i][0] = -1/(2*pow(pi,2))*func(i,0,dx,dy);
              p_new[i][COL-1] = -1/(2*pow(pi,2))*func(i,COL-1,dx,dy);
          }
        }

        //--------------------------------------------------------------------//
        //                        Convergence Criteria                        //
        //--------------------------------------------------------------------//
        for (i=2;i<mpi_info.nx_mpi;i++){
            for (j=2;j<mpi_info.ny_mpi;j++){
                SUM1_loc += fabs(p_new[i][j]);
                SUM2_loc += fabs(p_new[i+1][j] + p_new[i-1][j]
                                 + pow(beta,2)*(p_new[i][j+1] + p_new[i][j-1])
                                 - (2+2*pow(beta,2))*p_new[i][j]
                                 - dx*dx*func(i+ista-1,j+jsta-1,dx,dy));
            }
        }

        MPI_Allreduce(&SUM1_loc,&SUM1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&SUM2_loc,&SUM2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        if ( SUM2/SUM1 < tol ){

            *iter = it;
            end_t = clock();
            *tot_time = (double)(end_t - start_t)/(CLOCKS_PER_SEC);

            //----------------------------------------------------------------//
            //                      Writing functions                         //
            //----------------------------------------------------------------//
            sprintf(loc_name,"%d.%s",mpi_info.myrank,file_name);
            sprintf(file_path,"%s%s",dir_name,loc_name);
            printf("%s %s\n",loc_name,file_path);

            stream=fopen(file_path,"w");
            fprintf(stream,"ZONE I=%d J=%d \n",mpi_info.nx_mpi,mpi_info.ny_mpi);
            for (i=1;i<=mpi_info.nx_mpi;i++){
                for (j=1;j<=mpi_info.ny_mpi;j++){
                    fprintf(stream,"%f %f %f \n",(i+ista-1)*dx,(j+jsta-1)*dy,
                                                 p_new[i][j]);
                }
            }
            fclose(stream);

            free(p_tmp);
            free(p_new);
            break;
        }
        if(mpi_info.myrank==0)
         printf("Iteration : %d, SUM1 : %f, SUM2 : %f, Ratio : %f \n",
                                                       it,SUM1,SUM2,SUM2/SUM1);

        //--------------------------------------------------------------------//
        //                               Update                               //
        //--------------------------------------------------------------------//
        for (i=1;i<=mpi_info.nx_mpi;i++){
            for (j=1;j<=mpi_info.ny_mpi;j++){
              p_loc[i][j] = p_new[i][j];
            }
        }
    }
}
