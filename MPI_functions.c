#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "def.h"

typedef struct mympi{
    int rank_sur[4];
    int nprocs;
    int myrank;
    int nx_mpi,ny_mpi;
    int mpisize_x,mpisize_y,mpirank_x,mpirank_y;
} MYMPI;

//----------------------------------------------------------------------------//
//                            MPI Setup function                              //
//----------------------------------------------------------------------------//
void mpi_setup(int nx, int ny, int mpi_xsize, int mpi_ysize, MYMPI *mpi_info)
{
    mpi_info->mpisize_x = mpi_xsize;
    mpi_info->mpisize_y = mpi_ysize;
    mpi_info->nx_mpi = nx/mpi_xsize;
    mpi_info->ny_mpi = ny/mpi_ysize;
    mpi_info->mpirank_x = mpi_info->myrank/mpi_ysize;
    mpi_info->mpirank_y = mpi_info->myrank%mpi_ysize;

    if (mpi_info->mpirank_x == 0) mpi_info->rank_sur[0] = -1;
    else mpi_info->rank_sur[0] = mpi_info->myrank - mpi_info->mpisize_y;

    if (mpi_info->mpirank_x == mpi_info->mpisize_x-1) mpi_info->rank_sur[1] = -1;
    else mpi_info->rank_sur[1] = mpi_info->myrank + mpi_info->mpisize_y;

    if (mpi_info->mpirank_y == 0) mpi_info->rank_sur[2] = -1;
    else mpi_info->rank_sur[2] = mpi_info->myrank - 1;

    if (mpi_info->mpirank_y == mpi_info->mpisize_y-1) mpi_info->rank_sur[3] = -1;
    else mpi_info->rank_sur[3] = mpi_info->myrank + 1;
}

//----------------------------------------------------------------------------//
//                        Communication subroutines                           //
//                                              - u is the splitted domain    //
//----------------------------------------------------------------------------//
void send_north(double **u, int nx_mpi, int ny_mpi, MYMPI *mpi_info)
{
    int i;
    double *sendbuf, *recvbuf ;
    MPI_Request req1, req2;
    MPI_Status status1, status2;

    sendbuf = (double *) malloc(nx_mpi*sizeof(double));
    recvbuf = (double *) malloc(nx_mpi*sizeof(double));

    for (i=1;i<=nx_mpi;i++){
      sendbuf[i-1] = u[i][ny_mpi];
    }

    MPI_Isend(sendbuf,nx_mpi,MPI_DOUBLE,mpi_info->rank_sur[3],101,MPI_COMM_WORLD,&req1);
    MPI_Irecv(recvbuf,nx_mpi,MPI_DOUBLE,mpi_info->rank_sur[2],101,MPI_COMM_WORLD,&req2);
    MPI_Wait(&req1,&status1);
    MPI_Wait(&req2,&status2);

    for (i=1;i<=nx_mpi;i++){
      u[i][0] = recvbuf[i-1];
    }

    free(sendbuf);
    free(recvbuf);
}

void send_south(double **u, int nx_mpi, int ny_mpi, MYMPI *mpi_info)
{
    int i;
    double *sendbuf, *recvbuf ;
    MPI_Request req1, req2;
    MPI_Status status1, status2;

    sendbuf = (double *) malloc(nx_mpi*sizeof(double));
    recvbuf = (double *) malloc(nx_mpi*sizeof(double));

    for (i=1;i<=nx_mpi;i++){
      sendbuf[i-1] = u[i][1];
    }

    MPI_Isend(sendbuf,nx_mpi,MPI_DOUBLE,mpi_info->rank_sur[2],102,MPI_COMM_WORLD,&req1);
    MPI_Irecv(recvbuf,nx_mpi,MPI_DOUBLE,mpi_info->rank_sur[3],102,MPI_COMM_WORLD,&req2);
    MPI_Wait(&req1,&status1);
    MPI_Wait(&req2,&status2);

    for (i=1;i<=nx_mpi;i++){
      u[i][ny_mpi+1] = recvbuf[i-1];
    }

    free(sendbuf);
    free(recvbuf);
}

void send_west(double **u, int nx_mpi, int ny_mpi, MYMPI *mpi_info)
{
    int j;
    double *sendbuf, *recvbuf ;
    MPI_Request req1, req2;
    MPI_Status status1, status2;

    sendbuf = (double *) malloc(ny_mpi*sizeof(double));
    recvbuf = (double *) malloc(ny_mpi*sizeof(double));

    for (j=1;j<=ny_mpi;j++){
      sendbuf[j-1] = u[nx_mpi][j];
    }

    MPI_Isend(sendbuf,ny_mpi,MPI_DOUBLE,mpi_info->rank_sur[1],103,MPI_COMM_WORLD,&req1);
    MPI_Irecv(recvbuf,ny_mpi,MPI_DOUBLE,mpi_info->rank_sur[0],103,MPI_COMM_WORLD,&req2);
    MPI_Wait(&req1,&status1);
    MPI_Wait(&req2,&status2);

    for (j=1;j<=ny_mpi;j++){
      u[0][j]= recvbuf[j-1];
    }

    free(sendbuf);
    free(recvbuf);
}

void send_east(double **u, int nx_mpi, int ny_mpi, MYMPI *mpi_info)
{
    int j;
    double *sendbuf, *recvbuf ;
    MPI_Request req1, req2;
    MPI_Status status1, status2;

    sendbuf = (double *) malloc(ny_mpi*sizeof(double));
    recvbuf = (double *) malloc(ny_mpi*sizeof(double));

    for (j=1;j<=ny_mpi;j++){
      sendbuf[j-1] = u[1][j];
    }

    MPI_Isend(sendbuf,ny_mpi,MPI_DOUBLE,mpi_info->rank_sur[0],104,MPI_COMM_WORLD,&req1);
    MPI_Irecv(recvbuf,ny_mpi,MPI_DOUBLE,mpi_info->rank_sur[1],104,MPI_COMM_WORLD,&req2);
    MPI_Wait(&req1,&status1);
    MPI_Wait(&req2,&status2);

    for (j=1;j<=ny_mpi;j++){
      u[nx_mpi+1][j] = recvbuf[j-1];
    }

    free(sendbuf);
    free(recvbuf);
}
