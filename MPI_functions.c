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
