#include <mpi.h>
#include <cstdio>
#include <chrono>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <time.h>

int main(int argc, char *argv[])
{
    srand((unsigned int)time(nullptr));

    int seed = rand() % 10000;

    int rank, p;

    int x0 = atoi(argv[1]);
    int A = atoi(argv[2]);
    int B = atoi(argv[3]);
    int P = atoi(argv[4]);

    MPI_Init(nullptr, nullptr); MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD,&p);

    MPI_Finalize();
}