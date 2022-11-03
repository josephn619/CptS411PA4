//#include <mpi.h>
#include <cstdio>
#include <chrono>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <time.h>

void Matrix_Mult_TwoXTwo(std::vector<int> A, std::vector<int>B, std::vector<int>RecvBuff, int P)
{
    RecvBuff[0] = ( (A[0] * B[0]) + (A[1] * B[3])) % P;
    RecvBuff[1] = ( (A[0] * B[1]) + (A[1] * B[3])) % P;
    RecvBuff[2] = ( (A[2] * B[0]) + (A[3] * B[2])) % P;
    RecvBuff[3] = ( (A[2] * B[1]) + (A[3] * B[3])) % P;
}

int main(int argc, char *argv[])
{
    srand((unsigned int)time(nullptr));

    int seed = rand() % 10000;

    int rank, p;

    int x0 = atoi(argv[1]);
    int A = atoi(argv[2]);
    int B = atoi(argv[3]);
    int P = atoi(argv[4]);

    std::vector<int> I = new std::vector<int>(4,0);
    I.emplace_back(1);
    I.emplace_back(0);
    I.emplace_back(0);
    I.emplace_back(1);

    std::vector<int> M = new std::vector<int>(4,0);
    M.emplace_back(A);
    M.emplace_back(0);
    M.emplace_back(B);
    M.emplace_back(1);

    MPI_Init(nullptr, nullptr); MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD,&p);

    MPI_Finalize();
}

