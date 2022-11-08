//#include <mpi.h>
#include <cstdio>
#include <chrono>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <time.h>

void serial_baseline(int n) {
    for (int i = 0; i < n; i++) {
        // x_i = (x_(i-1)*a +1*b) % p
    }
}

void serial_matrix() {
    continue;
}

void matrix_mult_2x2(std::vector<int> A, std::vector<int>B, std::vector<int>RecvBuff, int P) {
    RecvBuff[0] = ((A.at(0) * B.at(0)) + (A.at(1) * B.at(3))) % P;
    RecvBuff[1] = ((A.at(0) * B.at(1)) + (A.at(1) * B.at(3))) % P;
    RecvBuff[2] = ((A.at(2) * B.at(0)) + (A.at(3) * B.at(2))) % P;
    RecvBuff[3] = ((A.at(2) * B.at(1)) + (A.at(3) * B.at(3))) % P;
}

int main(int argc, char *argv[]) {
    srand((unsigned int)time(nullptr));

    int rank, p;

    //assert(p>=1); assert(n%p==0); assert(n>p);

    int x0 = atoi(argv[1]); //seed
    int n = atoi(argv[2]);
    int A = atoi(argv[3]); //constant
    int B = atoi(argv[4]); //constant
    int P = atoi(argv[5]); //constant (PRIME)

    std::vector<int> I = new std::vector<int>(4,0);
    I.emplace_back(1); I.emplace_back(0); I.emplace_back(0); I.emplace_back(1);

    std::vector<int> M = new std::vector<int>(4,0);
    M.emplace_back(A); M.emplace_back(0); M.emplace_back(B); M.emplace_back(1);

    MPI_Init(nullptr, nullptr); MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD,&p);

    MPI_Finalize();

    return 0;
}