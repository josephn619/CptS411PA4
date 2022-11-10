#include <mpi.h>
#include <cstdio>
#include <iostream>
#include <chrono>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <time.h>

void printMatrix(int A[], int rank)
{
    std::cout << "A " << rank << ":\t{" << A[0] << ", " << A[1] <<"}\n\t{" << A[2] << ", " << A[3] << "}" << std::endl;
}

//randoms must be of size n
void serial_baseline(int n, int A, int B, int P, int seed, int randoms[])
{
    randoms[0] = seed;
    
    for (int i = 1; i < n; i++)  
    {
        randoms[i] = ((A * randoms[i-1]) + B) % P;
    }
}

//randoms must be of size n
void serial_matrix(int n, int A, int B, int P, int seed, int randoms[]) 
{
    randoms[0] = seed;
    for (int i = 1; i < n; i++)  
    {
        randoms[i] = ((A * randoms[i-1]) + B) % P;
    }
}

int summation(int A, int n) 
{
    int sum = 0;
    for (int i = 1; i < n; i++) {
        sum += pow(A,i);
    }
    return sum;
}

void calc_for_i(int A, int B, int i, std::vector<int>RecvBuff) 
{
    RecvBuff[0] = pow(A,i);
    RecvBuff[1] = 0;
    RecvBuff[2] = B*summation(A, i);
    RecvBuff[3] = 1;
}

void matrix_mult_2x2(int A[], int B[], int P) 
{
    int RecvBuff[4] = {0,0,0,0};
    RecvBuff[0] = ((A[0] * B[0]) + (A[1] * B[2])) % P;
    RecvBuff[1] = ((A[0] * B[1]) + (A[1] * B[3])) % P;
    RecvBuff[2] = ((A[2] * B[0]) + (A[3] * B[2])) % P;
    RecvBuff[3] = ((A[2] * B[1]) + (A[3] * B[3])) % P;
    A[0] = RecvBuff[0];
    A[1] = RecvBuff[1];
    A[2] = RecvBuff[2];
    A[3] = RecvBuff[3];
}

void vector_transformation(int vector[], int matrix[], int P)
{
    vector[0] = (vector[0]*matrix[0] + vector[1]*matrix[2]) % P; 
    vector[1] = (vector[0]*matrix[1] + vector[1]*matrix[3]) % P;
}

//randoms is n/p where p is number of processors
void parallel_prefix(int n, int A, int B, int P, int seed, int process_matrix[], int rank, int p)
{
    int index = n/p * rank;    
    int matrix_power = 1;
    int I[4]= {1,0,0,1};
    int world_size;
    int SendM[] = {A,0,B,1};
    int RecvM[] = {0,0,0,0};
    int target = 0;
    int i = 1;
    MPI_Status status;
    for (int j = 2; j <= ((p-1)*n/p); j = j<<1)
    {
        target = (rank ^ i);

        MPI_Sendrecv(SendM, 4, MPI_INT, target, 0, RecvM, 4, MPI_INT, target, 0, MPI_COMM_WORLD, &status);
        matrix_mult_2x2(SendM, RecvM, P);

        if ((j) == index)
        {
            process_matrix[0] = SendM[0];
            process_matrix[1] = SendM[1];
            process_matrix[2] = SendM[2];
            process_matrix[3] = SendM[3];
        }

    }
    
    if (rank==0)
    {
        process_matrix[0] = I[0];
        process_matrix[1] = I[1];
        process_matrix[2] = I[2];
        process_matrix[3] = I[3];
    }
    //printMatrix(process_matrix, rank);
}

void calculateRandoms(int process_matrix[], int randoms[], int seed, int rank, int P, int n, int p, int A, int B)
{
    
    int vector[2] = {seed, 1};
    int base[4] = {A,0,B,1};

    for (int i = 0; i < n/p; i++)
    {
        //get vector[0] to have the right value;
        vector_transformation(vector, process_matrix, P);
        randoms[i] = vector[0];
        std::cout << "Rank: " << rank << "\t Random Value: " << randoms[i] << std::endl;
        matrix_mult_2x2(process_matrix, base, P);
    }
    
    
}

int main(int argc, char *argv[]) 
{
    srand((unsigned int)time(nullptr));
    int seed = 0;
    int P = 7919;
    int rank, p, world_size;

    //assert(p>=1); assert(n%p==0); assert(n>p);

    //int x0 = atoi(argv[1]); //seed
    //int n = atoi(argv[2]); //size
    //int A = atoi(argv[3]); //constant
    //int B = atoi(argv[4]); //constant
    //int P = atoi(argv[5]); //constant (PRIME)


    MPI_Init(nullptr, nullptr); MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD,&p);

    int randoms[1024];
    int process_matrix[4] = {0,0,0,0};

    parallel_prefix(1024, 1, 1, 7919, 0, process_matrix, rank, p);
    calculateRandoms(process_matrix, randoms, seed, rank, P, 1024, p, 1, 1);

    MPI_Finalize();

    return 0;
}