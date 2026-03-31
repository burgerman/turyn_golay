#include <iostream>
#include <string>
#include <mpi.h>
#include "utilities.h"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int n = (argc > 1) ? std::stoi(argv[1]) : 7;
    int m = n + 1;
    if (rank == 0) std::cout << "Searching n=" << n << std::endl;
    findQuadruple(n, m, rank, size);
    MPI_Finalize();
    return 0;
}
