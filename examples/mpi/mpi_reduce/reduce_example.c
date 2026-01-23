#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Each process contributes its rank value to the sum
    int local_value = world_rank;
    int global_sum;

    // print local values
    printf("Process %d contributing value %d\n", world_rank, local_value);
    
    // Perform the reduction operation
    // Root process will hold the reduced result
    MPI_Reduce(&local_value, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // The root process (rank 0) will print the global sum
    if (world_rank == 0) {
        printf("The sum of all ranks is %d\n", global_sum);
    }

    MPI_Finalize();
    return 0;
}
