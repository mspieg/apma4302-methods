#include <stdio.h>
#include <omp.h>

int main() {
    int sum = 0;
    int n = 10; // The upper limit to sum up to

    // Parallel region starts here
    #pragma omp parallel for reduction(+:sum)
    for (int i = 1; i <= n; ++i) {
        sum += i;
        // Optional: Print the thread number and iteration for demonstration
        printf("Thread %d adding %d\n", omp_get_thread_num(), i);
    }

    // After the parallel region
    printf("Total sum of numbers from 1 to %d is %d\n", n, sum);
    printf("Total number of threads used: %d\n", omp_get_max_threads());

    return 0;
}
