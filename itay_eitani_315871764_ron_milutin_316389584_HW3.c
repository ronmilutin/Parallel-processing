#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#define GRAINS_PER_RUN 100000
#define HALL_TO_PAPER_RATIO 0.0001
#define NUM_RUNS 128
#define n 10000

int main() {
    
	double sum = 0.0;
	double start_time = omp_get_wtime();
	
    // Perform the independent runs in parallel
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < NUM_RUNS; i++) {
		
        int count_fall = 0; // Number of grains fell in one iteration
		unsigned int thread_id_seed = omp_get_thread_num() + time(NULL); //(unsigned int)time(NULL) was removeeeeeedddd
		
		int Nn;
        for (int j = 0; j < n; j++) {
            // Simulate the dropping of grains of sand
			Nn = GRAINS_PER_RUN - count_fall;
			
            for (int k = 0; k < Nn; k++) {
                double rand_num = ((double)rand_r(&thread_id_seed) / (double)RAND_MAX);
                if (rand_num <= HALL_TO_PAPER_RATIO) {
                    count_fall++;
                }
            }
        }
		Nn = GRAINS_PER_RUN - count_fall;
        // Calculate the estimate for 'e' in this run
        double estimate = GRAINS_PER_RUN / (double)Nn;
        sum += estimate;
    }
	
	double end_time = omp_get_wtime();
	double	run_time = end_time - start_time;
    // Calculate the final estimate for 'e'
    double final_estimate = sum / NUM_RUNS;

    printf("Final estimate for 'e' is: %.12f\n", final_estimate);
	printf("RunTime: %.6f seconds\n", run_time);

    return 0;
}