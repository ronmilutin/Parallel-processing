#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

typedef struct {
    double pos_x;   // x axis location of the body
    double pos_y;	// y axis location of the body
    double vel_x;	// velocity in x axis direction of the body
    double vel_y;	// velocity in y axis direction of the body
} Obj;

int main(int argc, char** argv) {
	#define N 992
	#define T_MAX 1800
	double  start_time, end_time, exec_time=0;
	int i, j, k;
	int t;
	double G = 6.674 * pow(10, -11);  // Gravitational constant
	double M = 2 * pow(10, 30); 
	Obj *local_obj, *global_obj;
	
	double Fx = 0, Fy = 0;
	int rank;
	int size;
	double x_curr, y_curr, vx_curr, vy_curr;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int n_per_process = N / size;
	local_obj  = (Obj*)malloc(n_per_process * sizeof(Obj));
	global_obj = (Obj*)malloc(N * sizeof(Obj));
	
	srand(time(NULL) * rank / size);
	for (i = 0; i < n_per_process; i++) {
		local_obj[i].pos_x = (double)rand() / RAND_MAX;  // x position [0, 1]
		local_obj[i].pos_y = (double)rand() / RAND_MAX;  // y position [0, 1]
		local_obj[i].vel_x = 0.5 * 200 + ((double)rand() / RAND_MAX) * 1 * 200;  // vx [0.5v, 1.5v]
		local_obj[i].vel_y = 0.5 * 200 + ((double)rand() / RAND_MAX) * 1 * 200;  // vy [0.5v, 1.5v]
	}
	
	MPI_Allgather(local_obj, n_per_process * sizeof(Obj), MPI_BYTE, 
                  global_obj, n_per_process * sizeof(Obj), MPI_BYTE, MPI_COMM_WORLD);
	
	if (rank == 0) {
		printf("Inital Galaxy:\n");
		for (i = 0; i < N; i++) {
			printf("Star %d: pos x=%f, pos y=%f, vel x=%f, vel y=%f.\n",
				i, global_obj[i].pos_x, global_obj[i].pos_y, global_obj[i].vel_x, global_obj[i].vel_y);
		}
	}	
	
	start_time = MPI_Wtime();
	for (t = 0; t < T_MAX; t++) {
		for (i = 0; i < n_per_process; i ++) {
			x_curr = local_obj[i].pos_x * 100 * 9 * pow(10, 15);  // x position in m
			y_curr = local_obj[i].pos_y * 100 * 9 * pow(10, 15);  // y position in m
			vx_curr = local_obj[i].vel_x;  // x velocity in km/sec
			vy_curr = local_obj[i].vel_y;  // y velocity in km/sec
				
				Fx = 0, Fy = 0;
				for (j = 0; j < N; j++) {
					if (j != i + rank * n_per_process) {
						double x = global_obj[j].pos_x * 100 * 9 * pow(10, 15);  // x position of other star in m
						double y = global_obj[j].pos_y * 100 * 9 * pow(10, 15);  // y position of other star in m

						double r = sqrt(pow(x - x_curr, 2) + pow(y - y_curr, 2));  // distance between stars in km
						double F = (G * M * M) / (r * r);  // gravitational force between stars in N

						double angle = atan2(y - y_curr, x - x_curr);  // angle of force vector
						Fx += F * cos(angle);
						Fy += F * sin(angle);
					}
				}
			
			double ax = Fx / M;  // x acceleration in m/s^2
			double ay = Fy / M;  // y acceleration in m/s^2

			double dt = pow(10,11);  // time step
			double new_vx = vx_curr + ax * dt;  // updated x velocity in km/sec
			double new_vy = vy_curr + ay * dt;  // updated y velocity in km/sec


			double new_x = (x_curr + new_vx * dt) / (100 * 9 * pow(10, 15));  // new x position [0, 1]
			double new_y = (y_curr + new_vy * dt) / (100 * 9 * pow(10, 15));  // new y position [0, 1]

			if (new_x < 0) {
				new_x = (double)rand() / RAND_MAX;  // reflect back inward
			}
			if (new_x > 1) {
				new_x = (double)rand() / RAND_MAX; // reflect back inward
			}
			if (new_y < 0) {
				new_y = (double)rand() / RAND_MAX;  // reflect back inward
			}
			if (new_y > 1) {
				new_y = (double)rand() / RAND_MAX;  // reflect back inward
			}
		
		
			local_obj[i].pos_x = new_x;
			local_obj[i].pos_y = new_y;
			local_obj[i].vel_x = new_vx;
			local_obj[i].vel_y = new_vy;
		}
		
		MPI_Allgather(local_obj, n_per_process * sizeof(Obj), MPI_BYTE, 
                  global_obj, n_per_process * sizeof(Obj), MPI_BYTE, MPI_COMM_WORLD);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if (t == T_MAX / 2 && rank == 0){
			printf("Middle Galaxy:\n");
			for (i = 0; i < N; i++) {
				printf("Star %d: pos x=%f, pos y=%f, vel x=%f, vel y=%f.\n",
					i, global_obj[i].pos_x, global_obj[i].pos_y, global_obj[i].vel_x, global_obj[i].vel_y);
			}
		} 
		
	}
		
		
	
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0){
		printf("Final Galaxy:\n");
		for (i = 0; i < N; i++) {
			printf("Star %d: pos x=%f, pos y=%f, vel x=%f, vel y=%f.\n",
				i, global_obj[i].pos_x, global_obj[i].pos_y, global_obj[i].vel_x, global_obj[i].vel_y);
	   }
	}   
	if (rank == 0) {
		end_time = MPI_Wtime();
		exec_time= end_time-start_time;
		printf("Computation time: %f seconds\n", exec_time);
	}
	
	free(local_obj);
	free(global_obj);
	
	MPI_Finalize();
	return 0;
}