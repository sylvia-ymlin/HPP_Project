#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

/**
    We need to rebuild the quadtree in each time integration.
*/

// Define constants
#define EPSILON_O 1e-3
#define CHUNK_SIZE 8

// Particle_Node
typedef struct{
    double pos_x;
    double pos_y;
    double mass;
}PNode;

// QNode
typedef struct{
    bool is_leaf;
    PNode* particle;
    struct QNode* child[4];
}QNode;

// typedef struct{
//     double x_min;
//     double x_max;
//     double y_min;
//     double y_max;
// }region;


int main(int argc, char* argv[]) {
    double time_tol = omp_get_wtime();

    if (argc != 7) {
        printf("You should enter the following parameters in order:\n");
        printf("N filname nsteps delta_t graphics n_threads\n");
        return 1;
    }

    int N = atoi(argv[1]);
    char* filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int graphics = atoi(argv[5]);
    int n_threads = atoi(argv[6]);

    FILE* data_file = fopen(filename, "rb");
    if (data_file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    // malloc
    Particle_Node* 
    double* pos_x = (double*)malloc(N * sizeof(double));
    double* pos_y = (double*)malloc(N * sizeof(double));
    double* mass = (double*)malloc(N * sizeof(double));
    double* vx = (double*)malloc(N * sizeof(double));
    double* vy = (double*)malloc(N * sizeof(double));
    double* brightness = (double*)malloc(N * sizeof(double));

    double* fx = (double*)malloc(N * sizeof(double));
    double* fy = (double*)malloc(N * sizeof(double));
    double* mass_inver = (double*)malloc(N * sizeof(double));

    double* acc

    for (int i = 0; i < N; i++) {
        fread(&pos_x[i], sizeof(double), 1 , data_file);
        fread(&pos_y[i], sizeof(double), 1 , data_file);
        fread(&mass[i], sizeof(double), 1 , data_file);
        mass_inver[i] = 1.0 / mass[i];
        fread(&vx[i], sizeof(double), 1 , data_file);
        fread(&vy[i], sizeof(double), 1 , data_file);
        fread(&brightness[i], sizeof(double), 1, data_file);
    }

    fclose(data_file);
    
    // Time integration
    double G = 100.0 / N;
    for (int step = 0; step < nsteps; step++) {
        int i, j;
        // Reset forces
        #pragma omp parallel for num_threads(n_threads)
        for (i = 0; i < N; i++) {
            fx[i] = 0.0;
            fy[i] = 0.0;
        }

        #pragma omp parallel for reduction(+:fx[:N], fy[:N]) schedule(static, CHUNK_SIZE) num_threads(n_threads)
        for (i = N; i >= 0; i--) {
            double pos_x_i = pos_x[i];
            double pos_y_i = pos_y[i];
            double mass_i = mass[i];
            double sum_i_x = 0.0, sum_i_y= 0.0;

            for (int j = i+1; j < N; j++) {
                double r_x = pos_x_i - pos_x[j];
                double r_y = pos_y_i - pos_y[j];
                double r_squared = r_x * r_x + r_y * r_y;
                double r_plummer = sqrt(r_squared) + EPSILON_O;
                double force_factor = -G * mass_i * mass[j] / (r_plummer * r_plummer * r_plummer);
                double force_ij_x = force_factor * r_x;
                double force_ij_y = force_factor * r_y;

                sum_i_x += force_ij_x;
                sum_i_y += force_ij_y;
                
                fx[j] -= force_ij_x;
                fy[j] -= force_ij_y; 
                
            }
            
            fx[i] += sum_i_x;
            fy[i] += sum_i_y;

        }
        
        #pragma omp parallel for num_threads(n_threads)
        for (i = 0; i < N; i++) {
            double acc_x = fx[i] * mass_inver[i];
            double acc_y = fy[i] * mass_inver[i];

            pos_x[i] += delta_t * vx[i];
            pos_y[i] += delta_t * vy[i];

            vx[i] += delta_t * fx[i] * mass_inver[i];
            vy[i] += delta_t * fy[i] * mass_inver[i];

            

        }
    }

    FILE* rfile = fopen("result.gal", "w");
    if (rfile == NULL) {
        printf("Error opening file!\n");
        return 1;
    }
    
    for (int i = 0; i < N; i++) {
        fwrite(&pos_x[i], sizeof(double), 1 , data_file);
        printf("%lf", pos_x[i]);
        fwrite(&pos_y[i], sizeof(double), 1 , data_file);
        printf("%lf\n", pos_y[i]);
        fwrite(&mass[i], sizeof(double), 1 , data_file);
        fwrite(&vx[i], sizeof(double), 1 , data_file);
        fwrite(&vy[i], sizeof(double), 1 , data_file);
        fwrite(&brightness[i], sizeof(double), 1, data_file);
    }

    free(pos_x);
    free(pos_y);
    free(mass);
    free(mass_inver);
    free(vx);
    free(vy);
    free(fx);
    free(fy);
    free(brightness);

    printf("f_std tests took %7.8f wall seconds.\n", omp_get_wtime() - time_tol);

    return 0;
}


