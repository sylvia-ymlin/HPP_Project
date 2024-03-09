#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdbool.h>

/**
    We need to rebuild the quadtree in each time integration.
*/

// Define constants
#define EPSILON_O 1e-3
#define CHUNK_SIZE 8
#define INITIAL_LB 0.0
#define INITIAL_RB 1.0
#define INITIAL_DB 0.0
#define INITIAL_UB 1.0
#define THETA_MAX 0.45

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
    double LB, RB, DB, UB;
}QNode;

// Region


QNode* create_new_QNode(int index, double LB, double RB, double DB, double UB){
    QNode* new_QNode = malloc(sizeof(QNode));
    if (new_QNode != NULL) {
        new_QNode->is_leaf = true;
        new_QNode->particle = NULL;
        for(int i=0; i<4; i++){
            new_QNode->child[i] = NULL;
        }
        double mid_x = 0.5 * (LB + RB);
        double mid_y = 0.5 * (DB + UB);
        new_QNode->LB = (index == 0 || index == 1) * LB + (index == 2 || index == 3) * mid_x;
        new_QNode->RB = (index == 0 || index == 1) * mid_x + (index == 2 || index == 3) * RB;
        new_QNode->DB = (index == 0 || index == 2) * DB + (index == 1 || index == 3) * mid_y;
        new_QNode->RB = (index == 0 || index == 2) * mid_y + (index == 1 || index == 3) * UB;
    }else{
        printf("Memory allocation failed");
    }
    return new_QNode;
}

PNode* create_new_PNode(double pos_x, double pos_y, double mass){
    PNode* new_PNode = malloc(sizeof(PNode));
    if(new_PNode != NULL){
        new_PNode->pos_x = pos_x;
        new_PNode->pos_y = pos_y;
        new_PNode->mass = mass;
    }else{
        printf("Memory allocation failed");
    }
    return new_PNode;
}

void insert(QNode* qNode, PNode* particle){
    // A subdomain can only hold one particle.
    if(qNode->is_leaf){
        qNode->is_leaf = false;
        /**
            typedef struct{
                bool is_leaf;
                PNode* particle;
                struct QNode* child[4];
                double LB, RB, DB, UB;
            }QNode;
        */
        /**        UB 
            |-------|-------|
            |   1   |   3   |
         LB |-------|-------| RB
            |   0   |   2   |
            |-------|-------|  
                   DB
        */
        // Fission, dividing into new subdomains.
        int index = (qNode->pos_x > mid_x) + (qNode->pos_y > mid_y) +
                    (qNode->pos_x > mid_x && qNode->pos_y<mid_y) + 
                    (qNode->pos_x > mid_x && qNode->pos_y> mid_y);

        qNode->child[index] = create_new_QNode(index, qNode->LB, qNode->RB, qNode->DB, qNode->UB);
        qNode->child[index]->particle = qNode->particle;
        qNode->child[index]->parallel->parent = qNode;
        
        qNode->particle = create_new_PNode();
    }

    int index = (particle->pos_x > mid_x) + (particle->pos_y > mid_y) +
                    (particle->pos_x > mid_x && particle->pos_y<mid_y) + 
                    (particle->pos_x > mid_x && particle->pos_y> mid_y);

    if(qNode->child[id] == NULL){
        qNode->child[id] = create_new_QNode();
        qNode->child[id]-> particle = particle;
    }else{
        insert(qNode->child[id], particle);
    }
}

void postOderTraversal_calculate(QNode* qNode){
    if(qNode->is_leaf){
        return;
    }

    for (int i = 0; i < 4; i++){
        postOderTraversal_calculate(qNode->child[i]);
    }

    for (int i = 0; i < 4; ++i){
        if(qNode->child[i]!=NULL){
            qNode->parallel->pos_x += qNode->child[i]->particle->mass * qNode->child[i]->particle->pos_x;
            qNode->parallel->pos_y += qNode->child[i]->particle->mass * qNode->child[i]->particle->pos_y;
            qNode->particle->mass += qNode->child[i]->particle->mass;
            qNode->parallel->pos_x /= qNode->particle->mass;
            qNode->parallel->pos_y /= qNode->particle->mass;
        }
    }
}

void barnesHut(PNode* particle, QNode* qNode, double* fx, double* fy){
    double theta;
    while(qNode){
        if(qNode->is_leaf == false){
            double mid_x = 0.5*(qNode->LB+qNode->RB);
            double mid_y = 0.5*(qNode->DB+qNode->UB);
            double width = (qNode->RB - qNode->LB);
            theta = width / sqrt((particle->pos_x - mid_x) * (particle->pos_x - mid_x) + 
                                (particle->pos_y - mid_y) * (particle->pos_y - mid_y));
        }
        if(qNode->is_leaf || theta <= THETA_MAX){
            double r_x = particle->pos_x - qNode->particle->pos_x;
            double r_y = particle->pos_y - qNode->particle->pos_y;
            double r_squared = r_x * r_x + r_y * r_y;
            double r_plummer = sqrt(r_squared) + EPSILON_O;
            double force_factor = -G * mass_i * qNode->particle->mass / (r_plummer * r_plummer * r_plummer);
            *fx += force_factor * r_x;
            *fy += force_factor * r_y;
            return;
        }else{
            for(int i=0; i< 4; i++){
                barnesHut(particle, qNode->child[i], fx, fy);
            }
        }
    }
    return;
}

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
    PNode* particles = (PNode*)malloc(N*sizeof(PNode));
    double* vx = (double*)malloc(N * sizeof(double));
    double* vy = (double*)malloc(N * sizeof(double));
    double* brightness = (double*)malloc(N * sizeof(double));

    double* fx = (double*)malloc(N * sizeof(double));
    double* fy = (double*)malloc(N * sizeof(double));
    double* mass_inver = (double*)malloc(N * sizeof(double));

    // double* acc_x = (double*)malloc(N * sizeof(double));
    // double* acc_y = (double*)malloc(N * sizeof(double));

    for (int i = 0; i < N; i++) {
        fread(&particles[i], sizeof(PNode), 1 , data_file);
        mass_inver[i] = 1.0 / particles[i]->mass;
        fread(&vx[i], sizeof(double), 1 , data_file);
        fread(&vy[i], sizeof(double), 1 , data_file);
        fread(&brightness[i], sizeof(double), 1, data_file);
        // acc_x[i] = 0.0;
        // acc_y[i] = 0.0;
    }

    fclose(data_file);
    
    // Time integration
    double G = 100.0 / N;
    for (int step = 0; step < nsteps; step++) {
        // build the tree
        QNode* qTree = malloc(sizeof(QNode));
        qTree->LB = INITIAL_LB;
        qTree->RB = INITIAL_RB;
        qTree->DB = INITIAL_DB;
        qTree->UB = INITIAL_UB;
        qTree->is_leaf = true;
        qTree->particle = particles[0];
        for(int i=0; i<N; i++){
            insert(qTree, particles[i]);
        }

        // calculate the group mass and center position
        // Quadtree postorder traversal
        postOderTraversal_calculate(qTree);
        
        // Force calculate: Barnes-Hut Algorithm
        for(int i=0; i< N; i++){
            barnesHut(particles[i], qTree, &fx[i], &fy[i]);
        }

        // update
        for (i = 0; i < N; i++) {
            // double new_acc_x = fx[i] * mass_inver[i];
            // double new_acc_y = fy[i] * mass_inver[i];

            // pos_x[i] += delta_t * vx[i] + 0.5 * acc_x[i] * delta_t * delta_t;
            // pos_y[i] += delta_t * vy[i] + 0.5 * acc_y[i] * delta_t * delta_t;

            // vx[i] += 0.5 * (acc_x[i] + new_acc_x) * delta_t;
            // vy[i] += 0.5 * (acc_y[i] + new_acc_y) * delta_t;

            // acc_x[i] = new_acc_x;
            // acc_y[i] = new_acc_y;

            vx[i] += delta_t * fx[i] * mass_inver[i];
            vy[i] += delta_t * fy[i] * mass_inver[i];

            particles[i].pos_x += delta_t * vx[i];
            particles[i].pos_y += delta_t * vy[i];
        }
    }

    FILE* rfile = fopen("result.gal", "w");
    if (rfile == NULL) {
        printf("Error opening file!\n");
        return 1;
    }
    
    for (int i = 0; i < N; i++) {
        fwrite(&particles[i], sizeof(PNode), 1 , data_file);
        fwrite(&vx[i], sizeof(double), 1 , data_file);
        fwrite(&vy[i], sizeof(double), 1 , data_file);
        fwrite(&brightness[i], sizeof(double), 1, data_file);
    }

    free(particles);
    free(mass_inver);
    free(vx);
    free(vy);
    free(fx);
    free(fy);
    free(brightness);

    printf("f_std tests took %7.8f wall seconds.\n", omp_get_wtime() - time_tol);

    return 0;
}


