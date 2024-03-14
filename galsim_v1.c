#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define EPSILON_O 1e-3
#define CHUNK_SIZE 8
#define INITIAL_LB 0.0
#define INITIAL_RB 1.0
#define INITIAL_DB 0.0
#define INITIAL_UB 1.0
#define THETA_MAX 0.25

double G;
int N;

typedef struct PNode{
    double pos_x;
    double pos_y;
    double mass;
} PNode;

typedef struct QNode{ 
    unsigned int is_leaf;
    PNode* particle;
    struct QNode* child[4];
    /**        UB
        |-------|-------|
        |   1   |   3   |
     LB |-------|-------| RB
        |   0   |   2   |
        |-------|-------|
            DB
    */
    double LB, RB, DB, UB;
} QNode;


// Function declarations
QNode* create_new_QNode(int index, double LB, double RB, double DB, double UB); // Creates a new QNode with specified boundaries
int insert(QNode* qNode, PNode* particle); // Inserts a particle into the quadtree
void barnesHut(PNode* particle, QNode* qNode, double* fx, double* fy); // Applies the Barnes-Hut algorithm to calculate the forces on particles
void destroy(QNode* root); // Destroys the quadtree and frees memory


void preOrder(QNode* qNode){
    if(qNode == NULL || qNode -> is_leaf){
        return;
    }else{
        qNode->particle->pos_x /= qNode->particle->mass;
        qNode->particle->pos_y /= qNode->particle->mass;
        preOrder(qNode->child[0]);
        preOrder(qNode->child[1]);
        preOrder(qNode->child[2]);
        preOrder(qNode->child[3]);
    }
}

int main(int argc, char* argv[]) {
    double time_tol = omp_get_wtime();

    if (argc != 6) {
        printf("You should enter the following parameters in order:\n");
        printf("N filname nsteps delta_t graphics n_threads\n");
        return 1;
    }

    int N = atoi(argv[1]);
    char* filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int graphics = atoi(argv[5]);
    // int n_threads = atoi(argv[6]);
    int n_threads = 8;

    FILE* data_file = fopen(filename, "rb");
    if (data_file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    PNode particles[N];
    double vx[N];
    double vy[N];
    double brightness[N];

    double fx[N];
    double fy[N];
    double mass_inver[N];

    double acc_x[N];
    double acc_y[N]; 

    for (int i = 0; i < N; i++) {
        fread(&particles[i], sizeof(PNode), 1, data_file);
        mass_inver[i] = 1.0 / particles[i].mass;
        fread(&vx[i], sizeof(double), 1, data_file);
        fread(&vy[i], sizeof(double), 1, data_file);
        fread(&brightness[i], sizeof(double), 1, data_file);
        acc_x[i] = 0.0;
        acc_y[i] = 0.0;
    }

    fclose(data_file);

    // Time integration
    G = 100.0 / N;
    for (int step = 0; step < nsteps; step++) {
        #pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(n_threads)
        for(int i=0; i<N; i++){
            fx[i] = 0.0;
            fy[i] = 0.0;
        }
        // create the tree and initialize the first node
        QNode* qTree = create_new_QNode(-1, INITIAL_LB, INITIAL_RB, INITIAL_DB, INITIAL_UB);
        qTree->particle = &particles[0];

        // build the tree
        for (int i = 1; i < N; i++) {
            if(insert(qTree, &particles[i]) == 1){
                return;
            }
        }

        preOrder(qTree);

        // Force calculate: Barnes-Hut Algorithm
        #pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(n_threads)
        for (int i = 0; i < N; i++) {
            barnesHut(&particles[i], qTree, &fx[i], &fy[i]);
        }

        // update
        #pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(n_threads)
        for (int i = 0; i < N; i++) {
            vx[i] += delta_t * fx[i] * mass_inver[i];
            vy[i] += delta_t * fy[i] * mass_inver[i];

            particles[i].pos_x += delta_t * vx[i];
            particles[i].pos_y += delta_t * vy[i];

            // particles[i].pos_x += delta_t * vx[i] + 0.5 * delta_t * delta_t * acc_x[i];
            // particles[i].pos_y += delta_t * vy[i] + 0.5 * delta_t * delta_t * acc_y[i];
            // double new_acc_x = fx[i] * mass_inver[i];
            // double new_acc_y = fy[i] * mass_inver[i];
            // vx[i] += 0.5 * delta_t * (new_acc_x + acc_x[i]);
            // vy[i] += 0.5 * delta_t * (new_acc_y + acc_y[i]);
            // acc_x[i] = new_acc_x;
            // acc_y[i] = new_acc_y;

            if(particles[i].pos_x < INITIAL_LB || particles[i].pos_x > INITIAL_RB || particles[i].pos_y < INITIAL_DB || particles[i].pos_y > INITIAL_UB){
                printf("At least one particle is out of the region, and the simulation has been terminated.\n");
                exit(0);
            }
        }

        destroy(qTree);
    }

    FILE* rfile = fopen("result.gal", "w");

    if (rfile == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    for (int i = 0; i < N; i++) {
        fwrite(&particles[i], sizeof(PNode), 1, rfile);
        fwrite(&vx[i], sizeof(double), 1, rfile);
        fwrite(&vy[i], sizeof(double), 1, rfile);
        fwrite(&brightness[i], sizeof(double), 1, rfile);
    }

    printf("f_std tests took %7.8f wall seconds.\n", omp_get_wtime() - time_tol);

    return 0;
}

QNode* create_new_QNode(int index, double LB, double RB, double DB, double UB) {
    QNode* new_QNode = malloc(sizeof(QNode));
    if (new_QNode != NULL) {
        for (int i = 0; i < 4; i++) {
            new_QNode->child[i] = NULL;
        }
        new_QNode->is_leaf = 1;
        double mid_x = 0.5 * (LB + RB);
        double mid_y = 0.5 * (DB + UB);
        new_QNode->LB = (index == -1) * LB + (index == 0 || index == 1) * LB + (index == 2 || index == 3) * mid_x;
        new_QNode->RB = (index == -1) * RB + (index == 0 || index == 1) * mid_x + (index == 2 || index == 3) * RB;
        new_QNode->DB = (index == -1) * DB +(index == 0 || index == 2) * DB + (index == 1 || index == 3) * mid_y;
        new_QNode->UB = (index == -1) * UB +(index == 0 || index == 2) * mid_y + (index == 1 || index == 3) * UB;
    } else {
        printf("Memory allocation failed");
    }
    return new_QNode;
}

int insert(QNode* qNode, PNode* particle) {
    double mid_x = 0.5 * (qNode->LB + qNode->RB);
    double mid_y = 0.5 * (qNode->UB + qNode->DB);
    if (qNode->is_leaf) {
        if(particle->pos_x == qNode->particle->pos_x && particle->pos_y == qNode->particle->pos_y){
            printf("Two particles are detected at the same location and the simulation terminates.\n");
            return 1;
        }
        qNode->is_leaf = 0;

        // int index = (qNode->particle->pos_x > mid_x) + (qNode->particle->pos_y > mid_y) +
        //             (qNode->particle->pos_x > mid_x && qNode->particle->pos_y < mid_y) +
        //             (qNode->particle->pos_x > mid_x && qNode->particle->pos_y > mid_y);
                
        int index = (qNode->particle->pos_y > mid_y) + 2 * (qNode->particle->pos_x > mid_x); 
        // printf("%d %d", index, index_test);
        qNode->child[index] = create_new_QNode(index, qNode->LB, qNode->RB, qNode->DB, qNode->UB);
        qNode->child[index]->particle = qNode->particle;
        qNode->particle = malloc(sizeof(PNode));
    
        qNode->particle->mass = qNode->child[index]->particle->mass;
        qNode->particle->pos_x = qNode->child[index]->particle->mass * qNode->child[index]->particle->pos_x;
        qNode->particle->pos_y = qNode->child[index]->particle->mass * qNode->child[index]->particle->pos_y;
    }


    qNode->particle->pos_x += particle->mass * particle->pos_x;
    qNode->particle->pos_y += particle->mass * particle->pos_y;
    qNode->particle->mass += particle->mass;

    // int index = (particle->pos_x > mid_x) + (particle->pos_y > mid_y) +
    //             (particle->pos_x > mid_x && particle->pos_y < mid_y) +
    //             (particle->pos_x > mid_x && particle->pos_y > mid_y);

    int index = (particle->pos_y > mid_y) + 2 * (particle->pos_x > mid_x); 

    // printf(" %d %d\n", index, index_test);

    if (qNode->child[index] == NULL) {
        qNode->child[index] = create_new_QNode(index, qNode->LB, qNode->RB, qNode->DB, qNode->UB);
        qNode->child[index]->particle = particle;
    } else {
        insert(qNode->child[index], particle);
    }
    return 0;
}

void barnesHut(PNode* particle, QNode* qNode, double* fx, double* fy) {
    double theta;

    if (!qNode || qNode->particle == particle) {
        return;
    }

    if (!qNode->is_leaf) {
        double mid_x = 0.5 * (qNode->LB + qNode->RB);
        double mid_y = 0.5 * (qNode->DB + qNode->UB);
        double width = (qNode->RB - qNode->LB);
        theta = width / sqrt((particle->pos_x - mid_x) * (particle->pos_x - mid_x) +
                            (particle->pos_y - mid_y) * (particle->pos_y - mid_y));
    }

    if (qNode->is_leaf || theta <= THETA_MAX) {
        double r_x = particle->pos_x - qNode->particle->pos_x;
        double r_y = particle->pos_y - qNode->particle->pos_y;
        double r_squared = r_x * r_x + r_y * r_y;
        double r_plummer = sqrt(r_squared) + EPSILON_O;
        double force_factor = -G * particle->mass * qNode->particle->mass / (r_plummer * r_plummer * r_plummer);
        *fx += force_factor * r_x;
        *fy += force_factor * r_y;
        // *fx += -G * particle->mass * qNode->particle->mass / 
        //     pow(sqrt(pow(particle->pos_x - qNode->particle->pos_x, 2) + 
        //     pow(particle->pos_y - qNode->particle->pos_y, 2)) + EPSILON_O, 3) * 
        //     (particle->pos_x - qNode->particle->pos_x);
        // *fy += -G * particle->mass * qNode->particle->mass / 
        //     pow(sqrt(pow(particle->pos_x - qNode->particle->pos_x, 2) + 
        //     pow(particle->pos_y - qNode->particle->pos_y, 2)) + EPSILON_O, 3) * 
        //     (particle->pos_y - qNode->particle->pos_y);
    } else {
        barnesHut(particle, qNode->child[0], fx, fy);
        barnesHut(particle, qNode->child[1], fx, fy);
        barnesHut(particle, qNode->child[2], fx, fy);
        barnesHut(particle, qNode->child[3], fx, fy);
    }
}

void destroy(QNode* root){
    if(root == NULL){
        return;
    }else if(root->is_leaf){
        free(root);
    }else{
        destroy(root->child[0]);
        destroy(root->child[1]);
        destroy(root->child[2]);
        destroy(root->child[3]);
        free(root->particle);
        free(root);
    }
}