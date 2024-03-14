#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define EPSILON_O 1e-3
#define CHUNK_SIZE 8
#define THETA_MAX 0.25

double G;
int N;

typedef struct PNode {
    double pos_x;
    double pos_y;
    double mass;
} PNode;

typedef struct TNode {  // =>64

    float LB, RB, DB, UB;    //=> 16
    struct TNode* child[4];  // =>32
    PNode* particle;         // => 8
    bool is_leaf;            // => 1, padding => 7

    /**        UB
        |-------|-------|
        |   1   |   3   |
     LB |-------|-------| RB
        |   0   |   2   |
        |-------|-------|
            DB
    */
} TNode;

// Function declarations
typedef struct QNode QNode;  // Add the struct definition for QNode

TNode* create_new_TNode(int index, double LB, double RB, double DB, double UB);  // Creates a new QNode with specified boundaries
int insert(TNode* qNode, PNode* particle);                                       // Inserts a particle into the quadtree
void barnesHut(PNode* particle, TNode* qNode, double* fx, double* fy);           // Applies the Barnes-Hut algorithm to calculate the forces on particles
void destroy(TNode* root);                                                       // Destroys the quadtree and frees memory

void preOrder(TNode* tNode) {
    if (tNode == NULL || tNode->is_leaf) {
        return;
    } else {
        double inverse_mass = 1.0 / tNode->particle->mass;
        tNode->particle->pos_x *= inverse_mass;
        tNode->particle->pos_y *= inverse_mass;
        preOrder(tNode->child[0]);
        preOrder(tNode->child[1]);
        preOrder(tNode->child[2]);
        preOrder(tNode->child[3]);
    }
}

int main(int argc, char* argv[]) {
    // printf("%d\n", sizeof(TNode));
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
    // printf("%lf", particles[0].pos_x);

    // Time integration
    G = 100.0 / N;
    for (int step = 0; step < nsteps; step++) {
#pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(n_threads)
        for (int i = 0; i < N; i++) {
            fx[i] = 0.0;
            fy[i] = 0.0;
        }

        float LB = 0.0;
        float RB = 1.0;
        float DB = 0.0;
        float UB = 1.0;

        // build the tree
        float square1[4] = {LB, 0.5 * (LB + RB), DB, 0.5 * (DB + UB)};
        float square2[4] = {0.5 * (LB + RB), RB, DB, 0.5 * (DB + UB)};
        float square3[4] = {LB, 0.5 * (LB + RB), 0.5 * (DB + UB), UB};
        float square4[4] = {0.5 * (LB + RB), RB, 0.5 * (DB + UB), UB};
        int group1_ofParticles[N];
        int N1 = 0;
        int group2_ofParticles[N];
        int N2 = 0;
        int group3_ofParticles[N];
        int N3 = 0;
        int group4_ofParticles[N];
        int N4 = 0;
        for (int i = 0; i < N; i++) {
            if (particles[i].pos_x >= square1[0] && particles[i].pos_x <= square1[1] &&
                particles[i].pos_y >= square1[2] && particles[i].pos_y <= square1[3]) {
                group1_ofParticles[N1] = i;
                N1++;
            }
            if (particles[i].pos_x >= square2[0] && particles[i].pos_x <= square2[1] &&
                particles[i].pos_y >= square2[2] && particles[i].pos_y <= square2[3]) {
                group2_ofParticles[N2] = i;
                N2++;
            }
            if (particles[i].pos_x >= square3[0] && particles[i].pos_x <= square3[1] &&
                particles[i].pos_y >= square3[2] && particles[i].pos_y <= square3[3]) {
                group3_ofParticles[N3] = i;
                N3++;
            }
            if (particles[i].pos_x >= square4[0] && particles[i].pos_x <= square4[1] &&
                particles[i].pos_y >= square4[2] && particles[i].pos_y <= square4[3]) {
                group4_ofParticles[N4] = i;
                N4++;
            }
        }

        TNode* tTree1;
        TNode* tTree2;
        TNode* tTree3;
        TNode* tTree4;
#pragma omp parallel num_threads(n_threads)
        {
#pragma omp sections
            {
#pragma omp section
                {
                    tTree1 = create_new_TNode(-1, square1[0], square1[1], square1[2], square1[3]);
                    tTree1->particle = &particles[group1_ofParticles[0]];
                    for (int i = 1; i < N1; i++) {
                        insert(tTree1, &particles[group1_ofParticles[i]]);
                    }
                    preOrder(tTree1);
                }

#pragma omp section
                {
                    tTree2 = create_new_TNode(-1, square2[0], square2[1], square2[2], square2[3]);
                    tTree2->particle = &particles[group2_ofParticles[0]];
                    for (int i = 1; i < N2; i++) {
                        insert(tTree2, &particles[group2_ofParticles[i]]);
                    }
                    preOrder(tTree2);
                }

#pragma omp section
                {
                    tTree3 = create_new_TNode(-1, square3[0], square3[1], square3[2], square3[3]);
                    tTree3->particle = &particles[group3_ofParticles[0]];
                    for (int i = 1; i < N3; i++) {
                        insert(tTree3, &particles[group3_ofParticles[i]]);
                    }
                    preOrder(tTree3);
                }

#pragma omp section
                {
                    tTree4 = create_new_TNode(-1, square4[0], square4[1], square4[2], square4[3]);
                    tTree4->particle = &particles[group4_ofParticles[0]];
                    for (int i = 1; i < N4; i++) {
                        insert(tTree4, &particles[group4_ofParticles[i]]);
                    }
                    preOrder(tTree4);
                }
            }
        }

        TNode* tTree = create_new_TNode(-1, 0.0, 1.0, 0.0, 1.0);
        tTree->child[0] = tTree1;
        tTree->child[1] = tTree2;
        tTree->child[2] = tTree3;
        tTree->child[3] = tTree4;
        tTree->is_leaf = 0;
        tTree->particle = malloc(sizeof(PNode));
        tTree->particle->mass = tTree1->particle->mass + tTree2->particle->mass + tTree3->particle->mass + tTree4->particle->mass;
        tTree->particle->pos_x = tTree1->particle->mass * tTree1->particle->pos_x + tTree2->particle->mass * tTree2->particle->pos_x + tTree3->particle->mass * tTree3->particle->pos_x + tTree4->particle->mass * tTree4->particle->pos_x;
        tTree->particle->pos_y = tTree1->particle->mass * tTree1->particle->pos_y + tTree2->particle->mass * tTree2->particle->pos_y + tTree3->particle->mass * tTree3->particle->pos_y + tTree4->particle->mass * tTree4->particle->pos_y;
        tTree->particle->pos_x /= tTree->particle->mass;
        tTree->particle->pos_y /= tTree->particle->mass;

// Force calculate: Barnes-Hut Algorithm
#pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(n_threads)
        for (int i = 0; i < N; i++) {
            barnesHut(&particles[i], tTree, &fx[i], &fy[i]);
        }

// printf("%lu\n", sizeof(tTree->particle));

// update
#pragma omp parallel for schedule(guided, CHUNK_SIZE) num_threads(n_threads)
        // f_std tests took 0.27875203 wall seconds
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

            if (particles[i].pos_x < LB || particles[i].pos_x > RB || particles[i].pos_y < DB || particles[i].pos_y > UB) {
                printf("At least one particle is out of the region, and the simulation has been terminated.\n");
                exit(0);
            }
        }

        destroy(tTree);
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

TNode* create_new_TNode(int index, double LB, double RB, double DB, double UB) {
    TNode* new_TNode = malloc(sizeof(TNode));
    if (new_TNode != NULL) {
        for (int i = 0; i < 4; i++) {
            new_TNode->child[i] = NULL;
        }
        new_TNode->is_leaf = 1;
        double mid_x = 0.5 * (LB + RB);
        double mid_y = 0.5 * (DB + UB);
        new_TNode->LB = (index == -1) * LB + (index == 0 || index == 1) * LB + (index == 2 || index == 3) * mid_x;
        new_TNode->RB = (index == -1) * RB + (index == 0 || index == 1) * mid_x + (index == 2 || index == 3) * RB;
        new_TNode->DB = (index == -1) * DB + (index == 0 || index == 2) * DB + (index == 1 || index == 3) * mid_y;
        new_TNode->UB = (index == -1) * UB + (index == 0 || index == 2) * mid_y + (index == 1 || index == 3) * UB;
    } else {
        printf("Memory allocation failed");
    }
    return new_TNode;
}

int insert(TNode* tNode, PNode* particle) {
    double mid_x = 0.5 * (tNode->LB + tNode->RB);
    double mid_y = 0.5 * (tNode->UB + tNode->DB);
    if (tNode->is_leaf) {
        if (particle->pos_x == tNode->particle->pos_x && particle->pos_y == tNode->particle->pos_y) {
            printf("Two particles are detected at the same location and the simulation terminates.\n");
            return 1;
        }
        tNode->is_leaf = 0;

        // int index = (tNode->particle->pos_x > mid_x) + (tNode->particle->pos_y > mid_y) +
        //             (tNode->particle->pos_x > mid_x && tNode->particle->pos_y < mid_y) +
        //             (tNode->particle->pos_x > mid_x && tNode->particle->pos_y > mid_y);

        int index = (tNode->particle->pos_y > mid_y) + 2 * (tNode->particle->pos_x > mid_x);
        // printf("%d %d", index, index_test);
        tNode->child[index] = create_new_TNode(index, tNode->LB, tNode->RB, tNode->DB, tNode->UB);
        tNode->child[index]->particle = tNode->particle;
        tNode->particle = malloc(sizeof(PNode));

        tNode->particle->mass = tNode->child[index]->particle->mass;
        tNode->particle->pos_x = tNode->child[index]->particle->mass * tNode->child[index]->particle->pos_x;
        tNode->particle->pos_y = tNode->child[index]->particle->mass * tNode->child[index]->particle->pos_y;
    }

    tNode->particle->pos_x += particle->mass * particle->pos_x;
    tNode->particle->pos_y += particle->mass * particle->pos_y;
    tNode->particle->mass += particle->mass;

    // int index = (particle->pos_x > mid_x) + (particle->pos_y > mid_y) +
    //             (particle->pos_x > mid_x && particle->pos_y < mid_y) +
    //             (particle->pos_x > mid_x && particle->pos_y > mid_y);

    int index = (particle->pos_y > mid_y) + 2 * (particle->pos_x > mid_x);

    // printf(" %d %d\n", index, index_test);

    if (tNode->child[index] == NULL) {
        tNode->child[index] = create_new_TNode(index, tNode->LB, tNode->RB, tNode->DB, tNode->UB);
        tNode->child[index]->particle = particle;
    } else {
        insert(tNode->child[index], particle);
    }
    return 0;
}

void barnesHut(PNode* particle, TNode* tNode, double* fx, double* fy) {
    double theta;

    if (!tNode || tNode->particle == particle) {
        return;
    }

    if (!tNode->is_leaf) {
        double mid_x = 0.5 * (tNode->LB + tNode->RB);
        double mid_y = 0.5 * (tNode->DB + tNode->UB);
        double width = (tNode->RB - tNode->LB);
        theta = width / sqrt((particle->pos_x - mid_x) * (particle->pos_x - mid_x) +
                             (particle->pos_y - mid_y) * (particle->pos_y - mid_y));
    }

    if (tNode->is_leaf || theta <= THETA_MAX) {
        double r_x = particle->pos_x - tNode->particle->pos_x;
        double r_y = particle->pos_y - tNode->particle->pos_y;
        double r_squared = r_x * r_x + r_y * r_y;
        double r_plummer = sqrt(r_squared) + EPSILON_O;
        double force_factor = -G * particle->mass * tNode->particle->mass / (r_plummer * r_plummer * r_plummer);
        *fx += force_factor * r_x;
        *fy += force_factor * r_y;
        // *fx += -G * particle->mass * tNode->particle->mass /
        //     pow(sqrt(pow(particle->pos_x - tNode->particle->pos_x, 2) +
        //     pow(particle->pos_y - tNode->particle->pos_y, 2)) + EPSILON_O, 3) *
        //     (particle->pos_x - tNode->particle->pos_x);
        // *fy += -G * particle->mass * tNode->particle->mass /
        //     pow(sqrt(pow(particle->pos_x - tNode->particle->pos_x, 2) +
        //     pow(particle->pos_y - tNode->particle->pos_y, 2)) + EPSILON_O, 3) *
        //     (particle->pos_y - tNode->particle->pos_y);
    } else if (tNode->is_leaf == 0) {
        barnesHut(particle, tNode->child[0], fx, fy);
        barnesHut(particle, tNode->child[1], fx, fy);
        barnesHut(particle, tNode->child[2], fx, fy);
        barnesHut(particle, tNode->child[3], fx, fy);
    }
}

void destroy(TNode* root) {
    if (root == NULL) {
        return;
    } else if (root->is_leaf) {
        free(root);
    } else {
        destroy(root->child[0]);
        destroy(root->child[1]);
        destroy(root->child[2]);
        destroy(root->child[3]);
        free(root->particle);
        free(root);
    }
}