#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


/**
    We need to rebuild the quadtree in each time integration.
    ./a.out 2000 ./input_data/ellipse_N_02000.gal 200 1e-5 0
    ./compare_gal_files/compare_gal_files 2000 result.gal ./ref_output_data/ellipse_N_02000_after200steps.gal     
*/

// Define constants
#define EPSILON_O 1e-3
#define CHUNK_SIZE 8
#define INITIAL_LB 0.0
#define INITIAL_RB 1.0
#define INITIAL_DB 0.0
#define INITIAL_UB 1.0
#define THETA_MAX 0.25


// Particle_Node
typedef struct PNode{
    double pos_x;
    double pos_y;
    double mass;
} PNode;

// QNode
typedef struct QNode{
    bool is_leaf;
    PNode* particle;
    struct QNode* child[4];
    double LB, RB, DB, UB;
} QNode;

double G;
int N;
// double THETA_MAX;

// stack
typedef struct SNode
{
    QNode* qNode;
    bool is_visited;
    struct SNode* next;
} SNode;

SNode* init_Stack(){
    SNode* stack = (SNode*)malloc(sizeof(SNode));
    stack->qNode = NULL;
    stack->next = NULL;
    stack->is_visited = false;
    return stack;
}

// with a head to simply the push
void push(QNode* qNode, SNode* stack){
    SNode* new_SNode = (SNode*)malloc(sizeof(SNode));
    new_SNode->qNode = qNode;
    new_SNode->next = stack->next;
    new_SNode->is_visited = false;
    stack->next = new_SNode;
}

QNode* pop(SNode* stack){
    SNode* pop_node = stack->next;
    stack->next = pop_node->next;
    QNode* popped_qNode = pop_node->qNode;
    free(pop_node);
    return popped_qNode; 
}

QNode* create_new_QNode(int index, double LB, double RB, double DB, double UB) {
    QNode* new_QNode = malloc(sizeof(QNode));
    if (new_QNode != NULL) {
        new_QNode->is_leaf = true;
        new_QNode->particle = NULL;
        for (int i = 0; i < 4; i++) {
            new_QNode->child[i] = NULL;
        }
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

PNode* create_new_PNode(double pos_x, double pos_y, double mass) {
    PNode* new_PNode = malloc(sizeof(PNode));
    if (new_PNode != NULL) {
        new_PNode->pos_x = pos_x;
        new_PNode->pos_y = pos_y;
        new_PNode->mass = mass;
    } else {
        printf("Memory allocation failed");
    }
    return new_PNode;
}

// Recursion will cause stack overflow
void insert(QNode* qNode, PNode* particle) {
    while(true){
        // A subdomain can only hold one particle.
        double mid_x = 0.5 * (qNode->LB + qNode->RB);
        double mid_y = 0.5 * (qNode->UB + qNode->DB);
        if (qNode->is_leaf) {
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
            int index = (qNode->particle->pos_x > mid_x) + (qNode->particle->pos_y > mid_y) +
                        (qNode->particle->pos_x > mid_x && qNode->particle->pos_y < mid_y) +
                        (qNode->particle->pos_x > mid_x && qNode->particle->pos_y > mid_y);

            qNode->child[index] = create_new_QNode(index, qNode->LB, qNode->RB, qNode->DB, qNode->UB);
            qNode->child[index]->particle = qNode->particle;
            qNode->particle = NULL;
        }

        int index = (particle->pos_x > mid_x) + (particle->pos_y > mid_y) +
                    (particle->pos_x > mid_x && particle->pos_y < mid_y) +
                    (particle->pos_x > mid_x && particle->pos_y > mid_y);

        if (qNode->child[index] == NULL) {
            qNode->child[index] = create_new_QNode(index, qNode->LB, qNode->RB, qNode->DB, qNode->UB);
            qNode->child[index]->particle = particle;
            break;
        } else {
            qNode = qNode->child[index];
        }
    }

}

// since the deepth of the tree is unpredictable, it is unable to use resursion to implement the traversal

void postOderTraversal_calculate(QNode* qNode) {
    // if (qNode == NULL || qNode->is_leaf) {
    //     return;
    // }

    SNode* stack = init_Stack();
    SNode* current;
    QNode* current_qNode;
    push(qNode, stack);

    
    while (stack->next != NULL) {
        // fetch the top node
        current = stack->next;
        current_qNode = current->qNode;
        
        if(!current->is_visited){
            current->is_visited = true;
            for (int i = 0; i < 4; ++i) {
                if (current_qNode->child[i] && !current_qNode->child[i]->particle) {
                    // all_visited = 0;
                    // // 非叶子结点才能入栈
                    push(current_qNode->child[i], stack);
                }
            }
        }else{
            double pos_x = 0.0, pos_y = 0.0, mass = 0.0;
            for (int i = 0; i < 4; ++i) {
                if (current_qNode->child[i]) {
                    pos_x += current_qNode->child[i]->particle->mass * current_qNode->child[i]->particle->pos_x;
                    pos_y += current_qNode->child[i]->particle->mass * current_qNode->child[i]->particle->pos_y;
                    mass += current_qNode->child[i]->particle->mass;
                }
            }
            current_qNode->particle = create_new_PNode(pos_x / mass, pos_y / mass, mass);
            pop(stack);
        }
    }
    free(stack);
}

// pre-order
void barnesHut(PNode* particle, QNode* qNode, double* fx, double* fy) {
    double theta;

    SNode* stack = init_Stack();
    QNode* current;
    push(qNode, stack);

    while (stack->next != NULL){
        current = pop(stack);

        // group
        if(!current->is_leaf){
            double mid_x = 0.5 * (current->LB + current->RB);
            double mid_y = 0.5 * (current->DB + current->UB);
            double width = (current->RB - current->LB);
            theta = width / sqrt((particle->pos_x - mid_x) * (particle->pos_x - mid_x) +
                                    (particle->pos_y - mid_y) * (particle->pos_y - mid_y));
        }
        // compute the force
        if((current->is_leaf && (current->particle != particle) )|| theta <= THETA_MAX){
            double r_x = particle->pos_x - current->particle->pos_x;
            double r_y = particle->pos_y - current->particle->pos_y;
            double r_squared = r_x * r_x + r_y * r_y;
            double r_plummer = sqrt(r_squared) + EPSILON_O;
            double force_factor = -G * particle->mass * current->particle->mass / (r_plummer * r_plummer * r_plummer);
            // printf("%d <-> %d : %lf\t\t%lf\n", current->particle, particle, force_factor * r_x, force_factor * r_y);
            *fx += force_factor * r_x;
            *fy += force_factor * r_y;              
        }else{
            for(int i=0; i< 4; i++){
                if(current->child[i]){
                    push(current->child[i], stack);
                }
            }
        }
    }   
    free(stack);
}

void destroy(QNode* root){
    if(root == NULL){
        return;
    }else if(root->is_leaf){
        free(root);
    }else{
        for(int i=0; i<4; i++){
            destroy(root->child[i]);
        }
        free(root->particle);
        free(root);
    }
}

int main(int argc, char* argv[]) {
    // double time_tol = omp_get_wtime();

    // if (argc != 8) {
    //     printf("You should enter the following parameters in order:\n");
    //     printf("N filname nsteps delta_t graphics n_threads\n");
    //     return 1;
    // }

    // int N = atoi(argv[1]);
    // char* filename = argv[2];
    // int nsteps = atoi(argv[3]);
    // double delta_t = atof(argv[4]);
    // int graphics = atoi(argv[5]);
    // int n_threads = atoi(argv[6]);
    // double THETA_MAX = atof(argv[7]);

    N = 2000;
    char* filename = "./input_data/ellipse_N_02000.gal";
    int nsteps = 200;
    double delta_t = 1e-5;


    FILE* data_file = fopen(filename, "rb");
    if (data_file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    // malloc
    PNode* particles = (PNode*)malloc(N * sizeof(PNode));
    double* vx = (double*)malloc(N * sizeof(double));
    double* vy = (double*)malloc(N * sizeof(double));
    double* brightness = (double*)malloc(N * sizeof(double));

    double* fx = (double*)malloc(N * sizeof(double));
    double* fy = (double*)malloc(N * sizeof(double));
    double* mass_inver = (double*)malloc(N * sizeof(double));


    // double* acc_x = (double*)malloc(N * sizeof(double));
    // double* acc_y = (double*)malloc(N * sizeof(double));

    for (int i = 0; i < N; i++) {
        fread(&particles[i], sizeof(PNode), 1, data_file);
        mass_inver[i] = 1.0 / particles[i].mass;
        fread(&vx[i], sizeof(double), 1, data_file);
        fread(&vy[i], sizeof(double), 1, data_file);
        fread(&brightness[i], sizeof(double), 1, data_file);
        // acc_x[i] = 0.0;
        // acc_y[i] = 0.0;
    }

    fclose(data_file);

    // Time integration
    G = 100.0 / N;
    for (int step = 0; step < nsteps; step++) {
        for(int i=0; i<N; i++){
            fx[i] = 0.0;
            fy[i] = 0.0;
        }
        // create the tree and initialize the first node
        QNode* qTree = create_new_QNode(-1, INITIAL_LB, INITIAL_RB, INITIAL_DB, INITIAL_UB);
        qTree->particle = &particles[0];

        // build the tree
        for (int i = 1; i < N; i++) {
            insert(qTree, &particles[i]);
        }

        // calculate the group mass and center position
        // Quadtree postorder traversal
        postOderTraversal_calculate(qTree);

        // Force calculate: Barnes-Hut Algorithm
        // #pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(n_threads)
        for (int i = 0; i < N; i++) {
            barnesHut(&particles[i], qTree, &fx[i], &fy[i]);
        }

        // update
        // #pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(n_threads)
        for (int i = 0; i < N; i++) {
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

        destroy(qTree);
    }

    FILE* rfile = fopen("result.gal", "w");
    if (rfile == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    for (int i = 0; i < N; i++) {
        fwrite(&particles[i], sizeof(PNode), 1, rfile);
        // printf("%.16f\t\t\t%.16f\n", particles[i].pos_x, particles[i].pos_y);
        // printf("%.16f\t\t\t%.16f\n", fx[i], fy[i]);
        fwrite(&vx[i], sizeof(double), 1, rfile);
        fwrite(&vy[i], sizeof(double), 1, rfile);
        fwrite(&brightness[i], sizeof(double), 1, rfile);
    }

    // free(particles);
    free(mass_inver);
    free(vx);
    free(vy);
    free(fx);
    free(fy);
    free(brightness);

    // printf("f_std tests took %7.8f wall seconds.\n", omp_get_wtime() - time_tol);

    return 0;
}
