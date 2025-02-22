#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/********************MACRO*************************/

#define G 1.
#define N 10000
#define TMAX 1000
#define DT 0.1
#define SIZE 100000
#define THETA 0.5 
#define EPS 1.
#define NE 0
#define NO 1
#define SO 2
#define SE 3

/********************STRUCT*************************/

typedef struct vector {
    double x,y;
} vector;

typedef struct Body {
    vector r, v, a;
    double m;
} Body;

typedef struct Quad {
    vector r;
    double size;
} Quad;

typedef struct Node {
    Body *body;
    vector r_cdm;
    double M;
    Quad region;
    struct Node *children[4];
    int is_a_leaf;
} Node;

/********************FUNZIONI*************************/

void generateGalaxy(Body *bodies, vector center, vector velocity, int n, double radius, double mass);
Node *createNode(Quad region);
int chooseQuad(vector r, Body *body);
void insert(Node *node, Body *body);
void computeForce(Node *node, Body *Body);
void updatePosition(Body *body, FILE *out);
void freeTree(Node *root);
void simulate(Body *bodies);

/********************MAIN*************************/

int main() {
    srand48(time(NULL));

    Body planets[N];

    generateGalaxy(&planets[0], (vector){0,0}, (vector){0,0}, N, 200, N);

    simulate(planets);

    return 0;
}

/********************FUNZIONI*************************/

void generateGalaxy(Body *bodies, vector center, vector velocity, int n, double radius, double mass) {
    for (int i = 0; i < n; i++) {
        double r = radius * sqrt(drand48());  // Distribuzione uniforme in area
        double theta = 2 * M_PI * drand48();
        double v = sqrt(G*N/(radius*radius*radius))*r;  // Velocità orbitale kepleriana

        bodies[i].r.x = center.x + r * cos(theta),
        bodies[i].r.y = center.y + r * sin(theta),
        bodies[i].v.x = velocity.x - v * sin(theta),
        bodies[i].v.y = velocity.y + v * cos(theta),
        bodies[i].m = mass/n;
    }
}

/*void generateGalaxy(Body *bodies, vector center, vector velocity, int n, double radius, double mass) {
    for(int i = 0; i < n; i++) {
        bodies[i].r.x = center.x - radius + 2 * radius * drand48();
        bodies[i].r.y = center.y - radius + 2 * radius * drand48();
        bodies[i].v.x = velocity.x;
        bodies[i].v.y = velocity.y;
        bodies[i].m = mass/n;
    }
}*/

Node *createNode(Quad region) {
    Node *node = (Node *)malloc(sizeof(Node));

    node->body = NULL;

    node->r_cdm = (vector){0,0};
    node->M = 0;

    node->region=region;

    for(int i = 0; i < 4; i++) {
        node->children[i] = NULL;
    }

    node->is_a_leaf = 1;

    return node;
}

int chooseQuad(vector r, Body *body) {
    if(body->r.x >= r.x && body->r.y >= r.y) {
        return NE;
    } else if(body->r.x < r.x && body->r.y >= r.y) {
        return NO;
    } else if(body->r.x < r.x && body->r.y < r.y) {
        return SO;
    } else {
        return SE;
    }
}

void insert(Node *node, Body *body) {
    if(node->body == NULL && node->is_a_leaf) {
        node->body = body;
        node->M = body->m;
        node->r_cdm = body->r;
        return;
    }

    double halfsize = node->region.size * 0.5;

    if(node->body != NULL && node->is_a_leaf) {
        node->children[NE] = createNode((Quad){(vector){node->region.r.x + halfsize * 0.5,node->region.r.y + halfsize * 0.5}, halfsize});
        node->children[NO] = createNode((Quad){(vector){node->region.r.x - halfsize * 0.5,node->region.r.y + halfsize * 0.5}, halfsize});
        node->children[SO] = createNode((Quad){(vector){node->region.r.x - halfsize * 0.5,node->region.r.y - halfsize * 0.5}, halfsize});
        node->children[SE] = createNode((Quad){(vector){node->region.r.x + halfsize * 0.5,node->region.r.y - halfsize * 0.5}, halfsize});

        insert(node->children[chooseQuad(node->region.r, node->body)], node->body);
        node->body = NULL;
        node->is_a_leaf = 0;
    }

    insert(node->children[chooseQuad(node->region.r, body)], body);

    node->r_cdm.x = (node->r_cdm.x * node->M + body->r.x * body->m) / (node->M + body->m);
    node->r_cdm.y = (node->r_cdm.y * node->M + body->r.y * body->m) / (node->M + body->m);

    node->M += body->m;
}

void computeForce(Node *node, Body *body) {
    vector dr;
    double dist,k;
    if(node->is_a_leaf && node->body != body) {
        dr = (vector){node->r_cdm.x - body->r.x, node->r_cdm.y - body->r.y};
        dist = sqrt(dr.x * dr.x + dr.y * dr.y) + EPS;
        k = G * node->M / pow(dist,3);

        body->a.x += (k * dr.x);
        body->a.y += (k * dr.y);
    }
    if(!node->is_a_leaf) {
            dr = (vector){node->r_cdm.x - body->r.x, node->r_cdm.y - body->r.y};
            dist = sqrt(dr.x * dr.x + dr.y * dr.y) + EPS;
        if(node->region.size / dist < THETA) {
            
            k = G * node->M / pow(dist,3);

            body->a.x += (k * dr.x);
            body->a.y += (k * dr.y);
        } else {
            for(int i = 0; i < 4; i++) {
                computeForce(node->children[i], body);
            }
        }
    }
}

void updatePosition(Body *body, FILE *out) {
    body->v.x += body->a.x * DT;
    body->v.y += body->a.y * DT;
    body->r.x += body->v.x * DT;
    body->r.y += body->v.y * DT;

    fprintf(out, "%lf %lf\n", body->r.x, body->r.y);
}

void freeTree(Node *node) {
    if(node == NULL) {
        return;
    }

    for(int i = 0; i < 4; i ++) {
        freeTree(node->children[i]);
    }

    free(node);
}

void printQuadTree(Node *node, FILE *out) {
    if(node == NULL) {
        return;
    } else {
        if(!node->is_a_leaf) {
            double half = node->region.size / 2.0;
fprintf(out, "%lf %lf\n", node->region.r.x - half, node->region.r.y - half);
fprintf(out, "%lf %lf\n", node->region.r.x + half, node->region.r.y - half);
fprintf(out, "%lf %lf\n", node->region.r.x + half, node->region.r.y + half);
fprintf(out, "%lf %lf\n", node->region.r.x - half, node->region.r.y + half);
fprintf(out, "%lf %lf\n\n", node->region.r.x - half, node->region.r.y - half);  // Chiude il quadrato
            for(int i = 0; i < 4; i++) {
                printQuadTree(node->children[i], out);
            }
        } else {
            double half = node->region.size / 2.0;
fprintf(out, "%lf %lf\n", node->region.r.x - half, node->region.r.y - half);
fprintf(out, "%lf %lf\n", node->region.r.x + half, node->region.r.y - half);
fprintf(out, "%lf %lf\n", node->region.r.x + half, node->region.r.y + half);
fprintf(out, "%lf %lf\n", node->region.r.x - half, node->region.r.y + half);
fprintf(out, "%lf %lf\n\n", node->region.r.x - half, node->region.r.y - half);  // Chiude il quadrato
        }
    }
}

void simulate(Body *bodies) {
    vector origin = {0,0};
    FILE *out = fopen("dati.dat", "w");
    char filename[50];

    for (int i = 0; i < TMAX / DT; i++) {
        fprintf(stdout, "State:%.1lf%%\r", i*DT/TMAX*100);
        fflush(stdout);
        Node *root = createNode((Quad){origin, SIZE});

        for (int j = 0; j < N; j++) {
            insert(root, &bodies[j]);
        }

        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < N; j++) {
            bodies[j].a.x = bodies[j].a.y = 0;
            computeForce(root, &bodies[j]);
        }

        #pragma omp parallel for schedule(dynamic)
        for(int j = 0; j < N; j++) {
            updatePosition(&bodies[j], out);
        }

        /*sprintf(filename, "quadranti/quadrants_%d.dat", i);
        FILE *quad_out = fopen(filename, "w");
        printQuadTree(root, quad_out);
        fclose(quad_out);*/

        freeTree(root);
    }

    fclose(out);
}