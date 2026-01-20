#include "galactic_dynamics.h"
#include "config.h"
#include "models.h"

//***********************EDDINGTON***************

#define N_DF 10000 // discretizzazione DF

double psi(double r, double (*phi)(double), double r_max) {
    return -phi(r) + phi(r_max);
}
double eps(double r, double v, double (*phi)(double), double r_max) {
    return psi(r, phi, r_max) - 0.5 * v * v;
}

//eddington
void sorting(double *array, double *coupled_array, int dim) {

    for (int i = 0; i < dim - 1; i++) {
        int min_idx = i;
        for (int j = i + 1; j < dim; j++) {
            if (array[j] < array[min_idx]) {
                min_idx = j;
            }
        }

        double tmp = array[i];
        double tmp2 = coupled_array[i];

        array[i] = array[min_idx];
        coupled_array[i] = coupled_array[min_idx];

        array[min_idx] = tmp;
        coupled_array[min_idx] = tmp2;
        
    }
}
void gradient_nonuniform(double *y, double *x, int n, double *dy_dx) {
    if (n == 1) { dy_dx[0] = NAN; return; }
    // bordo sinistro: forward difference
    dy_dx[0] = (y[1] - y[0]) / (x[1] - x[0]);
    // interni: centered difference non-uniforme
    for (int i = 1; i < n - 1; i++) {
        double dx = x[i+1] - x[i-1];
        dy_dx[i] = (y[i+1] - y[i-1]) / dx;
    }
    // bordo destro: backward difference
    dy_dx[n-1] = (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]);
}
double trapezoid(const double *x, const double *y, int n) {
    double I = 0.0;
    for (int i = 0; i < n - 1; i++) {
        double dx = x[i+1] - x[i];
        I += 0.5 * (y[i+1] + y[i]) * dx;
    }
    return I;
}
void eddington(double *eps_array, double *f_array, double (*rho)(double), double (*phi)(double), double r_min, double r_max) {

    int dim = N_DF;
    
    double *r_array = (double *)malloc(dim * sizeof(double));
    double *psi_array = (double *)malloc(dim * sizeof(double));
    double *rho_array = (double *)malloc(dim * sizeof(double));

    double log_min = log(r_min);
    double log_max = log(r_max);
    double step = (log_max - log_min) / (dim - 1);
    for (int i = 0; i < dim; i++) {    
        r_array[i] = exp(log_min + i * step);
        rho_array[i] = rho(r_array[i]);
        psi_array[i] = psi(r_array[i], phi, r_max);
    }

    free(r_array);

    sorting(psi_array, rho_array, dim);

    double *drho_dpsi = (double *)malloc((dim) * sizeof(double));
    double *d2rho_dpsi2 = (double *)malloc((dim) * sizeof(double));

    gradient_nonuniform(rho_array, psi_array, dim, drho_dpsi);
    gradient_nonuniform(drho_dpsi, psi_array, dim, d2rho_dpsi2);
    
    free(rho_array);

    double C = 1.0 / (sqrt(8)*M_PI*M_PI);

    double *integrand = (double *)malloc((dim) * sizeof(double));
    //FILE *out = fopen("test.dat", "w");
    for(int i = 0; i < dim; i++) {
        eps_array[i] = psi_array[i];
        for(int j = 0; j < i; j++) {
            integrand[j] = d2rho_dpsi2[j] / sqrt(eps_array[i] - psi_array[j]);
        }
        f_array[i] = C * (trapezoid(psi_array, integrand, i) + drho_dpsi[0] / sqrt(eps_array[i]));
        if(f_array[i] < 0) f_array[i] = 0;
        //fprintf(out, "%lf %lf %lf\n", eps_array[i], f_array[i], f_th(eps_array[i]));
    }
    //fclose(out);

    free(integrand);
    free(psi_array);
    free(drho_dpsi);
    free(d2rho_dpsi2);
}

//sampling
double interpol(double target, double *x, double *y, int n) {

    if (target <= x[0]) {
        double m = (y[1] - y[0]) / (x[1] - x[0]);
        return y[0] + m * (target - x[0]);
    }

    if (target >= x[n-1]) {
        double m = (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]);
        return y[n-1] + m * (target - x[n-1]);
    }

    // --- Ricerca binaria del segmento ---
    int lo = 0;
    int hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (target < x[mid])
            hi = mid;
        else
            lo = mid;
    }

    // Ora target è tra x[lo] e x[hi]
    double m = (y[hi] - y[lo]) / (x[hi] - x[lo]);
    return y[lo] + m * (target - x[lo]);
}
double estimate_r(double (*rho)(double), double r_max) {
    const int N_GRID = 400;
    double M = 0.0;
    for (int i = 0; i <= N_GRID; ++i) {
        double r = r_max * (double)i / (double)N_GRID;
        double pr = r * r * rho(r);
        if (pr > M) M = pr;
    }
    return 1.05 * M;
}
double estimate_v(double r, double vmax, double *eps_array, double *f_array, double (*phi)(double), double r_max) {
    const int N_GRID = 400;
    double M = 0.0;
    for (int i = 0; i <= N_GRID; ++i) {
        double v = vmax * (double)i / (double)N_GRID;
        double pv = v * v * interpol(eps(r,v,phi,r_max), eps_array, f_array, N_DF);
        if (pv > M) M = pv;
    }
    return 1.05 * M;
}
vector randomUnitDirection() {
    double u = drand48();
    double v = drand48(); 

    double theta = 2.0 * M_PI * u;            
    double z = 2.0 * v - 1.0;               
    double r = sqrt(1.0 - z*z);

    vector dir = { r * cos(theta), r * sin(theta), z };
    return dir;
}
void generateGalaxy(Body *bodies, Galaxy galaxy) {
    double *f_array = (double *)malloc(N_DF * sizeof(double));
    double *eps_array = (double *)malloc(N_DF * sizeof(double));

    eddington(eps_array, f_array, galaxy.rho, galaxy.phi, galaxy.r_min, galaxy.r_max);

    for (int i = 0; i < galaxy.n; i++) {
        fprintf(stdout, "Generazione Galassia %d, N = %d/%d                 \r", galaxy.idx, i, galaxy.n);
        fflush(stdout);

        // sampling posizioni dall'inversa della cumulativa
        double M_r = estimate_r(galaxy.rho, galaxy.r_max);
        double u_r, pr, r;
        do {
            r = galaxy.r_max * drand48();
            u_r = drand48();
            pr = r * r * galaxy.rho(r);
        } while (u_r > pr / M_r);

        vector unit_vector_r = randomUnitDirection();
        bodies[i].r.x = r * unit_vector_r.x + galaxy.r.x;
        bodies[i].r.y = r * unit_vector_r.y + galaxy.r.y;
        bodies[i].r.z = r * unit_vector_r.z + galaxy.r.z;

        //rejection sampling
        double vmax = sqrt(2 * psi(r, galaxy.phi, galaxy.r_max));
        double M_v = estimate_v(r, vmax, eps_array, f_array, galaxy.phi, galaxy.r_max);
        double u_v, pv, v;
        do {
            v = vmax * drand48();
            u_v = drand48();
            pv = v * v * interpol(eps(r, v, galaxy.phi, galaxy.r_max), eps_array, f_array, N_DF);
        } while (u_v > pv / M_v);

        vector unit_vector_v = randomUnitDirection();
        bodies[i].v.x = v * unit_vector_v.x + galaxy.v.x;
        bodies[i].v.y = v * unit_vector_v.y + galaxy.v.y;
        bodies[i].v.z = v * unit_vector_v.z + galaxy.v.z;

        bodies[i].m = galaxy.m / galaxy.n;
    }
}


//************************BARNES-HUT*********************

// indici degli ottanti
#define UNE 0
#define UNO 1
#define USO 2
#define USE 3
#define DNE 4
#define DNO 5
#define DSO 6
#define DSE 7

#define DT   0.0001     // Gyr 
#define TMAX 1.0      // Gyr

#define THETA 0.6     // Barnes–Hut
#define SIZE  50.0     // pc 
#define EPS   0.3   // pc 

#define N_BIN 200
#define R_BIN 15.     // pc

typedef struct Octant {
    vector r;
    double size;
} Octant;

typedef struct Node {
    Body *body;
    vector r_cdm;
    double M;
    Octant region;
    struct Node *children[8];
    int is_a_leaf;
} Node;

typedef struct { long long n; double mean; double m2; } stats_t;

static inline void welford_update(stats_t* s, double x) {
    s->n++;
    double delta  = x - s->mean;
    s->mean      += delta / (double)s->n;
    double delta2 = x - s->mean;
    s->m2        += delta * delta2;
}

static inline stats_t welford_merge(stats_t a, stats_t b) {
    if (b.n == 0) return a;
    if (a.n == 0) return b;
    double delta = b.mean - a.mean;
    stats_t out;
    out.n    = a.n + b.n;
    out.mean = a.mean + delta * (double)b.n / (double)out.n;
    out.m2   = a.m2 + b.m2 + delta*delta * (double)a.n * (double)b.n / (double)out.n;
    return out;
}


//analisi
double kinetic_energy(Body *bodies) {
    double T = 0.0;

    #pragma omp parallel for reduction(+:T) // per evitare race-condition
    for (int i = 0; i < N; i++) {
        double vx = bodies[i].v.x;
        double vy = bodies[i].v.y;
        double vz = bodies[i].v.z;
        T += 0.5 * bodies[i].m * (vx*vx + vy*vy + vz*vz);
    }
    
    return T;
}
double virial_W(Body *bodies) {
    double W = 0.0;

    #pragma omp parallel for reduction(+:W)
    for (int i = 0; i < N; i++) {
        W += bodies[i].m * ( bodies[i].r.x * bodies[i].a.x + bodies[i].r.y * bodies[i].a.y + bodies[i].r.z * bodies[i].a.z);
    }
    return W;
}


Node *createNode(Octant region) {
    Node *node = (Node *)malloc(sizeof(Node));

    node->body = NULL;

    node->r_cdm = (vector){0.,0.,0.};
    node->M = 0.;

    node->region=region;

    for(int i = 0; i < 8; i++) {
        node->children[i] = NULL;
    }

    node->is_a_leaf = 1;

    return node;
}
int chooseOctant(vector r, Body *body) {
    //scelta del quadrante per body rispetto ad r (il centro del quadrante in cui era)

    if (body->r.z >= r.z) {  // Parte "Up"
        if (body->r.y >= r.y) { // Nord
            if (body->r.x >= r.x) return UNE;  // Up Nord Est
            else                  return UNO;  // Up Nord Ovest
        } else { // Sud
            if (body->r.x >= r.x) return USE;  // Up Sud Est
            else                  return USO;  // Up Sud Ovest
        }
    } else {  // Parte "Down"
        if (body->r.y >= r.y) { // Nord
            if (body->r.x >= r.x) return DNE;  // Down Nord Est
            else                  return DNO;  // Down NorD Est
        } else { // Sud
            if (body->r.x >= r.x) return DSE;  // Down Sud Est
            else                  return DSO;  // Down Sud Ovest
        }
    }
}
void insert(Node *node, Body *body) {
    //se il corpo esce dal reticolo, escludilo (per evitare crash); si potrebbe sennò scegliere SIZE dinamicamente
    if(body->r.x < -SIZE * 0.5 || body->r.x > SIZE * 0.5 || body->r.y < -SIZE * 0.5 || body->r.y > SIZE * 0.5 || body->r.z < - SIZE * 0.5 || body->r.z > SIZE * 0.5) return;

    //se il nodo non contiene un corpo ed è una foglia, inserisci il corpo
    if(node->body == NULL && node->is_a_leaf) {
        node->body = body;
        node->M = body->m;
        node->r_cdm = body->r;
        return;
    }

    double halfsize = node->region.size * 0.5; 
    double childHalf = halfsize * 0.5;

    //se il nodo contiene un corpo ed è una foglia, crea gli 8 sottonodi ed inserisci il corpo che era nel nodo nel sottoquadrante giusto
    if(node->body != NULL && node->is_a_leaf) {

        node->children[UNE] = createNode((Octant){ (vector){ node->region.r.x + childHalf, node->region.r.y + childHalf, node->region.r.z + childHalf }, halfsize });
        node->children[UNO] = createNode((Octant){ (vector){ node->region.r.x - childHalf, node->region.r.y + childHalf, node->region.r.z + childHalf }, halfsize });
        node->children[USO] = createNode((Octant){ (vector){ node->region.r.x - childHalf, node->region.r.y - childHalf, node->region.r.z + childHalf }, halfsize });
        node->children[USE] = createNode((Octant){ (vector){ node->region.r.x + childHalf, node->region.r.y - childHalf, node->region.r.z + childHalf }, halfsize });
        node->children[DNE] = createNode((Octant){ (vector){ node->region.r.x + childHalf, node->region.r.y + childHalf, node->region.r.z - childHalf }, halfsize });
        node->children[DNO] = createNode((Octant){ (vector){ node->region.r.x - childHalf, node->region.r.y + childHalf, node->region.r.z - childHalf }, halfsize });
        node->children[DSO] = createNode((Octant){ (vector){ node->region.r.x - childHalf, node->region.r.y - childHalf, node->region.r.z - childHalf }, halfsize });
        node->children[DSE] = createNode((Octant){ (vector){ node->region.r.x + childHalf, node->region.r.y - childHalf, node->region.r.z - childHalf }, halfsize });
        
        insert(node->children[chooseOctant(node->region.r, node->body)], node->body);
        node->body = NULL;
        node->is_a_leaf = 0;
    }

    //inserisci ora il corpo nel nodo corretto
    insert(node->children[chooseOctant(node->region.r, body)], body);


    //aggiorna il CDM e la sua massa
    node->r_cdm.x = (node->r_cdm.x * node->M + body->r.x * body->m) / (node->M + body->m);
    node->r_cdm.y = (node->r_cdm.y * node->M + body->r.y * body->m) / (node->M + body->m);
    node->r_cdm.z = (node->r_cdm.z * node->M + body->r.z * body->m) / (node->M + body->m);

    node->M += body->m;
}
void computeForce(Node *node, Body *body) {
    if (!node) return;

    //se è una foglia e contiene un corpo diverso da quello di cui sto calcolando la forza, calcola la forza 
    if (node->is_a_leaf) {
        if (node->body && node->body != body) {
            double dx = node->r_cdm.x - body->r.x;
            double dy = node->r_cdm.y - body->r.y;
            double dz = node->r_cdm.z - body->r.z;
            double r = sqrt(dx*dx + dy*dy + dz*dz + EPS*EPS);

            body->a.x += G * node->M * dx / (r * r * r);
            body->a.y += G * node->M * dy / (r * r * r);
            body->a.z += G * node->M * dz / (r * r * r);

            body->phi -= G * node->M / r;
        }
        return;
    }

    double dx = node->r_cdm.x - body->r.x;
    double dy = node->r_cdm.y - body->r.y;
    double dz = node->r_cdm.z - body->r.z;
    double r = sqrt(dx*dx + dy*dy + dz*dz + EPS*EPS);

    //se il criterio di BH è soddisfatto, calcola la forza con il CDM
    if (node->region.size / r < THETA) {
        body->a.x += G * node->M * dx / (r * r * r);
        body->a.y += G * node->M * dy / (r * r * r);
        body->a.z += G * node->M * dz / (r * r * r);

        body->phi -= G * node->M / r;

    //sennò scendi nei nodi figli
    } else {
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                computeForce(node->children[i], body);
            }
        }
    }
}
void updatePosition(Body *bodies, int i, FILE *energy, double *E0) {
    // velocity-verlet:
    // x(t_n+1) = x(t_n) + v(t_n)*DT + 0.5 * a(t_n) * DT^2
    // v(t_n+1) = v(t_n) + 0.5 * (a(t_n) + a(t_n+1)) * DT

    #pragma omp parallel for schedule(static)
    for (int j = 0; j < N; j++) {
        if(i != 0) {
            bodies[j].v.x += 0.5 * bodies[j].a.x * DT;
            bodies[j].v.y += 0.5 * bodies[j].a.y * DT;
            bodies[j].v.z += 0.5 * bodies[j].a.z * DT;
        }
    }

    double T = kinetic_energy(bodies);
    double W = virial_W(bodies);
    double U = 0.;

    #pragma omp parallel for reduction(+:U)
    for (int j = 0; j < N; j++) {
        U += 0.5 * bodies[j].phi * bodies[j].m;
    }

    if(i == 0) *E0 = T + U;

    fprintf(energy, "%lf %lf %lf\n", i * DT, (T + U - *E0)/ *E0, 2 * T / -W);

    #pragma omp parallel for schedule(static)
    for (int j = 0; j < N; j++) {
        bodies[j].r.x += bodies[j].v.x * DT + 0.5 * bodies[j].a.x * DT * DT;
        bodies[j].r.y += bodies[j].v.y * DT + 0.5 * bodies[j].a.y * DT * DT;
        bodies[j].r.z += bodies[j].v.z * DT + 0.5 * bodies[j].a.z * DT * DT;

        bodies[j].v.x += 0.5 * bodies[j].a.x * DT;
        bodies[j].v.y += 0.5 * bodies[j].a.y * DT;
        bodies[j].v.z += 0.5 * bodies[j].a.z * DT;
    }
}
void freeTree(Node *node) {
    if(node == NULL) {
        return;
    }

    for(int i = 0; i < 8; i ++) {
        freeTree(node->children[i]);
    }

    free(node);
}
void simulate(Body *bodies) {
    vector origin = {0.,0.,0.};
    omp_set_num_threads(omp_get_max_threads());  // setting del numero di threads al massimo

    FILE *energy = fopen("dati/energy.dat", "w");  // t, |DeltaE / E|, 2T/-U (viriale) 
    FILE *histo = fopen("dati/histo.dat", "w"); // r, N di corpi in r_i +- DR/2, 4pir^2rho (il tutto per ogni frame, quindi N_bin*N_frame colonne)
    FILE *hmr = fopen("dati/hm_radius.dat", "w");
    int r_bar[N_BIN] = {0};
    int r_dm[N_BIN] = {0};
    double E0;
    double std_bar;
    stats_t s_bar;
    stats_t local;
    

    for (int i = 0; i < TMAX / DT; i++) {

        fprintf(stdout, "TIME: %.4lf/%.3lf Gyr                                        \r", i*DT, TMAX);
        fflush(stdout);

        Node *root = createNode((Octant){origin, SIZE});  //radice dell'OctaTree

        //riempimento albero
        for (int j = 0; j < N; j++) {
            insert(root, &bodies[j]);
        }



        //istogramma
        for (int k = 0; k < N_BIN; ++k) { r_bar[k]=0; r_dm[k]=0; } 
        double hm_radius_bar = NAN, hm_radius_dm = NAN, hm_radius_tot = NAN;

        double dr = R_BIN / (double)N_BIN;

        // --- istogramma barioni ---

        #pragma omp parallel for 
        for (int j = 0; j < N1; ++j) {
            double dx = bodies[j].r.x - root->r_cdm.x;
            double dy = bodies[j].r.y - root->r_cdm.y;
            double dz = bodies[j].r.z - root->r_cdm.z;
            double d = sqrt(dx*dx + dy*dy + dz*dz);
            if (d < R_BIN) {
                int k = (int)floor(d / dr);
                if (k >= 0 && k < N_BIN) {
                    #pragma omp atomic
                    r_bar[k]++;
                }
            }
        }

        // --- istogramma DM ---
        #pragma omp parallel for
        for (int j = N1; j < N; ++j) {
            double dx = bodies[j].r.x - root->r_cdm.x;
            double dy = bodies[j].r.y - root->r_cdm.y;
            double dz = bodies[j].r.z - root->r_cdm.z;
            double d = sqrt(dx*dx + dy*dy + dz*dz);
            if (d < R_BIN) {
                int k = (int)floor(d / dr);
                if (k >= 0 && k < N_BIN) {
                    #pragma omp atomic
                    r_dm[k]++;
                }
            }
        }

        s_bar = (stats_t){0, 0.0, 0.0};
        local = (stats_t){0, 0.0, 0.0};

        #pragma omp for nowait
        for (int j = 0; j < N1; ++j) {
            double vx = bodies[j].v.x;
            double vy = bodies[j].v.y;
            double vz = bodies[j].v.z;
            double v  = sqrt(vx*vx + vy*vy + vz*vz); // dispersione del modulo
            welford_update(&local, v);
        }
        #pragma omp critical
        { s_bar = welford_merge(s_bar, local); }


        double var_bar = (s_bar.n > 1) ? s_bar.m2 / (double)(s_bar.n - 1) : NAN;   // varianza campionaria (consigliata)
        std_bar = (var_bar >= 0.0) ? sqrt(var_bar) : NAN;


        long long mbar = 0, mdm = 0;
        for (int k = 0; k < N_BIN; ++k) {
            mbar += r_bar[k];
            mdm  += r_dm[k];

            double r_center = (k + 0.5) * dr;

            double cum_bar = (double)mbar / (double)N1;
            double cum_dm  = (double)mdm  / (double)N2;
            double cum_tot = (double)(mbar * (double)(M1) / N1 + mdm * (double)(M2) / N2) / (double)(M1 + M2); // uguali masse-particella

            if (isnan(hm_radius_bar) && cum_bar >= 0.5) hm_radius_bar = r_center;
            if (isnan(hm_radius_dm)  && cum_dm  >= 0.5) hm_radius_dm  = r_center;
            if (isnan(hm_radius_tot) && cum_tot >= 0.5) hm_radius_tot = r_center;

            double y_bar = (double)(r_bar[k])/(N1*dr);
            double y_dm = (double)(r_dm[k])/(N2*dr);
            double y_tot = y_bar * M1 / (M1 + M2) + y_dm * M2 / (M1 + M2);

            fprintf(histo, "%lf %lf %lf %lf %lf %lf %lf\n", r_center, y_bar, y_dm, y_tot, cum_bar, cum_dm, cum_tot);
        }

        fprintf(hmr, "%lf %lf %lf %lf %lf\n", i*DT, hm_radius_bar, hm_radius_dm, hm_radius_tot, std_bar*9.78*0.0001);
        //







        #pragma omp parallel for schedule(static)
        for (int j = 0; j < N; j++) {
            bodies[j].a.x = bodies[j].a.y = bodies[j].a.z = 0;
            bodies[j].phi = 0.;
        }

        #pragma omp parallel for schedule(dynamic, 8)
        for (int j = 0; j < N; j++) {
            computeForce(root, &bodies[j]);
        }
        
        updatePosition(bodies, i, energy, &E0);

        freeTree(root);
    }

    fclose(energy);
    fclose(histo);
    fclose(hmr);
}
