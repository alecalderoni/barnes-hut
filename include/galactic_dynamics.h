#pragma once

#include "config.h"

typedef struct vector {
    double x, y, z;
} vector;

typedef struct Body {
    vector r, v, a;
    double m, phi;
} Body;

typedef struct Galaxy {
    int idx;
    vector r;
    vector v;
    double m;
    int n;
    double (*rho)(double);
    double (*phi)(double);
    double r_min;
    double r_max;
} Galaxy;

//eddington
void eddington(double *eps_array, double *f_array, double (*rho)(double), double (*phi)(double), double r_min, double r_max);
void generateGalaxy(Body *bodies, Galaxy galaxy);

//BarnesHut
void simulate(Body *bodies);
