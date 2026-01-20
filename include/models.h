#pragma once

#include "config.h"

#define B  (1.5) // pc  
#define RS 10. 
#define c 20


// Galaxy 1 parameters
#define N1 10000
#define M1 (1e3)   // M_sun 
#define R_MIN1 0.001  //pc
#define R_MAX1 15. //pc

// Galaxy 2 parameters
#define N2 0000
#define M2 (0.*1e8)      
#define R_MIN2 0.001  
#define R_MAX2 100.  

#define RHO (0.*M2 / (4 * M_PI * RS*RS*RS * (log(1 + c) - (c / (1 + c)))))  // M_sun / kpc^3


#define N (N1 + N2) // numero di corpi totali

static inline double rho_plummer(double r, double mass, double scale) {
    return 0.75 * mass * scale * scale / (M_PI * pow(scale * scale + r * r, 2.5));
}

static inline double phi_plummer(double r, double mass, double scale) {
    return -G * mass / sqrt(r * r + scale * scale);
}

static inline double rho_nfw(double r, double rho, double rs) {
    return rho * rs / (r * pow(1 + r/rs,2));
}

static inline double phi_nfw(double r, double rho, double rs) {
    return -4 * M_PI * G * rho * rs*rs*rs * log(1 + r/rs) / r;
}



static inline double rho1(double r) { return rho_plummer(r, M1, B); }
static inline double phi1(double r) { return phi_plummer(r, M1, B) + phi_nfw(r, RHO, RS); }

static inline double rho2(double r) { return rho_nfw(r, RHO, RS); }
static inline double phi2(double r) { return phi_plummer(r, M1, B) + phi_nfw(r, RHO, RS);}
