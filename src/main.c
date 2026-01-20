#include "galactic_dynamics.h"
#include "models.h"
#include "config.h"

int main() {
    srand48(time(NULL));

    Body *planets = (Body *)malloc(N * sizeof(Body));

    //
    Galaxy galaxy_stars;
    galaxy_stars.idx = 1;
    galaxy_stars.r = (vector){0.,0.,0.};
    galaxy_stars.v = (vector){0.,0.,0.};
    galaxy_stars.m = M1;
    galaxy_stars.n = N1;
    galaxy_stars.r_max = R_MAX1;
    galaxy_stars.r_min = R_MIN1;
    galaxy_stars.phi = phi1;
    galaxy_stars.rho = rho1;

    generateGalaxy(planets, galaxy_stars);
    //
    
    //
    Galaxy galaxy_dm;
    galaxy_dm.idx = 2;
    galaxy_dm.r = (vector){0.,0.,0.};
    galaxy_dm.v = (vector){0.,0.,0.};
    galaxy_dm.m = M2;
    galaxy_dm.n = N2;
    galaxy_dm.r_max = R_MAX2;
    galaxy_dm.r_min = R_MIN2;
    galaxy_dm.phi = phi2;
    galaxy_dm.rho = rho2;

    generateGalaxy(&planets[N1], galaxy_dm);
    //
    

    simulate(planets);

    free(planets);

    return 0;
}
