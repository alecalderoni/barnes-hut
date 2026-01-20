# Barnes–Hut 3D N-Body Simulation

This project implements a three-dimensional gravitational N-body simulator based on the Barnes–Hut algorithm, developed for the numerical study of astrophysical systems such as galaxies. The Barnes–Hut method reduces the computational complexity of the N-body problem from O(N²) to O(N log N) by hierarchically grouping distant particles and approximating their gravitational interaction through their center of mass. In three dimensions, this is achieved using an octree spatial decomposition.

The simulation is fully configurable through header files located in the `include` directory, allowing physical and numerical parameters to be modified without changing the core source code. In particular, galaxy models are defined and loaded in `include/models.h`, where it is possible to configure galaxy properties, load predefined models, or define new ones.

Galaxies are generated using the Eddington inversion method, which produces particle distributions consistent with an isotropic distribution function in dynamical equilibrium. This allows the construction of realistic initial conditions suitable for long-term gravitational evolution and stability studies.

The project also includes a data analysis stage aimed at studying the physical evolution of the system. Typical analyses include total energy and energy conservation, spatial and velocity distributions of particles, density profiles, and the time evolution of relevant physical observables. Simulation outputs can be exported for further visualization or external post-processing.

The code is compiled using the provided Makefile. After compilation, the simulation can be executed directly from the command line. The project is designed to support simulations with a large number of particles and serves as a foundation for further extensions such as parallelization, real-time visualization, or GPU acceleration.

Project structure:
.
├── src/             Source code  
├── include/         Headers and configuration  
│   └── models.h     Galaxy models  
├── data/            Simulation outputs  
├── analysis/        Data analysis tools  
├── Makefile  
├── README.md  

Build and run:
make  
./simulation  

(The executable name may vary depending on the Makefile configuration.)

This project was developed for educational and research purposes in computational physics, with a focus on numerical methods for gravitational dynamics and galaxy modeling.

References: J. Barnes and P. Hut, A Hierarchical O(N log N) Force-Calculation Algorithm; A. S. Eddington, The Distribution of Stars in Globular Clusters; J. Binney and S. Tremaine, Galactic Dynamics.

Author: Alessandro Calderoni, Physics Student.

License: MIT License.
