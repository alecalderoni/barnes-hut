# Barnes–Hut 3D N-Body Simulation

This project is a **three-dimensional N-body gravitational simulator** based on the **Barnes–Hut algorithm**, developed to study the dynamical evolution of astrophysical systems such as **galaxies**.  
It includes **realistic galaxy generation using the Eddington method** and a dedicated **data analysis pipeline**.

---

## 📌 Features

- 🌌 **3D gravitational N-body simulation**
- 🌲 **Barnes–Hut algorithm** with an **octree** spatial decomposition
- ⚙️ Simulation parameters configurable via header files
- 🌠 Custom **galaxy loading and generation**
- 📊 Built-in **data analysis**
- 🚀 Designed for simulations with a large number of particles

---

## 🧠 Barnes–Hut Algorithm

The Barnes–Hut algorithm reduces the computational complexity of the N-body problem from:

\[
\mathcal{O}(N^2) \rightarrow \mathcal{O}(N \log N)
\]

by grouping distant particles and approximating their gravitational influence through their **center of mass**.

In three dimensions, the algorithm relies on an **octree** to recursively subdivide space.

---

## 🌌 Galaxy Models

Galaxies are generated using the **Eddington inversion method**, which produces particle distributions consistent with an isotropic distribution function in dynamical equilibrium.

📁 **Key file:**
