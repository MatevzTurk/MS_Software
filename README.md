# 2D Metropolis Monte Carlo (MMC) and Molecular Dynamics (MD) Simulation Package

This repository contains simulation software for **two-dimensional particle systems** interacting through an **isotropic pair potential**. Both **MMC** and **MD** methods are implemented, with support for both **NVT** and **NPT** statistical ensembles.

---

## Features

### Metropolis Monte Carlo (MC)
- **NVT ensemble**  
- **NPT ensemble**

### Molecular Dynamics (MD)
- Velocity-Verlet integration  
- **NVT ensemble** (Berendsen and Nose-Hoover thermostats are implemented)  
- **NPT ensemble** (Berendsen and Nose-Hoover thermostats and barostats are implemented)

---

## Repository Structure

Each simulation type is provided as a standalone instance containing:

- A `.txt` file with input parameters  
- A `.f90` Fortran source file implementing the simulation
- An example of an output file `Output_example.dat`

Layout:

```
MC_NVT/
    MC_NVT_parametri.txt
    MC_NVT.f90
    Output_example.dat

MC_NPT/
    MC_NPT_parametri.txt
    MC_NPT.f90
    Output_example.dat

MD_NVT/
    MD_NVT_parametri.txt
    MD_NVT.f90
    Output_example.dat

MD_NPT/
    MD_NPT_parametri.txt
    MD_NPT.f90
    Output_example.dat

LICENSE
README.md
```

---

## Requirements

- Fortran compiler (e.g., `gfortran`, `ifort`)

Compile with:

```bash
gfortran simulation_name.f90 -o simulation_name
```

Run with:

```bash
./simulation_name
```

---

## Description

The simulation scripts implement standard algorithms used in statistical mechanics of 2D systems:

- MC simulations use the **Metropolis acceptance criterion**  
- MD simulations use **velocity-Verlet** integration  
- Thermostats or barostats maintain the target ensemble  
- Periodic boundary conditions and isotropic interaction potentials are supported  

---

## License

This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.

---

## Author

Developed by **Matevž Turk**  
Personal webpage: https://matevzturk.github.io

