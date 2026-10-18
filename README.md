# Charged Particle Drift Simulations

![C++](https://img.shields.io/badge/Language-C%2B%2B17-blue)
![Python](https://img.shields.io/badge/Analysis-Python%203-yellow)
![CMake](https://img.shields.io/badge/Build-CMake%203.20%2B-green)
![Integrator](https://img.shields.io/badge/Integrator-RK4-orange)
![Physics](https://img.shields.io/badge/Physics-Single--Particle%20Drifts-red)
![License](https://img.shields.io/badge/License-MIT-lightgrey)

> A compact full-orbit charged-particle simulator for reproducing and validating the main single-particle drifts discussed in introductory plasma physics: **E×B**, **grad-B**, **curvature**, **nonuniform-E / finite-Larmor-radius**, and **polarization** drift.

---

## Overview

This project integrates the normalized Lorentz equations

```math
\frac{d\mathbf r}{dt}=\mathbf v,
\qquad
\frac{d\mathbf v}{dt}=\frac{q}{m}\left(\mathbf E+\mathbf v\times\mathbf B\right)
```

with a classical fourth-order Runge-Kutta method (RK4). The C++ code computes the full particle orbit directly; no drift velocity is inserted into the equations of motion.

Python scripts then read the generated CSV trajectories and perform independent analytical checks, guiding-center diagnostics, conservation tests, curve fitting, and visualization.

The repository contains the five required drift cases plus a uniform-magnetic-field benchmark used to validate the integrator.

---

## Project Layout

```text
charged_particle_sim/
├─ CMakeLists.txt
├─ LICENSE
├─ README.md
│
├─ include/
│  ├─ Vec3.hpp             # 3D vector operations
│  ├─ Particle.hpp         # Particle state and kinetic energy
│  ├─ Fields.hpp           # FieldCase and field parameters
│  └─ Integrator.hpp       # Lorentz acceleration and RK4 interface
│
├─ src/
│  ├─ main.cpp             # Active benchmark, time loop, CSV output
│  ├─ Fields.cpp           # Electromagnetic field models
│  └─ Integrator.cpp       # Full-orbit Lorentz integration with RK4
│
├─ analysis/
│  ├─ plot_uniform_B.py
│  ├─ plot_exb.py
│  ├─ plot_gradB.py
│  ├─ plot_curvature.py
│  ├─ plot_nonuniform_E.py
│  └─ plot_polarization.py
│
└─ results/
   ├─ *.csv                         # Numerical trajectories
   ├─ Charged_Particle_Drift_Simulations.tex
   ├─ Charged_Particle_Drift_Simulations_CN.tex
   └─ figures/                      # Analysis and validation plots
```

---

## Requirements

### C++ simulation

- C++17 compiler (`g++`, Clang, or MSVC)
- CMake >= 3.20
- Ninja is optional

### Python analysis

- Python 3
- NumPy
- Matplotlib

Install the Python dependencies with:

```bash
python -m pip install numpy matplotlib
```

---

## Build and Run

Run all commands from the repository root.

### Standard CMake build

```bash
cmake -S . -B build
cmake --build build
```

On Windows with the Ninja + GCC setup used during development:

```powershell
cmake -S . -B build -G Ninja -DCMAKE_CXX_COMPILER=g++
cmake --build build
.\build\charged_particle_sim.exe
```

On Linux/macOS:

```bash
./build/charged_particle_sim
```

The executable writes a CSV file to `results/` with columns

```text
t,x,y,z,vx,vy,vz,kinetic_energy
```

### Selecting a benchmark

The field models for all cases are already implemented in `src/Fields.cpp`. The active case is selected in the **Benchmark configuration** block near the top of `src/main.cpp`.

For each run, set:

1. the particle initial position and velocity,
2. `fields.type`,
3. the relevant field parameters,
4. `final_time`, and
5. the output CSV path.

The parameter sets used in the report are listed below.

---

## Benchmark Parameters

All simulations use normalized units with

```text
m = 1
q = +1
Δt = 0.01
```

The simulated particle is therefore a **generic positively charged reference particle**, not a specific ion species in SI units.

| Case | `FieldCase` | `r(0)` | `v(0)` | Field parameters | `t_final` | Output | Analysis |
|---|---|---:|---:|---|---:|---|---|
| Uniform B | `UniformB` | `(0,0,0)` | `(1,0,0.2)` | `B0 = 1` | 20 | `uniform_B.csv` | `plot_uniform_B.py` |
| E×B | `ExB` | `(0,0,0)` | `(1,0,0.2)` | `B0 = 1`, `E0 = 0.1` | 20 | `exb.csv` | `plot_exb.py` |
| Grad-B | `GradB` | `(0,0,0)` | `(1,0,0.2)` | `B0 = 1`, `gradient = 0.05` | 100 | `gradB.csv` | `plot_gradB.py` |
| Curvature | `Curvature` | `(10,0,0)` | `(0.2,1,0)` | `B0 = 1` | 100 | `curvature.csv` | `plot_curvature.py` |
| Nonuniform E | `NonuniformE` | `(0,0,0)` | `(1,0,0.2)` | `B0 = 1`, `E0 = 0.02`, `wavenumber = 0.4` | 60 | `nonuniform_E.csv` | `plot_nonuniform_E.py` |
| Polarization | `Polarization` | `(0,0,0)` | `(0,-0.2,0.2)` | `B0 = 1`, `E0 = 0.2`, `omega = 0.05` | 250 | `polarization.csv` | `plot_polarization.py` |

The repository currently keeps one benchmark active at a time so that the simulation setup remains explicit and easy to inspect.

---

## Run the Analysis

After generating the corresponding CSV file, run its Python analysis script from the repository root. For example:

```bash
python analysis/plot_nonuniform_E.py
```

or

```bash
python analysis/plot_polarization.py
```

Each script:

- loads the C++ trajectory,
- evaluates the relevant analytical or guiding-center prediction,
- computes quantitative errors or conservation diagnostics,
- saves figures under `results/figures/`, and
- displays the plots interactively.

The analysis is intentionally separate from the C++ integrator so that theoretical drift expressions are used only for **validation**, not to generate the simulated motion.

---

## Implemented Physics

- **Uniform B benchmark** — validates cyclotron motion, Larmor radius, helical motion, and kinetic-energy conservation.
- **E×B drift** — reproduces the charge- and mass-independent drift `E × B / B²`.
- **Grad-B drift** — uses a divergence-free weakly nonuniform magnetic field and compares the full orbit with first-order guiding-center theory.
- **Curvature drift** — uses circular constant-magnitude magnetic-field lines to isolate curvature from grad-B effects.
- **Nonuniform electric field** — demonstrates the finite-Larmor-radius correction to the local E×B drift and compares it with the direct gyro-average.
- **Polarization drift** — uses a slowly time-varying perpendicular electric field and compares the full orbit with both low-frequency drift theory and the exact forced solution.

For derivations, normalization, numerical results, figures, and quantitative comparisons, see the report in `results/`.

---

## Numerical Notes

- RK4 is used for every benchmark with `Δt = 0.01`.
- Uniform-field and exactly solvable cases are compared directly with analytical full-orbit solutions.
- Grad-B and curvature results are compared with **first-order guiding-center theory**, so their percent-level differences are physical approximation errors rather than RK4 errors.
- The nonuniform-E case explicitly tests finite-Larmor-radius physics.
- Energy or work conservation is checked using the appropriate invariant for each field configuration.

---

## Author

**Juan Pablo Solís Ruiz**  
M.Sc. Physics — Plasma Physics  
University of Science and Technology of China (USTC)

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE).
