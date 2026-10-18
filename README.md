# Charged Particle Drift Simulator

C++17 full-orbit simulation of a charged particle in prescribed electromagnetic fields, with Python analysis and visualization.

Implemented cases:

- Uniform magnetic field benchmark
- E x B drift
- Grad-B drift
- Curvature drift
- Spatially nonuniform electric field drift
- Polarization drift

## Build

```powershell
cmake -S . -B build -G Ninja -DCMAKE_CXX_COMPILER=g++
cmake --build build
.\build\charged_particle_sim.exe
```

The active benchmark is configured in `src/main.cpp`. The parameter sets used for all reported cases are listed in the accompanying report.

## Analysis

The scripts in `analysis/` read CSV files from `results/` and generate validation diagnostics and figures.

Example:

```powershell
python .\analysis\plot_nonuniform_E.py
```

Python dependencies: NumPy and Matplotlib.

## Numerical model

The C++ code integrates

`dr/dt = v`

`dv/dt = (q/m) [E(r,t) + v x B(r,t)]`

using a fourth-order Runge-Kutta method. Drift velocities are not inserted explicitly; they emerge from the full Lorentz-force trajectory.
