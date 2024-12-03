# Lithium-Air Battery Simulation

This repository contains a C++ program that simulates the electrochemical behavior of a lithium-air battery system. The simulation uses numerical methods to solve differential equations that describe the state variables of the battery over time.

## Overview

The simulation models:
- Electron conduction within the battery.
- Lithium ion conductivity in the electrolyte.
- Diffusion of lithium ions.
- Oxygen diffusion.
- Changes in porosity due to chemical reactions.

## Prerequisites

Before running this simulation, ensure you have:
- C++ Compiler with support for C++11 or later (e.g., GCC, Clang).

## Installation

- Clone the Repository.
- Use a C++ compiler to compile `main.cpp`. For example:
```sh
g++ -o simulate_battery main.cpp
```

## Usage

To run the simulation:
```
./simulate_battery
```

## Output Files

The simulation will generate several output files:
- `output_voltage.txt`: Voltage evolution over time..
- `output_phi.txt`: Electric potential distribution.
- `output_phi_li.txt`: Lithium potential distribution.
- `output_c_li.txt`: Lithium concentration distribution.
- `output_c_o2.txt`: Oxygen concentration distribution.
- `output_eps.txt`: Porosity distribution.

## Code Structure

- **Constants and Parameters**: Defined at the beginning of `main.cpp` for easy adjustment of simulation parameters like battery dimensions, physical constants, etc.
- **Function Definitions**:
  - `Solve`: Gauss elimination for solving linear systems.
  - `RC`, `dRC_*`: Reaction rate for oxygen reduction and its derivatives.
  - `RA`, `dRA_dphi_li`: Reaction rate for lithium oxidation and its derivative.
  - `ComputeRHS`: Computes the right-hand side of the differential equations.
  - `ComputeJ`: Constructs the Jacobian matrix for the system.
  - `Norm`: Calculates the Euclidean norm of a vector.
  - `porosity`: Checks if porosity remains positive to prevent simulation from unphysical states.
- **Main Simulation Loop**:
  - Initializes state variables.
  - Iterates over time steps, solving the differential equations until porosity becomes non-positive, indicating battery degradation.

## Key Equations

- **Electron Conduction**:
```
sigma * d^2phi/dx^2 = -F * RC(phi, phi_li, c_o2, eps)
```

- **Lithium Conductivity**:
```
kappa * d^2phi_li/dx^2 - kappa_d * d^2ln(c_li)/dx^2 = F * (RC - RA)
```

- **Lithium and Oxygen Diffusion**:
  Uses Fick's second law with reaction terms for lithium and oxygen.
  
- **Porosity Change**:
```
d(eps)/dt = - (M_li2o2 / 2 * rho_li2o2) * RC
```
