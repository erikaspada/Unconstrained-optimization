# Unconstrained Optimization Project

This repository contains the code developed for the *Unconstrained Optimization* project, part of the course **Numerical Optimization for Large-Scale Problems**.

The goal of this project is to compare and test the performance of three unconstrained optimization algorithms on several objective functions without constraints.

---

The following optimization methods are implemented and tested:

- **Fletcher-Reeves Conjugate Gradient (FR-CG)** — `FR_CG_bcktrck.m`
- **Inexact Newton Method** — `innewton_general.m`
- **Nelder-Mead Method** — `nelder_mead.m`

---

The algorithms have been tested on the following well-known nonlinear functions:

- `rosenbrock.m` — Rosenbrock function  
- `chainedpowell.m` — Chained Powell function  
- `banded_trigonometric.m` — Banded trigonometric function  
- `generalized_brown.m` — Generalized Brown function

Each function is accompanied by its gradient implementation (e.g., `rosenbrock_grad.m`, `generalized_brown_grad.m`, etc.).

---

## Requirements
- MATLAB
- File `.mat` with datas: `generalized_brown.mat`, `forcing_terms.mat`

---
## Authors
The authors of the project are:
- Erika Spada s318375@studenti.polito.it
- Lisa Gamarro s318426@studenti.polito.it 
- Carlotta Agnese Trovati s318396@studenti.polito.it 




