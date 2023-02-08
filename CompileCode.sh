#!/bin/bash
# Compiles code for production runs (with optimisations)

# Run code for normal simulations with optimisation options enabled
ifort -heap-arrays -O2 -ipo -fp-model precise -real-size 64 -mkl GlobalVariables.f90 MatrixOperations.f90 Geometric.f90  PrimConTransformation.f90 Fluxes.f90 BCandIC.f90 CNSSubroutines.f90 Load.f90 SolverCNS.f90

# debug
#ifort -heap-arrays -O0 -g -traceback -ipo -fp-model precise -real-size 64 -mkl GlobalVariables.f90 MatrixOperations.f90 Geometric.f90  PrimConTransformation.f90 Fluxes.f90 BCandIC.f90 CNSSubroutines.f90 Load.f90 SolverCNS.f90



