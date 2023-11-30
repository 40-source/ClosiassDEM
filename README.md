# ClosiassDEM

Closiass is a nuclear simulation for the Monte Carlo resolution of the Boltzmann equation in neutronics.

## Compilation in C++:

For standard simulation:
```bash
g++ -Ofast -march=native V1.cpp -o DEM
```
For large simulations:
```bash
g++ -Ofast -fopenmp -march=native V2.cpp -o DEMV2
```
You can add -std=c++17 for better compatibility.

## Example of parameters:
```bash
./DEM 200.0 0.1 10.0 3.0 1.0 12.0 6.0 2.0
```
