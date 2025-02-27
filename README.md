# Quantum Monte Carlo simulations of the transverse field Ising model using the stochastic series expansion method

## How to run this programme

1. To generate an input file containing all bonds (*bonds.txt*), run
	```python genlattice.py chain {n}```
	to generate a chain of *n* spins with periodic boundary conditions (PBC). To generate
	a simple square lattice with PBC, run
	```python genlattice.py square {nx} {ny}```
2. Modify *setup.json* as needed. The input parameters are:
	* *eqsweeps*: number of equilibration Monte Carlo (MC) sweeps per spin
	* *bins*: number of bins to average into
	* *avsweeps*: number of averaging MC sweeps per bin per spin
	* *max_expansion_order*: initial maximum series expansion order cutoff (this value must be 1 or greater,
	but will be adapted automatically during equilibration)
3. Compile *qmc.cpp*:
	```g++ qmc.cpp -g -o qmc```
4. Run executable with inline parameters:
	```./qmc {output directoy} {temperature} {transverse field}```
