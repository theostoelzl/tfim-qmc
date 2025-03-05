# Quantum Monte Carlo simulations of the transverse field Ising model using the stochastic series expansion method

## What this programme is

This programme is am implementation of a quantum Monte Carlo algorithm (QMC) using stochastic series expansion (SSE) developed by Anders Sandvik.
It is intended to simulate transverse-field Ising models (TFIMs) on arbitrary lattices and with arbitrary interactions. The only input this
programme requires are a file containing all bonds in the model and their respective coupling constants (*bonds.txt*) as well as a
configuration file (*setup.txt*) containing details such as the number of sweeps, etc.

## How to run this programme

1. To generate an input file containing all bonds (*bonds.txt*), run
	```python genlattice.py chain {n}```
	to generate a chain of *n* spins with periodic boundary conditions (PBC). To generate
	a simple square lattice with PBC, run
	```python genlattice.py square {nx} {ny}```.
	An output file (*bonds_{lattice}.txt*) will be generated, which you will need to rename to *bonds.txt*
	and put into the same directory as the main programme.
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

## Equilibration issues

## Licenses

(C) 2025 Theo Stölzl

This software is published under the GNU GPLv3 license
(https://www.gnu.org/licenses/gpl-3.0.en.html).

This software contains a JSON interface class by
Niels Lohmann under the MIT license (https://opensource.org/licenses/MIT).
(C) 2013-2025 Niels Lohmann