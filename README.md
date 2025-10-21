# Quantum Monte Carlo for transverse-field Ising models

## What this programme is

This programme is am implementation of a quantum Monte Carlo algorithm (QMC) using stochastic series expansion (SSE) developed by Anders Sandvik.
It is intended to simulate transverse-field Ising models (TFIMs) on arbitrary lattices and with arbitrary interactions. The only input this
programme requires are a file containing all bonds in the model and their respective coupling constants (*bonds.txt*) as well as a
configuration file (*setup.txt*) containing details such as the number of sweeps, etc.

## How to run this programme

1. To generate an input file containing all bonds (*bonds.txt*), run
	```python genlattice.py {lattice} {n1} {n2} ...```
	where ```n1```, ```n2```, and so on are the lattice lengths along each spatial direction.
	Available lattices are (all with periodic boundary conditions): ```chain``` (1-D chain), ```square``` (simple square), 
	```inverse_square``` (see Sandvik 2018), ```square_crossings``` (J1-J2 square), ```shastry``` (Shastry-Sutherland).
	As of the latest version, coupling constants need to be modified in the code directly.
	An output file (*bonds_{lattice}.txt*) will be generated, which you will need to rename to *bonds.txt*
	and put into the same directory as the main programme.
2. Modify *setup.json* as needed. The input parameters are:
	* *eqsweeps*: number of equilibration Monte Carlo (MC) sweeps per spin
	* *bins*: number of bins to average into
	* *avsweeps*: number of averaging MC sweeps per bin per spin
	* *operator_mode*: type of spin operator; either *pauli* (Pauli matrices) or *one-half* (spin-1/2 operators)
3. Compile *qmc.cpp*:
	```g++ -O3 -march=native qmc.cpp -o qmc``` (or using any other C++ compiler of your choice)
4. Run executable with inline parameters:
	```./qmc {output directoy} {temperature} {transverse field} {longitudinal field}```
5. The programme will generate two files in its directory, *spins.txt* and *opstring.txt*, which will be used automatically as input if a subsequent simulation is ran in the same directory.

## Licenses

(C) 2025 Theo Stölzl

This software is published under the GNU GPLv3 license
(https://www.gnu.org/licenses/gpl-3.0.en.html).

This software contains a JSON interface class by
Niels Lohmann under the MIT license (https://opensource.org/licenses/MIT).
(C) 2013-2025 Niels Lohmann
