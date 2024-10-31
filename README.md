# Quantum Monte Carlo simulations of the transverse field Ising model

1. Run *genlattice.py* to generate file containing bonds (*bonds.txt*)
2. Modify *setup.json* as needed. The input parameters are:
	* *eqsweeps*: number of equilibration Monte Carlo (MC) sweeps
	* *avsweeps*: number of averaging MC sweeps
	* *max_expansion_order*: maximum series expansion order cutoff
	* *temperature*: temperature multiplied by the Boltzmann constant
	* *transverse_field*: transverse field in units of the coupling constant
3. Compile *qmc.cpp*:
	```g++ qmc.cpp -g -o qmc```
