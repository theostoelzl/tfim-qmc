# Quantum Monte Carlo simulations of the transverse field Ising model

1. Run *genlattice.py* to generate file containing bonds (*bonds.txt*)
2. Modify *setup.json* as needed. The input parameters are:
	* *eqsweeps*: number of equilibration Monte Carlo (MC) sweeps per spin
	* *bins*: number of bins to average into
	* *avsweeps*: number of averaging MC sweeps per bin per spin
	* *max_expansion_order*: maximum series expansion order cutoff
	* *temperature* (obsolete): *kT*, temperature multiplied by the Boltzmann constant
	* *transverse_field* (obsolete): transverse field *h* in units of the coupling constant *J*
3. Compile *qmc.cpp*:
	```g++ qmc.cpp -g -o qmc```
4. Run executable with inline parameters:
	```./qmc {output directoy} {temperature} {transverse field}```
