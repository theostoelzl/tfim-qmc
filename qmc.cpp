/*

Quantum Monte Carlo simulation of transverse-field
Ising models by stochastic series expansion.

(C) 2025 Theo Stölzl

This software is published under the GNU GPLv3 license
(https://www.gnu.org/licenses/gpl-3.0.en.html).

This software contains a JSON interface class by
Niels Lohmann under the MIT license (https://opensource.org/licenses/MIT).
(C) 2013-2025 Niels Lohmann

*/


#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include "json.hpp"

using namespace std;

int random_conf(int spins[], int nspins, mt19937 &rng);
int diagonal_updates(int spins[], int nspins, int bonds[][2], double couplings[], int nbonds,
	int opstring[][3], int effexporder, double temp, double transfield, double longfield, 
	bool use_pauli_ops, mt19937 &rng);
int flip_spin_clusters(int spins[], int nspins, int bonds[][2],
	int nbonds, int opstring[][3], int effexporder, mt19937 &rng);
int shift_update(int opstring[][3], int effexporder, int spins[], int nspins, 
	int bonds[][2], double couplings[], int nbonds, mt19937 &rng);
int cluster_updates(int spins[], int nspins, int bonds[][2],
	int nbonds, int opstring[][3], int effexporder, mt19937 &rng);
int adjust_maxexporder(int opstring[][3], int effexporder, mt19937 &rng, int scale_up_to=0);
int measure_observables(int spins[], int nspins, int bonds[][2], int nbonds,
	int opstring[][3], int effexporder, double temp, double transfield, 
	int obs_count, double obs_exporder[], double obs_exporder_sq[], double obs_magn[], 
	double obs_magn_sq[], double obs_magn_quad[], double obs_trans_magn[], 
	double obs_trans_magn_sq[], double obs_corr[500][500]);

int main(int argc, char** argv) {

	// ----- Initialise random number generator -----
	int rngseed = 178;
	mt19937 rng;
	rng.seed(time(NULL)+100000*rngseed);

	// ----- Read in bonds from input file -----

	// Open input file containing bonds
	ifstream bonds_infile;
	bonds_infile.open("bonds.txt");

	// Read in number of spins and bonds
	string in_nspins, in_nbonds;
	getline(bonds_infile, in_nspins, '\t');
	getline(bonds_infile, in_nbonds, '\n');
	cout << "spins: " << in_nspins << ", bonds: " << in_nbonds << "\n";

	// First, initialise number of spins and bonds
	int nspins = stoi(in_nspins);
	int nbonds = stoi(in_nbonds);

	// Then, initialise arrays for all spins and bonds
	int spins[nspins] =  {0};
	int bonds[nbonds][2] = {0};
	double couplings[nbonds] = {0};

	// Read in bonds
	string myline[4];
	for (int l = 0; l < nbonds; l++) {
		for (int i = 0; i < 3; i++) {
			getline(bonds_infile, myline[i], '\t');
		}
		getline(bonds_infile, myline[3]);

		// Store spin indices and coupling constants in bonds array
		bonds[l][0] = stoi(myline[1]);
		bonds[l][1] = stoi(myline[2]);
		couplings[l] = stod(myline[3]);

		if (l < 10) {
			cout << bonds[l][0] << "\t" << bonds[l][1] << "\t" << couplings[l] << "\n";
		}
	}

	cout << "...\n";

	// ----- Read in simulation parameters -----

	ifstream f("setup.json");
	nlohmann::json data = nlohmann::json::parse(f);
	int eqsweeps = data["eqsweeps"];
	int bins = data["bins"];
	int avsweeps = data["avsweeps"];
	string op_mode = data["operator_mode"];
	bool use_pauli_ops = (op_mode == "pauli") ? true : false;
	double temp = 0;
	double transfield = 0;
	double longfield = 0;
	// Overwrite using in-line arguments
	string out_path = argv[1];
	temp = stod(argv[2]);
	cout << "temp " << temp << "\n";
	transfield = stod(argv[3]);
	cout << "transverse field " << transfield << "\n";
	longfield = stod(argv[4]);
	cout << "longitudinal field " << longfield << "\n";
	
	f.close();

	// Take input numbers of sweeps as sweeps per spin
	eqsweeps *= nspins;
	avsweeps *= nspins;

	// ----- Get initial state or generate new one -----

	// Try opening input file containing spins
	ifstream spins_infile;
	spins_infile.open("spins.txt");
	if (spins_infile.good()) {
		// Read in spins
		string spin_in;
		for (int l = 0; l < nspins; l++) {
			getline(spins_infile, spin_in, '\n');

			// Set spin
			spins[l] = stoi(spin_in);
		}
	} else {
		random_conf(spins, stoi(in_nspins), rng);
	}
	spins_infile.close();

	double avgm = 0;
	for (int i = 0; i < nspins; i++) {
		//cout << spins[i] << "\t";
		avgm += spins[i];
	}
	avgm = avgm / nspins;
	cout << "average initial magnetisation: " << avgm << "\n---------\n";

	// ----- Start Monte Carlo sweeps -----

	// Set initial truncated (effective) expansion order
	int effexporder = 10;

	// Initialise operator string
	int opstring[100000][3] = {};

	for (int p = 0; p < 100000; p++) {
		// Type of operator at position p
		// 0 - identity op
		// 1 - diagonal op
		// 2 - off-diagonal op
		// 3 - longitudinal field op
		opstring[p][0] = 0;
		// If operator is acting bonds, this contains bond index
		// otherwise it's -1
		opstring[p][1] = -1;
		// If operator is acting on spins, this contains spin index
		// otherwise it's -1
		opstring[p][2] = -1;
	}

	// Try opening input file containing opstring
	ifstream ops_infile;
	ops_infile.open("opstring.txt");
	if (ops_infile.good()) {
		// Read in maximum expansion order
		string in_effexporder;
		getline(ops_infile, in_effexporder, '\n');
		effexporder = stoi(in_effexporder);

		// Read in opstring
		string myline[3];
		for (int l = 0; l < effexporder; l++) {
			for (int i = 0; i < 2; i++) {
				getline(ops_infile, myline[i], '\t');
			}
			getline(ops_infile, myline[2], '\n');

			// Set operator
			opstring[l][0] = stoi(myline[0]);
			opstring[l][1] = stoi(myline[1]);
			opstring[l][2] = stoi(myline[2]);
		}
	}
	ops_infile.close();

	cout << "Started equilibration ..." << "\n";
	cout << flush;

	// Mean maximum expansion order
	double avg_max_exporder = 0;
	int avg_from = eqsweeps - 1000*nspins;

	// Start equilibration
	for (int i = 0; i < eqsweeps; i++) {

		// Diagonal updates to insert / remove operators
		diagonal_updates(spins, nspins, bonds, couplings, nbonds, opstring, 
				effexporder, temp, transfield, longfield, use_pauli_ops, rng);
		
		if (longfield == 0) {
			flip_spin_clusters(spins, nspins, bonds, nbonds, opstring, effexporder, rng);
		}

		// Cluster updates to vary diagonal / off-diagonal ops
		cluster_updates(spins, nspins, bonds, nbonds, opstring, effexporder, rng);
		
		// Shift updates to move around bond and transverse-field operators
		// at little to no costf
		shift_update(opstring, effexporder, spins, nspins, bonds, couplings, nbonds, rng);

		// Adjust largest expansion order
		effexporder = adjust_maxexporder(opstring, effexporder, rng);

		if (i > avg_from) {
			avg_max_exporder += effexporder/(double) (1000*nspins);
		}

	}

	// Set new maximum expansion order and scale opstring
	int new_maxexporder = ceil(avg_max_exporder);
	effexporder = adjust_maxexporder(opstring, effexporder, rng, new_maxexporder);

	cout << "Finished equilibration." << "\n";
	cout << "Maximum expansion order: " << effexporder << "\n";
	cout << "Started averaging ..." << "\n";
	cout << flush;

	// Open output file
	ofstream out_file; 
	out_file.open(out_path+"/results_t_"+to_string(temp)+"_g_"+to_string(transfield)+"_h_"+to_string(longfield)+".csv");
	out_file << "bin,exporder,exporder_sq,magn,magn_sq,magn_quad,trans_magn,trans_magn_sq\n";
	
	// Setup for measuring observables
	double obs_exporder[bins] = {0};
	double obs_exporder_sq[bins] = {0};
	double obs_magn[bins] = {0};
	double obs_magn_sq[bins] = {0};
	double obs_magn_quad[bins] = {0};
	double obs_trans_magn[bins] = {0};
	double obs_trans_magn_sq[bins] = {0};
	double obs_corr[500][500] = {0};
	
	// Start averaging sweeps
	for (int n = 0; n < bins; n++) {
		
		double obs_exporder_bin[avsweeps] = {0};
		double obs_exporder_sq_bin[avsweeps] = {0};
		double obs_magn_bin[avsweeps] = {0};
		double obs_magn_sq_bin[avsweeps] = {0};
		double obs_magn_quad_bin[avsweeps] = {0};
		double obs_trans_magn_bin[avsweeps] = {0};
		double obs_trans_magn_sq_bin[avsweeps] = {0};

		for (int j = 0; j < avsweeps; j++) {

			// Diagonal updates to insert / remove operators
			diagonal_updates(spins, nspins, bonds, couplings, nbonds, opstring, 
					effexporder, temp, transfield, longfield, use_pauli_ops, rng);
			
			if (longfield == 0) {
				flip_spin_clusters(spins, nspins, bonds, nbonds, opstring, effexporder, rng);
			}

			// Cluster updates to vary diagonal / off-diagonal ops
			cluster_updates(spins, nspins, bonds, nbonds, opstring, effexporder, rng);

			// Shift updates to move around bond and transverse-field operators
			// at little to no cost
			shift_update(opstring, effexporder, spins, nspins, bonds, couplings, nbonds, rng);
			
			// Measure observables
			measure_observables(spins, nspins, bonds, nbonds, opstring,
					effexporder, temp, transfield, j, obs_exporder_bin, 
					obs_exporder_sq_bin, obs_magn_bin, obs_magn_sq_bin, 
					obs_magn_quad_bin, obs_trans_magn_bin, obs_trans_magn_sq_bin,
					obs_corr);

		}
		
		// Now average over bin and add to big arrays
		for (int j = 0; j < avsweeps; j++) {
			obs_exporder[n] += obs_exporder_bin[j]/(double)avsweeps;
			obs_exporder_sq[n] += obs_exporder_sq_bin[j]/(double)avsweeps;
			obs_magn[n] += obs_magn_bin[j]/(double)avsweeps;
			obs_magn_sq[n] += obs_magn_sq_bin[j]/(double)avsweeps;
			obs_magn_quad[n] += obs_magn_quad_bin[j]/(double)avsweeps;
			obs_trans_magn[n] += obs_trans_magn_bin[j]/(double)avsweeps;
			obs_trans_magn_sq[n] += obs_trans_magn_sq_bin[j]/(double)avsweeps;
		}

		// Write observables to output file
		out_file << to_string(n)+","+to_string(obs_exporder[n])+","+to_string(obs_exporder_sq[n])+","
				+to_string(obs_magn[n])+","+to_string(obs_magn_sq[n])+","+to_string(obs_magn_quad[n]) +
				+","+to_string(obs_trans_magn[n])+","+to_string(obs_trans_magn_sq[n])+"\n" << flush;
	}

	out_file.close();

	// Write spin-spin correlations to dedicated file
	ofstream corrs_file;
	corrs_file.open(out_path+"/corrs_t_"+to_string(temp)+"_g_"+to_string(transfield)+"_h_"+to_string(longfield)+".dat");
	corrs_file << nspins << "\n";
	for (int si = 0; si < nspins; si++) {
		for (int sj = 0; sj < nspins; sj++) {
			corrs_file << obs_corr[si][sj]/(double)(avsweeps*bins) << "\t";
		}
		corrs_file << "\n";
	}

	cout << "Finished averaging." << "\n";
	cout << flush;
	
	// Find average observables
	double avgn = 0;
	avgm = 0;
	for (int n = 0; n < bins; n++) {
		avgn += obs_exporder[n]/(double)(bins);
		//cout << obs_magn[n] << "\n";

		avgm += obs_magn[n]/(double)(bins);
	}
	
	cout << "average expansion order: " << avgn << "\n";
	cout << "average magnetisation: " << avgm << "\n";
	
	// Write spins to output file
	ofstream spins_out;
	spins_out.open("spins.txt");
	for (int s = 0; s < nspins; s++) {
		spins_out << spins[s] << "\n";
	}
	spins_out.close();

	// Write opstring to output file
	ofstream ops_out;
	ops_out.open("opstring.txt");
	ops_out << effexporder << "\n";
	for (int p = 0; p < effexporder; p++) {
		ops_out << opstring[p][0] << "\t"
			<< opstring[p][1] << "\t"
			<< opstring[p][2] << "\n";
	}
	ops_out.close();

	return 0;
	
}

int random_conf(int spins[], int nspins, mt19937 &rng) {
	
	// Random stuff
	uniform_real_distribution<double> uni_dist(0,1);

	// Iterate through lattice and assign random spin directions 
	for (int i = 0; i < nspins; i++) {
		spins[i] = (uni_dist(rng) < 0.5) ? -1 : 1;
	}

	return 0;
	
}

int diagonal_updates(int spins[], int nspins, int bonds[][2], double couplings[], int nbonds,
	int opstring[][3], int effexporder, double temp, double transfield, double longfield,
	bool use_pauli_ops, mt19937 &rng) {

	// Initialise random number generator
	uniform_real_distribution<double> uni_dist(0,1);
	uniform_int_distribution<int> rand_spin(0,nspins-1);
	uniform_int_distribution<int> rand_bond(0,nbonds-1);

	// Evaluate current expansion order
	int exporder = effexporder;
	for (int i = 0; i < effexporder; i++) {
		if (opstring[i][0] == 0) {
			exporder = exporder - 1;
		}
	}

	// Initialise variables
	int op_type = 0, randb = 0, rands = 0, s1i = 0, s1 = 0, s2i = 0, s2 = 0, si = 0;
	double coin = 0, paccept = 0, bj = 0;

	// Set up acceptance probabilities
	double paccept_insert_bond = 0, paccept_insert_trans = 0, paccept_insert_long = 0,
		paccept_remove_bond = 0, paccept_remove_trans = 0, paccept_remove_long = 0;
	
	// Iterate over all positions p in operator string
	for (int p = 0; p < effexporder; p++) {
		op_type = opstring[p][0];
		
		// Evaluate the type of operator present at p
		if (op_type == 0) {
			// Identity operator
			
			// Flip a coin and decide if to proceed with 
			// insertion of bond or spin operator
			coin = uni_dist(rng);
			if (coin < 1/(double)3) {
				// Try to insert a diagonal bond operator
				
				// Choose random bond to insert into
				randb = rand_bond(rng);
				
				// Get corresponding spins
				s1i = bonds[randb][0];
				s1 = spins[s1i];
				s2i = bonds[randb][1];
				s2 = spins[s2i];

				// Get coupling constant of bond
				bj = couplings[randb];
				
				// Check if spins are parallel
				// if bj < 0, anti-ferromagnetic bond
				// if bj > 0, ferromagnetic bond
				if (bj*s1*s2 > 0) {
					// Spins are arranged correctly so try to insert diagonal operator

					// Evaluate acceptance probability
					if (use_pauli_ops) {
						paccept = (1/temp) * 3 * nbonds * 2 * abs(bj) /(double) (effexporder - exporder); // Pauli matrices
					} else {
						paccept = (1/temp) * 3 * nbonds * abs(bj) /(double) (2*(effexporder - exporder)); // spin 1/2
					}

					// Attempt move
					coin = uni_dist(rng);
					if (coin < paccept) {
						// Insert diagonal bond operator
						opstring[p][0] = 1;
						opstring[p][1] = randb;
						opstring[p][2] = -1;
						exporder++;
					}
				}
			} else if (coin < 2/(double)3) {
				// Try to insert a diagonal spin transverse operator
				
				// Get random position
				rands = rand_spin(rng);
				
				// Evaluate acceptance probability				
				if (use_pauli_ops) {
					paccept = (1/temp) * 3 * nspins * transfield /(double) (effexporder - exporder); // Pauli matrices
				} else {
					paccept = (1/temp) * 3 * nspins * transfield /(double) (2*(effexporder - exporder)); // spin 1/2
				}
				
				// Attempt move
				if (uni_dist(rng) < paccept) {
					// Insert diagonal spin operator
					opstring[p][0] = 1;
					opstring[p][1] = -1;
					opstring[p][2] = rands;
					exporder++;
				}
			} else {
				// Try to insert a diagonal spin longitudinal operator
				
				// Get random position
				rands = rand_spin(rng);
				
				// Check if spin is aligned with field
				if (spins[rands] == 1) {
					// Evaluate acceptance probability
					if (use_pauli_ops) {
						// Pauli matrices
						paccept = (1/temp) * 3 * nspins * 2 * longfield /(double) (effexporder - exporder);
					} else {
						// spin 1/2
						paccept = (1/temp) * 3 * nspins * longfield /(double) (effexporder - exporder);
					}
					
					// Attempt move
					if (uni_dist(rng) < paccept) {
						// Insert diagonal spin operator
						opstring[p][0] = 3;
						opstring[p][1] = -1;
						opstring[p][2] = rands;
						exporder++;
					}
				}
			}
		} else if (op_type == 1) {
			// Diagonal operator
			
			// Check for error
			if (opstring[p][1] > -1 && opstring[p][2] > -1) {
				cout << "Something's wrong!" << "\n";
				cout << flush;
			}
			
			// Check if bond or spin operator
			if (opstring[p][1] > -1) {
				// Operator is acting on bonds
				
				// Get coupling constant of bond
				bj = couplings[opstring[p][1]];

				// Evaluate acceptance probability
				if (use_pauli_ops) {
					// Pauli matrices
					paccept = (effexporder - exporder + 1) /(double) ((1/temp) * 3 * nbonds * 2 * abs(bj));
				} else {
					// spin 1/2
					paccept = 2 * (effexporder - exporder + 1) /(double) ((1/temp) * 3 * nbonds * abs(bj));
				}
				
				// Attempt move
				if (uni_dist(rng) < paccept) {
					// Remove diagonal bond operator
					opstring[p][0] = 0;
					opstring[p][1] = -1;
					opstring[p][2] = -1;
					exporder--;
				}
			} else if (opstring[p][2] > -1) {
				// Operator is transverse field op
				
				// Evaluate acceptance probability
				if (use_pauli_ops) {
					// Pauli matrices
					paccept = (effexporder - exporder + 1) /(double) ((1/temp) * 3 * nspins * transfield);
				} else {
					// spin 1/2
					paccept = 2 * (effexporder - exporder + 1) /(double) ((1/temp) * 3 * nspins * transfield);
				}
				
				// Attempt move
				if (uni_dist(rng) < paccept) {
					// Remove diagonal spin operator
					opstring[p][0] = 0;
					opstring[p][1] = -1;
					opstring[p][2] = -1;
					exporder--;
				}
			}
		} else if (op_type == 3) {
			// Longitudinal field operator
			
			// Evaluate acceptance probability
			if (use_pauli_ops) {
				// Pauli matrices
				paccept = (effexporder - exporder + 1) /(double) ((1/temp) * 3 * nspins * 2 * longfield);
			} else {
				// spin 1/2
				paccept = (effexporder - exporder + 1) /(double) ((1/temp) * 3 * nspins * longfield);
			}
			
			// Attempt move
			if (uni_dist(rng) < paccept) {
				// Remove longitudinal field operator
				opstring[p][0] = 0;
				opstring[p][1] = -1;
				opstring[p][2] = -1;
				exporder--;
			}
		} else {
			// Off-diagonal operator

			// Modify spin according to opstring
			si = opstring[p][2];
			spins[si] *= -1;
		}
	}

	return 0;
	
}

int flip_spin_clusters(int spins[], int nspins, int bonds[][2],
	int nbonds, int opstring[][3], int effexporder, mt19937 &rng) {
	
	// Random stuff
	uniform_real_distribution<double> uni_dist(0,1);

	// ----- Construct vertex link list -----

	// List of first and last vertex legs of each spin
	int firsts[nspins] = {0};
	int lasts[nspins] = {0};
	for (int i = 0; i < nspins; i++) {
		firsts[i] = -1;
		lasts[i] = -1;
	}
	// List of vertex links
	int links[4*effexporder] = {0};
	for (int i = 0; i < 4*effexporder; i++) {
		links[i] = -1;
	}
	// List of operator types each vertex leg belongs to
	int linkops[4*effexporder] = {0};
	for (int i = 0; i < 4*effexporder; i++) {
		// This encodes the type of operator each
		// vertex leg belongs too
		// 0 - identity op
		// 1 - bond op
		// 2 - transverse field op
		// 3 - longitudinal field op
		linkops[i] = 0;
	}

	// Initialise variables
	int v0 = 0, v1 = 0, v2 = 0, s1i = 0, s2i = 0, bond = 0, sp = 0, si = 0;

	// First, deal with links within bounds of opstring
	// To do this, iterate over opstring
	for (int p = 0; p < effexporder; p++) {
		// Check if operator at p is acting on bonds or spins
		if (opstring[p][0] > 0 && opstring[p][1] > -1) {
			// Operator is acting on bonds

			// Set vertex index
			v0 = 4*p;
			// Get bond index the op is acting on
			bond = opstring[p][1];
			// Get indices of spins belonging to bond
			s1i = bonds[bond][0];
			s2i = bonds[bond][1];
			// Store op type
			linkops[v0] = 1;
			linkops[v0+1] = 1;
			linkops[v0+2] = 1;
			linkops[v0+3] = 1;
			// Get last vertex indices of each spin
			v1 = lasts[s1i];
			v2 = lasts[s2i];
			// Check if spins have appeared in vertex before
			if (v1 > -1) {
				// First spin has appeared before, so set link
				links[v1] = v0;
				links[v0] = v1;
			} else {
				// Set first appearance of spin as current leg
				firsts[s1i] = v0;
			}
			if (v2 > -1) {
				// Second spin has appeared before, so set link
				links[v2] = v0+1; // ?? are we sure this is supposed to be v0 and not v0+1?
				links[v0+1] = v2;
			} else {
				// Set first appearance of spin as current p
				firsts[s2i] = v0+1;
			}
			// Update latest vertex appearance
			lasts[s1i] = v0+2;
			lasts[s2i] = v0+3;
		}
	}

	// ----- Then trace all vertex loops and do flip updates -----

	// Initialise variables
	bool finished = false, fliploop = true, allspinops = false,
		cont = false, bondops_loop[effexporder] = {false},
        spins_loop[nspins] = {false}, visited[4*effexporder] = {false};
	for (int i = 0; i < 4*effexporder; i++) {
		visited[i] = false;
	}
	for (int i = 0; i < effexporder; i++) {
		bondops_loop[i] = false;
	}
	for (int i = 0; i < nspins; i++) {
		spins_loop[i] = false;
	}
	int currentv = 0, linkedv = 0, lastv = 0, i = 0, p = 0, 
		linkedp = 0, change_counter = 0, bi = 0;
	double coin = 0;

	// Iterate over all vertex legs
	for (int v = 0; v < 4*effexporder; v++) {
		// Check if vertex leg is connected to any other legs
		if (links[v] > -1 && !visited[v]) {
            // Set up variables for loop
			finished = false;
            fliploop = true;

            // Set current leg as initial leg
			currentv = v;

            // See if initial leg belongs to bond or spin op
            if (linkops[currentv] == 1) {
                // Initial leg belongs to bond op

				visited[currentv] = true;

                // Get opstring index
                p = floor(currentv/(double)4);
                
                // Operator is part of loop
                bondops_loop[p] = true;

				// Get spins belonging to bond
				bi = opstring[p][1];
				s1i = bonds[bi][0];
				s2i = bonds[bi][1];

				// Set spins as part of cluster
				spins_loop[s1i] = true;
				spins_loop[s2i] = true;
            } else if (linkops[currentv] == 2) {
                cout << "nope!";
            }
            
			// Traverse loop
			while (!finished) {
                // Reset counter of how many ops have been added to loop
                change_counter = 0;

                // Iterate over all operators
                for (int pi = 0; pi < effexporder; pi++) {
                    if (bondops_loop[pi]) {
                        // Iterate over legs of op
                        // Go through legs belonging to bond op
                        for (int i = 0; i < 4; i++) {
                            // Get current vertex leg
                            currentv = 4*pi+i;
                            visited[currentv] = true;
                            
                            linkedv = links[currentv];
							
							// Check if leg is linked to anything
							if (linkedv > -1) {
								// Leg is linked to something, so get linked op
								linkedp = floor(linkedv/(double)4);
								// Check if this linked op is a bond op
								if (linkops[linkedv] == 1) {
									// Check if bond op is already part of loop
									if (!bondops_loop[linkedp]) {
										// Bond op hasn't been added yet, so add it
										bondops_loop[linkedp] = true;
										change_counter++;

										// Get spins belonging to bond
										bi = opstring[linkedp][1];
										s1i = bonds[bi][0];
										s2i = bonds[bi][1];

										// Set spins as part of cluster
										spins_loop[s1i] = true;
										spins_loop[s2i] = true;
									}
								}
							}
                        }
                    }
                }

				// If loop hasn't grown in this iteration, finish the loop
                if (change_counter == 0) {
                    finished = true;
                }
	        }

            // Check if loop should be flipped
            if (uni_dist(rng) < 0.5) {
                // Loop is flippable

                // Iterate over all spin ops and change type
                for (int si = 0; si < nspins; si++) {
                    // Check if spin op is part of loop
                    if (spins_loop[si]) {
                        // Flip spin
						spins[si] *= -1;
                    }
                }
            }
        }

		// Empty list of bond and spin ops in loop
		for (int i = 0; i < effexporder; i++) {
			bondops_loop[i] = false;
		}
		for (int si = 0; si < nspins; si++) {
			spins_loop[si] = false;
        }
	}

	// Finally, flip all spins without associated bond ops
	for (int si = 0; si < nspins; si++) {
		// Check if spin op has no associated bond ops
		if (firsts[si] == -1 && uni_dist(rng) < 0.5) {
			// Flip spin with 1/2 prob
			spins[si] *= -1;
		}
	}

	return 0;

}

int shift_update(int opstring[][3], int effexporder, int spins[], int nspins, 
	int bonds[][2], double couplings[], int nbonds, mt19937 &rng) {

	// Random stuff
	uniform_real_distribution<double> uni_dist(0,1);
	uniform_int_distribution<int> rand_spin(0,nspins-1);

	// ----- First, swap exchange operators and identity operators -----
	
	// List of opstring positions of transverse and id ops
	int tf_op[effexporder] = {0};
	int id_op[effexporder] = {0};
	for (int i = 0; i < effexporder; i++) {
		tf_op[i] = -1;
		id_op[i] = -1;
	}
	// Iterate over opstring and count number of 
	// transverse and id ops in the opstring
	int n_tf = 0;
	int n_id = 0;
	for (int p = 0; p < effexporder; p++) {
		if (opstring[p][0] == 0) {
			// Identity operator
			id_op[n_id] = p;
			n_id++;
		} else if (opstring[p][0] == 1 && opstring[p][2] > -1) {
			// Diagonal transverse-field operator
			tf_op[n_tf] = p;
			n_tf++;
		}
	}

	int remove_tf = -1, insert_id = -1, insert_spin = -1, offset = 0,
		loop_until = 0, rand = -1;
	
	// More transverse ops than id ops?
	bool more_tf = (n_tf > n_id) ? true : false;
	if (more_tf) {
		loop_until = n_id;
	} else {
		loop_until = n_tf;
	}

	// Iterate over all ops of which we have more
	for (int i = 0; i < loop_until; i++) {
		// Choose if op will be moved
		if (uni_dist(rng) < 0.5) {
			// Proceed with moving op

			if (more_tf) {
				// Choose random transverse-field op
				uniform_int_distribution<int> rand_op(0,n_tf-1);
				rand = rand_op(rng);
				remove_tf = tf_op[rand];
				insert_id = id_op[i];
			} else {
				// Choose random identity op
				uniform_int_distribution<int> rand_op(0,n_id-1);
				rand = rand_op(rng);
				remove_tf = tf_op[i];
				insert_id = id_op[rand];
			}
			// DEBUG
			if (remove_tf == -1 || insert_id == -1) {
				cout << "HALT" << flush;
			}
			
			// Swap ops
			opstring[insert_id][0] = 1;
			insert_spin = rand_spin(rng);
			opstring[insert_id][2] = insert_spin;
			opstring[remove_tf][0] = 0;
			opstring[remove_tf][2] = -1;

			if (more_tf) {
				// Drop operator we just moved from list
				for (int j = 0; j < n_tf-1; j++) {
					// Is this identity operator the one we just replaced?
					if (tf_op[j] == remove_tf) {
						// Yes, it is
						offset = 1;
					}
					tf_op[j] = tf_op[j+offset];
				}
				tf_op[n_tf-1] = -1;
				n_tf--;
			} else {
				// Update number of identity operators and list of positions
				for (int j = 0; j < n_id-1; j++) {
					// Is this identity operator the one we just replaced?
					if (id_op[j] == insert_id) {
						// Yes, it is
						offset = 1;
					}
					id_op[j] = id_op[j+offset];
				}
				id_op[n_id-1] = -1;
				n_id--;
			}
		}
		offset = 0;
	}

	// ----- Now shift bond ops -----

	int randb = -1;
	int s1i = 0, s2i = 0, s1 = 0, s2 = 0;
	double bj_old = 0, bj_new = 0, paccept = 0;

	uniform_int_distribution<int> rand_bond(0,nbonds-1);

	for (int p = 0; p < effexporder; p++) {
		if (opstring[p][0] == 1 && opstring[p][1] > -1) {
			// Bond operator
			
			// Decide whether to attempt shift
			if (uni_dist(rng) < 0.5) {
				// Get current coupling constant
				bj_old = couplings[opstring[p][1]];

				// Choose random bond to shift to
				randb = rand_bond(rng);

				// Get corresponding spins
				s1i = bonds[randb][0];
				s1 = spins[s1i];
				s2i = bonds[randb][1];
				s2 = spins[s2i];

				// Get coupling constant of bond
				bj_new = couplings[randb];
				paccept = abs(bj_new/(double) bj_old);
				
				// Check if spins are parallel
				// if bj < 0, anti-ferromagnetic bond
				// if bj > 0, ferromagnetic bond
				if (bj_new*s1*s2 > 0 && uni_dist(rng) < paccept) {
					// Shift operator to other bond
					opstring[p][1] = randb;
				}
			}
		} else if (opstring[p][0] == 2) {
			spins[opstring[p][2]] *= -1;
		}
	}

	return 0;
}

int cluster_updates(int spins[], int nspins, int bonds[][2],
	int nbonds, int opstring[][3], int effexporder, mt19937 &rng) {
	
	// Random stuff
	uniform_real_distribution<double> uni_dist(0,1);

	// ----- Construct vertex link list -----

	// List of first and last vertex legs of each spin
	int firsts[nspins] = {0};
	int lasts[nspins] = {0};
	for (int i = 0; i < nspins; i++) {
		firsts[i] = -1;
		lasts[i] = -1;
	}
	// List of vertex links
	int links[4*effexporder] = {0};
	for (int i = 0; i < 4*effexporder; i++) {
		links[i] = -1;
	}
	// List of operator types each vertex leg belongs to
	int linkops[4*effexporder] = {0};
	for (int i = 0; i < 4*effexporder; i++) {
		// This encodes the type of operator each
		// vertex leg belongs too
		// 0 - identity op
		// 1 - bond op
		// 2 - transverse field op
		// 3 - longitudinal field op
		linkops[i] = 0;
	}
	bool cross_boundary_link[4*effexporder] = {false};

	// Initialise variables
	int v0 = 0, v1 = 0, v2 = 0, s1i = 0, s2i = 0, bond = 0, 
		sp = 0, si = 0;

	// First, deal with links within bounds of opstring
	// To do this, iterate over opstring
	for (int p = 0; p < effexporder; p++) {
		// Check if operator at p is acting on bonds or spins
		if (opstring[p][0] > 0 && opstring[p][1] > -1) {
			// Operator is acting on bonds

			// Set vertex index
			v0 = 4*p;
			// Get bond index the op is acting on
			bond = opstring[p][1];
			// Get indices of spins belonging to bond
			s1i = bonds[bond][0];
			s2i = bonds[bond][1];
			// Store op type
			linkops[v0] = 1;
			linkops[v0+1] = 1;
			linkops[v0+2] = 1;
			linkops[v0+3] = 1;
			// Get last vertex indices of each spin
			v1 = lasts[s1i];
			v2 = lasts[s2i];
			// Check if spins have appeared in vertex before
			if (v1 > -1) {
				// First spin has appeared before, so set link
				links[v1] = v0;
				links[v0] = v1;
			} else {
				// Set first appearance of spin as current leg
				firsts[s1i] = v0;
			}
			if (v2 > -1) {
				// Second spin has appeared before, so set link
				links[v2] = v0+1; // ?? are we sure this is supposed to be v0 and not v0+1?
				links[v0+1] = v2;
			} else {
				// Set first appearance of spin as current p
				firsts[s2i] = v0+1;
			}
			// Update latest vertex appearance
			lasts[s1i] = v0+2;
			lasts[s2i] = v0+3;
		} else if (opstring[p][0] > 0 && opstring[p][2] > -1) {
			// Operator is acting on spins

			// Set vertex index
			v0 = 4*p;
			// Get bond index the op is acting on
			sp = opstring[p][2];
			// Store op type
			if (opstring[p][0] == 3) {
				// Operator is longitudinal field
				linkops[v0] = 3;
				linkops[v0+2] = 3;
			} else {
				linkops[v0] = 2;
				linkops[v0+2] = 2;
			}
			// Get last vertex indices of each spin
			v1 = lasts[sp];
			// Check if spins have appeared in vertex before
			if (v1 > -1) {
				// First spin has appeared before
				links[v1] = v0;
				links[v0] = v1;
			} else {
				// Set first appearance of spin as current p
				firsts[sp] = v0;
			}
			// Update latest vertex appearance
			lasts[sp] = v0+2;
		}
	}

    // Finally, construct links across boundary
	int first = 0, last = 0, pfirst = 0, plast = 0;
	for (int i = 0; i < nspins; i++) {
		first = firsts[i];
		if (first > -1) {
			last = lasts[i];

			// Get operators first and last legs belong to
			pfirst = floor(first/(double)4);
			plast = floor(last/(double)4);

			links[first] = last;
			links[last] = first;
			cross_boundary_link[first] = true;
			cross_boundary_link[last] = true;
		}
	}

	// ----- Then trace all vertex loops and do flip updates -----

	// Initialise variables
	bool finished = false, fliploop = true, allspinops = false, cont = true, bondops_loop[effexporder] = {false},
        spinops_loop[effexporder] = {false}, visited[4*effexporder] = {false};
	int currentv = 0, linkedv = 0, lastv = 0, i = 0, p = 0, linkedp = 0, change_counter = 0,
        free_spins[nspins] = {false};
	for (int i = 0; i < nspins; i++) {
		// This encodes, whether a spin is free or part of
		// a loop that has been flipped
		// 0 - spin is free and no op is acting on it
		// 1 - spin isn't free but it isn't part of any
		//     flipped loop
		// 2 - spin isn't free and it's part of a flipped loop
		// (currently only items with value 2 are used, so might
		// convert this to a bool array)
		free_spins[i] = 0;
	}
	double coin = 0;
	bool flip_spin[nspins] = {false};

	// Iterate over all vertex legs
	for (int v = 0; v < 4*effexporder; v++) {
		// Check if vertex leg is connected to any other legs
		if (links[v] > -1 && !visited[v]) {
            // Set up variables for loop
            cont = true;
			finished = false;
            fliploop = true;

            // Set current leg as initial leg
			currentv = v;

            // See if initial leg belongs to bond or spin op
            if (linkops[currentv] == 1) {
                // Initial leg belongs to bond op

				visited[currentv] = true;

                // Get opstring index
                p = floor(currentv/(double)4);
                
                // Operator is part of loop
                bondops_loop[p] = true;
            } else if (linkops[currentv] == 2) {
                // Initial leg belongs to transverse field op

                visited[currentv] = true;

                // Get linked leg
                linkedv = links[currentv];
                
				// Check if this operator is linked to anything
				if (linkedv > -1) {
					// Check type of operator this leg is linked to
					if (linkops[linkedv] == 1) {
						// Linked to bond op
						currentv = linkedv;

						// Get opstring index
						p = floor(currentv/(double)4);
						
						// Operator is part of loop
						bondops_loop[p] = true;
					} else if (linkops[linkedv] == 2) {
						// Linked to spin op

						visited[linkedv] = true;

						// Get op indices belonging to this cluster
						p = floor(currentv/(double)4);
						linkedp = floor(linkedv/(double)4);
						
						// Check if the two ops are different
						if (p != linkedp) {
							// Set ops flippable
							spinops_loop[p] = true;
							spinops_loop[linkedp] = true;

							// Check if cluster reaches across time boundary
							if (cross_boundary_link[linkedv]) {
								// Cluster is across boundary, so get spin index and
								// set "to be flipped"
								si = opstring[linkedp][2];
								flip_spin[si] = true;
							}
						} else {
							// Same operator, so cannot flip
							fliploop = false;
						}

						// Stop further loop from happening
						cont = false;
					}
				}
            }
            
			// Traverse loop
			while (cont && !finished) {
                // Reset counter of how many ops have been added to loop
                change_counter = 0;

                // Iterate over all operators
                for (int pi = 0; pi < effexporder; pi++) {
                    if (bondops_loop[pi]) {
                        // Iterate over legs of op
                        // Go through legs belonging to bond op
                        for (int i = 0; i < 4; i++) {
                            // Get current vertex leg
                            currentv = 4*pi+i;
                            visited[currentv] = true;
                            
                            linkedv = links[currentv];
							
							// Check if leg is linked to anything
							if (linkedv > -1) {
								// Leg is linked to something, so get linked op
								linkedp = floor(linkedv/(double)4);
								// Check if this linked op is a bond op
								if (linkops[linkedv] == 1) {
									// Check if bond op is already part of loop
									if (!bondops_loop[linkedp]) {
										// Bond op hasn't been added yet, so add it
										bondops_loop[linkedp] = true;
										change_counter++;
									}
								}
							}
                        }
                    }
                }

				// If loop hasn't grown in this iteration, finish the loop
                if (change_counter == 0) {
                    finished = true;
                }
	        }

            // Now see if all legs of the bond ops in the cluster
            // belong to spin ops or other bond ops in the cluster
            for (int pi = 0; pi < effexporder; pi++) {
				// Check if operator is part of current loop
                if (bondops_loop[pi]) {
					// Operator is part of current loop
					
					// Get spins associated with operator and set as "part of loop"
					bond = opstring[pi][1];
					s1i = bonds[bond][0];
					s2i = bonds[bond][1];
					free_spins[s1i] = 1;
					free_spins[s2i] = 1;

                    // Go through legs
                    for (int i = 0; i < 4; i++) {
                        // Get current vertex leg
                        currentv = 4*pi+i;
                        visited[currentv] = true;
                        
                        linkedv = links[currentv];
                        linkedp = floor(linkedv/(double)4);

						// Only proceed if this leg is linked to anything
						if (linkedv > -1) {
							// Check if this linked op is a bond op or spin op
							if (linkops[linkedv] == 1) {
								// !! why do what's below ??

								// Check if cluster reaches across time boundary
								if (cross_boundary_link[linkedv]) {
									// Cluster is across boundary, so get spin index and
									// set "to be flipped"
									bond = opstring[linkedp][1];
									s1i = bonds[bond][0];
									s2i = bonds[bond][1];
									if ((linkedv-4*linkedp) == 3 || (linkedv-4*linkedp) == 1) {
										flip_spin[s2i] = true;
									} else if ((linkedv-4*linkedp) == 2 || (linkedv-4*linkedp) == 0) {
										flip_spin[s1i] = true;
									}
								}

								// Check if linked bond op is already part of loop
								if (!bondops_loop[linkedp]) {
									// This SHOULD NOT happen, a mere safety measure
									
									// !!!! ?? fliploop = false;
								}
							} else if (linkops[linkedv] == 2) {
								// This leg is linked to a spin op

								visited[linkedv] = true;
								
								// Check if spin op is already part of loop
								if (spinops_loop[linkedp]) {
									// Spin op has already been added on different leg
									fliploop = false;
								} else {
									// Spin op hasn't been added before
								
									// Add spin op to list of ops that should be flipped
									
									spinops_loop[linkedp] = true;
									
									// Check if link reaches across time boundary
									if (cross_boundary_link[linkedv]) {
										// Cluster is across boundary, so get spin index and
										// set "to be flipped"
										si = opstring[linkedp][2];
										flip_spin[si] = true;
									}
								}
							} else {
								// This leg isn't connected to anything
								fliploop = false;
							}
						} else {
							// This leg isn't connected to anything
							fliploop = false;
						}
                    }
                }
			}

            // Check if loop should be flipped
            if (fliploop && uni_dist(rng) < 0.5) {
                // Loop is flippable

                // Iterate over all spin ops and change type
                for (int pi = 0; pi < effexporder; pi++) {
                    // Check if spin op is part of loop
                    if (spinops_loop[pi]) {
                        // Flip spin op
                        opstring[pi][0] = (opstring[pi][0] == 1) ? 2 : 1;

                        // Set spin as part of flipped cluster
                        sp = opstring[pi][2];
                        free_spins[sp] = 2;
                    }
                }
				// Cluster has been flipped, so add flippable spins to list
				for (int s = 0; s < nspins; s++) {
					if (flip_spin[s]) {
						firsts[s] = -2;
					}
				}
            } else {
				// Loop isn't flipped
				// ...
				// Nothing happens here, so this is redundant
            }
        }

		// Empty list of bond and spin ops in loop
		for (int i = 0; i < effexporder; i++) {
			bondops_loop[i] = false;
			spinops_loop[i] = false;
		}
		// ? Necessary ?
		for (int s = 0; s < nspins; s++) {
			flip_spin[s] = false;
		}
    }

	// ----- Finally, take care of flipping spins in the periodic state -----

	// Iterate over spins i
	for (int i = 0; i < nspins; i++) {
		// Check if spin is free
		//if (free_spins[i] == 0) {
		if (firsts[i] == -1) {
			// Spin is free, so flip with 1/2 probability
			if (uni_dist(rng) < 0.5) {
				spins[i] *= -1;
			}
		} else if (firsts[i] == -2) {
			// Spin isn't free but is part of a cluster that has been flipped
			// across the periodic time boundary
			spins[i] *= -1;
		}

		// Reset free-ness of spin
		free_spins[i] = 0;
	}

	return 0;

}

int adjust_maxexporder(int opstring[][3], int effexporder, mt19937 &rng, int scale_up_to) {

	// Initialise random number generator
	uniform_real_distribution<double> uni_dist(0,1);

	// Evaluate current expansion order
	int exporder = effexporder;
	for (int i = 0; i < effexporder; i++) {
		if (opstring[i][0] == 0) {
			exporder = exporder - 1;
		}
	}

	// Initialise new expansion order
	int newexporder = 1;
	//cout << exporder << " " << newexporder << "\n";

	// To distribute operators in new opstring evenly, define
	// insertion probability
	int ops_inserted, old_op = 0;
	double insert_prob = 0;
	
	// Is there a target max expansion order?
	if (scale_up_to == 0) {
		// Target max expansion order
		if (exporder < 30) {
			// This might be necessary to account for fluctuations
			// in expansion order at high T and low expansion orders
			insert_prob = 0.2;
		} else {
			insert_prob = 0.67;
		}
	} else {
		insert_prob = exporder/(double) scale_up_to;
	}

	// Initialise variables
	int newopstring[100000][3];
	for (int p = 0; p < 100000; p++) {
		newopstring[p][0] = 0;
		newopstring[p][1] = -1;
		newopstring[p][2] = -1;
	}
	int opcount = 0;
	bool op_inserted;

	for (int p = 0; p < effexporder; p++) {
		// Check if operator isn't identity
		if (opstring[p][0] > 0) {
			// Try to insert operator p until it worked
			op_inserted = false;

			while (!op_inserted) {
				// Check if probability is fulfilled
				if (uni_dist(rng) < insert_prob) {
					newopstring[opcount][0] = opstring[p][0];
					newopstring[opcount][1] = opstring[p][1];
					newopstring[opcount][2] = opstring[p][2];
					opcount++;

					// Move on to next operator
					op_inserted = true;
					old_op++;
				} else {
					// Insert identity operator
					newopstring[opcount][0] = 0;
					newopstring[opcount][1] = -1;
					newopstring[opcount][2] = -1;
					opcount++;
				}
			}
		}
	}

	// If no operators are in the opstring, set the new max exporder as 5
	newexporder = (opcount > 0) ? opcount : 5;

	for (int p = 0; p < newexporder; p++) {
		opstring[p][0] = newopstring[p][0];
		opstring[p][1] = newopstring[p][1];
		opstring[p][2] = newopstring[p][2];
	}

	return newexporder;
	
}

int measure_observables(int spins[], int nspins, int bonds[][2], int nbonds,
	int opstring[][3], int effexporder, double temp, double transfield, int obs_count,
	double obs_exporder[], double obs_exporder_sq[], double obs_magn[], 
	double obs_magn_sq[], double obs_magn_quad[], double obs_trans_magn[], 
	double obs_trans_magn_sq[], double obs_corr[500][500]) {

	// ----- Measure internal energy (prop. to expansion order) -----

	// Initialise expansion order
	int exporder = 0;

	// Evaluate current expansion order
	exporder = effexporder;
	for (int i = 0; i < effexporder; i++) {
		if (opstring[i][0] == 0) {
			exporder = exporder - 1;
		}
	}
	obs_exporder[obs_count] = exporder;
	obs_exporder_sq[obs_count] = pow(exporder, 2);

	// Copy spins to modify later on
	int spins_mod[nspins] = {};
	for (int i = 0; i < nspins; i++) {
		spins_mod[i] = spins[i];
	}
	// Index of spin operator is acting on
	int si = -1;

	// Initialise variables
	int spin = 0;
	double avg_magn = 0, avg_magn_sq = 0, avg_magn_quad = 0;
	int magn = 0;

	// Initialise variables
	int s1i = 0, s1 = 0, s2i = 0, s2 = 0, total_corr = 0;
	int corrs[nspins][nspins] = {0};
	double avg_corr = 0;
	bool zerops = false;
	int iterations = effexporder;
	int divideby = exporder;

	// Check if expansion order is zero to allow averaging
	// loop to work if it is
	if (exporder == 0) {
		zerops = true;
		iterations = 1;
		divideby = 1;
	}

	for (int p = 0; p < iterations; p++) {
					//cout << "loop1";
		// Check if current spin state is about to be changed,
		// (or if there the opstring is empty)
		// so to avoid including states between which there is
		// only an identity
		if (opstring[p][0] != 0 || zerops) {
			// Reset observables
			magn = 0;

			// Get observables of current state

			// ----- Measure magnetisation -----

			// Iterate over all spins
			for (int i = 0; i < nspins; i++) {
				magn += spins_mod[i];
			}

			avg_magn += abs(magn)/(double)divideby;
			avg_magn_sq += pow(magn, 2)/(double)divideby;
			avg_magn_quad += pow(magn, 4)/(double)divideby;
			
			// ----- Find spin-spin correlations -----

			// Iterate over all spin pairs
			for (int s1i = 0; s1i < nspins; s1i++) {
				for (int s2i = 0; s2i < s1i; s2i++) {
					// Get spin states
					s1 = spins_mod[s1i];
					s2 = spins_mod[s2i];

					// Add product to list
					obs_corr[s1i][s2i] += s1*s2/(double)divideby;
					obs_corr[s2i][s1i] += s1*s2/(double)divideby;
				}
			}
		}

		// Apply operators to state (propagate state)
		if (opstring[p][0] == 2) {
			// Operator acting is off-diagonal

			if (opstring[p][2] > -1) {
				// Operator is acting on spin

				si = opstring[p][2];

				// Flip spin
				spins_mod[si] *= -1;
			}
		}
	}

	
	// Get averages and add to lists
	//avg_magn = avg_magn/(double)(nspins); //*effexporder);
	obs_magn[obs_count] = avg_magn/(double)nspins; // abs(avg_magn);
	obs_magn_sq[obs_count] = avg_magn_sq/(double) pow(nspins, 2);
	obs_magn_quad[obs_count] = avg_magn_quad/(double) pow(nspins, 4);

	// Evaluate transverse magnetisation (see Sandvik book chapter)
	// To do so, count all transverse field operators
	int count_trans_ops = 0;
	for (int p = 0; p < effexporder; p++) {
		// Check if operator at position p is trans field op
		if (opstring[p][0] == 2 && opstring[p][2] > -1) {
			count_trans_ops++;
		}
	}
	double trans_magn = 0;
	trans_magn = count_trans_ops*temp/(double) transfield;
	obs_trans_magn[obs_count] = trans_magn/(double) nspins;
	obs_trans_magn_sq[obs_count] = pow(trans_magn, 2)/(double) pow(nspins, 2);

	// Need to implement more observables ?
	// ...

	return 0;
	
}