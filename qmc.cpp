/*

Quantum Monte Carlo simulation of transverse-field
Ising models by stochastic series expansion.

(C) 2025 Theo Stölzl

This software is licensed under the CC-BY-SA 4.0 license
(https://creativecommons.org/licenses/by-sa/4.0/).

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

// TO IMPLEMENT:
// 1 - flipping of single spins after flip operator updates
//     and flipping all spins pairwise who have bond ops in between (DONE)
// 2 - flipping of cluster spins (equiv to Fig 59, Sandvik) ? (DONE)
// 3 - flip updates of clusters of 2 spin ops sandwiching
//     bond ops (1 or several) (DONE)
// 4 - expansion cut off (DONE)
// 5 - input file (DONE)
// 6 - J-coupling for individual bonds (DONE)
// 7 - measuring observables (heat capacity !)
// 8 - look into why max cutoff needs to be about 3x the mean expansion order?
// 9 - averaging into bins (DONE)
// 10 - set up random functions for whole programme, not in each
//      function separately (DONE)

int random_conf(int spins[], int nspins, mt19937 &rng);
int diagonal_updates(int spins[], int nspins, int bonds[][2], double couplings[], int nbonds,
	int opstring[][3], int effexporder, double temp, double hfield, mt19937 &rng, ofstream &cont_out);
int cluster_updates(int spins[], int nspins, int bonds[][2],
	int nbonds, int opstring[][3], int effexporder, mt19937 &rng, ofstream &cont_out);
int adjust_maxexporder(int opstring[][3], int effexporder, mt19937 &rng);
int measure_observables(int spins[], int nspins, int bonds[][2], int nbonds,
	int opstring[][3], int effexporder, int obs_count, double obs_exporder[], 
	double obs_exporder_sq[], double obs_magn[], double obs_magn_sq[], 
	double obs_magn_quad[], double obs_nn_corr[]);

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
	int spins[nspins];
	int bonds[nbonds][2];
	double couplings[nbonds];

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

		cout << bonds[l][0] << "\t" << bonds[l][1] << "\t" << couplings[l] << "\n";
	}

	// ----- Read in simulation parameters -----

	ifstream f("setup.json");
	nlohmann::json data = nlohmann::json::parse(f);
	int eqsweeps = data["eqsweeps"];
	int bins = data["bins"];
	int avsweeps = data["avsweeps"];
	int maxexporder = data["max_expansion_order"];
	double temp = data["temperature"];
	double hfield = data["transverse_field"];
	// Overwrite using in-line arguments
	temp = stod(argv[2]);
	cout << "temp " << temp << "\n";
	hfield = stod(argv[3]);
	cout << "field " << hfield << "\n";
	string out_path = argv[1];
	
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
	int effexporder = maxexporder;

	// Initialise operator string
	int opstring[100000][3] = {};
	/*for (int i = 0; i < effexporder; i++) {
		//opstring.push_back(vector<int> {0, -1, -1});
		opstring[i][0] = 0;
		opstring[i][1] = -1;
		opstring[i][2] = -1;
	}*/

	for (int p = 0; p < 100000; p++) {
		// Type of operator at position p
		// 0 - identity op
		// 1 - diagonal op
		// 2 - off-diagonal op
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

			// Set opstring
			opstring[l][0] = stoi(myline[0]);
			opstring[l][1] = stoi(myline[1]);
			opstring[l][2] = stoi(myline[2]);
		}
	}
	ops_infile.close();

	cout << "Started equilibration ..." << "\n";
	cout << flush;

	// Start equilibration
	for (int i = 0; i < eqsweeps; i++) {

		ofstream cont_out; 
		cont_out.open("sweeps/cont_out_"+to_string(i)+".txt");
		
		cont_out << "===== BEGIN SWEEP =====\n\n";

		// Debug spins
		cont_out << "----- BEFORE SPINS -----\n";
		for (int s = 0; s < nspins; s++) {
			cont_out << spins[s] << "\n";
		}
		cont_out << "----------\n";

		cont_out << "----- BEFORE OPSTRING -----\n";
		for (int p = 0; p < effexporder; p++) {
			cont_out << p << "\t" << opstring[p][0] << "\t" << opstring[p][1] << "\t" << opstring[p][2] << "\n";
		}
		cont_out << "----------\n";

		// Diagonal updates to insert / remove operators
		diagonal_updates(spins, nspins, bonds, couplings, nbonds, opstring, 
				effexporder, temp, hfield, rng, cont_out);
		
		cont_out << "----- AFTER DIAGONAL SPINS -----\n";
		for (int s = 0; s < nspins; s++) {
			cont_out << spins[s] << "\n";
		}
		cont_out << "----------\n";
		
		// Cluster updates to vary diagonal / off-diagonal ops
		cluster_updates(spins, nspins, bonds, nbonds, opstring, effexporder, rng, cont_out);
		
		// Adjust largest expansion order
		effexporder = adjust_maxexporder(opstring, effexporder, rng);
		//cout << "max exp order: " << effexporder << "\n";
		//effexporder = maxexporder;

	/*
				cout << "\n";
			for (int p = 0; p < effexporder; p++) {
				cout << opstring[p][0] << "\t" << opstring[p][1] << "\t" << opstring[p][2] << "\n";
			}

			cout << "\n\n";
	*/

		// Debug spins
		cont_out << "----- AFTER SPINS -----\n";
		for (int s = 0; s < nspins; s++) {
			cont_out << spins[s] << "\n";
		}
		cont_out << "----------\n";

		cont_out << "----- AFTER OPSTRING -----\n";
		for (int p = 0; p < effexporder; p++) {
			cont_out << p << "\t" << opstring[p][0] << "\t" << opstring[p][1] << "\t" << opstring[p][2] << "\n";
		}
		cont_out << "----------\n";
		cont_out.close();

		/*
		// DEBUG
		int counter_spin_flips[nspins] = {0};
		int spini = 0;

		// DEBUG
		for (int x = 0; x < effexporder; x++) {
			if (opstring[x][0] == 2) {
				if (opstring[x][2] > -1) {
					counter_spin_flips[opstring[x][2]]++;
				} else {
					cout << "ERROR\n";
				}
			}
		}

		for (int x = 0; x < nspins; x++) {
			if (counter_spin_flips[x] % 2 != 0) {
				cout << "ABORT at sweep " << i;
				cout << " spin " << x << " " << counter_spin_flips[x]++;
				exit(0);
			}
		}
		*/
	}

	cout << "Finished equilibration." << "\n";
	cout << "Maximum expansion order: " << effexporder << "\n";
	cout << "Started averaging ..." << "\n";
	cout << flush;

	// Open output file
	ofstream out_file; 
	out_file.open(out_path+"/results_t_"+to_string(temp)+"_h_"+to_string(hfield)+".csv");
	out_file << "bin,exporder,exporder_sq,magn,magn_sq,magn_quad,nn_corr\n";
	
	// Setup for measuring observables
	double obs_exporder[bins];
	double obs_exporder_sq[bins];
	double obs_magn[bins];
	double obs_magn_sq[bins];
	double obs_magn_quad[bins];
	double obs_nn_corr[bins];
	for (int n = 0; n < bins; n++) {
		obs_exporder[n] = 0;
		obs_exporder_sq[n] = 0;
		obs_magn[n] = 0;
		obs_magn_sq[n] = 0;
		obs_magn_quad[n] = 0;
		obs_nn_corr[n] = 0;
	}
	
	// Start averaging sweeps
	for (int n = 0; n < bins; n++) {
		
		double obs_exporder_bin[avsweeps];
		double obs_exporder_sq_bin[avsweeps];
		double obs_magn_bin[avsweeps];
		double obs_magn_sq_bin[avsweeps];
		double obs_magn_quad_bin[avsweeps];
		double obs_nn_corr_bin[avsweeps];
		for (int k = 0; k < avsweeps; k++) {
			obs_exporder_bin[k] = 0;
			obs_exporder_sq_bin[k] = 0; 
			obs_magn_bin[k] = 0;
			obs_magn_sq_bin[k] = 0;
			obs_magn_quad_bin[k] = 0;
			obs_nn_corr_bin[k] = 0;
		}

		for (int j = 0; j < avsweeps; j++) {

			ofstream cont_out; 
			cont_out.open("sweeps/cont_out_avg_"+to_string(j)+".txt");

			// Diagonal updates to insert / remove operators
			diagonal_updates(spins, nspins, bonds, couplings, nbonds, opstring, 
					effexporder, temp, hfield, rng, cont_out);

			// Cluster updates to vary diagonal / off-diagonal ops
			cluster_updates(spins, nspins, bonds, nbonds, opstring, effexporder, rng, cont_out);
			
			// Measure observables
			measure_observables(spins, nspins, bonds, nbonds, opstring,
					effexporder, j, obs_exporder_bin, obs_exporder_sq_bin,
					obs_magn_bin, obs_magn_sq_bin, obs_magn_quad_bin,
					obs_nn_corr_bin);

			cont_out.close();

			/*
			// DEBUG
			int counter_spin_flips[nspins] = {0};
			int spini = 0;		

			// DEBUG
			for (int x = 0; x < effexporder; x++) {
				if (opstring[x][0] == 2) {
					if (opstring[x][2] > -1) {
						counter_spin_flips[opstring[x][2]]++;
					} else {
						cout << "ERROR\n";
					}
				}
			}

			for (int x = 0; x < nspins; x++) {
				if (counter_spin_flips[x] % 2 != 0) {
					cout << "ABORT at sweep e" << n << " " << j;
					exit(0);
				}
			}
			*/

		}

		// Now average over bin and add to big arrays
		for (int j = 0; j < avsweeps; j++) {
			obs_exporder[n] += obs_exporder_bin[j]/(double)avsweeps;
			obs_exporder_sq[n] += obs_exporder_sq_bin[j]/(double)avsweeps;
			obs_magn[n] += obs_magn_bin[j]/(double)avsweeps;
			obs_magn_sq[n] += obs_magn_sq_bin[j]/(double)avsweeps;
			obs_magn_quad[n] += obs_magn_quad_bin[j]/(double)avsweeps;
			obs_nn_corr[n] += obs_nn_corr_bin[j]/(double)avsweeps;
		}
		
		out_file << to_string(n)+","+to_string(obs_exporder[n])+","+to_string(obs_exporder_sq[n])+","
				+to_string(obs_magn[n])+","+to_string(obs_magn_sq[n])+","+to_string(obs_magn_quad[n])+","
				+to_string(obs_nn_corr[n])+"\n" << flush;
	
	}

	out_file.close();

	cout << "Finished averaging." << "\n";
	cout << flush;
	
	// Find average observables
	double avgn = 0, avgcorr;
	avgm = 0;
	for (int n = 0; n < bins; n++) {
		avgn += obs_exporder[n]/(double)(bins);
		//cout << obs_magn[n] << "\n";

		avgm += obs_magn[n]/(double)(bins);
		avgcorr += obs_nn_corr[n]/(double)(bins);
	}
	
	cout << "average expansion order: " << avgn << "\n";
	cout << "average magnetisation: " << avgm << "\n";
	cout << "average nn correlation: " << avgcorr << "\n";
	
	/*
	for (int i = 0; i < nspins; i++) {
		cout << spins[i] << "\t"; 
		if (i%16 == 0) {
			cout << "\n"; 
		}
	}
	*/

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
	int opstring[][3], int effexporder, double temp, double hfield, mt19937 &rng, ofstream &cont_out) {

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
	int op_type, randb, rands, s1i, s1, s2i, s2, si;
	double coin, paccept, bj;

	// Iterate over all positions p in operator string
	for (int p = 0; p < effexporder; p++) {
		op_type = opstring[p][0];
		
		// Evaluate the type of operator present at p
		if (op_type == 0) {
			// Identity operator
			
			// Flip a coin and decide if to proceed with 
			// insertion of bond or spin operator
			coin = uni_dist(rng);
			if (coin < 0.5) {
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
					// paccept = (1/temp) * nbonds * bj /(double) (2*(effexporder - exporder)); // original (wrong!)
					// paccept = (1/temp) * nbonds * bj /(double) (effexporder - exporder); // fixed (spin 1/2)
					paccept = (1/temp) * 2 * nbonds * 2 * bj /(double) (effexporder - exporder); // spin 1
					
						//cout << "in" << exporder << "\t" << effexporder << "\t" << 1/temp 
								//<< "\t" << nbonds << "\t" << bj << "\t" << paccept << "\n";

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
			} else {
				// Try to insert a diagonal spin operator
				
				// Get random position
				rands = rand_spin(rng);
				
				// Evaluate acceptance probability
				//paccept = (1/temp) * nspins * hfield /(double) (effexporder - exporder);
				paccept = (1/temp) * 2 * nspins * hfield /(double) (effexporder - exporder);
				
				// Attempt move
				if (uni_dist(rng) < paccept) {
					// Insert diagonal spin operator
					opstring[p][0] = 1;
					opstring[p][1] = -1;
					opstring[p][2] = rands;
					exporder++;
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
				//paccept = 2*(effexporder - exporder + 1) /(double) ((1/temp) * nbonds * bj); // original (wrong!)
				// paccept = (effexporder - exporder + 1) /(double) ((1/temp) * nbonds * bj); // spin 1/2
				paccept = (effexporder - exporder + 1) /(double) ((1/temp) * 2 * nbonds * 2 * bj); // spin 1

				//cout << "out" << exporder << "\t" << effexporder << "\t" << 1/temp 
						//<< "\t" << nbonds << "\t" << bj << "\t" << paccept << "\n";
				
				// Attempt move
				if (uni_dist(rng) < paccept) {
					// Remove diagonal bond operator
					opstring[p][0] = 0;
					opstring[p][1] = -1;
					opstring[p][2] = -1;
					exporder--;
				}
			} else if (opstring[p][2] > -1) {
				// Operator is acting on spins
				
				// Evaluate acceptance probability
				//paccept = (effexporder - exporder + 1) /(double) ((1/temp) * nspins * hfield);
				paccept = (effexporder - exporder + 1) /(double) ((1/temp) * 2 * nspins * hfield);
				
				// Attempt move
				if (uni_dist(rng) < paccept) {
					// Remove diagonal spin operator
					opstring[p][0] = 0;
					opstring[p][1] = -1;
					opstring[p][2] = -1;
					exporder--;
				}
			}
		} else {
			// Off-diagonal operator
			
			// Modify spin according to opstring
			si = opstring[p][2];
			spins[si] *= -1;
			cont_out << "flipped spin " << si << "\n";
		}
	}

	/* DEBUG
	// How many bond ops?
	int count_bond = 0;
	int count_sp = 0;
	for (int i = 0; i < effexporder; i++) {
		if (opstring[i][0] > 0 && opstring[i][1] > -1) {
			count_bond++;
		}
		if (opstring[i][0] > 1 && opstring[i][2] > -1) {
			count_sp++;
		}
	}
	// cout << count_bond << " bond ops\t" << count_sp << " spin ops\n";
	*/

	return 0;
	
}

int cluster_updates(int spins[], int nspins, int bonds[][2],
	int nbonds, int opstring[][3], int effexporder, mt19937 &rng, ofstream &cont_out) {
	
	// Random stuff
	uniform_real_distribution<double> uni_dist(0,1);

	// ----- Construct vertex link list -----

	// List of first and last vertex legs of each spin
	int firsts[nspins];
	int lasts[nspins];
	for (int i = 0; i < nspins; i++) {
		firsts[i] = -1;
		lasts[i] = -1;
	}
	// List of vertex links
	int links[4*effexporder];
	for (int i = 0; i < 4*effexporder; i++) {
		links[i] = -1;
	}
	// List of operator types each vertex leg belongs to
	int linkops[4*effexporder];
	for (int i = 0; i < 4*effexporder; i++) {
		// This encodes the type of operator each 
		// vertex leg belongs too
		// 0 - identity op
		// 1 - bond op
		// 2 - spin op
		linkops[i] = 0;
	}

	// Initialise variables
	int v0, v1, v2, s1i, s2i, bond, sp, si;

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
			linkops[v0] = 2;
			linkops[v0+2] = 2;
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

	/* LEGACY - need to figure this out !
	
	Deactivate cross-boundary links completely to avoid illegal moves,
	where two operators are updates and consequently the periodic state
	would have to be updated, which might lead to clashes with other
	operators that aren't involved in the cluster update.

    // Finally, construct links across boundary
	int first, last, pfirst, plast;
	for (int i = 0; i < nspins; i++) {
		first = firsts[i];
		if (first > -1) {
			last = lasts[i];

			// Get operators first and last legs belong to
			pfirst = floor(first/(double)4);
			plast = floor(last/(double)4);

			// To avoid bond ops linking to themselves and
			// make our lives easier later when tracing loops,
			// ignore links between bond ops across the periodic
			// time boundary
			if ((linkops[first] != 1 || linkops[last] != 1) && pfirst != plast) {
				links[first] = last;
				links[last] = first;
			}
		}
	}

	*/

	// ----- Then trace all vertex loops and do flip updates -----

	// Initialise variables
	bool finished, fliploop, allspinops, cont, bondops_loop[effexporder],
        spinops_loop[effexporder], visited[4*effexporder];
	for (int i = 0; i < 4*effexporder; i++) {
		visited[i] = false;
	}
	for (int i = 0; i < effexporder; i++) {
		bondops_loop[i] = false;
        spinops_loop[i] = false;
	}
	int currentv, linkedv, lastv, i, p, linkedp, change_counter,
        free_spins[nspins];
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
	double coin;

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
                // Initial leg belongs to spin op

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
						
						// Set ops flippable
						spinops_loop[p] = true;
						spinops_loop[linkedp] = true;

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
				
				// DEBUG
				int counter_spin_flips[nspins] = {0};
				int spini = 0;
				
				string output = "";

				cont_out << "-- FLIPPABLE OPSTRING --\n";
				for (int p = 0; p < effexporder; p++) {
					cont_out << p << "\t" << opstring[p][0] << "\t" << opstring[p][1] << "\t" << opstring[p][2] << "\n";
				}
				cont_out << "----\n";

				output += "-- FLIPPABLE OPSTRING --\n";
				for (int p = 0; p < effexporder; p++) {
					output += to_string(p) + "\t" + to_string(opstring[p][0]) + "\t" 
						+ to_string(opstring[p][1]) + "\t" + to_string(opstring[p][2]) + "\n";
				}
				output += "----\n";

                // Loop is flippable
				//cout << "-----\n";
                // Iterate over all spin ops and change type
                for (int pi = 0; pi < effexporder; pi++) {
                    // Check if spin op is part of loop
                    if (spinops_loop[pi]) {
                        // Flip spin op
						// DEBUG
						cont_out << pi << "\t" << opstring[pi][0] << "\t"
							<< opstring[pi][1] << "\t"
							<< opstring[pi][2] << "\t";

                        opstring[pi][0] = (opstring[pi][0] == 1) ? 2 : 1;

						// DEBUG
						spini = opstring[pi][2];

						counter_spin_flips[spini]++;
						// DEBUG
						cont_out << opstring[pi][0] << "\t";
						cont_out << "flip op\n";

						output += "flip op " + to_string(pi) + "\n";

					// cout << "flip op\n";
                        // Set spin as part of flipped cluster
                        sp = opstring[pi][2];
                        free_spins[sp] = 2;
                    }
                }

				cont_out << "-- FLIPPED OPSTRING --\n";
				for (int p = 0; p < effexporder; p++) {
					cont_out << p << "\t" << opstring[p][0] << "\t" << opstring[p][1] << "\t" << opstring[p][2] << "\n";
				}
				cont_out << "----\n";
				//cout << "-----\n";

				/*
				// DEBUG
				for (int x = 0; x < nspins; x++) {
					if (counter_spin_flips[x] % 2 != 0) {
						ofstream outputf;
						outputf.open("errors_"+to_string(effexporder)+".txt");
						outputf << output;
						outputf << "-- FLIPPED OPSTRING --\n";
						for (int p = 0; p < effexporder; p++) {
							outputf << p << "\t" << opstring[p][0] << "\t" << opstring[p][1] << "\t" << opstring[p][2] << "\n";
						}
						outputf << "----\n";
						outputf.close();

						cout << "ABORT";
						exit(0);
					}
				}
				*/

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
		} else if (free_spins[i] == 2) {
			// Spin isn't free but is part of a cluster that has been flipped
			//spins[i] *= -1;
			//cout << i << " flipped spin" << spins[i] << "\n";
		}

		// Reset free-ness of spin
		free_spins[i] = 0;
	}
	
	//cout << "sweep\n";

	return 0;

}

int adjust_maxexporder(int opstring[][3], int effexporder, mt19937 &rng) {

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
	int newexporder = 0;
	//cout << exporder << " " << newexporder << "\n";

	// To distribute operators in new opstring evenly, define
	// insertion probability
	int ops_inserted, old_op = 0;
	//double insert_prob = (newexporder-exporder)/(double)(newexporder+1);
	double insert_prob = 0.6;
	//double insert_prob = 0.1;
	//cout << insert_prob << "\n";

	// Initialise variables
	int newopstring[100000][3];
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

	newexporder = opcount;
	//cout << newexporder << "\n";

	for (int p = 0; p < newexporder; p++) {
		opstring[p][0] = newopstring[p][0];
		opstring[p][1] = newopstring[p][1];
		opstring[p][2] = newopstring[p][2];
	}

	return newexporder;
	
}

int measure_observables(int spins[], int nspins, int bonds[][2], int nbonds,
	int opstring[][3], int effexporder, int obs_count, double obs_exporder[], 
	double obs_exporder_sq[], double obs_magn[], double obs_magn_sq[], 
	double obs_magn_quad[], double obs_nn_corr[]) {

	// ----- Measure internal energy (prop. to expansion order) -----

	// Initialise expansion order
	int exporder;

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
	double avg_magn, avg_magn_sq, avg_magn_quad = 0;
	int magn = 0;

	// Initialise variables
	int s1i, s1, s2i, s2, total_corr = 0;
	int nn_corrs[nspins] = {0};
	double avg_corr;

	// Check if expansion order is zero to allow averaging
	if (effexporder == 0) {
	// if (true) {
		// Get average of just periodic basis state

		// Reset observables
		magn = 0;

		// Get observables of current state

		// ----- Measure magnetisation -----

		// Iterate over all spins
		for (int i = 0; i < nspins; i++) {
			// Initial spin state
			//spin = spins_mod[i];
			//cout << spin;
			
			/*
			// For each spin, iterate over opstring
			for (int p = 0; p < effexporder; p++) {
				avg_magn += double(spin);
				// Check if an off-diagonal operator is acting on this spin
				if (opstring[p][0] == 2 && opstring[p][2] == i) {
					spin = spin*(double)(-1);
				}
			}
			*/
			
			// DEBUG
			/*
			if (spin != spin) {
				cout << "spin is nan" << "\n";
			}*/

			magn += spins_mod[i];
		}

		avg_magn += abs(magn);

		/*
		// DEBUG
		if (avg_magn != avg_magn) {
			cout << "avg_magn is nan" << "\n";
		}
		*/

		avg_magn_sq += pow(magn, 2);
		avg_magn_quad += pow(magn, 4);
		
		// ----- Find nearest neighbour spin-spin correlations -----

		// Iterate over all bonds
		for (int b = 0; b < nbonds; b++) {
			// Get spins belonging to bond b
			s1i = bonds[b][0];
			s2i = bonds[b][1];

			// Get spin states
			s1 = spins_mod[s1i];
			s2 = spins_mod[s2i];

			// Add product to list
			nn_corrs[s1i] = s1*s2;
			nn_corrs[s2i] = s1*s2;
		}

		// Now iterate over all spins and get average correlation
		total_corr = 0;
		for (int s = 0; s < nspins; s++) {
			total_corr += nn_corrs[s];
		}
	} else {
		// Iterate over all propagated states
		for (int p = 0; p < effexporder; p++) {
			// Check if current spin state is about to be changed,
			// so to avoid including states between which there is
			// only an identity
			if (opstring[p][0] != 0) {
				// Reset observables
				magn = 0;

				// Get observables of current state

				// ----- Measure magnetisation -----

				// Iterate over all spins
				for (int i = 0; i < nspins; i++) {
					// Initial spin state
					//spin = spins_mod[i];
					//cout << spin;
					
					/*
					// For each spin, iterate over opstring
					for (int p = 0; p < effexporder; p++) {
						avg_magn += double(spin);
						// Check if an off-diagonal operator is acting on this spin
						if (opstring[p][0] == 2 && opstring[p][2] == i) {
							spin = spin*(double)(-1);
						}
					}
					*/
					
					// DEBUG
					/*
					if (spin != spin) {
						cout << "spin is nan" << "\n";
					}
					*/

					magn += spins_mod[i];
				}

				avg_magn += abs(magn)/(double)exporder;

				// DEBUG
				if (avg_magn != avg_magn) {
					cout << "avg_magn is nan" << "\n";
				}

				avg_magn_sq += pow(magn, 2)/(double)exporder;
				avg_magn_quad += pow(magn, 4)/(double)exporder;
				
				// ----- Find nearest neighbour spin-spin correlations -----

				// Iterate over all bonds
				for (int b = 0; b < nbonds; b++) {
					// Get spins belonging to bond b
					s1i = bonds[b][0];
					s2i = bonds[b][1];

					// Get spin states
					s1 = spins_mod[s1i];
					s2 = spins_mod[s2i];

					// Add product to list
					nn_corrs[s1i] = s1*s2;
					nn_corrs[s2i] = s1*s2;
				}

				// Now iterate over all spins and get average correlation
				total_corr = 0;
				for (int s = 0; s < nspins; s++) {
					total_corr += nn_corrs[s];
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
	}

	
	// Get averages and add to lists
	//avg_magn = avg_magn/(double)(nspins); //*effexporder);
	obs_magn[obs_count] = avg_magn/(double)nspins; // abs(avg_magn);
	obs_magn_sq[obs_count] = avg_magn_sq/(double)nspins;
	obs_magn_quad[obs_count] = avg_magn_quad/(double)nspins;

	// Normalise correlation and add to list
	avg_corr = total_corr/(double)(exporder*nspins);
	obs_nn_corr[obs_count] = avg_corr;

	// Need to implement more observables
	// ...

	return 0;
	
}