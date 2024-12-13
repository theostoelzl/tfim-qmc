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
// 2 - flipping of cluster spins (equiv to Fig 59, Sandvik) ? (DONE ?)
// 3 - flip updates of clusters of 2 spin ops sandwiching
//     bond ops (1 or several) (DONE)
// 4 - expansion cut off
// 5 - input file (DONE)
// 6 - J-coupling for individual bonds (DONE)
// 7 - measuring observables
// 8 - look into why max cutoff needs to be about 3x the mean expansion order?
// 9 - averaging into bins
// 10 - set up random functions for whole programme, not in each
//      function separately

int random_conf(int spins[], int nspins);
int diagonal_updates(int spins[], int nspins, int bonds[][3], int nbonds,
	int opstring[][3], int effexporder, double temp, double hfield);
int cluster_updates(int spins[], int nspins, int bonds[][3],
	int nbonds, int opstring[][3], int effexporder);
int adjust_maxexporder(int opstring[][3], int effexporder);
int measure_observables(int spins[], int nspins, int bonds[][3], int nbonds,
	int opstring[][3], int effexporder, int obs_count, int obs_exporder[], 
	double obs_magn[], double obs_magn_sq[], double obs_magn_quad[],
	double obs_nn_corr[]);

int main() {

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
	int bonds[nbonds][3];

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
		bonds[l][2] = stoi(myline[3]);

		cout << bonds[l][0] << "\t" << bonds[l][1] << "\t" << bonds[l][2] << "\n";
	}

	// ----- Read in simulation parameters -----

	ifstream f("setup.json");
	nlohmann::json data = nlohmann::json::parse(f);
	int eqsweeps = data["eqsweeps"];
	int avsweeps = data["avsweeps"];
	int maxexporder = data["max_expansion_order"];
	double temp = data["temperature"];
	double hfield = data["transverse_field"];
	
	// Take input numbers of sweeps as sweeps per spin
	eqsweeps *= nspins;
	avsweeps *= nspins;

	// ----- Randomise initial state -----
	random_conf(spins, stoi(in_nspins));

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

	for (int p = 0; p < effexporder; p++) {
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

	cout << "Started equilibration ..." << "\n";
	
	// Start equilibration
	for (int i = 0; i < eqsweeps; i++) {
		
		// Diagonal updates to insert / remove operators
		diagonal_updates(spins, nspins, bonds, nbonds, opstring, 
				effexporder, temp, hfield);
		
		// Cluster updates to vary diagonal / off-diagonal ops
		cluster_updates(spins, nspins, bonds, nbonds, opstring, effexporder);
		
		// Adjust largest expansion order
		effexporder = adjust_maxexporder(opstring, effexporder);
		//cout << "max exp order: " << effexporder << "\n";
		//effexporder = maxexporder;
	}

	cout << "Finished equilibration." << "\n";
	cout << "Maximum expansion order: " << effexporder << "\n";
	cout << "Started averaging ..." << "\n";
	
	// Setup for measuring observables
	int obs_exporder[avsweeps/100];
	double obs_magn[avsweeps/100];
	double obs_magn_sq[avsweeps/100];
	double obs_magn_quad[avsweeps/100];
	double obs_nn_corr[avsweeps/100];
	for (int n = 0; n < avsweeps/100; n++) {
		obs_exporder[n] = -1;
		obs_magn[n] = -1;
		obs_nn_corr[n] = -2;
	}
	int obs_count = 0;
	
	// Start averaging sweeps
	for (int i = 0; i < avsweeps; i++) {
		
		// Diagonal updates to insert / remove operators
		diagonal_updates(spins, nspins, bonds, nbonds, opstring, 
				effexporder, temp, hfield);
		
		// Cluster updates to vary diagonal / off-diagonal ops
		cluster_updates(spins, nspins, bonds, nbonds, opstring, effexporder);

		// Measure observables (currently only expansion order)
		if (i%100 == 0) {
			measure_observables(spins, nspins, bonds, nbonds,
				opstring, effexporder, obs_count, obs_exporder, 
				obs_magn, obs_magn_sq, obs_magn_quad,
				obs_nn_corr);
			obs_count++;
		}
		
	}

	cout << "Finished averaging." << "\n";
	
	// Find average observables
	double avgn = 0, avgcorr;
	avgm = 0;
	for (int n = 0; n < obs_count; n++) {
		avgn += obs_exporder[n]/(double)(avsweeps/100);
		//cout << obs_magn[n] << "\n";

		avgm += obs_magn[n]/(double)(avsweeps/100);
		avgcorr += obs_nn_corr[n]/(double)(avsweeps/100);
	}
	
	cout << "average expansion order: " << avgn << "\n";
	cout << "average magnetisation: " << avgm << "\n";
	cout << "average nn correlation: " << avgcorr << "\n";

	
	for (int i = 0; i < nspins; i++) {
		cout << spins[i] << "\t"; 
		if (i%16 == 0) {
			cout << "\n"; 
		}
	}
	
	return 0;
	
}

int random_conf(int spins[], int nspins) {
	
	// Initialise random number generator
	int rngseed = 178;
	mt19937 rng;
	rng.seed(time(NULL)+100000*rngseed);
	uniform_real_distribution<double> uni_dist(0,1);

	// Iterate through lattice and assign random spin directions 
	for (int i = 0; i < nspins; i++) {
		spins[i] = (uni_dist(rng) < 0.5) ? -1 : 1;
	}
	return 0;
	
}

int diagonal_updates(int spins[], int nspins, int bonds[][3], int nbonds,
	int opstring[][3], int effexporder, double temp, double hfield) {

	// Initialise random number generator
	int rngseed = 178;
	mt19937 rng;
	rng.seed(time(NULL)+100000*rngseed);
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
				bj = bonds[randb][2];	
				
				// Check if spins are parallel
				// if bj < 0, anti-ferromagnetic bond
				// if bj > 0, ferromagnetic bond
				if (bj*s1*s2 > 0) {
					// Spins are parallel so try to insert diagonal operator
					paccept = (1/temp) * nbonds * bj / (2*(effexporder - exporder));
				
					// Attempt move
					if (uni_dist(rng) < paccept) {
						// Insert diagonal bond operator
						opstring[p][0] = 1;
						opstring[p][1] = randb;
						exporder++;
					}
				}
			} else {
				// Try to insert a diagonal spin operator
				
				// Get random position
				rands = rand_spin(rng);
				
				// Evaluate acceptance probability
				paccept = (1/temp) * nspins * hfield /(double) (effexporder - exporder);
				
				// Attempt move
				if (uni_dist(rng) < paccept) {
					// Insert diagonal spin operator
					opstring[p][0] = 1;
					opstring[p][2] = rands;
					exporder++;
				}
			}
		} else if (op_type == 1) {
			// Diagonal operator
			
			// Check for error
			if (opstring[p][1] > -1 && opstring[p][2] > -1) {
				cout << "Something's wrong!" << "\n";
			}
			
			// Check if bond or spin operator
			if (opstring[p][1] > -1) {
				// Operator is acting on bonds
				
				// Get coupling constant of bond
				bj = bonds[opstring[p][1]][2];

				// Evaluate acceptance probability
				paccept = 2*(effexporder - exporder + 1) /(double) ((1/temp) * nbonds * bj);
				
				// Attempt move
				if (uni_dist(rng) < paccept) {
					// Remove diagonal bond operator
					opstring[p][0] = 0;
					opstring[p][1] = -1;
					exporder--;
				}
			} else if (opstring[p][2] > -1) {
				// Operator is acting on spins
				
				// Evaluate acceptance probability
				paccept = (effexporder - exporder + 1) / ((1/temp) * nspins * hfield);
				
				// Attempt move
				if (uni_dist(rng) < paccept) {
					// Remove diagonal spin operator
					opstring[p][0] = 0;
					opstring[p][2] = -1;
					exporder--;
				}
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

int cluster_updates(int spins[], int nspins, int bonds[][3],
	int nbonds, int opstring[][3], int effexporder) {
		
	int rngseed = 178;
	mt19937 rng;
	rng.seed(time(NULL)+100000*rngseed);
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

    // Finally, construct links across boundary
	int first, last;
	for (int i = 0; i < nspins; i++) {
		first = firsts[i];
		if (first > -1) {
			last = lasts[i];
			// To avoid bond ops linking to themselves and
			// make our lives easier later when tracing loops,
			// ignore links between bond ops across the periodic
			// time boundary
			if (linkops[first] != 1 || linkops[last] != 1) {
				links[first] = last;
				links[last] = first;
			}
		}
	}

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

                // Get opstring index
                p = floor(currentv/(double)4);
                
                // Operator is part of loop
                bondops_loop[p] = true;
            } else {
                // Initial leg belongs to spin op

                visited[v] = true;

                // Get linked leg
                linkedv = links[currentv];
                
                // Check type of operator this leg is linked to
                if (linkops[linkedv] == 1) {
                    // Linked to bond op
                    currentv = linkedv;

                    // Get opstring index
                    p = floor(currentv/(double)4);
                    
                    // Operator is part of loop
                    bondops_loop[p] = true;
                } else {
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
            
			// Traverse loop
			while (cont && !finished) {
                // Reset counter of how many ops have been added to loop
                change_counter = 0;

                // Iterate over all operators
                for (int p = 0; p < effexporder; p++) {
                    if (bondops_loop[p]) {
                        // Iterate over legs of op
                        // Go through legs belonging to bond op
                        for (int i = 0; i < 4; i++) {
                            // Get current vertex leg
                            currentv = 4*p+i;
                            visited[currentv] = true;
                            
                            linkedv = links[currentv];
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

				// If loop hasn't grown in this iteration, finish the loop
                if (change_counter == 0) {
                    finished = true;
                }
	        }

            // Now see if all legs of the bond ops in the cluster
            // belong to spin ops or other bond ops in the cluster
            for (int p = 0; p < effexporder; p++) {
				// Check if operator is part of current loop
                if (bondops_loop[p]) {
					// Operator is part of current loop
					
					// Get spins associated with operator and set as "part of loop"
					bond = opstring[p][1];
					s1i = bonds[bond][0];
					s2i = bonds[bond][1];
					free_spins[s1i] = 1;
					free_spins[s2i] = 1;

                    // Go through legs
                    for (int i = 0; i < 4; i++) {
                        // Get current vertex leg
                        currentv = 4*p+i;
                        visited[currentv] = true;
                        
                        linkedv = links[currentv];
                        linkedp = floor(linkedv/(double)4);

                        // Check if this linked op is a bond op
                        if (linkops[linkedv] == 1) {
                            // Check if bond op is already part of loop
                            if (!bondops_loop[linkedp]) {
                                // This SHOULD NOT happen, a mere safety measure
                                fliploop = false;
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
                    }
                }
            }

            // Check if loop is flippable
            if (fliploop && uni_dist(rng) < 0.5) {
                // Loop is flippable
                
                // Iterate over all spin ops and change type
                for (int p = 0; p < effexporder; p++) {
                    // Check if spin op is part of loop
                    if (spinops_loop[p]) {
                        // Flip spin op
                        opstring[p][0] = (opstring[p][0] == 1) ? 2 : 1;

                        // Set spin as part of flipped cluster
                        sp = opstring[p][2];
                        free_spins[sp] = 2;
                    }
                }
            } else {
				// Loop isn't flipped
				// ...
				// Nothing happens here, so this is redundant
            }

            // Empty list of bond and spin ops in loop
            for (int i = 0; i < effexporder; i++) {
                bondops_loop[i] = false;
                spinops_loop[i] = false;
            }
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
			spins[i] *= -1;
		}

		// Reset free-ness of spin
		free_spins[i] = 0;
	}

	return 0;

}

int adjust_maxexporder(int opstring[][3], int effexporder) {

	// Initialise random number generator
	int rngseed = 178;
	mt19937 rng;
	rng.seed(time(NULL)+100000*rngseed);
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
	double insert_prob = 0.33;
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

int measure_observables(int spins[], int nspins, int bonds[][3], int nbonds,
	int opstring[][3], int effexporder, int obs_count, int obs_exporder[], 
	double obs_magn[], double obs_magn_sq[], double obs_magn_quad[],
	double obs_nn_corr[]) {

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

	// ----- Measure magnetisation -----
	
	// Initialise variables
	int spin = 0;
	double avg_magn = 0;

	// Iterate over all spins
	for (int i = 0; i < nspins; i++) {
		// Initial spin state
		spin = spins[i];
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
		avg_magn += spin;
	}


	// Get average and add to list
	avg_magn = avg_magn/(double)(nspins); //*effexporder);
	obs_magn[obs_count] = abs(avg_magn);
	
	// ----- Find nearest neighbour spin-spin correlations -----

	// Initialise variables
	int s1i, s1, s2i, s2, total_corr = 0;
	int nn_corrs[nspins] = {0};
	double avg_corr;

	// Iterate over all bonds
	for (int b = 0; b < nbonds; b++) {
		// Get spins belonging to bond b
		s1i = bonds[b][0];
		s2i = bonds[b][1];

		// Get spin states
		s1 = spins[s1i];
		s2 = spins[s2i];

		// Add product to list
		nn_corrs[s1i] = s1*s2;
		nn_corrs[s2i] = s1*s2;
	}

	// Now iterate over all spins and get average correlation
	total_corr = 0;
	for (int s = 0; s < nspins; s++) {
		total_corr += nn_corrs[s];
	}
	// Normalise correlation and add to list
	avg_corr = total_corr/(double)nspins;
	obs_nn_corr[obs_count] = avg_corr;

	// Need to implement more observables
	// ...

	return 0;
	
}