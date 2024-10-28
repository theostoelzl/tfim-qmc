#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include "json.hpp"


// TO IMPLEMENT:
// 1 - flipping of single spins after flip operator updates
//     and flipping all spins pairwise who have bond ops in between (DONE)
// 2 - flipping of cluster spins (equiv to Fig 59, Sandvik) ?
// 3 - flip updates of clusters of 2 spin ops sandwiching
//     bond ops (1 or several)
// 4 - expansion cut off
// 5 - input file (DONE)
// 6 - J-coupling for individual bonds (DONE)



const int nx = 10;
const int ny = 10;

int random_conf(int spins[], int nspins);
int get_bonds(int lattice[nx][ny], int nx, int ny, int lattice_bonds[][2][2], int nbonds);
int diagonal_updates(int spins[], int nspins, int bonds[][3], int nbonds,
	int opstring[][3], int effexporder, double temp, double hfield);
int cluster_updates(int spins[], int nspins, int bonds[][3],
	int nbonds, int opstring[][3], int effexporder);
int adjust_maxexporder(int opstring[][3], int effexporder);

int main() {
	
	/*
	// Find all bonds and get number of bonds nbonds
	// For the simple square lattice nbonds = nx*ny
	//int nbonds = 2*nx*ny;
	int lattice_bonds[nbonds][2][2];
	get_bonds(lattice, nx, ny, lattice_bonds, nbonds);
	
	// Map bonds from 2D lattice to 1D string
	int bonds[nbonds][3];
	for (int b = 0; b < nbonds; b++) {
		// Get coords of sites connected by bonds
		int s1x = lattice_bonds[b][0][0];
		int s1y = lattice_bonds[b][0][1];
		int s2x = lattice_bonds[b][1][0];
		int s2y = lattice_bonds[b][1][1];
		
		// Find spin indices for both sites
		int s1i = spins_map[s1x][s1y];
		int s2i = spins_map[s2x][s2y];
		
		// Store 1D spin indices in fresh bonds array
		bonds[b][0] = s1i;
		bonds[b][1] = s2i;
		bonds[b][2] = 1;
	}
	std::cout << std::size(bonds) << "\n";
	*/

	// ----- Read in bonds from input file -----

	// Open input file containing bonds
	std::ifstream bonds_infile; 
	bonds_infile.open("bonds.txt");

	// Read in number of spins and bonds
	std::string in_nspins, in_nbonds;
	std::getline(bonds_infile, in_nspins, '\t');
	std::getline(bonds_infile, in_nbonds, '\n');
	std::cout << "spins: " << in_nspins << ", bonds: " << in_nbonds << "\n";

	// First, initialise number of spins and bonds
	int nspins = stoi(in_nspins);
	int nbonds = stoi(in_nbonds);

	// Then, initialise arrays for all spins and bonds
	int spins[nspins];
	int bonds[nbonds][3];

	// Read in bonds
	std::string myline[4];
	for (int l = 0; l < nbonds; l++) {
		for (int i = 0; i < 3; i++) {
			std::getline(bonds_infile, myline[i], '\t');
		}
		std::getline(bonds_infile, myline[3]);
		for (int i = 0; i < 4; i++) {
			std::cout << myline[i] << '\t';
		}
		std::cout << '\n';

		// Store spin indices and coupling constants in bonds array
		bonds[l][0] = stoi(myline[1]);
		bonds[l][1] = stoi(myline[2]);
		bonds[l][2] = stoi(myline[3]);
	}

	// ----- Read in simulation parameters -----

	std::ifstream f("setup.json");
	nlohmann::json data = nlohmann::json::parse(f);
	int eqsweeps = data["eqsweeps"];
	int avsweeps = data["avsweeps"];
	int maxexporder = data["max_expansion_order"];
	double temp = data["temperature"];
	double hfield = data["transverse_field"];
	// Change this !
	double jcoupling = 1.0;

	std::cout << data["eqsweeps"] << "\n";
	
	// ----- Randomise initial state -----
	random_conf(spins, stoi(in_nspins));

	double avgm = 0;
	for (int i = 0; i < nspins; i++) {
		std::cout << spins[i] << "\t";
		avgm += spins[i];
	}
	avgm = avgm / nspins;
	std::cout << avgm << "\n";

	// ----- Start Monte Carlo sweeps -----

	// Set initial truncated (effective) expansion order
	int effexporder = maxexporder;

	// Initialise operator string
	int opstring[maxexporder][3];
	for (int i = 0; i < maxexporder; i++) {
		opstring[i][0] = 0;
		opstring[i][1] = -1;
		opstring[i][2] = -1;
	}
	
	// Start equilibration
	for (int i = 0; i < eqsweeps; i++) {
		
		// Diagonal updates to insert / remove operators
		diagonal_updates(spins, nspins, bonds, nbonds, opstring, 
				effexporder, temp, hfield);
		
		// Cluster updates to vary diagonal / off-diagonal ops
		cluster_updates(spins, nspins, bonds, nbonds, opstring, effexporder);
		
		// Adjust largest expansion order
		//effexporder = adjust_maxexporder(opstring);
		//effexporder = maxexporder;
		
	}
	
	
	int nn[avsweeps/20];
	int nnn = 0;
	
	// Start averaging sweeps
	for (int i = 0; i < avsweeps; i++) {
		
		// Diagonal updates to insert / remove operators
		nnn = diagonal_updates(spins, nspins, bonds, nbonds, opstring, 
				effexporder, temp, hfield);
		
		if (i%20 == 0) {
			nn[i/20] = nnn;
		}
		
		// Cluster updates to vary diagonal / off-diagonal ops
		cluster_updates(spins, nspins, bonds, nbonds, opstring, effexporder);

		// Measure observables
		//measure_observables(lattice, nx, ny, opstring);
		
	}
	
	for (int i = 0; i < nspins; i++) {
		std::cout << spins[i] << "\t";
	}
	
	double avgn = 0;
	for (int n = 0; n < avsweeps/20; n++) {
		avgn += nn[n]/(double)(avsweeps/20);
	}
	
	std::cout << avgn;
	
	return 0;
	
}

int random_conf(int spins[], int nspins) {
	
	// Initialise random number generator
	int rngseed = 178;
	std::mt19937 rng;
	rng.seed(time(NULL)+100000*rngseed);
	std::uniform_real_distribution<double> uni_dist(0,1);

	// Iterate through lattice and assign random spin directions 
	for (int i = 0; i < nspins; i++) {
		spins[i] = (uni_dist(rng) < 0.5) ? -1 : 1;
	}
	return 0;
	
}

int get_bonds(int lattice[nx][ny], int nx, int ny, int lattice_bonds[][2][2], int nbonds) {
	
	// Empty bonds array
	for (int b = 0; b < nbonds; b++) {
		lattice_bonds[b][0][0] = -1;
		lattice_bonds[b][0][1] = -1;
		lattice_bonds[b][1][0] = -1;
		lattice_bonds[b][1][1] = -1;
	}
	
	// !! SIMPLE SQUARE LATTICE
	
	// Set up bond index for loop
	int bi = 0;	
	// Iterate over lattice
	for (int x = 0; x < nx; x++) {
		for (int y = 0; y < ny; y++) {
			// List of nearest neighbours
			int nn[4][2];
			nn[0][0] = (x+1 == nx) ? 0 : (x+1);
			nn[0][1] = y;
			nn[1][0] = (x-1 == -1) ? (nx-1) : (x-1)%(nx);
			nn[1][1] = y;
			nn[2][0] = x;
			nn[2][1] = (y+1 == ny) ? 0 : (y+1);
			nn[3][0] = x;
			nn[3][1] = (y-1 == -1) ? (ny-1) : (y-1)%(ny);
			
			// Go through nearest neighbours
			for (int nni = 0; nni < 4; nni++) {
				int nnx = nn[nni][0];
				int nny = nn[nni][1];
				
				// Check if bond exists
				bool exists = false;
				for (int b = 0; b < nbonds+1; b++) {
					if ( (lattice_bonds[b][0][0] == x && lattice_bonds[b][0][1] == y
						&& lattice_bonds[b][1][0] == nnx && lattice_bonds[b][1][1] == nny) ||
						(lattice_bonds[b][1][0] == x && lattice_bonds[b][1][1] == y
						&& lattice_bonds[b][0][0] == nnx && lattice_bonds[b][0][1] == nny) ) {
						// Bond exists already
						exists = true;
						break; 
					}
				}
				
				// If bond doesn't exist, save it
				if (!exists) {
					lattice_bonds[bi][0][0] = x;
					lattice_bonds[bi][0][1] = y;
					lattice_bonds[bi][1][0] = nnx;
					lattice_bonds[bi][1][1] = nny;
					bi++;
				}
			}
		}
	}
	
	return 0;
}

int diagonal_updates(int spins[], int nspins, int bonds[][3], int nbonds,
	int opstring[][3], int effexporder, double temp, double hfield) {

	// Initialise random number generator
	int rngseed = 178;
	std::mt19937 rng;
	rng.seed(time(NULL)+100000*rngseed);
	std::uniform_real_distribution<double> uni_dist(0,1);
	std::uniform_int_distribution<int> rand_spin(0,nspins-1);
	std::uniform_int_distribution<int> rand_bond(0,nbonds-1);

	// Evaluate current expansion order
	int exporder = effexporder;
	for (int i = 0; i < effexporder; i++) {
		if (opstring[i][0] == 0) {
			exporder = exporder - 1;
		}
	}
	//std::cout << "exporder" << exporder << "\n";

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
				
				// Check if spins are anti-parallel
				if (s1*s2 < 0) {
					// Spins are anti-parallel so try to insert diagonal operator
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
				// !! Make sure to change acceptance prob when using different lattice
				paccept = (1/temp) * nspins * hfield / (effexporder - exporder);
				
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
				std::cout << "Something's wrong!" << "\n";
			}
			
			// Check if bond or spin operator
			if (opstring[p][1] > -1) {
				// Operator is acting on bonds
				
				// Get coupling constant of bond
				bj = bonds[opstring[p][1]][2];

				// Evaluate acceptance probability
				// !! Make sure to change acceptance prob when using different lattice
				paccept = (effexporder - exporder + 1) / ((1/temp) * nbonds * bj);
				
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
				// !! Make sure to change acceptance prob when using different lattice
				paccept = (effexporder - exporder + 1) / (2 * (1/temp) * nspins * hfield);
				
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
			
			/*
			for (int i = 0; i < 9; i++) {
				std::cout << spins[i] << "\t";
			}
			std::cout << "\n";
			*/
		
		}
	}
	
	//std::cout<<exporder<<"\n";
	return exporder;
	
}

int cluster_updates(int spins[], int nspins, int bonds[][3],
	int nbonds, int opstring[][3], int effexporder) {
		
	int rngseed = 178;
	std::mt19937 rng;
	rng.seed(time(NULL)+100000*rngseed);
	std::uniform_real_distribution<double> uni_dist(0,1);

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
		linkops[i] = 0;
	}

	// Initialise variables
	int v0, v1, v2, s1i, s2i, bond, sp, si;

	// First, deal with links within bounds of opstring
	// To do this, iterate over opstring
	for (int p = 0; p < effexporder; p++) {
		// Check if operator at p is acting on bonds or spins
		if (opstring[p][0] == 0 && opstring[p][1] > -1) {
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
			// Get last vertex indices of each spin
			v1 = lasts[s1i];
			v2 = lasts[s2i];
			// Check if spins have appeared in vertex before
			if (v1 > -1) {
				// First spin has appeared before
				links[v1] = v0;
				links[v0] = v1;
			} else {
				// Set first appearance of spin as current p
				firsts[s1i] = v0;
			}
			if (v2 > -1) {
				// Second spin has appeared before
				links[v2] = v0+1; // ?? are we sure this is supposed to be v0 and not v0+1?
				links[v0+1] = v2;
			} else {
				// Set first appearance of spin as current p
				firsts[s2i] = v0+1;
			}
			// Update latest vertex appearance
			lasts[s1i] = v0+2;
			lasts[s2i] = v0+3;
		} else if (opstring[p][0] == 0 && opstring[p][2] > -1) {
			// Operator is acting on spins

			// Set vertex index
			v0 = 4*p;
			// Get bond index the op is acting on
			sp = opstring[p][2];
			// Store op type
			linkops[v0] = 2;
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
			last = firsts[i];
			links[first] = last;
			links[last] = first;
		}
	}

	// ----- Then trace all vertex loops and do flip updates -----

	// Initialise variables
	bool finished, fliploop, allspinops, flipops[effexporder], visited[4*effexporder];
	for (int i = 0; i < 4*effexporder; i++) {
		visited[i] = false;
	}
	for (int i = 0; i < effexporder; i++) {
		flipops[i] = false;
	}
	int currentv, linkedv, lastv, i, p, free_spins[nspins];
	for (int i = 0; i < nspins; i++) {
		free_spins[i] = 0;
	}
	double coin;

	// Iterate over all vertex legs
	for (int v = 0; v < 4*effexporder; v++) {
		// Check if vertex leg is connected to any other legs
		// !! should also check if leg has been visited before
		if (links[v] > -1 && !visited[v]) {
			finished = false;
			fliploop = true;
			currentv = v;
			while (!finished) {
				visited[currentv] = true;

				// Only for spin ops
				if (linkops[currentv] == 2) {
					// Get the operator p this vertex leg belongs to
					i = 0;
					while ( (currentv-i) % 4 != 0) {
						p = currentv/4;
						i++;
					}
					// Add the op this vertex leg belongs to as "to be flipped"
					flipops[p] = true;
				}

				// Check type of operator this leg is connected to
				if (linkops[currentv] == 1) {
					// Bond op

					// Check if leg is linked to bond op
					linkedv = links[currentv];
					if (linkedv > -1) {
						if (linkops[linkedv] == 1) {
							// Leg is linked to a leg on bond op
							
							// Change to linked leg on other bond op
							currentv = linkedv;
						} else if (linkops[linkedv] == 2) {
							// Leg is linked to a leg on spin op
							// Go to next leg on this bond op (or return to first leg)
							currentv += ( (currentv+1)%4 == 0) ? -3 : 1;
						} else {
							// Leg isn't linked to anything

							// Abandon this loop cluster
							finished = true;
							fliploop = false;
						}
					} else {
						// Leg isn't linked to anything

						// Abandon this loop cluster
						finished = true;
						fliploop = false;
					}
				} else if (linkops[currentv] == 2) {
					// Spin op
					lastv = currentv;
					currentv = links[currentv];
				}

				// Finish loop if we have reached original leg is reached again
				if (currentv == v) {
					finished = true;
				}
			}

			// Check if we should flip this loop
			coin = uni_dist(rng);
			if (fliploop && coin < 0.5) {
				// Loop is to be flipped

				// Iterate over all ops and see if they are part of the loop
				for (int p = 0; p < effexporder; p++) {
					if (flipops[p]) {
						// Operator is part of loop, so flip
						opstring[p][0] = (opstring[p][0] == 1) ? 2 : 1;
						// Get spins belonging to operator that is being flipped
						free_spins[opstring[p][2]] = 2;
					}
				}
			}

			// Reset list of ops to be flipped
			for (int i = 0; i < effexporder; i++) {
				flipops[i] = false;
			}
		}
	}

	// ----- Finally, take care of flipping spins in the periodic state -----

	// Iterate over spins i
	for (int i = 0; i < nspins; i++) {
		// Check if spin is free
		if (free_spins[i] == 0) {
			// Spin is free, so flip with 1/2 probability
			if (uni_dist(rng) < 0.5) {
				spins[i] *= -1;
			}
		} else if (free_spins[i] == 2) {
			// Spin isn't free but is part of a cluster that has been flipped
			spins[i] *= -1;
		}
	}

	return 0;

}

int adjust_maxexporder(int opstring[][3]) {
		
	return 0;

}

int adjust_maxexporder(int opstring[][3], int effexporder) {
	
	int newopstring[][3] = {};
	int newopcounter = 0;
	
	// Copy all non-unity operators into new opstring
	for (int p = 0; p < effexporder; p++) {
		if (opstring[p][0] > 0) {
			newopstring[newopcounter][0] = opstring[p][0];
			newopstring[newopcounter][1] = opstring[p][1];
			newopstring[newopcounter][2] = opstring[p][2];
			
		}
	}	
	
	return 0;
	
}