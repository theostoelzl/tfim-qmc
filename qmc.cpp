#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <random>
#include <string>
#include <vector>

// TO IMPLEMENT:
// 1 - flipping of single spins after flip operator updates
//     and flipping all spins pairwise who have bond ops in between
// 2 - flipping of cluster spins (equiv to Fig 59, Sandvik) ?
// 3 - flip updates of clusters of 2 spin ops sandwiching
//     bond ops (1 or several)
// 4 - expansion cut off
// 5 - input file

const int nx = 16;
const int ny = 16;
// !! SQUARE LATTICE
const int nspins = nx*ny;
// !! SQUARE LATTICE
const int nbonds = 2*nx*ny;
const double temp = 1.0;
const double hfield = 1.0;
const double jcoupling = 1.0;

const int maxexporder = 1000;

const int eqsweeps = 10000;
const int avsweeps = 100000;

int random_conf(int lattice[nx][ny], int nx, int ny);
int get_bonds(int lattice[nx][ny], int nx, int ny, int lattice_bonds[nbonds][2][2], int nbonds);
int diagonal_updates(int spins[nspins], int nx, int ny, int bonds[nbonds][2],
	int nbonds, int opstring[maxexporder][3], int effexporder);
int flip_updates(int spins[nspins], int bonds[nbonds][2],
	int nbonds, int opstring[maxexporder][3], int effexporder);
int cluster_updates(int lattice[nx][ny], int nx, int ny, int opstring[maxexporder][3]);
int adjust_maxexporder(int opstring[maxexporder][3]);

int main() {
	
	// Initialise lattice
	int lattice[nx][ny];
	random_conf(lattice, nx, ny);
	int spins[nspins];
	int spins_map[nx][ny];
	int i = 0;
	for (int x = 0; x < nx; x++) {
		for (int y = 0; y < ny; y++) {
			spins[i] = lattice[x][y];
			spins_map[x][y] = i;
			i++;
		}
	}
	
	// Find all bonds and get number of bonds nbonds
	// For the simple square lattice nbonds = nx*ny
	//int nbonds = 2*nx*ny;
	int lattice_bonds[nbonds][2][2];
	get_bonds(lattice, nx, ny, lattice_bonds, nbonds);
	
	// Map bonds from 2D lattice to 1D string
	int bonds[nbonds][2];
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
	}
	std::cout << std::size(bonds) << "\n";
	
	// Initialise operator string
	int opstring[maxexporder][3];
	for (int i = 0; i < maxexporder; i++) {
		opstring[i][0] = 0;
		opstring[i][1] = -1;
		opstring[i][2] = -1;
	}
	
	double avgm = 0;
	for (int i = 0; i < nspins; i++) {
		std::cout << spins[i] << "\t";
		avgm += spins[i];
	}
	avgm = avgm / nspins;
	std::cout << avgm << "\n";
	
	int effexporder = maxexporder;	
	
	// Start equilibration
	for (int i = 0; i < eqsweeps; i++) {
		
		// Diagonal updates to insert / remove operators
		diagonal_updates(spins, nx, ny, bonds, nbonds, opstring, effexporder);
		
		// Cluster updates to vary diagonal / off-diagonal ops
		//cluster_updates(lattice, nx, ny, opstring);
		// CAN ONLY FLIP TWO AT ONCE !
		flip_updates(spins, bonds, nbonds, opstring, effexporder);
		
		// Adjust largest expansion order
		//effexporder = adjust_maxexporder(opstring);
		//effexporder = maxexporder;
		
	}
	
	
	int nn[avsweeps/20];
	int nnn = 0;
	
	// Start averaging sweeps
	for (int i = 0; i < avsweeps; i++) {
		
		// Diagonal updates to insert / remove operators
		nnn = diagonal_updates(spins, nx, ny, bonds, nbonds, opstring, effexporder);
		
		if (i%20 == 0) {
			nn[i/20] = nnn;
		}
		
		// Cluster updates to vary diagonal / off-diagonal ops
		//cluster_updates(lattice, nx, ny, opstring);
		flip_updates(spins, bonds, nbonds, opstring, effexporder);
		
		// Measure observables
		//measure_observables(lattice, nx, ny, opstring);
		
	}
	
	for (int i = 0; i < nspins; i++) {
		std::cout << spins[i] << "\t";
	}
	std::cout << spins[i] << "\n";
	
	double avgn = 0;
	for (int n = 0; n < avsweeps/20; n++) {
		avgn += nn[n]/(double)(avsweeps/20);
	}
	
	std::cout << avgn;
	
	return 0;
	
}

int random_conf(int lattice[nx][ny], int nx, int ny) {
	
	// Initialise random number generator
	int rngseed = 178;
	std::mt19937 rng;
	rng.seed(time(NULL)+100000*rngseed);
	std::uniform_real_distribution<double> uni_dist(0,1);

	// Iterate through lattice and assign random spin directions 
	for (int x = 0; x < nx; x++) {
		for (int y = 0; y < ny; y++) {
			lattice[x][y] = (uni_dist(rng) < 0.5) ? -1 : 1;
		}
	}
	return 0;
	
}

int get_bonds(int lattice[nx][ny], int nx, int ny, int lattice_bonds[nbonds][2][2], int nbonds) {
	
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

int diagonal_updates(int spins[nspins], int nx, int ny, int bonds[nbonds][2], 
	int nbonds, int opstring[maxexporder][3], int effexporder) {

	// Initialise random number generator
	int rngseed = 178;
	std::mt19937 rng;
	rng.seed(time(NULL)+100000*rngseed);
	std::uniform_real_distribution<double> uni_dist(0,1);
	std::uniform_int_distribution<int> rand_x(0,nx);
	std::uniform_int_distribution<int> rand_y(0,ny);
	// !! Make sure to change rand_bond bound when using different lattice
	std::uniform_int_distribution<int> rand_spin(0,nspins-1);
	std::uniform_int_distribution<int> rand_bond(0,nbonds-1);

	// Evaluate current expansion order
	int exporder = maxexporder;
	for (int i = 0; i < maxexporder; i++) {
		if (opstring[i][0] == 0) {
			exporder = exporder - 1;
		}
	}
	//std::cout << "exporder" << exporder << "\n";

	// Iterate over all positions p in operator string
	for (int p = 0; p < maxexporder; p++) {
		int op_type = opstring[p][0];
		
		// Evaluate the type of operator present at p
		if (op_type == 0) {
			// Identity operator
			
			// Flip a coin and decide if to proceed with 
			// insertion of bond or spin operator
			double coin = uni_dist(rng);
			if (coin < 0.5) {
				// Try to insert a diagonal bond operator
				
				// Choose random bond to insert into
				int randb = rand_bond(rng);
				
				// Get corresponding spins
				int s1i = bonds[randb][0];
				int s1 = spins[s1i];
				int s2i = bonds[randb][1];
				int s2 = spins[s2i];		
				
				// Check if spins are anti-parallel
				if (s1*s2 < 0) {
					// Spins are anti-parallel so try to insert diagonal operator
					double paccept = (1/temp) * nbonds / (2*(maxexporder - exporder));
				
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
				int rands = rand_spin(rng);
				
				// Evaluate acceptance probability
				// !! Make sure to change acceptance prob when using different lattice
				double paccept = (1/temp) * nspins * hfield / (maxexporder - exporder);
				
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
				
				// Evaluate acceptance probability
				// !! Make sure to change acceptance prob when using different lattice
				double paccept = (maxexporder - exporder + 1) / ((1/temp) * nbonds);
				
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
				double paccept = (maxexporder - exporder + 1) / (2 * (1/temp) * nspins * hfield);
				
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
			int si = opstring[p][2];
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

int flip_updates(int spins[nspins], int bonds[nbonds][2],
	int nbonds, int opstring[maxexporder][3], int effexporder) {
	
	// Initialise random number generator
	int rngseed = 178;
	std::mt19937 rng;
	rng.seed(time(NULL)+100000*rngseed);
	std::uniform_real_distribution<double> uni_dist(0,1);
	std::uniform_int_distribution<int> rand_x(0,nx);
	std::uniform_int_distribution<int> rand_y(0,ny);
	// !! Make sure to change rand_bond bound when using different lattice
	std::uniform_int_distribution<int> rand_spin(0,nspins-1);
	std::uniform_int_distribution<int> rand_bond(0,nbonds-1);

	// Evaluate current expansion order
	int exporder = maxexporder;
	for (int i = 0; i < maxexporder; i++) {
		if (opstring[i][0] == 0) {
			exporder = exporder - 1;
		}
	}
	//std::cout << "exporder" << exporder << "\n";

	// List of spins indicating whether each is free (no ops)
	// 0 - free, 1 - not free, 2 - not free but flipped cluster
	int free_spins[nspins];
	for (int i = 0; i < nspins; i++) {
		free_spins[nspins] = 0;
	}	

	// Deal with diagonal -> off-diagonal flips first
	
	// List of whether we will just have flipped op or not
	bool just_flipped[maxexporder];
	for (int i = 0; i < maxexporder; i++) {
		just_flipped[i] = false;
	}
	
	// List for opstring index p at which diagonal spin op appears first
	int clusters_start[nspins];
	for (int n = 0; n < nspins; n++) {
		clusters_start[n] = -1;
	}
	// List for opstring index p at which diagonal bond op appears
	int clusters_bond[nspins];
	for (int n = 0; n < nspins; n++) {
		clusters_bond[n] = -1;
	}
	
	// Initialise variables
	int bd, sp, prevp, s1, s2;
	
	// Iterate over opstrings
	for (int p = 0; p < maxexporder; p++) {
		// Evaluate the type of operator present at p
		if (opstring[p][0] == 1) {
			// Diagonal operator
			
			// Check for error
			if (opstring[p][1] > -1 && opstring[p][2] > -1) {
				std::cout << "Something's wrong!" << "\n";
			}
			
			// Indices of bond / spin this op is acting on (one of which is -1)
			bd = opstring[p][1];
			sp = opstring[p][2];
			
			// Check if spin operator
			if (sp > -1) {
				// Operator is acting on spins
				
				// Set spin as not free
				free_spins[sp] = (free_spins[sp] == 0) ? 1 : free_spins[sp];

				// Check if there is another diagonal op acting on this spin
				if (clusters_start[sp] > -1) {
					// We already found a spin op acting on this spin
					
					// opstring index of earlier spin op acting on this spin
					prevp = clusters_start[sp];
					
					// Check if the two spin ops have a bond op in between
					if (clusters_bond[sp] > prevp && clusters_bond[sp] < p) {
						// If bond op found, set current p as new clusters_start
						
						clusters_start[sp] = p;
						clusters_bond[sp] = -1;
					} else {
						// If no bond op found, flip both spin ops with 1/2 probability
						if (uni_dist(rng) < 0.5) {
							opstring[prevp][0] = 2;
							opstring[p][0] = 2;
							
							// Save information which operators we just flipped
							just_flipped[prevp] = true;
							just_flipped[p] = true;

							// Set spin as part of flipped cluster
							free_spins[sp] = 2;
							
							// Reset stored indices
							clusters_start[sp] = -1;
							clusters_bond[sp] = -1;
						} else {
							// Reset stored indices
							clusters_start[sp] = p;
							clusters_bond[sp] = -1;
						}
						
					}
				} else {
					// This is the first spin op we've found
					clusters_start[sp] = p;
				}
			} else {
				// Operator is acting on bonds
				
				// Get spins belonging to this bond
				s1 = bonds[bd][0];
				s2 = bonds[bd][1];
				
				// We only care about bond ops if we've already found a
				// diagonal spin op before !
				// So check if spin op has been found for either spin and if yes
				// add bond op to list
				if (clusters_start[s1] > -1) {
					clusters_bond[s1] = p;
				}
				if (clusters_start[s2] > -1) {
					clusters_bond[s2] = p;
				}
			}
		}
	}
	
	// Then deal with off-diagonal -> diagonal

	// Empty list for opstring index p at which diagonal spin op appears first
	for (int n = 0; n < nspins; n++) {
		clusters_start[n] = -1;
	}
	// Empty lIst for opstring index p at which diagonal bond op appears
	for (int n = 0; n < nspins; n++) {
		clusters_bond[n] = -1;
	}
	
	// Iterate over opstrings
	for (int p = 0; p < maxexporder; p++) {
		// Only proceed if operator at p hasn't just been flipped
		if (!just_flipped[p]) {
			// Operator hasn't been flipped in this sweep
			
			// Evaluate the type of operator present at p
			if (opstring[p][0] == 2) {
				// Diagonal operator
				
				// Check for error
				if (opstring[p][1] > -1 && opstring[p][2] > -1) {
					std::cout << "Something's wrong!" << "\n";
				}
				
				// Indices of bond / spin this op is acting on (one of which is -1)
				bd = opstring[p][1];
				sp = opstring[p][2];
				
				// Check if spin operator
				if (sp > -1) {
					// Operator is acting on spins
					
					// Set spin as not free
					free_spins[sp] = (free_spins[sp] == 0) ? 1 : free_spins[sp];

					// Check if there is another off-diagonal op acting on this spin
					if (clusters_start[sp] > -1) {
						// We already found a spin op acting on this spin
						
						// opstring index of earlier spin op acting on this spin
						prevp = clusters_start[sp];
						
						// Check if the two spin ops have a bond op in between
						if (clusters_bond[sp] > prevp && clusters_bond[sp] < p) {
							// If bond op found, set current p as new clusters_start
							
							clusters_start[sp] = p;
							clusters_bond[sp] = -1;
						} else {
							// No bond op found in between

							// Now check if cluster has been flipped already
							// from diagonal to off-diagonal
							if (!just_flipped[prevp] && !just_flipped[p]) {
								// Cluster hasn't been flipped,
								// flip both spin ops with 1/2 probability
								if (uni_dist(rng) < 0.5) {
									opstring[prevp][0] = 1;
									opstring[p][0] = 1;

									// Set spin as part of flipped cluster
									free_spins[sp] = 2;
									
									// Reset stored indices
									clusters_start[sp] = -1;
									clusters_bond[sp] = -1;
								} else {
									// Reset stored indices
									clusters_start[sp] = p;
									clusters_bond[sp] = -1;
								}
							}
						}
					} else {
						// This is the first spin op we've found
						clusters_start[sp] = p;
					}
				} else {
					// Operator is acting on bonds
					
					// Get spins belonging to this bond
					s1 = bonds[bd][0];
					s2 = bonds[bd][1];
					
					// We only care about bond ops if we've already found a
					// diagonal spin op before !
					// So check if spin op has been found for either spin and if yes
					// add bond op to list
					if (clusters_start[s1] > -1) {
						clusters_bond[s1] = p;
					}
					if (clusters_start[s2] > -1) {
						clusters_bond[s2] = p;
					}
				}
			}
		}
	}

	// ---- Then check for clusters of two bonded spins with diagonal
	//      spin ops sandwiching a bond op -----

	// Empty list for opstring index p at which diagonal spin op appears first
	for (int n = 0; n < nspins; n++) {
		clusters_start[n] = -1;
	}
	// Empty list for opstring index p at which diagonal bond op appears
	for (int n = 0; n < nspins; n++) {
		clusters_bond[n] = -1;
	}
	// Empty list of just flipped operators
	for (int i = 0; i < maxexporder; i++) {
		just_flipped[i] = false;
	}
	// List for ends of clusters
	int clusters_end[nspins];
	for (int n = 0; n < nspins; n++) {
		clusters_end[n] = -1;
	}

	// Init variable
	int othersp = -1;

	// Iterate over opstrings
	for (int p = 0; p < maxexporder; p++) {
		// Evaluate the type of operator present at p
		if (opstring[p][0] == 1) {
			// Diagonal operator
			
			// Check for error
			if (opstring[p][1] > -1 && opstring[p][2] > -1) {
				std::cout << "Something's wrong!" << "\n";
			}
			
			// Indices of bond / spin this op is acting on (one of which is -1)
			bd = opstring[p][1];
			sp = opstring[p][2];
			
			// Check if spin operator
			if (sp > -1) {
				// Operator is acting on spins
				
				// Set spin as not free
				free_spins[sp] = (free_spins[sp] == 0) ? 1 : free_spins[sp];

				// Check if there is another diagonal op acting on this spin
				if (clusters_start[sp] > -1) {
					// We already found a spin op acting on this spin
					
					// opstring index of earlier spin op acting on this spin
					prevp = clusters_start[sp];
					
					// Check if the two spin ops have a bond op in between
					if (clusters_bond[sp] > prevp && clusters_bond[sp] < p) {
						
						// Bond op between, now check if other spin in bond
						// is part of a cluster
						
						othersp = (opstring[clusters_bond[sp]][0] == sp) 
										? opstring[clusters_bond[sp]][1] 
										: opstring[clusters_bond[sp]][0];

						// Has other spin an associated first spin op?
						if (clusters_start[othersp] > -1 && clusters_start[othersp] < p) {
							// Has other spin an associated second spin op?
							if (clusters_end[othersp] > clusters_bond[sp]) {
								// Cluster is complete, flip with 1/2 probability
								if (uni_dist(rng) < 0.5) {
									// Change op type to off-diagonal of all ops in cluster
									opstring[clusters_start[sp]][0] = 2;
									opstring[p][0] = 2;
									opstring[clusters_start[othersp]][0] = 2;
									opstring[clusters_end[othersp]][0] = 2;

									// Reset all counters
									clusters_start[sp] = -1;
									clusters_end[sp] = -1;
									clusters_bond[sp] = -1;
									clusters_start[othersp] = -1;
									clusters_start[othersp] = -1;
									clusters_bond[othersp] = -1;

									// Add ops to list of just flipped clusters
									just_flipped[clusters_start[sp]] = true;
									just_flipped[p] = true;
									just_flipped[clusters_start[othersp]] = true;
									just_flipped[clusters_end[othersp]] = true;

									// Set spins as part of flipped clusters
									free_spins[sp] = 2;
									free_spins[othersp] = 2;
								} else {
									// Cluster isn't flipped, so start new cluster from here
									clusters_start[sp] = p;
									clusters_end[sp] = -1;
									clusters_bond[sp] = -1;
									clusters_start[othersp] = clusters_end[othersp];
									clusters_end[othersp] = -1;
									clusters_bond[othersp] = -1;
								}
							} else {
								// Cluster is incomplete, so add current spin to cluster ends
								clusters_end[sp] = p;
							}
						}
					} else {
						// If no bond op found in between, reset counters
						clusters_start[sp] = p;
						clusters_end[sp] = -1;
						clusters_bond[sp] = -1;
					}
				} else {
					// This is the first spin op we've found
					clusters_start[sp] = p;
				}
			} else {
				// Operator is acting on bonds
				
				// Get spins belonging to this bond
				s1 = bonds[bd][0];
				s2 = bonds[bd][1];
				
				// We only care about bond ops if we've already found a
				// diagonal spin op before !
				// So check if spin op has been found for both spins and if yes
				// add bond op to list
				if (clusters_start[s1] > -1 && clusters_start[s2] > -1) {
					clusters_bond[s1] = p;
					clusters_bond[s2] = p;
				}
			}
		}
	}

	// ---- Then check for clusters of two bonded spins with off-diagonal
	//      spin ops sandwiching a bond op -----

	// Empty list for opstring index p at which diagonal spin op appears first
	for (int n = 0; n < nspins; n++) {
		clusters_start[n] = -1;
	}
	// Empty list for opstring index p at which diagonal bond op appears
	for (int n = 0; n < nspins; n++) {
		clusters_bond[n] = -1;
	}
	// Empty list for ends of clusters
	for (int n = 0; n < nspins; n++) {
		clusters_end[n] = -1;
	}

	// Reset variable
	othersp = -1;

	// Iterate over opstrings
	for (int p = 0; p < maxexporder; p++) {
		// Evaluate the type of operator present at p
		if (opstring[p][0] == 2) {
			// Diagonal operator
			
			// Check for error
			if (opstring[p][1] > -1 && opstring[p][2] > -1) {
				std::cout << "Something's wrong!" << "\n";
			}
			
			// Indices of bond / spin this op is acting on (one of which is -1)
			bd = opstring[p][1];
			sp = opstring[p][2];
			
			// Check if spin operator
			if (sp > -1) {
				// Operator is acting on spins
				
				// Set spin as not free
				free_spins[sp] = (free_spins[sp] == 0) ? 1 : free_spins[sp];

				// Check if there is another diagonal op acting on this spin
				if (clusters_start[sp] > -1) {
					// We already found a spin op acting on this spin
					
					// opstring index of earlier spin op acting on this spin
					prevp = clusters_start[sp];
					
					// Check if the two spin ops have a bond op in between
					if (clusters_bond[sp] > prevp && clusters_bond[sp] < p) {
						
						// Bond op between, now check if other spin in bond
						// is part of a cluster
						
						othersp = (opstring[clusters_bond[sp]][0] == sp) 
										? opstring[clusters_bond[sp]][1] 
										: opstring[clusters_bond[sp]][0];

						// Has other spin an associated first spin op?
						if (clusters_start[othersp] > -1 && clusters_start[othersp] < p) {
							// Has other spin an associated second spin op?
							if (clusters_end[othersp] > clusters_bond[sp]) {
								// Cluster is complete, now check if cluster 
								// hasn't just been flipped earlier
								if (!just_flipped[clusters_start[sp]] 
								&& !just_flipped[p]
								&& !just_flipped[clusters_start[othersp]] 
								&& !just_flipped[clusters_end[othersp]]) {
									// Cluster hasn't just been flipped,
									// so flip cluster with 1/2 probability
									if (uni_dist(rng) < 0.5) {
										// Change op type to diagonal of all ops in cluster
										opstring[clusters_start[sp]][0] = 1;
										opstring[p][0] = 1;
										opstring[clusters_start[othersp]][0] = 1;
										opstring[clusters_end[othersp]][0] = 1;

										// Reset all counters
										clusters_start[sp] = -1;
										clusters_end[sp] = -1;
										clusters_bond[sp] = -1;
										clusters_start[othersp] = -1;
										clusters_start[othersp] = -1;
										clusters_bond[othersp] = -1;

										// Set spins as part of flipped clusters
										free_spins[sp] = 2;
										free_spins[othersp] = 2;
									} else {
										// Cluster isn't flipped, so start new cluster from here
										clusters_start[sp] = p;
										clusters_end[sp] = -1;
										clusters_bond[sp] = -1;
										clusters_start[othersp] = clusters_end[othersp];
										clusters_end[othersp] = -1;
										clusters_bond[othersp] = -1;
									}
								}
							} else {
								// Cluster is incomplete, so add current spin to cluster ends
								clusters_end[sp] = p;
							}
						}
					} else {
						// If no bond op found in between, reset counters
						clusters_start[sp] = p;
						clusters_end[sp] = -1;
						clusters_bond[sp] = -1;
					}
				} else {
					// This is the first spin op we've found
					clusters_start[sp] = p;
				}
			} else {
				// Operator is acting on bonds
				
				// Get spins belonging to this bond
				s1 = bonds[bd][0];
				s2 = bonds[bd][1];
				
				// We only care about bond ops if we've already found a
				// diagonal spin op before !
				// So check if spin op has been found for both spins and if yes
				// add bond op to list
				if (clusters_start[s1] > -1 && clusters_start[s2] > -1) {
					clusters_bond[s1] = p;
					clusters_bond[s2] = p;
				}
			}
		}
	}
	
	// Now flip spins to change initial spin state

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

int cluster_updates(int lattice[nx][ny], int nx, int ny, int opstring[maxexporder][3]) {
	
	// To determine all clusters, go through opstring
	for (int p = 0; p < maxexporder; p++) {
		
	}
	
	return 0;
	
}

int adjust_maxexporder(int opstring[maxexporder][3]) {
		
	return 0;

}

int adjust_maxexporder(int opstring[maxexporder][3], int effexporder) {
	
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