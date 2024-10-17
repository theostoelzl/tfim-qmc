import sys
import numpy as np

def main():
    # Lattice size
    nx = int(sys.argv[1])
    ny = int(sys.argv[2])

    # Simple square lattice
    lattice = 2*np.random.randint(0, 2, size=(nx, ny), dtype=int) - 1

    print(lattice)
    print(np.mean(lattice))
    
    # Flatten lattice to 1D list of spins
    spins = []
    spins_map = np.ndarray((nx, ny), dtype=int)
    spins_reverse_map = []
    i = 0
    for x in range(nx):
        for y in range(ny):
            spins.append(lattice[x,y])
            spins_map[x,y] = i
            spins_reverse_map.append([x,y])
            i += 1

    # List of bonds
    bonds = []

    # Bond index
    bi = 0

    # Find nearest neighbours
    for x in range(nx):
        for y in range(ny):
            # Get 1D index of this spin
            sp = spins_map[x,y]

            # Find nearest neighbours
            coupling = 1
            nn = []
            nn.append(spins_map[ (x+1)%nx, y ])
            nn.append(spins_map[ (x-1)%nx, y ])
            nn.append(spins_map[ x, (y-1)%ny ])
            nn.append(spins_map[ x, (y-1)%ny ])

            # Go through neighbours
            for n in nn:
                # Check if bond exists
                exists = False
                for b in bonds:
                    if not exists:
                        if (b[0] == sp and b[1] == n) or (b[0] == n and b[1] == sp):
                            exists = True
                
                if not exists:
                    bonds.append([sp, n, coupling])
    
    print(bonds)
    print( spins_reverse_map[248] )
    print(len(bonds))

    # Save bonds to file
    with open("bonds.txt", "w") as out_file:
        out_file.write(f"{len(spins)}\t{len(bonds)}")
        for i in range(len(bonds)):
            out_file.write(f"\n{i}\t{bonds[i][0]}\t{bonds[i][1]}\t{bonds[i][2]}")
                

if __name__ == "__main__":
    main()