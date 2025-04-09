import sys
import numpy as np

def main():

    # Get type of lattice to generate
    lattype = sys.argv[1]

    # Generate lattice
    if lattype == "chain":
        # 1-D chain
        n = int(sys.argv[2])
        gen_lattice_chain(n)
    elif lattype == "inverse_square":
        n = int(sys.argv[2])
        gen_lattice_inv_square(n)
    elif lattype == "square":
        # 2-D simple square
        nx = int(sys.argv[2])
        ny = int(sys.argv[3])
        gen_lattice_square(nx, ny)
    elif lattype == "square_crossings":
        # 2-D square with diagonal crossings
        nx = int(sys.argv[2])
        ny = int(sys.argv[3])
        gen_lattice_sq_crossings(nx, ny)
    elif lattype == "shastry":
        # Shastry-Sutherland lattice
        nx = int(sys.argv[2])
        ny = int(sys.argv[3])
        gen_lattice_shastry(nx, ny)
    elif lattype == "shastry_j3":
        # Shastry-Sutherland lattice
        nx = int(sys.argv[2])
        ny = int(sys.argv[3])
        gen_lattice_shastry(nx, ny, j3=True)



def gen_lattice_chain(n):

    # Simple chain with periodic boundaries
    spins = 2*np.random.randint(0, 2, size=(n), dtype=int) - 1

    # List of bonds
    bonds = []

    # Bond index
    bi = 0

    # Find nearest neighbours
    for x in range(n):
        # Find nearest neighbours
        coupling = 1
        nn = []
        nn.append( (x+1)%n )
        nn.append( (x-1)%n )

        # Go through neighbours
        for ni in nn:
            # Check if bond exists
            exists = False
            for b in bonds:
                if not exists:
                    if (b[0] == x and b[1] == ni) or (b[0] == ni and b[1] == x):
                        exists = True
            
            if not exists:
                bonds.append([x, ni, coupling])

    print(bonds)
    #print( spins_reverse_map[248] )
    print(len(bonds))

    # Save bonds to file
    with open("bonds_chain.txt", "w") as out_file:
        out_file.write(f"{len(spins)}\t{len(bonds)}")
        for i in range(len(bonds)):
            out_file.write(f"\n{i}\t{bonds[i][0]}\t{bonds[i][1]}\t{bonds[i][2]}")


def gen_lattice_inv_square(n):

    # Simple chain with periodic boundaries
    spins = 2*np.random.randint(0, 2, size=(n), dtype=int) - 1

    # List of bonds
    bonds = []

    # Bond index
    bi = 0

    # Construct bonds among all spins
    # (like in Sandvik arXiv 2018 we're counting all pairs twice)
    for i in range(n):
        for j in range(n):
            if i != j:
                # Determine coupling
                coupling = 1/2*( 1/np.abs(i-j)**2 + 1/(n-np.abs(i-j))**2 )
                # Add bond to list
                bonds.append([ i, j, coupling])

    print(bonds)
    #print( spins_reverse_map[248] )
    print(len(bonds))

    # Save bonds to file
    with open("bonds_inverse_square.txt", "w") as out_file:
        out_file.write(f"{len(spins)}\t{len(bonds)}")
        for i in range(len(bonds)):
            out_file.write(f"\n{i}\t{bonds[i][0]}\t{bonds[i][1]}\t{bonds[i][2]}")
            

def gen_lattice_square(nx, ny):

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
            nn.append(spins_map[ x, (y+1)%ny ])
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
    #print( spins_reverse_map[248] )
    print(len(bonds))

    # Save bonds to file
    with open("bonds_square.txt", "w") as out_file:
        out_file.write(f"{len(spins)}\t{len(bonds)}")
        for i in range(len(bonds)):
            out_file.write(f"\n{i}\t{bonds[i][0]}\t{bonds[i][1]}\t{bonds[i][2]}")


def gen_lattice_sq_crossings(nx, ny):

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
            coupling = -1
            nn = []
            nn.append(spins_map[ (x+1)%nx, y ])
            nn.append(spins_map[ (x-1)%nx, y ])
            nn.append(spins_map[ x, (y+1)%ny ])
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

            # Find next-nearest neighbours
            coupling = -0.8
            nnn = []
            nnn.append(spins_map[ (x+1)%nx, (y+1)%ny ])
            nnn.append(spins_map[ (x+1)%nx, (y-1)%ny ])
            nnn.append(spins_map[ (x-1)%nx, (y+1)%ny ])
            nnn.append(spins_map[ (x-1)%nx, (y-1)%ny ])

            # Go through neighbours
            for n in nnn:
                # Check if bond exists
                exists = False
                for b in bonds:
                    if not exists:
                        if (b[0] == sp and b[1] == n) or (b[0] == n and b[1] == sp):
                            exists = True
                
                if not exists:
                    bonds.append([sp, n, coupling])
            
    
    print(bonds)
    #print( spins_reverse_map[248] )
    print(len(bonds))

    # Save bonds to file
    with open("bonds_square_crossings.txt", "w") as out_file:
        out_file.write(f"{len(spins)}\t{len(bonds)}")
        for i in range(len(bonds)):
            out_file.write(f"\n{i}\t{bonds[i][0]}\t{bonds[i][1]}\t{bonds[i][2]}")



def gen_lattice_shastry(nx, ny, j3 = 0):

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
            coupling = -1
            nn = []
            nn.append(spins_map[ (x+1)%nx, y ])
            nn.append(spins_map[ (x-1)%nx, y ])
            nn.append(spins_map[ x, (y+1)%ny ])
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

            # Find dimer bonds
            coupling_even = -7.0769
            nnn = []

            # Alternate dimer bond direction for each row
            if x%2 == 0:
                # To the left
                if y%2 != 0:
                    # Every other spin bonds to the left
                    nnn.append(spins_map[ (x+1)%nx, (y-1)%ny ])
            

            # Go through neighbours
            for n in nnn:
                # Check if bond exists
                exists = False
                for b in bonds:
                    if not exists:
                        if (b[0] == sp and b[1] == n) or (b[0] == n and b[1] == sp):
                            exists = True
                
                if not exists:
                    bonds.append([sp, n, coupling_even])

            # Other dimers
            nnn = []
            coupling_odd = -5.0769
            if x%2 != 0:
                # To the right
                if y%2 != 0:
                    # Every other spin bonds to the left
                    nnn.append(spins_map[ (x+1)%nx, (y+1)%ny ])

            # Go through neighbours
            for n in nnn:
                # Check if bond exists
                exists = False
                for b in bonds:
                    if not exists:
                        if (b[0] == sp and b[1] == n) or (b[0] == n and b[1] == sp):
                            exists = True
                
                if not exists:
                    bonds.append([sp, n, coupling_odd])

            # Add J3 couplings if wanted
            if j3:
                # Find dimer bonds
                coupling = -0.461538
                nnn = []

                # Alternate dimer bond direction for each row
                if x%2 == 0:
                    # To the left
                    if y%2 == 0:
                        # Every even index spin bonds to the bottom left
                        nnn.append(spins_map[ (x+1)%nx, (y-1)%ny ])
                    else:
                        # Every odd index spin bonds to the bottom right
                        nnn.append(spins_map[ (x+1)%nx, (y+1)%ny ])
                else:
                    # To the left
                    if y%2 != 0:
                        # Every odd index spin bonds to the bottom left
                        nnn.append(spins_map[ (x+1)%nx, (y-1)%ny ])
                    else:
                        # Every even index spin bonds to the bottom right
                        nnn.append(spins_map[ (x+1)%nx, (y+1)%ny ])

                # Go through neighbours
                for n in nnn:
                    # Check if bond exists
                    exists = False
                    for b in bonds:
                        if not exists:
                            if (b[0] == sp and b[1] == n) or (b[0] == n and b[1] == sp):
                                exists = True
                    
                    if not exists:
                        bonds.append([sp, n, coupling])
    
    print(bonds)
    #print( spins_reverse_map[248] )
    print(len(bonds))

    # Save bonds to file
    if j3:
        filename = "bonds_shastry_j3.txt"
    else:
        filename = "bonds_shastry.txt"
    with open(filename, "w") as out_file:
        out_file.write(f"{len(spins)}\t{len(bonds)}")
        for i in range(len(bonds)):
            out_file.write(f"\n{i}\t{bonds[i][0]}\t{bonds[i][1]}\t{bonds[i][2]}")




if __name__ == "__main__":
    main()
