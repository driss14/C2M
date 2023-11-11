class Atom:
    def __init__(self, aid, symbol, x, y, z):
        self.aid           = aid
        self.symbol        = symbol
        self.coordinates   = np.array([x, y, z])
        self.atomic_number = self.get_atomic_number()
        self.atomic_radii  = self.get_atomic_radii()
        self.covalent_radii  = self.get_covalent_radii()
        self.atomic_mass   = self.get_atomic_mass()
        self.neighbors = []  # Initialize an empty list for storing neighbors
        
    def get_atomic_number(self):
        return symbols.index(self.symbol)

    def get_atomic_radii(self):
        return radii[self.atomic_number]

    def get_covalent_radii(self):
        return covalent_radii[self.atomic_number]
    
    def get_atomic_mass(self):
        return masses[self.atomic_number]
    
    def calculate_bond_lengths(self, other):
        return distance(self.coordinates, other.coordinates)
    
    def calculate_angles(self, Atom1, Atom2):
        return angle(self.coordinates, Atom1.coordinates, Atom2.coordinates)
    
    def calculate_dihedral(self, Atom1, Atom2, Atom3):
        return dihedral(self.coordinates, Atom1.coordinates, Atom2.coordinates, Atom3.coordinates)
    
    def is_neighbor(self, other):
        bond_length = self.calculate_bond_lengths(other)
        return bond_length <= (self.covalent_radii + other.covalent_radii)
    
    def add_neighbor(self, other):
        if self.is_neighbor(other):
            self.neighbors.append(other.aid)

    def __str__(self):
        return f"{self.symbol} {self.coordinates[0]:.8f} {self.coordinates[1]:.8f} {self.coordinates[2]:.8f}"
