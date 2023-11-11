class Molecule:
    def __init__(self, xyz_file_path, empty = False):
        
        self.empty = empty 
        if self.empty == False:
            self.xyz_file_path = xyz_file_path
            self.atoms = self.read_xyz_file(xyz_file_path)
        else:
            self.atoms = []
    @staticmethod
    def read_xyz_file(file_path):
        atoms = []
        
        with open(file_path, 'r') as file:
            lines = file.readlines()
                
        #read number of atoms
        num_atoms = int(lines[0])
            
        #check if the file is a valid XYZ file
        if len(lines) < 2:
            raise ValueError("Invalid XYZ file format. It must contain at least 2 lines.")

        aid = 0
        for line in lines[2:num_atoms+2]:
            parts = line.split()
            aid  += 1
            symbol = parts[0]                
            coordinates = [float(parts[1]), float(parts[2]), float(parts[3])]
            atoms.append(Atom(aid, symbol, *coordinates))

        return atoms
        

    
    def get_atoms_formula(self):
        formula = ''
        atoms_symbol = [atom.symbol for atom in self.atoms]
        for sym in sorted(set(atoms_symbol)):
            count = atoms_symbol.count(sym)
            if count > 1:
                formula += sym + str(atoms_symbol.count(sym))
            else:
                formula += sym 
        return formula
    
    def get_atoms_number(self):
        return len(self.atoms)
    
    def get_total_mass(self):
        return sum([atom.get_atomic_mass() for atom in self.atoms])
    
    def get_number_electrons(self):
        return sum([atom.get_atomic_number() for atom in self.atoms])
    
    def get_center_mass(self):
        """
        Calculates the center of mass of the molecule.
        """
        center = np.zeros(3)
        for atom in self.atoms:
            center += atom.get_atomic_mass()*atom.coordinates
        if self.empty:
            return center
        else:
            return center/self.get_total_mass()
    
    def get_centroid(self):
        """
        Calculates the centroid of the molecule.
        """
        center = np.zeros(3)
        for atom in self.atoms:
            center += atom.coordinates
        if self.empty:
            return center
        else:
            return center/self.get_atoms_number()
    
    def get_moment_of_inertia(self):
        """
        Calculates the moment of inertia tensor of the molecule.
        """
        # Center of mass
        center = self.get_center_mass()
        # moment of inertia
        moment_of_inertia = np.zeros((3, 3))
        

        for atom in self.atoms:
            r = atom.coordinates - center

            # Update the moment of inertia tensor
            moment_of_inertia[0, 0] += atom.get_atomic_mass() * (r[1]**2 + r[2]**2) 
            moment_of_inertia[1, 1] += atom.get_atomic_mass() * (r[0]**2 + r[2]**2)
            moment_of_inertia[2, 2] += atom.get_atomic_mass() * (r[0]**2 + r[1]**2)
            moment_of_inertia[0, 1] -= atom.get_atomic_mass() * r[0] * r[1]
            moment_of_inertia[0, 2] -= atom.get_atomic_mass() * r[0] * r[2]
            moment_of_inertia[1, 2] -= atom.get_atomic_mass() * r[1] * r[2]

        # Fill in the symmetric elements of the tensor
        moment_of_inertia[1, 0] = moment_of_inertia[0, 1]
        moment_of_inertia[2, 0] = moment_of_inertia[0, 2]
        moment_of_inertia[2, 1] = moment_of_inertia[1, 2]

        return moment_of_inertia
    
    def get_principal_axis(self):
        """
        Calculates the principal axis (eigenvector) of the molecule's moment of inertia tensor.
        """
        eigenvalues, eigenvectors = np.linalg.eigh(self.get_moment_of_inertia())
        # Sort eigenverctors from highest to lowest
        idsorted = eigenvalues.argsort()
        principal_axis = eigenvectors[:, idsorted]

        return principal_axis
    
    def reorient(self):  
        """
        Center the molecule and align it with the Z-axis
        """
        center_mass = self.get_center_mass()
        self.translate(-center_mass[0], -center_mass[1], -center_mass[2])
        self.align_principal_axis_with_z_axis()
        
    def align_principal_axis_with_z_axis(self):
        """
        Align the molecule with the Z-axis
        """    
        principal_axis = self.get_principal_axis().T[0]
        
        rotation_matrix = rotmat_from_vec( principal_axis, np.array([0, 0, 1]))
        
        for atom in self.atoms:
            atom.coordinates = np.dot(atom.coordinates, rotation_matrix.T)
            
    def translate(self, dx, dy, dz):
        """
        Translates all atoms in the molecule by the specified amounts in the X, Y, and Z axes.
        """
        translation_vector = np.array([dx, dy, dz])
        for atom in self.atoms:
            atom.coordinates += translation_vector
    def rotate(self, axis, angle_degrees):
        # Convert angle from degrees to radians
        rotation_matrix = rotate(axis, angle_degrees)
        # Apply the rotation to all atoms in the molecule
        for atom in self.atoms:
            atom.coordinates = np.dot(rotation_matrix, atom.coordinates)
            
    def calculate_rmsd(self, other_molecule):
        """
        calculate the RMSD of two molecules with the same atomes and bonds.
        """
        if self.empty or other_molecule.empty:
            raise ValueError("One or both molecules are empty. Cannot calculate RMSD.")
        
        if len(self.atoms) != len(other_molecule.atoms):
            raise ValueError("Both molecules must have the same number of atoms for RMSD calculation.")

        squared_diff_sum = 0.0

        for atom1, atom2 in zip(self.atoms, other_molecule.atoms):
            diff = atom1.coordinates - atom2.coordinates
            squared_diff = np.sum(diff**2)
            squared_diff_sum += squared_diff

        rmsd = np.sqrt(squared_diff_sum / len(self.atoms))
        return rmsd
    
    def concatenate(self, other_molecule):
        if self.empty:
            raise ValueError("Cannot concatenate an empty molecule.")
        
        # Translate all atoms in Molecule 2
        for atom in other_molecule.atoms:
            self.atoms.append(atom)

    
    @staticmethod
    def calculate_bounding_box(molecule):
        if not molecule.atoms:
            raise ValueError("Cannot calculate bounding box for an empty molecule.")
        
        # Initialize min and max coordinates
        min_x, min_y, min_z = float('inf'), float('inf'), float('inf')
        max_x, max_y, max_z = float('-inf'), float('-inf'), float('-inf')
        
        # Iterate through the atoms to find min and max coordinates
        for atom in molecule.atoms:
            x, y, z = atom.coordinates
            min_x = min(min_x, x)
            min_y = min(min_y, y)
            min_z = min(min_z, z)
            max_x = max(max_x, x)
            max_y = max(max_y, y)
            max_z = max(max_z, z)
        
        # Calculate the dimensions of the bounding box
        width = max_x - min_x
        height = max_y - min_y
        depth = max_z - min_z
        
        # Calculate the center of the bounding box
        center_x = (max_x + min_x) / 2.0
        center_y = (max_y + min_y) / 2.0
        center_z = (max_z + min_z) / 2.0
        
        # Return the bounding box as a dictionary
        bounding_box = {
            "min_x": min_x,
            "min_y": min_y,
            "min_z": min_z,
            "max_x": max_x,
            "max_y": max_y,
            "max_z": max_z,
            "width": width,
            "height": height,
            "depth": depth,
            "center_x": center_x,
            "center_y": center_y,
            "center_z": center_z
        }
        
        return bounding_box

    def save_to_xyz(self, output_file_path):
        with open(output_file_path, 'w') as file:
            file.write(f"{len(self.atoms)}\n{self.get_atoms_formula()}\n")
            for atom in self.atoms:
                file.write(str(atom) + "\n")
    def __str__(self):
        return f"{self.get_atoms_formula()}"
