def read_xyz_file(file_path):
    """
    Reads an XYZ file and returns atom coordinates as a NumPy array.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    num_atoms = int(lines[0])
    atom_data = []
    
    for line in lines[2:2+num_atoms]:
        Axyz = line.split()
        atom_data.append([elements.index(Axyz[0]), float(Axyz[1]), float(Axyz[2]), float(Axyz[3])])
    
    return np.array(atom_data)
