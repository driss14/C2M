def normalize_vector(a):
    """ Returns the unit vector of a  """
    return a / np.linalg.norm(a)

def distance(a, b):
    return np.linalg.norm(b - a)
    
def angle(a, b, c):
    u1 = b - a
    u2 = c - b
    ang_rad = np.arccos(np.dot(u1, u2)/( np.linalg.norm(u1) * np.linalg.norm(u2)))
    return np.degrees(ang_rad)

def dihedral(a, b, c, d):
    u1 = b - a
    u2 = c - b
    u3 = d - c
    v1 = np.cross(u1, u2)
    v1 = v1 / np.linalg.norm(v1)
    v2 = np.cross(u2, u3)
    v2 = v2 / np.linalg.norm(v2)
    ang_rad = np.arccos(np.dot(v1, v2)/(np.linalg.norm(v1) * np.linalg.norm(v2)))
    sign =  np.sign(v1 * u3).sum()
    if sign != 0:
        ang_rad = ang_rad * sign
    return np.degrees(ang_rad)

def rotmat_from_vec(v1,v2):
    # Using Rodrigues' rotation formula
    # https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    # v=a×b
    # s=∥v∥  sin(o)
    # c=a⋅b  cos(o)
    v1[v1 == 0] = 1e-12
    
    # Normalize vectors
    v1 = normalize_vector(v1)
    v2 = normalize_vector(v2)
    
    # Calculation of dot and cross products
    v = np.cross(v1, v2) 
    s = np.linalg.norm(v)
    c = np.dot(v1, v2)
    if s == 0:
        rotation_matrix = np.eye(3) 
    else:
        v = v / s
    
    # R = I*c + K + K^2*(1-c)/(s^2),
    k = np.zeros((3, 3))
    k[0, 1] = -v[2]
    k[0, 2] =  v[1]
    k[1, 2] = -v[0]
    k[1, 0] = -k[0, 1]
    k[2, 0] = -k[0, 2]
    k[2, 1] = -k[1, 2]

    k2 = np.outer(v, v)

    rotation_matrix =  np.eye(3) * c + k * s + k2 / (1 + c)
    
    return rotation_matrix

#get rotation matrix from angles
def rotate(axis, angle_degrees):
    # Convert angle from degrees to radians
    angle_radians = np.radians(angle_degrees)

    # Normalize the axis vector
    axis = np.array(axis) / np.linalg.norm(axis)

    # Create the rotation matrix using the Rodrigues' rotation formula
    c = np.cos(angle_radians)
    s = np.sin(angle_radians)
    t = 1 - c
    x, y, z = axis
    rotation_matrix = np.array([
            [t*x*x + c, t*x*y - s*z, t*x*z + s*y],
            [t*x*y + s*z, t*y*y + c, t*y*z - s*x],
            [t*x*z - s*y, t*y*z + s*x, t*z*z + c]
        ])
    return rotation_matrix
