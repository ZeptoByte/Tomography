
# Compute Length of Interaction for a Muon Given a object
import numpy as np

def spherical_to_cartesian(azimuth, zenith):
    x = np.sin(zenith) * np.cos(azimuth)
    y = np.sin(zenith) * np.sin(azimuth)
    z = np.cos(zenith)
    return x, y, z

def compute_length(discriminant, p_x, p_y, p_z, dx, dy, dz, a, b):
    t1 = (-b+np.sqrt(discriminant))/(2*a)
    t2 = (-b-np.sqrt(discriminant))/(2*a)
    p1 = [p_x+t1*dx, p_y+t1*dy, p_z+t1*dz]
    p2 = [p_x+t2*dx, p_y+t2*dy, p_z+t2*dz]

    return np.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2)


def compute_interaction_length(c_z, r, p_x, p_y, p_z, zenith, azimuth):

    if r ==0:
        return c_z*0 + 0
    
    dx, dy, dz = spherical_to_cartesian(azimuth, zenith)


    a = dx**2 + dy**2 + dz**2
    b = 2*(dx*p_x + dy*p_y + dz*(p_z - c_z))
    c = p_x**2 + p_y**2 + (p_z -c_z)**2 - r**2
    
    # Solve the quadratic equation
    discriminant = b**2 - 4*a*c

    # Return Length that the muon travels in the sphere
    # Also Need to Return the Intersection Point that is at the lowest Z postion
    return np.where(discriminant>0, compute_length(discriminant, p_x, p_y, p_z, dx, dy, dz, a, b),0) 


def compute_point_ex(discriminant, p_x, p_y, p_z, dx, dy, dz, a, b):

    t1 = (-b+np.sqrt(discriminant))/(2*a)
    t2 = (-b-np.sqrt(discriminant))/(2*a)
    p1 = [p_x+t1*dx, p_y+t1*dy, p_z+t1*dz]
    p2 = [p_x+t2*dx, p_y+t2*dy, p_z+t2*dz]

    z1 = p_z+t1*dz
    z2 = p_z+t2*dz

    return np.where(z1<z2, p1, p2)

# Find the exit point of the muon from the sphere
def compute_exit_point(c_z, r, p_x, p_y, p_z, zenith, azimuth):

    dx, dy, dz = spherical_to_cartesian(azimuth, zenith)


    a = dx**2 + dy**2 + dz**2
    b = 2*(dx*p_x + dy*p_y + dz*(p_z - c_z))
    c = p_x**2 + p_y**2 + (p_z -c_z)**2 - r**2
    
    # Solve the quadratic equation
    discriminant = b**2 - 4*a*c

    return compute_point_ex(discriminant, p_x, p_y, p_z, dx, dy, dz, a, b)

def compute_point_ent(discriminant, p_x, p_y, p_z, dx, dy, dz, a, b):

    t1 = (-b+np.sqrt(discriminant))/(2*a)
    t2 = (-b-np.sqrt(discriminant))/(2*a)
    p1 = [p_x+t1*dx, p_y+t1*dy, p_z+t1*dz]
    p2 = [p_x+t2*dx, p_y+t2*dy, p_z+t2*dz]

    z1 = p_z+t1*dz
    z2 = p_z+t2*dz

    return np.where(z1>z2, p1, p2)

def compute_entry_point(c_z, r, p_x, p_y, p_z, zenith, azimuth):

    dx, dy, dz = spherical_to_cartesian(azimuth, zenith)


    a = dx**2 + dy**2 + dz**2
    b = 2*(dx*p_x + dy*p_y + dz*(p_z - c_z))
    c = p_x**2 + p_y**2 + (p_z -c_z)**2 - r**2
    
    # Solve the quadratic equation
    discriminant = b**2 - 4*a*c

    return compute_point_ent(discriminant, p_x, p_y, p_z, dx, dy, dz, a, b)





# Compute the correct scattering angle for the muon based on the energy and other parametrs
def cartesian_to_spherical(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    return r, theta, phi


# From the second intersection point, see if the muon passes through the detector
def hits_in_detector(p_x, p_y, p_z, zenith, azimuth, z, size):

    dx, dy, dz = spherical_to_cartesian(azimuth, zenith)

    t = (z - p_z)/dz

    x = p_x + t*dx
    y = p_y + t*dy

    return np.where(((np.abs(x) <= size)& (np.abs(y) <= size)), True, False) 

def scatter_muon(zenith, azimuth, scat_angle):
    
    scat_azi = np.random.uniform(-np.pi,np.pi, len(zenith))

    a = np.cos(zenith)*np.cos(azimuth)
    b = -np.sin(azimuth)
    c = np.sin(zenith)*np.cos(azimuth)
    d = np.cos(zenith)*np.sin(azimuth)
    e = np.cos(azimuth)
    f = np.sin(zenith)*np.sin(azimuth)
    g = -np.sin(zenith)
    h = 0
    i = np.cos(zenith)
    rotate_x = np.sin(scat_angle)*np.cos(scat_azi)
    rotate_y = np.sin(scat_angle)*np.sin(scat_azi)
    rotate_z = np.cos(scat_angle)
    d_x = a*rotate_x+b*rotate_y+c*rotate_z
    d_y = d*rotate_x+e*rotate_y+f*rotate_z
    d_z = g*rotate_x+h*rotate_y+i*rotate_z

    # Direction Vector of Final Particle Direction 
    r, theta, phi = cartesian_to_spherical(d_x, d_y, d_z)
    return theta, phi

if __name__ == '__main__':
    pass