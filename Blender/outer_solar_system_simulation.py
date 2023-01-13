import numpy as np
#import bpy

# Gravitational constant when the mass of the sun is 1.
G = 2.95912208286e-4

# Planet names and order
planets = ('Sun','Jupiter','Saturn','Uranus','Neptune','Pluto')

# The data below is obtained from here: https://ssd.jpl.nasa.gov/horizons.cgi

# Masses relative to the sun (the increased sun mass is to account for the inner planets)
masses = np.array([1.00000597682, 
                   0.000954786104043, 
                   0.000285583733151, 
                   0.0000437273164546, 
                   0.0000517759138449, 
                   6.571141277023631e-09])

# Positions of the planets in astronomical units (au) on September 5, 1994, at 0h00 GST.
positions = np.array([[0., 0., 0.],
                    [-3.502576677887171E+00, -4.111754751605156E+00,  9.546986420486078E-02],
                    [9.075323064717326E+00, -3.443060859273154E+00, -3.008002285860299E-01],
                    [8.309900066449559E+00, -1.782348877489204E+01, -1.738826162402036E-01],
                    [1.147049510166812E+01, -2.790203169301273E+01,  3.102324955757055E-01],
                    [-1.553841709421204E+01, -2.440295115792555E+01,  7.105854443660053E+00]])

# Velocities of the planets relative to the sun in au/day on September 5, 1994, at 0h00 GST.
velocities = np.array([[0., 0., 0.],
                    [5.647185685991568E-03, -4.540768024044625E-03, -1.077097723549840E-04],
                    [1.677252496875353E-03,  5.205044578906008E-03, -1.577215019146763E-04],
                    [3.535508197097127E-03,  1.479452678720917E-03, -4.019422185567764E-05],
                    [2.882592399188369E-03,  1.211095412047072E-03, -9.118527716949448E-05],
                    [2.754640676017983E-03, -2.105690992946069E-03, -5.607958889969929E-04]])

# Compute total linear momentum
ptot = (masses[:,np.newaxis]*velocities).sum(axis=0)

# Recompute velocities relative to the center of mass
velocities -= ptot/masses.sum()

# Linear momenta of the planets: p = m*v
momenta = masses[:,np.newaxis]*velocities

# Function for Newtonian acceleration field
def acc(x, masses = masses, G = G):
    N = masses.shape[0]
    d = x.shape[-1]
    dx_pairs = x[:, np.newaxis] - x[np.newaxis, :]
    msq_pairs = masses[:, np.newaxis]*masses[np.newaxis, :]
    
    # Remove self-self interactions
    dx_pairs = np.delete(dx_pairs.reshape((N*N,d)),slice(None,None,N+1), axis = 0).reshape((N,N-1,d))
    msq_pairs = np.delete(msq_pairs.reshape((N*N)),slice(None,None,N+1), axis = 0).reshape((N,N-1))
    
    # Compute pairwise distances
    dist_pairs = np.sqrt((dx_pairs**2).sum(axis=-1))
    
    # Compute the gravitational force using Newton's law
    forces = -G*(dx_pairs*msq_pairs[:,:,np.newaxis]/dist_pairs[:,:,np.newaxis]**3).sum(axis=1)
    
    # Return accelerations
    return forces/masses[:,np.newaxis]

# Select time step and total integration time (measured in days)
h = 100 # Time stepsize in days
totaltime = 100*365 # Total simulation time in days

# Preallocate output vectors at each step
t_out = np.arange(0.,totaltime,h)
x_out = np.zeros(t_out.shape + positions.shape, dtype=float)
x_out[0,:,:] = positions
v_out = np.zeros_like(x_out)
v_out[0,:,:] = velocities

# Use Symplectic Euler method for integration
for x0, x1, v0, v1 in zip(x_out[:-1],x_out[1:],v_out[:-1],v_out[1:]):
    x1[:,:] = x0 + h*v0
    v1[:,:] = v0 + h*acc(x1)
    
for p in range(len(planets)):
    print(p)
"""
# -------------------------
# Add the Blender code here
sun = bpy.data.objects['Sun']
jupiter = bpy.data.objects['Jupiter']
saturn = bpy.data.objects['Saturn']
uranus = bpy.data.objects['Uranus']
neptune = bpy.data.objects['Neptune']
pluto = bpy.data.objects['Pluto']
myobj = np.array([sun, jupiter, saturn, uranus, neptune, pluto])
bpy.context.scene.frame_start = 0
bpy.context.scene.frame_end = 365
# loop of frames and insert keyframes at every frame
nlast = bpy.context.scene.frame_end

for n in range(nlast):
    # Set frame like this
    bpy.context.scene.frame_set(n)
        for obj in myobj:
            obj.location.x = x_out[n, myobj.index(obj), 0]
            obj.location.y = x_out[n, myobj.index(obj), 1]
            obj.location.z = x_out[n, myobj.index(obj), 2]
            # Insert new keyframe for "location" like this
            obj.keyframe_insert(data_path="location")
 """