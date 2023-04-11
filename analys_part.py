from vpython import *
import random
import math
import json
import numpy as np
import matplotlib.pyplot as plt

# Define the range of t values for the simulation
t_min = 0
t_max = 0.2269
dt = 0.01
n = 1000

with open("param.json") as FileOfParam:
    Parametrs = json.load(FileOfParam)

Parametrs = Parametrs["Params1"]

# Define other variables
A = float(Parametrs[0]["A"]) # area
K = float(Parametrs[0]["K"])
b = float(Parametrs[0]["b"])
v = float(Parametrs[0]["v"])

# Define the particle's initial positions and velocities
x0 = [random.gauss(0, 2) for i in range(n)]
y0 = [random.gauss(0, 2) for i in range(n)]
z0 = [random.gauss(0, 5) for i in range(n)]
vx0 = [random.gauss(0, 1) for i in range(n)]
vy0 = [random.gauss(0, 1) for i in range(n)]
vz0 = [random.gauss(0, 1) for i in range(n)]

# Define the particle's position functions
def x(t):
    return ( A/2 * (math.sin((2 * K - 1) * b(t)) / (2 * K - 1) - math.sin((2 * K + 1) * b(t)) / (2 * K + 1)) )

def y(t):
    return ( A/2 * (math.cos((2 * K - 1) * b(t)) / (2 * K - 1) - math.cos( (2 * K + 1) * b(t) ) / (2 * K + 1) + 2 / (1 - 4 * K)) )

def z(t):
    return ( A * (math.sin(b(t)) - 1) )

def b(t):
    return ( v * t / A)

# Create a list of sphere objects to represent the particles
particles = [sphere(pos=vector(x0[i], y0[i], z0[i]), radius=0.1, color=color.red) for i in range(n)]

# Main simulation loop
t = t_min
while t < t_max:
    for i in range(n):
        # Update the particle's position based on the formulas x(t), y(t), and z(t)
        particles[i].pos = vector(x(t), y(t), z(t)) + vector(x0[i], y0[i], z0[i])
        
        # Update the particle's velocity based on the time step and the particle's acceleration
        vx = (x(t + dt) - x(t)) / dt
        vy = (y(t + dt) - y(t)) / dt
        vz = (z(t + dt) - z(t)) / dt
        particles[i].velocity = vector(vx, vy, vz) + vector(vx0[i], vy0[i], vz0[i])

    # Increment the time step
    t += dt

    # Wait for a short time to slow down the simulation
    rate(100)


# Set the particle positions and velocities at time t
for i in range(n):
    particles[i].pos = vector(x(t_min), y(t_min), z(t_min)) + vector(x0[i], y0[i], z0[i])
        
    # Update the particle's velocity based on the time step and the particle's acceleration
    vx = (x(t_min + dt) - x(t_min)) / dt
    vy = (y(t_min + dt) - y(t_min)) / dt
    vz = (z(t_min + dt) - z(t_min)) / dt
    particles[i].velocity = vector(vx, vy, vz) + vector(vx0[i], vy0[i], vz0[i])
    
    # Get the particle positions and velocities
x_data_min = [p.pos.x for p in particles]
vx_data_min = [p.velocity.x for p in particles]
y_data_min = [p.pos.y for p in particles]
vy_data_min = [p.velocity.y for p in particles]
z_data_min = [p.pos.z for p in particles]
vz_data_min = [p.velocity.z for p in particles]
    
# Plot the V_x(x) and V_y(y) graphs
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 5))
ax1.scatter(x_data_min, vx_data_min, s=1)
ax1.set_xlabel('x')
ax1.set_ylabel('V_x')
ax2.scatter(y_data_min, vy_data_min, s=1)
ax2.set_xlabel('y')
ax2.set_ylabel('V_y')
ax3.scatter(z_data_min, vz_data_min, s=1)
ax3.set_xlabel('z')
ax3.set_ylabel('V_z')
plt.show()

x_data_min = np.array(x_data_min)
vx_data_min = np.array(vx_data_min)
y_data_min = np.array(y_data_min)
vy_data_min = np.array(vy_data_min)
z_data_min = np.array(z_data_min)
vz_data_min = np.array(vz_data_min)
np.savez('data_min.npz', x_data_min=x_data_min, vx_data_min=vx_data_min,
         y_data_min=y_data_min, vy_data_min=vy_data_min, z_data_min=z_data_min,
         vz_data_min=vz_data_min)

print("Done MIN")

# Set the particle positions and velocities at time t
for i in range(n):
    particles[i].pos = vector(x(t_max), y(t_max), z(t_max)) + vector(x0[i], y0[i], z0[i])
        
    # Update the particle's velocity based on the time step and the particle's acceleration
    vx = (x(t_max + dt) - x(t_max)) / dt
    vy = (y(t_max + dt) - y(t_max)) / dt
    vz = (z(t_max + dt) - z(t_max)) / dt
    particles[i].velocity = vector(vx, vy, vz) + vector(vx0[i], vy0[i], vz0[i])
    
    # Get the particle positions and velocities
x_data_max = [p.pos.x for p in particles]
vx_data_max = [p.velocity.x for p in particles]
y_data_max = [p.pos.y for p in particles]
vy_data_max = [p.velocity.y for p in particles]
z_data_max = [p.pos.z for p in particles]
vz_data_max = [p.velocity.z for p in particles]
    
# Plot the V_x(x) and V_y(y) graphs
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 5))
ax1.scatter(x_data_max, vx_data_max, s=1)
ax1.set_xlabel('x')
ax1.set_ylabel('V_x')
ax2.scatter(y_data_max, vy_data_max, s=1)
ax2.set_xlabel('y')
ax2.set_ylabel('V_y')
ax3.scatter(z_data_max, vz_data_max, s=1)
ax3.set_xlabel('z')
ax3.set_ylabel('V_z')
plt.show()

x_data_max = np.array(x_data_max)
vx_data_max = np.array(vx_data_max)
y_data_max = np.array(y_data_max)
vy_data_max = np.array(vy_data_max)
z_data_max = np.array(z_data_max)
vz_data_max = np.array(vz_data_max)

np.savez('data_max.npz', x_data_max=x_data_max, vx_data_max=vx_data_max,
         y_data_max=y_data_max, vy_data_max=vy_data_max, z_data_max=z_data_max,
         vz_data_max=vz_data_max)

print("Done MAX")
