import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 1.0                     # Box length
N = 5                     # Number of points in space (x, y)
h = L / (N - 1)             # Step size
x = np.linspace(0, L, N)    # Mesh in x
y = np.linspace(0, L, N)    # Mesh in y
U = 1.0                     # Plate velocity
Re = 10.0                   # Reynolds number
alpha = 1.5                 # Relaxation parameter (SOR)
tn = 0                      # Iteration
psi = np.zeros((N, N))      # Streamfunction variable
omega = np.zeros((N, N))    # Vorticity variable

# Convergence criteria
l1_norm_omega = 1.0
l1_norm_psi = 1.0
l1_target = 1e-6

# Color map type
colorMap = plt.cm.jet

def L1norm ( new , old ) :
    norm = np.sum(np.abs(new - old))
    return norm


# Temporal loop
while l1_norm_omega > l1_target or l1_norm_psi > l1_target:
    # Calculate psi at internal grid points (BCs for psi are automatically enforced)
    psiold = psi.copy()
    for i in range(1, N - 1):
        for j in range(1, N - 1):
            psi[i, j] = 0.25 * ((psi[i + 1, j] + psi[i - 1, j] + psi[i, j + 1] + psi[i, j - 1]) + (h ** 2) * omega[i, j])
    
    # Calculate Vorticity at the wall; omit corner points due to singularity
    omega[1:-1, 0] = (-2. / (h ** 2)) * psi[1:-1, 1]  # bottom
    omega[1:-1, -1] = (-2. / (h ** 2)) * psi[1:-1, -2] - 2. * U / h  # top
    omega[0, 1:-1] = (-2. / (h ** 2)) * psi[1, 1:-1]  # left
    omega[-1, 1:-1] = (-2. / (h ** 2)) * psi[-2, 1:-1]  # right
    
    # Calculate Vorticity at internal grid points
    omegaold = omega.copy()
    for i in range(1, N - 1):
        for j in range(1, N - 1):
            omega[i, j] = alpha * ((omega[i - 1, j] + omega[i + 1, j] + omega[i, j - 1] + omega[i, j + 1]
             - 0.25 * Re * (psi[i, j + 1] - psi[i, j - 1]) * (omega[i + 1, j] - omega[i - 1, j])
             + 0.25 * Re * (psi[i + 1, j] - psi[i - 1, j]) * (omega[i, j + 1] - omega[i, j - 1])) / 4) + (1 - alpha) * omega[i, j]

    l1_norm_psi = L1norm(psi, psiold)
    l1_norm_omega = L1norm(omega, omegaold)
    tn = tn + 1

print("Iteration = ", tn)
print( "l1_norm_psi = " , l1_norm_psi)
print( "l1_norm_omega = " , l1_norm_omega)



# Plot stream-function
fig, ax = plt.subplots(figsize=(5.5, 5))
plt.contourf(x, y, psi.transpose(), 20, cmap=colorMap, vmin=-0.08, vmax=-0.01)
plt.title('Stream function for Re = ' + str(Re), fontsize=14)
plt.xlabel('x', fontsize=14)
plt.ylabel('y', fontsize=14)
plt.subplots_adjust(left=0.18)
plt.subplots_adjust(bottom=0.18)
plt.show()

# Plot vorticity
fig, ax = plt.subplots(figsize=(5.5, 5))
plt.contourf(x, y, omega.transpose(), 20, cmap=colorMap, vmin=-100, vmax=20)
plt.title('Vorticity for Re = ' + str(Re), fontsize=14)
plt.xlabel('x', fontsize=14)
plt.ylabel('y', fontsize=14)
plt.subplots_adjust(left=0.18)
plt.subplots_adjust(bottom=0.18)
plt.show()


# Compute velocities
u = np.zeros((N, N))
v = np.zeros((N, N))

for i in range(1, N - 1):
    for j in range(1, N - 1):
        u[i, j] = (psi[i, j + 1] - psi[i, j - 1]) / (2 * h)
        v[i, j] = -(psi[i + 1, j] - psi[i - 1, j]) / (2 * h)

u[:, -1] = U
u[:, 0] = 0
u[0, :] = 0
u[-1, :] = 0

v[:, -1] = 0
v[:, 0] = 0
v[0, :] = 0
v[-1, :] = 0

# Plot velocity field
fig, ax = plt.subplots(figsize=(5.5, 5))
plt.quiver(x, y, u.T, v.T)
plt.xlabel('x', fontsize=14)
plt.ylabel('y', fontsize=14)
plt.subplots_adjust(left=0.18)
plt.subplots_adjust(bottom=0.18)
plt.show()






