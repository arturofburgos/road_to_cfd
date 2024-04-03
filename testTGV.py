import math
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, sin, cos, acos, pi

# Mesh parameters
nx = 32
ny = 32
nz = 32
L = 2 * pi
dx = L / nx
dy = L / ny
dz = L / nz
x = np.arange(-pi, pi, dx)
y = np.arange(-pi, pi, dy)
z = np.arange(-pi, pi, dz)

# Physical parameters
Re = 100.  # Reynolds number
Pr = 0.723  # Prandtl number
gamma = 1.4  # Heat capacity ratio
rg = 287.06  # Gas constant
Chp = rg * gamma / (gamma - 1.)  # Heat capacity at constant pressure
Pref = 100.  # Reference pressure
Rref = 1.  # Reference density
Uref = 1.  # Reference velocity

# Temporal parameters
t = 0.
nmax = 20000
if nx == 32:
    dt = 0.0005
if nx == 64:
    dt = 0.00025

scheme_order = 10

if scheme_order == 2:
    # derivative and filter coeffs . o2
    a = [0., 1./2., 0., 0., 0., 0.]
    d = [0.5, -0.25, 0., 0., 0., 0.]
elif scheme_order == 4:
    # derivative and filter coeffs . o4
    a = [0., 2./3., -1./12., 0., 0., 0.]
    d = [0.375, -0.25, 0.0625, 0., 0., 0.]
elif scheme_order == 6:
    # derivative and filter coeffs . o6
    a = [0., 3./4., -3./20., 1./60., 0., 0.]
    d = [0.3125, -0.234375, 0.09375, -0.015625, 0., 0.]
elif scheme_order == 8:
    # derivative and filter coeffs . o8
    a = [0., 4./5., -1./5., 4./105., -1./280., 0.]
    d = [35./128., -7./32., 7./64., -1./32., 1./256., 0.]
elif scheme_order == 10:
    # derivative and filter coeffs . o10
    a = [0., 5./6, -5./21, 5./84, -5./504, 1./1260]
    d = [63./256, -105./512, 15./128, -45./1024, 5./512, -1./1024]
else:
    print("Error"); exit()

# Filter parameters
xi = 0.2
rho0 = Rref
rhou0 = 0
rhov0 = 0
rhow0 = 0
rhoe0 = Pref / (gamma - 1.)

def initialisation(Rref, Pref, Uref):
    u = np.zeros((nx, ny, nz))
    v = np.zeros((nx, ny, nz))
    p = np.zeros((nx, ny, nz))
    for i in range(0, nx):
        for j in range(0, ny):
            for k in range(0, nz):
                u[i, j, k] = Uref * sin(x[i]) * cos(y[j]) * cos(z[k])
                v[i, j, k] = -Uref * cos(x[i]) * sin(y[j]) * cos(z[k])
                p[i, j, k] = Pref + Rref * Uref ** 2 / 16. * (cos(2 * z[k]) + 2) * (cos(2 * x[i]) + cos(2 * y[j]))
    rho = np.full((nx, ny, nz), Rref)
    rhou = rho * u
    rhov = rho * v
    rhow = np.zeros((nx, ny, nz))
    rhoe = p / (gamma - 1.) + 0.5 * (rhou ** 2 + rhov ** 2 + rhow ** 2) / rho
    return rho, rhou, rhov, rhow, rhoe


def compute_cons2prim(rho, rhou, rhov, rhow, rhoe):
    u = rhou / rho
    v = rhov / rho
    w = rhow / rho
    p = (gamma - 1.) * (rhoe - 0.5 * rho * (u ** 2 + v ** 2 + w ** 2))
    return u, v, w, p

def calc_xderiv(U):
    Um1 = np.roll(U, 1, axis=0)
    Up1 = np.roll(U, -1, axis=0)
    Um2 = np.roll(U, 2, axis=0)
    Up2 = np.roll(U, -2, axis=0)
    Um3 = np.roll(U, 3, axis=0)
    Up3 = np.roll(U, -3, axis=0)
    Um4 = np.roll(U, 4, axis=0)
    Up4 = np.roll(U, -4, axis=0)
    Um5 = np.roll(U, 5, axis=0)
    Up5 = np.roll(U, -5, axis=0)
    Ux = (a[1] * (Up1 - Um1) + a[2] * (Up2 - Um2) + a[3] * (Up3 - Um3) + a[4] * (Up4 - Um4) + a[5] * (Up5 - Um5)) / dx
    return Ux

def calc_yderiv(U):
    Um1 = np.roll(U, 1, axis=1)
    Up1 = np.roll(U, -1, axis=1)
    Um2 = np.roll(U, 2, axis=1)
    Up2 = np.roll(U, -2, axis=1)
    Um3 = np.roll(U, 3, axis=1)
    Up3 = np.roll(U, -3, axis=1)
    Um4 = np.roll(U, 4, axis=1)
    Up4 = np.roll(U, -4, axis=1)
    Um5 = np.roll(U, 5, axis=1)
    Up5 = np.roll(U, -5, axis=1)
    Uy = (a[1] * (Up1 - Um1) + a[2] * (Up2 - Um2) + a[3] * (Up3 - Um3) + a[4] * (Up4 - Um4) + a[5] * (Up5 - Um5)) / dy
    return Uy

def calc_zderiv(U):
    Um1 = np.roll(U, 1, axis=2)
    Up1 = np.roll(U, -1, axis=2)
    Um2 = np.roll(U, 2, axis=2)
    Up2 = np.roll(U, -2, axis=2)
    Um3 = np.roll(U, 3, axis=2)
    Up3 = np.roll(U, -3, axis=2)
    Um4 = np.roll(U, 4, axis=2)
    Up4 = np.roll(U, -4, axis=2)
    Um5 = np.roll(U, 5, axis=2)
    Up5 = np.roll(U, -5, axis=2)
    Uz = (a[1] * (Up1 - Um1) + a[2] * (Up2 - Um2) + a[3] * (Up3 - Um3) + a[4] * (Up4 - Um4) + a[5] * (Up5 - Um5)) / dz
    return Uz

def calc_inviscid_fluxes(rho, rhou, rhov, rhow, rhoe):
    u = rhou / rho
    v = rhov / rho
    w = rhow / rho
    p = (gamma - 1.) * (rhoe - 0.5 * rho * (u ** 2 + v ** 2 + w ** 2))
    Frho = rhou
    Frhou = p + rhou * u
    Frhov = rhou * v
    Frhow = rhou * w
    Frhoe = (rhoe + p) * u
    Grho = rhov
    Grhou = rhov * u
    Grhov = p + rhov * v
    Grhow = rhov * w
    Grhoe = (rhoe + p) * v
    Hrho = rhow
    Hrhou = rhow * u
    Hrhov = rhow * v
    Hrhow = p + rhow * w
    Hrhoe = (rhoe + p) * w
    # Non-split form: unstable => filtering procedure to stabilize the scheme
    Krho = calc_xderiv(Frho) + calc_yderiv(Grho) + calc_zderiv(Hrho)
    Krhou = calc_xderiv(Frhou) + calc_yderiv(Grhou) + calc_zderiv(Hrhou)
    Krhov = calc_xderiv(Frhov) + calc_yderiv(Grhov) + calc_zderiv(Hrhov)
    Krhow = calc_xderiv(Frhow) + calc_yderiv(Grhow) + calc_zderiv(Hrhow)
    Krhoe = calc_xderiv(Frhoe) + calc_yderiv(Grhoe) + calc_zderiv(Hrhoe)
    return Krho, Krhou, Krhov, Krhow, Krhoe

def calc_viscous_fluxes(rho, rhou, rhov, rhow, rhoe):
    u = rhou / rho
    v = rhov / rho
    w = rhow / rho
    p = (gamma - 1.) * (rhoe - 0.5 * rho * (u ** 2 + v ** 2 + w ** 2))
    T = p / (rho * rg)
    visco = 1. / Re  # viscosity
    kappa = Chp * visco / Pr  # thermal coefficient k = mu * Cp / Pr
    S11 = calc_xderiv(u)
    S22 = calc_yderiv(v)
    S33 = calc_zderiv(w)
    S12 = 0.5 * (calc_xderiv(v) + calc_yderiv(u))
    S13 = 0.5 * (calc_xderiv(w) + calc_zderiv(u))
    S23 = 0.5 * (calc_yderiv(w) + calc_zderiv(v))
    dTdx = calc_xderiv(T)
    dTdy = calc_yderiv(T)
    dTdz = calc_zderiv(T)
    Frhou = - visco * 2. * (S11 - (S11 + S22 + S33) / 3.)
    Frhov = - visco * 2. * S12
    Frhow = - visco * 2. * S13
    Frhoe = - visco * 2. * (rhou * (S11 - (S11 + S22 + S33) / 3.) + rhov * S12 + rhow * S13) / rho - kappa * dTdx
    Grhou = - visco * 2. * S12
    Grhov = - visco * 2. * (S22 - (S11 + S22 + S33) / 3.)
    Grhow = - visco * 2. * S23
    Grhoe = - visco * 2. * (rhou * S12 + rhov * (S22 - (S11 + S22 + S33) / 3.) + rhow * S23) / rho - kappa * dTdy
    Hrhou = - visco * 2. * S13
    Hrhov = - visco * 2. * S23
    Hrhow = - visco * 2. * (S33 - (S11 + S22 + S33) / 3.)
    Hrhoe = - visco * 2. * (rhou * S13 + rhov * S23 + rhow * (S33 - (S11 + S22 + S33) / 3.)) / rho - kappa * dTdz
    Krhou = calc_xderiv(Frhou) + calc_yderiv(Grhou) + calc_zderiv(Hrhou)
    Krhov = calc_xderiv(Frhov) + calc_yderiv(Grhov) + calc_zderiv(Hrhov)
    Krhow = calc_xderiv(Frhow) + calc_yderiv(Grhow) + calc_zderiv(Hrhow)
    Krhoe = calc_xderiv(Frhoe) + calc_yderiv(Grhoe) + calc_zderiv(Hrhoe)
    return Krhou, Krhov, Krhow, Krhoe

def calc_filter(rho, rhou, rhov, rhow, rhoe):
    rhof = explicit_filtering(rho, rho0)
    rhouf = explicit_filtering(rhou, rhou0)
    rhovf = explicit_filtering(rhov, rhov0)
    rhowf = explicit_filtering(rhow, rhow0)
    rhoef = explicit_filtering(rhoe, rhoe0)
    return rhof, rhouf, rhovf, rhowf, rhoef

def explicit_filtering(U, U0):
    Uf = U - np.full((nx, ny, nz), U0)  # Separate ref field

    # Filtering in x
    Um1 = np.roll(Uf, 1, axis=0)
    Up1 = np.roll(Uf, -1, axis=0)
    Um2 = np.roll(Uf, 2, axis=0)
    Up2 = np.roll(Uf, -2, axis=0)
    Um3 = np.roll(Uf, 3, axis=0)
    Up3 = np.roll(Uf, -3, axis=0)
    Um4 = np.roll(Uf, 4, axis=0)
    Up4 = np.roll(Uf, -4, axis=0)
    Um5 = np.roll(Uf, 5, axis=0)
    Up5 = np.roll(Uf, -5, axis=0)
    Dx = d[0] * Uf + d[1] * (Up1 + Um1) + d[2] * (Up2 + Um2) + d[3] * (Up3 + Um3) + d[4] * (Up4 + Um4) + d[5] * (Up5 + Um5)

    # Filtering in y
    Um1 = np.roll(Uf, 1, axis=1)
    Up1 = np.roll(Uf, -1, axis=1)
    Um2 = np.roll(Uf, 2, axis=1)
    Up2 = np.roll(Uf, -2, axis=1)
    Um3 = np.roll(Uf, 3, axis=1)
    Up3 = np.roll(Uf, -3, axis=1)
    Um4 = np.roll(Uf, 4, axis=1)
    Up4 = np.roll(Uf, -4, axis=1)
    Um5 = np.roll(Uf, 5, axis=1)
    Up5 = np.roll(Uf, -5, axis=1)
    Dy = d[0] * Uf + d[1] * (Up1 + Um1) + d[2] * (Up2 + Um2) + d[3] * (Up3 + Um3) + d[4] * (Up4 + Um4) + d[5] * (Up5 + Um5)

    # Filtering in z
    Um1 = np.roll(Uf, 1, axis=2)
    Up1 = np.roll(Uf, -1, axis=2)
    Um2 = np.roll(Uf, 2, axis=2)
    Up2 = np.roll(Uf, -2, axis=2)
    Um3 = np.roll(Uf, 3, axis=2)
    Up3 = np.roll(Uf, -3, axis=2)
    Um4 = np.roll(Uf, 4, axis=2)
    Up4 = np.roll(Uf, -4, axis=2)
    Um5 = np.roll(Uf, 5, axis=2)
    Up5 = np.roll(Uf, -5, axis=2)
    Dz = d[0] * Uf + d[1] * (Up1 + Um1) + d[2] * (Up2 + Um2) + d[3] * (Up3 + Um3) + d[4] * (Up4 + Um4) + d[5] * (Up5 + Um5)

    Uf = U - xi * Dx - xi * Dy - xi * Dz
    return Uf

def plot_image(rho, rhou, rhov, rhow, rhoe):
    u, v, w, p = compute_cons2prim(rho, rhou, rhov, rhow, rhoe)
    
    # Subplot
    colormap = 'jet'  # 'RdBu'
    fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
    
    axes[0][0].set_title(r"$u$")
    im = axes[0][0].pcolor(x, y, u[:, :, 0].transpose(), shading='auto', cmap=colormap, vmin=-1, vmax=1)
    axes[0][0].set_ylabel(r"$y$")
    fig.colorbar(im, ax=axes[0][0])
    
    axes[0][1].set_title(r"$v$")
    im = axes[0][1].pcolor(x, y, v[:, :, 0].transpose(), shading='auto', cmap=colormap, vmin=-1, vmax=1)
    fig.colorbar(im, ax=axes[0][1])
    
    axes[1][0].set_title(r"$rhoe$")
    im = axes[1][0].pcolor(x, y, rhoe[:, :, 0].transpose(), shading='auto', cmap=colormap)
    axes[1][0].set_xlabel(r"$x$")
    axes[1][0].set_ylabel(r"$y$")
    fig.colorbar(im, ax=axes[1][0])
    
    axes[1][1].set_title(r"$p$")
    im = axes[1][1].pcolor(x, y, p[:, :, 0].transpose(), shading='auto', cmap=colormap)
    axes[1][1].set_xlabel(r"$x$")
    fig.colorbar(im, ax=axes[1][1])
    
    plt.draw()
    plt.pause(0.001)


# Initialization of conservative variables
rho, rhou, rhov, rhow, rhoe = initialisation(Rref, Pref, Uref)

# Temporal loop
file = open("./data.dat", "w")
for n in range(1, nmax):
    Krho, Krhou, Krhov, Krhow, Krhoe = calc_inviscid_fluxes(rho, rhou, rhov, rhow, rhoe)
    Lrhou, Lrhov, Lrhow, Lrhoe = calc_viscous_fluxes(rho, rhou, rhov, rhow, rhoe)
    
    rho = rho - dt * Krho
    rhou = rhou - dt * Krhou - dt * Lrhou
    rhov = rhov - dt * Krhov - dt * Lrhov
    rhow = rhow - dt * Krhow - dt * Lrhow
    rhoe = rhoe - dt * Krhoe - dt * Lrhoe
    
    rho, rhou, rhov, rhow, rhoe = calc_filter(rho, rhou, rhov, rhow, rhoe)
    
    n = n + 1
    t = t + dt
    u, v, w, p = compute_cons2prim(rho, rhou, rhov, rhow, rhoe)
    kk = 0.5 * (u ** 2 + v ** 2 + w ** 2)
    kk = np.sum(kk) / (nx * ny * nz)
    
    if n % 50 == 0:
        print(n, t, kk)
    
    file.write('{:>20}'.format(t))
    file.write('{:>20}'.format(kk))
    file.write('\n')
    
    if n % 5000 == 0:
        plot_image(rho, rhou, rhov, rhow, rhoe)

file.close()
plt.show()