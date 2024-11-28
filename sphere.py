import matplotlib.pyplot as plt
import numpy as np

G = 6.67e-11

def calc_mass(radius, rho):
    vol = np.pi * radius**3 * 4 / 3
    mass = rho * vol

    return mass

def distance(obs_x, obs_y, obs_z, sph_x, sph_y, sph_z):
    rx = obs_x - sph_x
    ry = obs_y - sph_y
    rz = obs_z - sph_z

    r = np.sqrt(rx**2 + ry**2 + rz**2)
    r3 = r**3

    return rx, ry, rz, r3

def grav_component(mass, rx, ry, rz, r3):
    grav_x = -G * mass * rx / r3
    grav_y = -G * mass * ry / r3
    grav_z = -G * mass * rz / r3

    grav_x = grav_x * 1e8
    grav_y = grav_y * 1e8
    grav_z = grav_z * 1e8

    return grav_x, grav_y, grav_z

def observator(length):
    obs_x = np.arange(start=0, stop=length+1, step=10, dtype=float)
    obs_y = np.arange(start=0, stop=length+1, step=10, dtype=float)

    obs_x, obs_y = np.meshgrid(obs_x, obs_y)

    obs_z = np.zeros_like(obs_x)

    return obs_x, obs_y, obs_z

def sphere_card(radius, sph_x, sph_y, sph_z):
    theta = np.linspace(0, 2 * np.pi, 100)
    phi = np.linspace(0, np.pi, 50)
    phi, theta = np.meshgrid(theta, phi)

    card_x = radius * np.sin(theta) * np.cos(phi)
    card_y = radius * np.sin(theta) * np.sin(phi)
    card_z = radius * np.cos(theta)

    card_x += sph_x
    card_y += sph_y
    card_z += sph_z

    return card_x, card_y, card_z

def plots(obs_x, obs_y, grav_z, sph_x, sph_y, sph_z, size):
    lim = [0, size/1000]

    fig = plt.figure(figsize=plt.figaspect(0.5))

    delta_gz = fig.add_subplot(1, 2, 1, projection="3d",
    xlabel="x (km)",
    ylabel="y (km)",
    zlabel="Gz (mGal)",
    title="Gz in function of x and y")
    delta_gz.plot_surface(obs_x/1000, obs_y/1000, grav_z, cmap="viridis")

    sp_sph = fig.add_subplot(1, 2, 2, projection="3d",
    xlim=lim,
    ylim=lim,
    zlim=lim,
    xlabel="x (km)",
    ylabel="y (km)",
    zlabel="z (km)",
    aspect="equal",
    title="Sphere Projection")
    sp_sph.plot_surface(sph_x/1000, sph_y/1000, sph_z/1000, cmap="viridis")

    plt.show()

def user_input():
    size = 1000 * float(input("Grid size (km): "))
    x = 1000 * float(input("Sphere's x (km): "))
    y = 1000 * float(input("Sphere's y (km): "))
    z = 1000 * float(input("Sphere's z (km): "))
    rad = float(input("Sphere's radius (m): "))
    rho = float(input("Sphere's density (kg/m³): "))
    
    return size, x, y, z, rad, rho

def main(debug=False):
    if debug:
        size = 12000
        x = 7000
        y = 6000
        z = 5000
        rad = 777
        rho = 16
    else:
        size, x, y, z, rad, rho = user_input()
    
    obs_x, obs_y, obs_z = observator(size)

    rx, ry, rz, r3 = distance(obs_x, obs_y, obs_z, x, y, z)

    mass = calc_mass(rad, rho)

    grav_x, grav_y, grav_z = grav_component(mass, rx, ry, rz, r3)

    sph_x, sph_y, sph_z = sphere_card(rad, x, y, z)

    plots(obs_x, obs_y, grav_z, sph_x, sph_y, sph_z, size)
    
main(True)