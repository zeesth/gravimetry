import numpy as np
import matplotlib.pyplot as plt

ISIGN = np.array([-1, 1])
GAMMA = 6.67e-11

def distance(length, P):
    obs_x = np.arange(start=0, stop=(length+1), step=10, dtype=float)
    obs_y = np.zeros(shape=1)
    obs_z = np.zeros(shape=1)

    rx1 = obs_x - P[0]
    ry1 = obs_y - P[1]
    rz1 = obs_z - P[2]
    rx2 = obs_x - P[3]
    ry2 = obs_y - P[4]
    rz2 = obs_z - P[5]

    x = [rx1, rx2]
    y = [ry1, ry2]
    z = [rz1, rz2]

    return x, y, z, obs_x

def grav_component(rho, x, y, z):
    t_sum = 0

    for i in range(2):
        for j in range(2):
            for k in range(2):
                vertex_r = np.sqrt(x[i]**2 + y[j]**2 + z[k]**2)

                ijk = ISIGN[i] * ISIGN[j] * ISIGN[k]

                arg1 = np.arctan2((x[i]*y[j]), (z[k]*vertex_r))
                arg2 = vertex_r + y[j]
                arg3 = vertex_r + x[i]
                arg2 = np.log(arg2)
                arg3 = np.log(arg3)

                t_sum += ijk * (z[k] * arg1 - x[i] * arg2 - y[j] * arg3)
    
    g = rho * GAMMA * t_sum * 1e8

    return g

def main():
    size, P1, rho, P2, rho2 = user_input()

    P1x, P1y, P1z, obs1_x = distance(size, P1)
    P2x, P2y, P2z, obs2_x = distance(size, P2)

    P1grav = grav_component(rho, P1x, P1y, P1z)
    P2grav = grav_component(rho2, P2x, P2y, P2z)

    gravsum = P1grav + P2grav
    plots(obs1_x, gravsum, P1grav, P2grav, size)

def plots(obs_x, grav, grav1, grav2, size):
    plt.plot(obs_x/1000, grav, label="Sum")
    plt.plot(obs_x/1000, grav1, label="Prism 1")
    plt.plot(obs_x/1000, grav2, label="Prism 2")
    plt.xlabel("x (km)")
    plt.ylabel("Gz (mGal)")
    plt.title("Gz in function of x")
    plt.legend()
    plt.xlim(0, size/1000)
    plt.grid()
    plt.show()

def user_input():
    size = 1000 * float(input("Grid size (km): "))
    x1 = 1000 * float(input("First prism's x1 (km): "))
    y1 = 1000 * float(input("First prism's y1 (km): "))
    z1 = 1000 * float(input("First prism's z1 (km): "))
    x2 = 1000 * float(input("First prism's x2 (km): "))
    y2 = 1000 * float(input("First prism's y2 (km): "))
    z2 = 1000 * float(input("First prism's z2 (km): "))
    rho = float(input("First prism's density (kg/m³): "))

    x3 = 1000 * float(input("Second prism's x3 (km): "))
    y3 = 1000 * float(input("Second prism's y3 (km): "))
    z3 = 1000 * float(input("Second prism's z3 (km): "))
    x4 = 1000 * float(input("Second prism's x4 (km): "))
    y4 = 1000 * float(input("Second prism's y4 (km): "))
    z4 = 1000 * float(input("Second prism's z4 (km): "))
    rho2 = float(input("Second prism's density (kg/m³): "))

    P1 = [x1, y1, z1, x2, y2, z2]
    P2 = [x3, y3, z3, x4, y4, z4]
    
    return size, P1, rho, P2, rho2



main()
