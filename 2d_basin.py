import matplotlib.pyplot as plt
import numpy as np

ISIGN = np.array([-1, 1])
G = 6.67e-11

def distance(obs_x, params):
    obs_y = np.zeros_like(obs_x)
    obs_z = np.zeros_like(obs_x)

    rx1 = obs_x - params[0]
    ry1 = obs_y - 1
    rz1 = obs_z - 1
    rx2 = obs_x - params[1]
    ry2 = obs_y - 1000
    rz2 = obs_z - params[2]

    x = [rx1, rx2]
    y = [ry1, ry2]
    z = [rz1, rz2]

    return x, y, z

def grav_component(rho, x, y, z):
    t_sum = 0

    for i in range(2):
        for j in range(2):
            for k in range(2):
                vertex_r = np.sqrt(x[i]**2 + y[j]**2 + z[k]**2)

                ijk = ISIGN[i] * ISIGN[j] * ISIGN[k]

                arg1 = np.arctan2((x[i]*y[j]), (z[k]*vertex_r))
                if np.any(arg1 < 0):
                    arg1[arg1 < 0] += 2*(np.pi)

                arg2 = vertex_r + y[j]
                arg3 = vertex_r + x[i]
                arg2 = np.log(arg2)
                arg3 = np.log(arg3)

                t_sum += ijk * (z[k] * arg1 - x[i] * arg2 - y[j] * arg3)
    
    gravz = rho * G * t_sum * 1e8

    return gravz

def main():
    size = 1000 * float(input("Grid size (km): "))
    bodies = int(input("How many prisms do you want to plot? "))

    obs_x, param_t, rho_t = init(bodies, size)

    grav_t = multi(bodies, param_t, rho_t, obs_x)

    grav_sum = sum(grav_t)
    plots(bodies, obs_x, grav_sum, grav_t)

def init(bodies, size):
    param_t = []
    rho_t = []

    for n in range(bodies):
        temp_xyz, temp_rho = get_params(n+1)
        param_t.append(temp_xyz)
        rho_t.append(temp_rho)
    
    obs_x = np.arange(start=0, stop=(size+1), step=10, dtype=float)

    return obs_x, param_t, rho_t

def multi(bodies, param_t, rho_t, obs_x):
    param_x = []
    param_y = []
    param_z = []
    grav_t = []

    for n in range(bodies):
        temp_x, temp_y, temp_z = distance(obs_x, param_t[n])
        param_x.append(temp_x)
        param_y.append(temp_y)
        param_z.append(temp_z)
        grav_t.append(grav_component(rho_t[n], param_x[n], param_y[n], param_z[n]))
    
    return grav_t

def plots(bodies, obs_x, grav_sum, grav_t):

    fig = plt.figure(figsize=[12.8, 4.8])

    total = fig.add_subplot(1, 2, 1,
    xlabel="x (km)",
    ylabel="Gz (mGal)",
    title="Total Gz")
    total.plot(obs_x/1000, grav_sum)
    total.grid()
    total.set_xlim(0, obs_x[-1]/1000)

    indiv = fig.add_subplot(1, 2, 2,
    xlabel="x (km)",
    ylabel="Gz (mGal)",
    title="Individual Gz")
    for n in range(bodies):
        indiv.plot(obs_x/1000, grav_t[n], label=f"Prism {n+1}")
    
    indiv.grid()
    indiv.set_xlim(0, obs_x[-1]/1000)
    indiv.legend()
    
    plt.show()


def get_params(n=1):
    x1 = (n-1) * 2000
    x2 = x1 + 2000
    z = 1000 * float(input(f"Prism number {n} depth (km): "))
    rho = float(input(f"Prism number {n} density (kg/m³): "))

    params = [x1, x2, z]
    
    return params, rho

main()