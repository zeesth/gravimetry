import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

ISIGN = np.array([-1, 1])
G = 6.67e-11

def distance(obs_x, obs_y, params):
    obs_z = np.zeros_like(obs_x)

    ## Calculates the distance between the observator and each vertex
    rx1 = obs_x - params[0]
    ry1 = obs_y - 1
    rz1 = obs_z - 0
    rx2 = obs_x - params[1]
    ry2 = obs_y - 1
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
                ## Calculates the vector module of the distance from the observer to a prism vertex
                vertex_r = np.sqrt(x[i]**2 + y[j]**2 + z[k]**2)

                ## Adjusts the sign for each vertex of the prism
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

    obs_x, obs_y, param_t, rho_t = init(bodies, size)

    grav_t = multi(bodies, param_t, rho_t, obs_x, obs_y)
    grav_sum = sum(grav_t)

    plots(bodies, size, param_t, obs_x, obs_y, grav_sum, grav_t)

def init(bodies, size):
    param_t = []
    rho_t = []

    for n in range(bodies):
        temp_xyz, temp_rho = user_input(n+1, size, bodies)
        param_t.append(temp_xyz)
        rho_t.append(temp_rho)

    obs_x = np.arange(start=0, stop=(size+1), step=10, dtype=float)
    obs_y = np.arange(start=0, stop=(1001), step=10, dtype=float)
    obs_x, obs_y = np.meshgrid(obs_x, obs_y)
    param_t = np.array(param_t)

    return obs_x, obs_y, param_t, rho_t

def multi(bodies, param_t, rho_t, obs_x, obs_y):
    param_x = []
    param_y = []
    param_z = []
    grav_t = []

    for n in range(bodies):
        temp_x, temp_y, temp_z = distance(obs_x, obs_y, param_t[n])
        param_x.append(temp_x)
        param_y.append(temp_y)
        param_z.append(temp_z)
        grav_t.append(grav_component(rho_t[n], param_x[n], param_y[n], param_z[n]))
    
    return grav_t

def plots(bodies, size, param_t, obs_x, obs_y, grav_sum, grav_t):
    fig = plt.figure(figsize=[13, 10])
    y_filter = np.where(obs_y==0)[0][0]
    prop = size/(bodies*1000)

    ## Plots total Gz
    total = fig.add_subplot(2, 2, 1,
    xlabel="x (km)",
    ylabel="Gz (mGal)",
    title="Total Gz in function of x")
    total.plot(obs_x[0]/1000, grav_sum[y_filter, :], label="Sum", color="k")
    total.grid()
    total.set_xlim(0, obs_x[0][-1]/1000)

    ## Plots each prism's Gz
    solo = fig.add_subplot(2, 2, 2,
    xlabel="x (km)",
    ylabel="Gz (mGal)",
    title="Gz in function of x")
    for n in range(bodies):
        solo.plot(obs_x[0]/1000, grav_t[n][y_filter, :], label=f"Prism {n+1}")
    solo.grid()
    solo.set_xlim(0, obs_x[0][-1]/1000)

    ## Plots body position
    body1 = fig.add_subplot(2, 2, 3,
    xlabel = "x (km)",
    ylabel = "z (km)")
    for i in range(bodies):
        body1.add_patch(Rectangle(((i) * prop, 0), prop, param_t[i][2]/1000, fc="y", ec="k"))
    body1.set_xlim(0, obs_x[0][-1]/1000)
    body1.set_ylim(0, (np.max(param_t[:, 2]))/800)
    body1.invert_yaxis()
    body1.set_aspect('equal')

    body2 = fig.add_subplot(2, 2, 4,
    xlabel = "x (km)",
    ylabel = "z (km)")
    for i in range(bodies):
        body2.add_patch(Rectangle(((i) * prop, 0), prop, param_t[i][2]/1000, fc="y", ec="k"))
    body2.set_xlim(0, obs_x[0][-1]/1000)
    body2.set_ylim(0, (np.max(param_t[:, 2]))/800)
    body2.invert_yaxis()
    body2.set_aspect('equal')
    
    plt.show()

def user_input(n=1, size=3, bodies=10):
    prop = size/bodies
    x1 = (n-1) * prop
    x2 = x1 + prop
    z = 1000 * float(input(f"Prism number {n} depth (km): "))
    rho = 2

    params = [x1, x2, z]
    return params, rho

main()
