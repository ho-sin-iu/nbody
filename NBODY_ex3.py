#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Project:     N-Body Simulations with REBOUND
Descrption:  Exercise 3: Stability of Saturn's rings
Author:      Sin-iu Ho (sin-iu.ho@student.uni-tuebingen.de)
Version:     2.1
Updated:     April 28, 2023
Annotation:
* n = 100 is super slow! (Try a better integration scheme, a different integrator?)
* No objects are removed during simulation.
* results:    (first unstable case occurs at)
                n = 8    -->   1.02 m_c
                n = 10   -->   1.03 m_c
                n = 36   -->   1.03 m_c
                n = 100  -->   1.03 m_c
"""
import os
import numpy as np
import rebound as rb
from matplotlib import pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
plt.rcParams["figure.dpi"] = 120

plot_orbit_or_not = True
save_orbit_or_not = False
newfoldername = "test" # save the plots in this folder
cycles = 10

#==========================================================================
# print customized orbit
def plot_orbit(sim):
    if plot_orbit_or_not == True:
        sim.move_to_hel()
        op = rb.OrbitPlot(sim, unitlabel="[AU]")
        op.ax.set_title(f"$n = {n:d}$,  $m = {mratio:.3f} m_c$\n{sim.t/period:.3f} cycles",
                        size="large")
        op.primary.set_sizes([200])
        size_li, color_li = [], []
        for i in range(sim.N-1):
            if sim.calculate_orbits()[i].e >= 1: # not-orbiting
                op.orbits[i].set_color("r")  
                op.orbits[i].set_linewidth(1)
                color_li.append("r")
                size_li.append(60)
            else: # orbiting
                op.orbits[i].set_color("y")  
                op.orbits[i].set_linewidth(0.5)
                color_li.append("orange")
                size_li.append(30)
        op.particles.set_color(color_li)
        op.particles.set_sizes(size_li)
        if save_orbit_or_not == True and sim.t != 0:
            re = f"_{removed:d}" if removed > 0 else ""
            os.chdir(newpath)
            fn = f"{n:d}_{mratio:.3f}_{re}.png"
            op.fig.savefig(fn, bbox_inches='tight', dpi=200)
            plt.close(op.fig)
    return None
#==========================================================================

print(f"cycles = {cycles}")
a0 = 1.0 # AU
mSun = 1.0
mSaturn = mSun * 2.85716656e-4
gamma_dict = {8: 2.412, 10: 2.375, 36: 2.306, 100: 2.300}  # Table 1 in Ref. [1]

# lists to be iterated; reduce the numbner of these lists to do a quick test
n_list = [8, 10, 36, 100] # number of small objects
m_list = np.linspace(0.9, 1.1, 20, endpoint=False) 
# m_list = mass of small objects, {0.9, 0.91, 0.92, ..., 1.09} * m_c
#==========================================================================
# Make directory
if save_orbit_or_not == True:
    cwd = os.getcwd()
    newpath = cwd + "/" + newfoldername
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    print(f"path: {newpath}")

# Main program
for n in n_list:
    # Calculate initial parameters
    print("\n"+"="*40 + f"\nn      = {n:d}")
    gamma = gamma_dict[n]
    print(f"gamma  = {gamma:.3f}")
    mCritical = gamma * mSaturn / n**3 # critical mass, Eq. 51 in Ref. [1], 
    print(f"m_c    = {mCritical:.3e} M_Sun\n") 
    
    k = np.arange(1, n, 1) # k = 1, 2, ..., n - 1
    In = np.sum( 1/(4 * np.sin(k*np.pi/n)) )  # definded in Ref. [1]
    omega = (1.0 * (mSaturn + mCritical * In) / a0**3)**0.5 # ang. velocity, Eq. 50 in Ref. [1]
    period = 2 * np.pi / omega 
    v0 = a0 * omega
    
    step = period * 1e-5
    duration = period * cycles
    #times = np.arange(0, duration, step)
    print(f"duration = {duration:.3f}")

    for mratio in m_list:
        m = mratio * mCritical
        removed = 0

        # Initializing
        sim = rb.Simulation()
        sim.integrator = "leapfrog" 
        sim.add(m = mSaturn) # Primary 
        rd = 0 # phase of the first small object
        for i in range(n): # i = 0, 1, ..., n-1
            phase = 2 * np.pi * (rd + i / n) # phase of every small objects
            x0 = a0 * np.cos(phase)
            y0 = a0 * np.sin(phase)
            vx0 = v0 * -np.sin(phase)
            vy0 = v0 * np.cos(phase)
            sim.add(m = m, x = x0, y=y0, vx=vx0, vy=vy0)
        sim.move_to_com()
        sim.dt = step
        
        # Integration
        for h in range(10):
            sim.steps(int(duration/step/10))
            achieved = sim.t/duration * 100
            print(f"\rtime = {sim.t:7.3f}   ({achieved:5.2f}%)", end="") # heartbeat 

        # Stability diagnosis
        # e_li: eccentricities of every small bodies
        e_li = np.array([orbits.e for orbits in sim.calculate_orbits()]) 
        gone = len(np.nonzero(e_li >= 1)[0])
        
        # Print result
        stability = "stable" if gone == 0 else "unstable"
        gonetxt = " "*10 if gone == 0 else f" ({gone:d} gone)"
        print(f"\r{mratio:.3f} x m_c | " + stability + gonetxt)
        
        plot_orbit(sim)
        sim = None

