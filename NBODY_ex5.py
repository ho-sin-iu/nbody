#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N-Body Simulations with REBOUND
Exercise 5: Resonant capture of planet
@author: Sin-iu Ho
based on example from 
https://github.com/dtamayo/reboundx/blob/master/ipython_examples/ExponentialRC.ipynb
"""
import rebound
import reboundx
import numpy as np
import matplotlib.pyplot as plt

duration = 1e5            # total duration of simulation, default=1e6
efolding = 1e4            # e-holding time (migration rate), default=1e5
plot_evolution = True       # plot the evolution of orbital elements?
plot_ratio_evolution = True # plot the evolution of the ratios of orbital elements?
plot_inifin = True          # plot the snapshot of orbits at the beginning and the end?
sample = 1000             # number of samples used to plot the evolution

####################

def printinfo(text):
    a1, a2 = sim.particles[1].a, sim.particles[2].a
    e1, e2 = sim.particles[1].e, sim.particles[2].e
    P1, P2 = sim.particles[1].P, sim.particles[2].P
    eratio = e1/e2 if np.abs(e2) > 1e-4 else np.NaN # if denominator too small
    vari, pl1, pl2, rati = "", "Planet 1", "Planet 2", "Ratio"
    print(text + "\n" + "-"*35)
    print(f"{vari:3s} {pl1:10s} {pl2:13s} {rati:s}")
    print("-"*35)
    print(f"{ylabels[0][1:-1]:s} {a1:10.3f} {a2:10.3f} {a1/a2:10.3f}")
    print(f"{ylabels[1][1:-1]:s} {e1:10.3f} {e2:10.3f} {eratio:10.3f}")
    print(f"{ylabels[2][1:-1]:s} {P1:10.3f} {P2:10.3f} {P1/P2:10.3f}")
    print("-"*35)
    return None

def plotorbit(text): 
    op = rebound.OrbitPlot(sim, unitlabel="[AU]", lw=2, xlim=[-30,30], ylim=[-30,30])
    op.primary.set_sizes([200])
    op.particles.set_color(['r', 'b'])
    op.orbits[0].set_color('r')
    op.orbits[1].set_color('b')
    op.orbits[0].set_linestyle("-")
    op.orbits[1].set_linestyle(":")
    op.fig.savefig(f"RC_{efolding:.1e}_{duration:.1e}_{text:s}.png", figsize=(6.6), dpi=150)
    return None

####################
# INITIALIZATION
times = np.linspace(0, duration, sample)
varray = np.zeros((2,3,len(times))) # 2 planets x 3 variables x len(times) sample points
tarray = np.zeros_like(times)
axtitles = ["Semimajor axis", "Eccentricity", "Period"]
ylabels = ["$a$", "$e$", "$P$"]
y2labels = ["$a_1/a_2$", "$P_1/P_2$"]

sim = rebound.Simulation() # Initiate rebound simulation
sim.units = ("yr", "AU", "Msun")
sim.add(m=1)
sim.add(m=5.1e-5, a=24., e=0.01, hash="neptune") # Add Neptune (pre-migration) at 24 AU
sim.add(m=6e-6, a=10., e=0.0, hash="planet") # Add an Earth-mass planet
sim.move_to_com()

rebx = reboundx.Extras(sim) # Initiate reboundx
mod_effect = rebx.load_force("exponential_migration") # Add the migration force
rebx.add_force(mod_effect) # Add the migration force

sim.particles["neptune"].params["em_aini"] = 24. # parameter 1: Neptune's initial semimajor axis
sim.particles["neptune"].params["em_afin"] = 10.0 # parameter 2: Neptune's final semimajor axis
sim.particles[1].params["em_tau_a"] = efolding # parameter 3: the migration e-folding time
# this parameter is to play around with
# using 1e5 the resonance will be 3:2, using 1e6 it's 2:1

# SIMULATION
print(f"e-folding = {efolding:.2e}")
print(f"duration = {duration:.2e}\n")
printinfo("[Initial]  " + "t = " + str(sim.t))
if plot_inifin == True:
    plotorbit("Initial")
for t, time in enumerate(times):   # Integrate the system for 1e6 yr
    sim.integrate(time)
    print(f"\rt = {time:.1f} yr", end="")
    for i in range(2):
        varray[i,0,t] = sim.particles[i+1].a
        varray[i,1,t] = sim.particles[i+1].e
        varray[i,2,t] = sim.particles[i+1].P
    tarray[t] = sim.t
print("\r              ")
printinfo("\n[Final]  " + "t = " + str(sim.t))
if plot_inifin == True:
    plotorbit("Final")  
#sim.save(f"RC_{efolding:.1e}_{duration:.1e}.bin")

####################
# PLOTS

if plot_evolution == True:
    # PLOT a, e & P over TIME
    figs, axs = plt.subplots(3,1, sharex=True, figsize=(6,7), dpi=150)
    for n in range(3):
        ax = axs[n]
        ax.plot(tarray, varray[0,n], 'r-', label="Planet 1")
        ax.plot(tarray, varray[1,n], 'b--', label="Planet 2")
        ax.set_xlabel("$t$ [yr]", size="large")
        ax.set_ylabel(ylabels[n], size="large", rotation=0, labelpad=10)
        ax.set_xscale("log")
        #ax.set_yscale("log")
        #ax.set_title(axtitles[n], size="large")
        ax.label_outer()
        ax.grid(True, which="both", color="0.8")
    handles, labels = ax.get_legend_handles_labels()
    figs.legend(handles, labels,
                bbox_to_anchor=(0.5, -0), loc="lower center",
                bbox_transform=figs.transFigure, ncol=2)
    title = f"Resonance capture\ne-folding time = {efolding:.1e},  duration = {duration:.1e}"
    figs.suptitle(title)
    #figs.subplots_adjust(hspace = 0.5)
    #figs.savefig(f"RC_{efolding:.1e}_{duration:.1e}.png")

if plot_ratio_evolution == True:
    # PLOT THE RATIOS OF a & P
    figr, axr = plt.subplots(2,1, sharex=True, figsize=(6,5), dpi=150)
    ratioplot = {0: 0, 1: 2}
    for m in ratioplot:
        n = ratioplot[m]
        ax = axr[m]
        ax.plot(tarray, varray[0,n]/varray[1,n], 'm-')
        ax.set_xlabel("$t$ [yr]", size="large")
        ax.set_ylabel(y2labels[m], size="large", labelpad=10)
        ax.set_xscale("log")
        #ax.set_yscale("log")
        ax.label_outer()
        ax.grid(b=True, which='major', color='0.7', ls='-')
        ax.grid(b=True, which='minor', color='0.8', ls='-')
    figr.suptitle(title)
    #figr.savefig(f"RC_{efolding:.1e}_{duration:.1e}_ratio.png")





