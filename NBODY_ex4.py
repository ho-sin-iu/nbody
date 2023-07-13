#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N-Body Simulations with REBOUND
Exercise 4: Jupiter and Kirkwood gaps
https://ssd.jpl.nasa.gov/horizons/app.html#/
@author: Sin-iu Ho
"""
import numpy as np
import rebound
import matplotlib.pyplot as plt
import time

N = 100 # numer of particles

# masses
mM = 3.21e-7 # MSun
mJ = 9.55e-4  #MSun
# semi-major axis
aM = 1.5273 # AU
aJ = 5.2028 # AU
# eccentricity
eM = 0.093
eJ = 0.048

##########################
def resona(Pli): 
    Pratio = Pli[0]/Pli[1] # Period of asteroid / Period of Jupiter
    return (Pratio)**(2/3) * aJ
res = [[1,3], [2,5], [3,7], [1,2]]

def plot_hist(save=False, final=False):
    # recording the final semi-major axes
    a_li = np.array([sim.particles[i].a for i in range(3, N+3, 1)])
    # plot the histogram of semi-major axes
    binwidth = 0.02 if final == False else 0.005 # AU
    xmin = Mars.a - 1 if final == False else 1.8
    xmax = Jupiter.a + 1 if final == False else 3.6
    color = "teal" if final == False else "lightblue"
    binssetting = np.arange(xmin, xmax, binwidth)
    fighi, axhi = plt.subplots()
    axhi.hist(a_li, bins=binssetting, color=color)
    # plot vertical lines that indicate semi-major axes of Mars and Jupiter
    ymin, ymax = axhi.get_ylim()
    if final == True:
        axhi.set_xlim((1.8, 3.6))
        axhi.set_ylim((ymin, ymax*1.25))
    ycen = (ymin + ymax) * 0.5
    if final == False:
        axhi.axvline(Mars.a, c ="crimson", ls="--")
        axhi.axvline(Jupiter.a, c ="sienna", ls="--")
        axhi.text(Mars.a+0.1, ycen, "Mars", rotation= "vertical")
        axhi.text(Jupiter.a+0.1, ycen, "Jupiter", rotation = "vertical")
    else:
        for r in res:
            ares = resona(r)
            axhi.axvline(ares, c ="teal", ls="--", lw=2.5)
            axhi.text(ares+0.04, ymax*1.1, str(r[1])+":"+str(r[0]), 
                      rotation= "vertical", size="x-large")
    axhi.set_xlabel("Semi-major axis $a$ [A.U.]")
    axhi.set_ylabel("Number")
    axhi.set_title(f"Distribution of semi-major axis\n$N$ = {N:d},   $t$ = {sim.t:.0f} yr")
    plt.show()
    #print("\nI plot at t = {:.3f} y.".format(year))
    if save == True:
        if final == False:
            fighi.savefig(f"JK_{N:d}_{sim.t:.0f}_hist.png", figsize=(6,5), dpi=200)
        else:
            fighi.savefig(f"JK_{N:d}_{sim.t:.0f}_hist_detailed.png", figsize=(6,5), dpi=200)

def get_e():  # get the eccentricities of all small particles 
    eli_ = [sim.particles[i].e for i in rangePar]
    return eli_

def get_a():  # get the semi-major axes of all small particles 
    ali_ = [sim.particles[i].a for i in rangePar]
    return ali_

def plot_ea(save=False):
    ali, eli = get_a(), get_e() # e & a for all small particles 
    figea, axea = plt.subplots()
    axea.scatter(ali, eli, s=2, c="teal")
    axea.plot(Mars.a, Mars.e, marker="p", ms=10, mfc="crimson", mew=0)
    axea.plot(Jupiter.a, Jupiter.e, marker="H", ms=15, mfc="sienna", mew=0)
    axea.set_xlabel("Semi-major axis $a$ [A.U.]")
    axea.set_ylabel("Eccentricity $e$")
    axea.set_xlim(Mars.a-1, Jupiter.a+1)
    axea.set_ylim(0,0.6)
    axea.set_title(f"Eccentricity vs. Semi-major axis\n$N$ = {N:d},   $t$ = {sim.t:.0f} yr")
    plt.show()
    if save == True:
        figea.savefig(f"JK_{N:d}_{sim.t:.0f}_ea.png", figsize=(6,5), dpi=200)
    
def plot_orbit(): # not recommended when N is large
    if N < 200:
        lcli = ["crimson", "sienna"]
        mcli = ["crimson", "sienna"]
        sli = [150, 200]
        lwli = [1.5, 3]
        for i in rangePar:
            lcli.append("teal")
            mcli.append("teal")
            sli.append(5)
            lwli.append(0.8)
        op = rebound.OrbitPlot(sim, unitlabel="[AU]", xlim=[-6,6], ylim=[-6,6])
        for i in range(N+2):
            op.orbits[i].set_linewidth(lwli[i])
            op.orbits[i].set_color(lcli[i])
        op.particles.set_color(mcli)
        op.particles.set_sizes(sli)
        op.primary.set_sizes([500])
    else:
        print("N >= 200")

##########################
# INITIALIZATION (SUN, MARS, JUPITER)

sim = rebound.Simulation()
sim.units = ("yr", "AU", "Msun")
sim.add(m=1.)
sim.add(m=mM, a=aM, e=eM)
sim.add(m=mJ, a=aJ, e=eJ)

sim.move_to_com()
sim.integrator = "leapfrog" 


# INITIALIZATION (SMALL OBJECTS)
for a in np.linspace(2, 4, N): 
    sim.add(a = a, # initally equally distributed semi-major axes
            f = np.random.rand() * 2*np.pi,  # random true anomalies
            omega = np.random.rand() * 2*np.pi, # random arguments of perihelion
            e = 0.5 * np.random.rand() # random eccentricities
            ) # mass set to 0 by default
rangePar = range(3, N+3, 1)
Sun, Mars, Jupiter = sim.particles[0], sim.particles[1], sim.particles[2]
sim.N_active = 3 # Only the first three objects (Sun, Mars, Jupiter) interact mutually.

##########################
# SIMULATION

du = 1000
st = 100
print(f"Number of small particles N = {N:d}")
print(f"Simulation duration t = {du:.3f} yr")
#plot_ea()
#plot_hist()
#plot_orbit()

sim.dt = 0.01
times = np.arange(0, du+st, st)
tick0 = time.time()

for t in times:
    tick = time.time()
    sim.integrate(t)
    #plot_ea()
    if t % 100 == 0:
        plot_hist(save=False)
        plot_ea(save=False)
    tock = time.time()
    dtsys = tock - tick
    tsys = tock - tick0
    print(f"\rt = {sim.t:.3f} yr,  clock = {tsys:.3f} s,  ticktock = {dtsys:.3f} s", end="")

plot_hist(save=False, final=True)
#plot_orbit()
