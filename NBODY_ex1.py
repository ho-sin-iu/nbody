#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N-Body Simulations with REBOUND
Exercise 1: The two-body problem
@author: Sin-iu Ho
"""
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.ticker as tck
import rebound

def error(array):
    errarray = np.ones_like(array)
    refindex = 1 if name == "Argument of pericenter" else 0
    for r, row in enumerate(array):
        ref = row[refindex]
        errarray[r] = np.abs(row - ref)/np.abs(ref)
    return errarray

# create a sim object
sim = rebound.Simulation()
# set the integrator to type REB_INTEGRATOR_LEAPFROG
sim.integrator = "leapfrog" 
# choice: leapfrog, IAS15, WHFast, Mercurius, TES, etc.
# https://rebound.readthedocs.io/en/latest/integrators/

# here add your particles ....
sim.add(m=1.0) 
sim.add(m=1e-3, a=1.0, e=0.3, hash = "planet")
# do not forget to move to the center of mass 
sim.move_to_com()
# create time array, let's say 1 orbit, plot 250 times per orbit 
Norbits = 1000
stepperorbit = 1
Nsteps = Norbits * stepperorbit
times = np.linspace(0, Norbits*2*np.pi, Nsteps)
x = np.zeros((sim.N, Nsteps)) # coordinates for both particles
y = np.zeros_like(x)

# use fixed time steps
poweri = 1
powerf = poweri + 4
dtlipower = np.arange(poweri, powerf, 1)
dtli = 1/(10**(dtlipower))
#dtli = np.append(dtli, 5e-6)

plotnames = ["Orbit", "Energy", "Angular momentum", "Eccentricity", "Semimajor axis", "Argument of pericenter"]
varnum = len(plotnames) - 1
xlabels = ["Time $t$" for h in range(varnum+1)]
xlabels[0] = "$x$ cooridnate"
notation = ["", "$E$", "$L_z$", "$e$", "$a$", "$\omega$"]
ylabels = ["$y$ cooridnate","$\log(\Delta E)$", "$\log(\Delta L)$", "$\log(\Delta e)$", "$\log(\Delta a)$", "$\log(\Delta \omega)$"]
colors = ['k', 'r', 'g', 'c', 'b', 'm']
var = np.zeros((varnum, Nsteps))

row, col = 2, 2
fh, fw = 5, 6

figset = [plt.subplots(row, col, sharey=True) for name in plotnames]

errli = np.zeros((len(dtli),varnum,2))

print("sim.integrator")
for n, name in enumerate(plotnames):
    fig = figset[n][0]
    (fh, fw) = (4, 6) if name == "Orbit" else (5, 6) 
    #if n != 0:
        #fig.subplots_adjust(left = 0.0) # 0.2 for ias15, 0.0 for leapfrog
    fig.set_figheight(fh)
    fig.set_figwidth(fw)
    fig.set_dpi(120)
    
for k, dt in enumerate(dtli):
    sim.dt = dt
    print(f"step size = {sim.dt:.1e}")
    # now integrate
    for i, t in enumerate(times):
        print(t, end="\r")
        sim.integrate(t, exact_finish_time=0)
        var[:,i] = [sim.energy(),
                         sim.angular_momentum()[2],
                         sim.particles["planet"].e,
                         sim.particles["planet"].a,
                         sim.particles["planet"].omega]
        for j in range(sim.N):
            x[j,i] = sim.particles[j].x
            y[j,i] = sim.particles[j].y 
    err = np.log10(error(var))
    
    for n, name in enumerate(plotnames):
        ax = figset[n][1][k//col,k%col]
        fig = figset[n][0]
        ax.grid(True)
        if name == "Orbit":
            tpad = 20
            for j in range(sim.N):
                ax.scatter(x[j], y[j], color=colors[n], s=.8)
            ax.set_aspect("equal")
        else:
            tpad = 25
            ax.scatter(times, err[n-1], color=colors[n], s=.8)
            #ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g$\pi$')) # for 1 cycle
            #ax.xaxis.set_major_locator(tck.MultipleLocator(base=1.0)) # for 1 cycle
            if sim.integrator == 'leapfrog':
                ybase = 3.0 if name != "Angular momentum" else 1.0 # for leapfrog
            else:
                ybase = 0.5 if name != "Angular momentum" else 0.2
            ax.yaxis.set_major_locator(tck.MultipleLocator(base=ybase)) 
        ax.set_title(f"$\Delta t = ${dtli[k]:.1g}") 
        
        print(f"... {name:s}: done")
        ax.label_outer()
        
        ax.tick_params(top=True, right=True, direction='in')
        title = f"{name:s} {notation[n]:s},  {sim.integrator:s},  {Norbits:d} x {stepperorbit:d} steps"
        axn = fig.add_subplot(111, frameon=False)
        axn.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        axn.set_xlabel(xlabels[n], size='large')
        axn.set_ylabel(ylabels[n], size='large', labelpad=12.0)
        axn.set_title(title, pad=tpad)
        
    for r, row in enumerate(error(var)):
        errli[k][r][0] = np.average(row)
        errli[k][r][1] = np.std(row)
        


for n, name in enumerate(plotnames):
    nametosave = name.replace(' ', '_')
    figset[n][0].savefig(f"{Norbits:d}x{stepperorbit:d}_{sim.integrator:s}_{nametosave:s}.png", dpi=500)

if sim.integrator == 'ias15':
    varnum -= 1
figee, axee = plt.subplots(nrows=varnum, ncols=1, figsize=(6,11),dpi=100)
for n in range(varnum):
    axe = axee[n]
    for r, dtrow in enumerate(errli):
        axe.errorbar(dtli[r], dtrow[n,0], yerr=dtrow[n,1], 
                     fmt='o-', color=colors[n+1], capsize=4,
                     markersize=4)
        axe.yaxis.set_major_locator(tck.MultipleLocator(5)) 
        axe.set_xlabel("Time step $\Delta t$")
        axe.set_xscale("log")
        axe.set_yscale("log")
        axe.set_title(f"{plotnames[n+1]:s} {notation[n+1]:s}", fontsize=10, y=1.0) 
        axe.label_outer()
        axe.grid(True, which="both", color="0.9")
figee.subplots_adjust(hspace = 0.45)
axnn = figee.add_subplot(111, frameon=False)
axnn.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
axnn.set_ylabel("Error", labelpad=15.0)
axnn.set_title(f"{sim.integrator:s},  {Norbits:d} x {stepperorbit:d} steps", size="x-large", pad=25)
figee.savefig(f"{Norbits:d}x{stepperorbit:d}_{sim.integrator:s}.png", dpi=500)