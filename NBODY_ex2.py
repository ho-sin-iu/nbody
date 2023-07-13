#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Project:     N-Body Simulations with REBOUND
Descrption:  Exercise 2: Stability of a planetary system
Author:      Sin-iu Ho (sin-iu.ho@student.uni-tuebingen.de)
Version:     2.6
Updated:     April 28, 2023
Annotation:
* Basic version
"""
import numpy as np
import matplotlib.pyplot as plt 
import rebound
plt.rcParams.update({'figure.max_open_warning': 0})

Ncycles = 100 # number of cycles, default = 10000
setups = {1: 1e-5, 2: 1e-7} # masses of Planet 2, for Setups 1 and 2
percents = np.array([10, 50, 100, 150, 1000])
# a list of inital orbital separations in percentage of Delta_c (see Line 58)
colors = ["r","g"]
names = ["a", "e"]
units = [" [AU]", ""]

#==========================================================================
# PARAMETERS FOR INITIAL SET-UP
m1 = 1.0        # mass of primary
m2 = 1e-5       # mass of Planet 1
a2 = 1.0        # semi-major axis of Planet 1
e2 = e3 = 0.0   # eccentricies of planets
omega2  = 0.0   # argument of pericenter of Planet 1
omega3 = np.pi  # argument of pericenter of Planet 2 (see Fig.2 of Ref.[1]), default = np.pi
period = 2*np.pi * (a2**3/m1)**0.5 # expected period of Planet 1

#==========================================================================
# CUSTOM FUNCTIONS
def rhill(a, m): # calculate Hill radius --> planets' radii
    return a*(m/(3.*m1))**(1/3) 

def collision_print_only(sim_pointer, collision): # modified according to Ref.[3]
    global enc_t
    sim = sim_pointer.contents # get simulation object from pointer
    enc_t.append(sim.t)        # record the time of close encounter           
    return 0                   # 0 means "Don't remove either particle"

#==========================================================================
# PARAMETERS FOR INTEGRATION
Nsteps = Ncycles * 100 # 100 steps for each orbit
duration = Ncycles * period / (2*np.pi) # in yr
times = np.linspace(0, duration, Nsteps) # time array, in yr
print(f"Simulation duration = {Ncycles:d} cycles = {duration:.1f} yr ")

#==========================================================================
# MAIN PROGRAM
for setup in setups:
    m3 = setups[setup]       # mass of Planet 2
    print("-"*27 + " [SETUP " + str(setup) + "] " + "-"*27)
    Nencs, Nenct = [], []
    
    # critical orbital separation of planets; formula given in Ref.[2]
    Delta_c = 2.40 * (m2/m1 + m3/m1)**(1/3) 
    print(f"critical separation (Delta_c) = {Delta_c:.4f} AU")

    # Iterate over the five values of "percents"
    for percent in percents:
        Delta = percent/100. * Delta_c # initial orbital separation
        a3 = a2 + Delta # semi-major axis of Planet 2
        r2, r3 = rhill(a2,m2), rhill(a3,m3) # radii of planets
    
        #=====================================
        # Initialization
        sim = rebound.Simulation() # create a simulation
        sim.integrator = "IAS15"   # set the integrator
        sim.collision = "direct"   # collision detection
        # call "collision_print_only" if a collision is found 
        sim.collision_resolve = collision_print_only
        sim.add(m=m1)              # add Primary
        sim.add(m=m2, a=a2, e=e2, r=r2, omega=omega2) # add Planet 1
        sim.add(m=m3, a=a3, e=e3, r=r3, omega=omega3) # add Planet 2
        sim.move_to_com()   # move to the center of mass
        ps = sim.particles  # an array of pointers that update as the simulation runs
    
        #=====================================
        # Integration and data recording
        a_li = np.zeros((sim.N-1, Nsteps))  # semi-major axis
        e_li = np.zeros_like(a_li)          # eccentricity
        d_li = np.zeros(Nsteps)             # distance between 2 planets
        enc_t = []                          # time of close encounter
        
        for i, t in enumerate(times):
            if i % 1000 == 0: # heartbeat: print simulator time every 1000 steps
                achieved = t/duration * 100
                print(f"\rtime = {t:.1f} yr ({achieved:.2f}%)"+" "*40, end="") 
            sim.integrate(t, exact_finish_time=0)
            dp = ps[2] - ps[1] 
            d_li[i] = np.sqrt(dp.x**2 + dp.y**2 + dp.z**2)
            for j in [0,1]: 
                a_li[j,i] = ps[j+1].a
                e_li[j,i] = ps[j+1].e
        
        #=====================================
        # Print orbits
        op = rebound.OrbitPlot(sim, unitlabel="[AU]",color=True)
        op.ax.set_title(f"Setup {setup}:  $\Delta$ = {percent:d}% x $\Delta_c$\n$t = {sim.t:.1f}$ yr")
        op.particles.set_color(colors)
        for j in [0,1]:
            op.orbits[j].set_color(colors[j])
            op.orbits[j].set_label(f"Planet {j+1:d}")
        op.ax.legend(loc="upper right")
        
        #=====================================
        # Diagnose the stability
        if enc_t != []:
            # group the continuous time in enc_t, 
            # criterion of "continuity" : time difference of adjecent entry < 0.1
            groups = np.split(enc_t, np.where((np.diff(enc_t, n=2) > 0.1))[0]+2)
            
            # collect the first and the last entries of each group
            endpts = [[gp[0], gp[-1]]  for gp in groups]
            Nenc = len(endpts)
        else:
            Nenc = 0
            
        #=====================================
        # Print the number of encounter
        stability = "unstable" if Nenc > 0 else "stable"
        Nencs.append(len(enc_t))    # in steps
        Nenct.append(Nenc)          # in "times"
        enctxt = f"encounter : {len(enc_t):6d} steps; {Nenc:3d} times"
        print(f"\r{percent:4d}% x Delta_c | {enctxt:s} | {stability:s}")
  
        #=====================================
        # Plot the evolution
        for q, data in enumerate([a_li, e_li]):
            fig, axe = plt.subplots(2, 1, dpi=100)
            for j in [0,1]:
                ax = axe[j]
                ax.plot(times[:], data[j][:], c=colors[j], lw = 1., ls="-",
                        label=f"Planet {j+1:d}", zorder=5)
                if enc_t != []:
                    for endpt in endpts:
                        ax.axvspan(endpt[0], endpt[1], color="c", alpha=.6, zorder=2)
                ax.set_xlabel("Time $t$ [yr]", size="large")
                ax.set_ylabel(f"${names[q]:s}_{j+1:d}$ {units[q]:s}", size="large")
                ax.label_outer()
                leg = ax.legend(loc="upper right")
                leg.set_zorder(6)
            fig.suptitle(f"Setup {setup}:  $\Delta$ = {percent:d}% x $\Delta_c$", size="large")

    print("\r"+"-"*65+"\n") 