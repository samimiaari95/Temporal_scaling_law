import os
from random_functions import *


def plot_pressure_profile():
    plt.rcParams.update({'font.size': 22})
    # plt.figure(figsize=(16, 9))
    fig, ax = plt.subplots(figsize=(16, 9))
    name='/p/project1/cslts/miaari1/python_scripts/Temporal_scaling_law/pressureprofileexample/infiltration'
    pressures = {}
    #pressures["z"] = [z/100 for z in range(5,400, 10)]
    pressures["z"] = [-1*z/10 for z in range(0,40, 1)]
    #pressures["z"] = list(reversed(pressures["z"]))

    for t in range(0, 210, 10):
        print(name + '.out.satur.'+('{:05d}'.format(t))+'.pfb')
        data = SLOTH.sloth.IO.read_pfb(name + '.out.press.'+ ('{:05d}'.format(t)) + '.pfb')

        plt.plot(data[:,0,0], list(reversed(pressures['z'])), color="black")
        #plt.plot(data[:,0,0], list(pressures['z']), color="black")
        #if t==0:
        #    plt.annotate(f"t={t} hr", xy=(min(data[:,0,0]+0.15), max(pressures['z'])-0.15), color="black")
        #elif t==205:
        #    plt.annotate(r"$t=t_{SS}$", xy=(min(data[:,0,0]+0.01), max(pressures['z'])-0.15), color="black")
        
        pressures[f"time={t}"] = data[:,0,0]

    # ax.gca().invert_yaxis()
    #ax.invert_yaxis()
    # your plotting code...
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel("Pressure head (m)")
    ax.set_ylabel("Depth (m)")
    fig.savefig(os.path.join("/p/project1/cslts/miaari1/python_scripts/Temporal_scaling_law", "pressrue_profile.png"))
