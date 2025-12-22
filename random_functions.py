import pandas as pd
import numpy as np
import SLOTH.sloth.IO
import os
import re
import random
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from utils import powerlaw_func, linear_law
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
from utils import make_dir

plt.rcParams.update({'font.size': 22})

DIRPATH = os.path.dirname(os.path.realpath(__file__))
SCRATCHPATH = "/p/scratch/cslts/miaari1"

def plot_ss_profiles():
    plt.rcParams.update({'font.size': 22})
    folderpath=os.path.join(DIRPATH, 'ss_profilesplot')
    z = list(range(41, 1, -1))
    z = [x*0.1 for x in z]

    files = [os.path.join(folderpath, x) for x in os.listdir(folderpath)]
    for file in files:
        name = os.path.basename(file)
        name = name.replace(".pfb", "")
        namelist = name.split("_")
        q = namelist[0].replace("-",".")
        k = namelist[1].replace("-",".")
        t = namelist[2]
        qk = float(q)/float(k)

        data = SLOTH.sloth.IO.read_pfb(file)

        plt.plot(data[:,0,0],z, color="black")
        plt.annotate(f"q/k={qk}", xy=(data[:,0,0][-1], z[-1]), color="black")
    
    plt.gca().invert_yaxis()
    plt.xlim([-4, 0.2])
    plt.xlabel("Pressure (m)")
    plt.ylabel("Soil depth (m)")
    plt.savefig(os.path.join(DIRPATH, 'ss_profilesplot', 'ss_profiles.png'))

def simulation_bash():
    filepath = os.path.join(DIRPATH, "parflow_simulations.sh")
    f = open(filepath, "r")
    sbatch = f.readlines()
    f.close()
    sbatch = sbatch[:19]
    run_command = "tclsh infiltration.tcl"
    run_command_2 = "tclsh exfiltration.tcl"    
    for i in range(0, 100):
        case_path = f"cd /p/scratch/cslts/miaari1/infexfcases/test_case{i}"
        sbatch.append(case_path)
        sbatch.append(run_command)
        sbatch.append(run_command_2)
    
    with open(os.path.join(DIRPATH, "parflow_simulations.sh"), 'w') as sbatchfile:
        sbatchfile.write('\n'.join(sbatch))
    sbatchfile.close()

def get_toplayer_pressure(pfd_path):
    data = SLOTH.sloth.IO.read_pfb(pfd_path)
    surface_pressure = data[-1,0,0]
    return surface_pressure

def steadystate_kinsol():
    # infiltration
    dir_path = os.path.join(SCRATCHPATH, "infexfcases")
    data = {"case_index": [], "q": [], "k": [], "d": [], "alfa": [], "n": [], "theta_r": [], "theta_s": [], "time": [], "toplayer_pressure": []}
    for case_index in range(len(os.listdir(dir_path))):
        print(case_index)
        reached = False
        case_dir = os.path.join(dir_path, f"test_case{case_index}")
        f = open(os.path.join(case_dir, "settings.txt"), "r")
        settings = f.readlines()
        f.close()

        kinsolfile = os.path.join(case_dir, "infiltration.out.kinsol.log")
        if not os.path.exists(kinsolfile):
            print("doesn't exit")
            continue
        filedata = open(kinsolfile)
        lines = filedata.readlines()
        time_kinsol = [x for x in lines if "KINSOL starting step for time" in x or "KINSolInit nni=    0  fnorm=" in x]
        for i in range(1, len(time_kinsol)):
            if "KINSolInit nni=    0  fnorm=" not in time_kinsol[i] and "KINSolInit nni=    0  fnorm=" not in time_kinsol[i-1]:
                index = i-1
                time = re.findall("\d+\.\d+",time_kinsol[index])
                #print(time)
                time = float(time[0])
                data["time"].append(time)
                data["k"].append(float(settings[1].replace("\n","")))
                data["q"].append(float(settings[2].replace("\n","")))
                data["alfa"].append(float(settings[3].replace("\n","")))
                data["n"].append(float(settings[4].replace("\n","")))
                data["theta_r"].append(float(settings[5].replace("\n","")))
                data["theta_s"].append(float(settings[6].replace("\n","")))
                data["d"].append(float(settings[7].replace("\n","")))
                data["case_index"].append(case_index)
                data["toplayer_pressure"].append(get_toplayer_pressure(os.path.join(case_dir, "infiltration.out.press.00001.pfb")))
                reached = True
                break
            if i == len(time_kinsol)-1 and not reached:
                print(f"not reached: {case_index}")
    print(len(data["time"]))
    df = pd.DataFrame(data)
    print(df)
    df.to_csv(os.path.join(os.path.dirname(os.path.realpath(__file__)), "outputs", "infiltration_toplayerpressure_simplescenarios.csv"), index=False)

def steadystate_kinsol_drain_inf(case_index):
    # infiltration related to drainage script
    dir_path = os.path.join(SCRATCHPATH, "infexfcases")

    reached = False
    case_dir = os.path.join(dir_path, f"test_case{case_index}")
    f = open(os.path.join(case_dir, "settings.txt"), "r")
    settings = f.readlines()
    f.close()
    kinsolfile = os.path.join(case_dir, "infiltration.out.kinsol.log")
    if not os.path.exists(kinsolfile):
        print(f" infiltration doesn't exit for {case_index}")
        return None
    filedata = open(kinsolfile)
    lines = filedata.readlines()
    time_kinsol = [x for x in lines if "KINSOL starting step for time" in x or "KINSolInit nni=    0  fnorm=" in x]
    for i in range(1, len(time_kinsol)):
        if "KINSolInit nni=    0  fnorm=" not in time_kinsol[i] and "KINSolInit nni=    0  fnorm=" not in time_kinsol[i-1]:
            index = i-1
            time = re.findall("\d+\.\d+",time_kinsol[index])

            time = float(time[0])
            reached = True
            break
        if i == len(time_kinsol)-1 and not reached:
            print(f" infiltration not reached: {case_index}")
            return None
    return time

def steadystate_kinsol_drainage():
    # exfiltration
    dir_path = os.path.join(SCRATCHPATH, "infexfcases")
    data = {"case_index": [], "q": [], "k": [], "d": [], "alfa": [], "n": [], "theta_r": [], "theta_s": [], "inf_time": [], "exf_time": [], "toplayer_pressure": []}
    for case_index in range(len(os.listdir(dir_path))):
        print(case_index)
        reached = False
        case_dir = os.path.join(dir_path, f"test_case{case_index}")
        f = open(os.path.join(case_dir, "settings.txt"), "r")
        settings = f.readlines()
        f.close()

        inf_time = steadystate_kinsol_drain_inf(case_index)
        if inf_time:
            kinsolfile = os.path.join(case_dir, "exfiltration.out.kinsol.log")
            if not os.path.exists(kinsolfile):
                print("doesn't exit")
                continue
            filedata = open(kinsolfile)
            lines = filedata.readlines()
            time_kinsol = [x for x in lines if "KINSOL starting step for time" in x or "KINSolInit nni=    0  fnorm=" in x]
            for i in range(1, len(time_kinsol)):
                if "KINSolInit nni=    0  fnorm=" not in time_kinsol[i] and "KINSolInit nni=    0  fnorm=" not in time_kinsol[i-1]:
                    index = i-1
                    time = re.findall("\d+\.\d+",time_kinsol[index])
                    #print(time)
                    time = float(time[0])
                    data["inf_time"].append(inf_time)
                    data["exf_time"].append(time)
                    data["k"].append(float(settings[1].replace("\n","")))
                    data["q"].append(float(settings[2].replace("\n","")))
                    data["alfa"].append(float(settings[3].replace("\n","")))
                    data["n"].append(float(settings[4].replace("\n","")))
                    data["theta_r"].append(float(settings[5].replace("\n","")))
                    data["theta_s"].append(float(settings[6].replace("\n","")))
                    data["d"].append(float(settings[7].replace("\n","")))
                    data["case_index"].append(case_index)
                    data["toplayer_pressure"].append(get_toplayer_pressure(os.path.join(case_dir, "exfiltration.out.press.00000.pfb")))
                    reached = True
                    break
                if i == len(time_kinsol)-1 and not reached:
                    print(f"not reached: {case_index}")
    print(len(data["exf_time"]))
    df = pd.DataFrame(data)
    print(df)
    df.to_csv(os.path.join(DIRPATH, "outputs", "inf_exf_specialcases_toplayer.csv"), index=False)



def parflow_namelist_infiltration():
    dir_path = os.path.join(SCRATCHPATH, "testcases")
    filepath = os.path.join(DIRPATH, "inputs", "infiltration.tcl")
    VG_params_path = os.path.join(DIRPATH, "inputs", "WaterFlowParameters.csv")
    f = open(filepath, "r")
    namelist = f.readlines()
    f.close()

    d_range = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    VG_params = pd.read_csv(VG_params_path)
    VG_indexes = [0,1,2,3,4,5,6,7,8,9,10,11]

    k_namelist = namelist[64]
    q_namelist = namelist[214]
    d_namelist = namelist[54]
    Nz_namelist = namelist[31]
    alfa_namelist1 = namelist[142]
    alfa_namelist2 = namelist[151]
    n_namelist1 = namelist[143]
    n_namelist2 = namelist[152]
    porosity_namelist = namelist[129]
    Sres_namelist = namelist[153]

    for i in range(1000):
        new_namelist = namelist
        case_path = os.path.join(dir_path, f"test_case{i}")
        os.mkdir(case_path)

        VG_index = random.choice(VG_indexes)
        
        alfa = VG_params["Alpha"].iloc[VG_index]
        n = VG_params["n"].iloc[VG_index]
        theta_r = VG_params["Sr"].iloc[VG_index]
        theta_s = VG_params["Ss"].iloc[VG_index]
        k = VG_params["Ks"].iloc[VG_index]

        q_k = random.uniform(0.0001, 1)
        q = q_k * k
        
        new_namelist[64] = k_namelist.replace("0.01",f"{k}")
        new_namelist[214] = q_namelist.replace("0.001", f"{q}")
        new_namelist[142] = alfa_namelist1.replace("1.0", f"{alfa}")
        new_namelist[151] = alfa_namelist2.replace("1.", f"{alfa}")
        new_namelist[143] = n_namelist1.replace("2.0", f"{n}")
        new_namelist[152] = n_namelist2.replace("2.0", f"{n}")
        new_namelist[129] = porosity_namelist.replace("0.5", f"{theta_s}") # porosity = theta_s
        #new_namelist[154] = Ssat_namelist.replace("1.0", f"{theta_s/theta_s}") # porosity/theta_s
        new_namelist[153] = Sres_namelist.replace("0.2", f"{theta_s/theta_r}") # porosity/theta_r
        d = random.choice(d_range)
        new_namelist[54] = d_namelist.replace("4.0", f"{d}")

        new_namelist[31] = Nz_namelist.replace("40", f"{int(d/0.1)}")

        settings = ["k,q,alfa,n,theta_r,theta_s,d",f"{k}",f"{q}",f"{alfa}",f"{n}",f"{theta_r}", f"{theta_s}",f"{d}"]
        with open(os.path.join(dir_path, f"test_case{i}", "settings.txt"),'w') as settingsfile:
            settingsfile.write('\n'.join(settings))
        settingsfile.close()

        with open(os.path.join(dir_path, f"test_case{i}", "infiltration.tcl"), 'w') as tclfile:
            tclfile.write('\n'.join(new_namelist))
        tclfile.close()
    return

def parflow_namelist_inex():
    dir_path = os.path.join(SCRATCHPATH, "testcases")
    filepath = os.path.join(DIRPATH, "inputs", "infiltration.tcl")
    exfilepath = os.path.join(DIRPATH, "inputs", "exfiltration.tcl")
    VG_params_path = os.path.join(DIRPATH, "inputs", "WaterFlowParameters.csv")
    
    f = open(filepath, "r")
    namelist = f.readlines()
    f.close()

    exf = open(exfilepath, "r")
    exnamelist = exf.readlines()
    exf.close()

    d_range = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    VG_params = pd.read_csv(VG_params_path)
    VG_indexes = [0,1,2,3,4,5,6,7,8,9,10,11]

    k_namelist = namelist[64]
    q_namelist = namelist[214]
    d_namelist = namelist[54]
    Nz_namelist = namelist[31]
    alfa_namelist1 = namelist[142]
    alfa_namelist2 = namelist[151]
    n_namelist1 = namelist[143]
    n_namelist2 = namelist[152]
    porosity_namelist = namelist[129]
    Ssat_namelist = namelist[154]
    Sres_namelist = namelist[153]

    for i in range(1400, 1600):
        new_namelist = namelist
        exnew_namelist = exnamelist
        case_path = os.path.join(dir_path, f"test_case{i}")
        os.mkdir(case_path)

        VG_index = random.choice(VG_indexes)
        
        alfa = VG_params["Alpha"].iloc[VG_index]
        n = VG_params["n"].iloc[VG_index]
        theta_r = VG_params["Sr"].iloc[VG_index]
        theta_s = VG_params["Ss"].iloc[VG_index]
        k = VG_params["Ks"].iloc[VG_index]

        q_k = random.uniform(0.0001, 0.001)
        q = q_k * k

        # infiltration        
        new_namelist[64] = k_namelist.replace("0.01",f"{k}")
        new_namelist[214] = q_namelist.replace("0.001", f"{q}")
        new_namelist[142] = alfa_namelist1.replace("1.0", f"{alfa}")
        new_namelist[151] = alfa_namelist2.replace("1.", f"{alfa}")
        new_namelist[143] = n_namelist1.replace("2.0", f"{n}")
        new_namelist[152] = n_namelist2.replace("2.0", f"{n}")
        new_namelist[129] = porosity_namelist.replace("0.5", f"{theta_s}") # porosity = theta_s
        new_namelist[153] = Sres_namelist.replace("0.2", f"{theta_r/theta_s}") # theta_r/porosity
        d = random.choice(d_range)
        new_namelist[54] = d_namelist.replace("4.0", f"{d}")
        new_namelist[31] = Nz_namelist.replace("40", f"{int(d/0.1)}")

        # exfiltration
        exnew_namelist[64] = k_namelist.replace("0.01",f"{k}")
        #exnew_namelist[214] = q_namelist.replace("0.001", f"{q}")
        exnew_namelist[142] = alfa_namelist1.replace("1.0", f"{alfa}")
        exnew_namelist[151] = alfa_namelist2.replace("1.", f"{alfa}")
        exnew_namelist[143] = n_namelist1.replace("2.0", f"{n}")
        exnew_namelist[152] = n_namelist2.replace("2.0", f"{n}")
        exnew_namelist[129] = porosity_namelist.replace("0.5", f"{theta_s}") # porosity = theta_s
        exnew_namelist[153] = Sres_namelist.replace("0.2", f"{theta_r/theta_s}") # theta_r/porosity
        #d = random.choice(d_range)
        exnew_namelist[54] = d_namelist.replace("4.0", f"{d}")
        exnew_namelist[31] = Nz_namelist.replace("40", f"{int(d/0.1)}")


        settings = ["k,q,alfa,n,theta_r,theta_s,d",f"{k}",f"{q}",f"{alfa}",f"{n}",f"{theta_r}", f"{theta_s}",f"{d}"]
        with open(os.path.join(dir_path, f"test_case{i}", "settings.txt"),'w') as settingsfile:
            settingsfile.write('\n'.join(settings))
        settingsfile.close()

        with open(os.path.join(dir_path, f"test_case{i}", "infiltration.tcl"), 'w') as tclfile:
            tclfile.write('\n'.join(new_namelist))
        tclfile.close()

        with open(os.path.join(dir_path, f"test_case{i}", "exfiltration.tcl"), 'w') as extclfile:
            extclfile.write('\n'.join(exnew_namelist))
        extclfile.close()
    return

def parflow_namelist_inex_fixed():
    dir_path = os.path.join(SCRATCHPATH, "infexfcases")
    filepath = os.path.join(DIRPATH, "inputs", "infiltration.tcl")
    exfilepath = os.path.join(DIRPATH, "inputs", "exfiltration.tcl")
    VG_params_path = os.path.join(DIRPATH, "inputs", "WaterFlowParameters.csv")
    
    f = open(filepath, "r")
    infnamelist = f.readlines()
    f.close()

    exf = open(exfilepath, "r")
    exnamelist = exf.readlines()
    exf.close()

    VG_params = pd.read_csv(VG_params_path)
    VG_indexes = [0,1,2,3,4,5,6,7,8,9,10,11]

    k_namelist = infnamelist[64]
    q_namelist = infnamelist[214]
    d_namelist = infnamelist[54]
    Nz_namelist = infnamelist[31]
    alfa_namelist1 = infnamelist[142]
    alfa_namelist2 = infnamelist[151]
    n_namelist1 = infnamelist[143]
    n_namelist2 = infnamelist[152]
    porosity_namelist = infnamelist[129]
    Ssat_namelist = infnamelist[154]
    Sres_namelist = infnamelist[153]
    #q_vg_list = {"0.1": [0,1], "0.01":[2,3,6], "0.001":[4,5,9]}
    q_vg_list = {"0.001":VG_indexes, "0.0001":VG_indexes, "0.00015":VG_indexes}
    dlist = [1, 2, 3]
    i = 0
    for q,v in q_vg_list.items():
        for ind in v:
            for d in dlist:
                alfa = VG_params["Alpha"].iloc[ind]
                n = VG_params["n"].iloc[ind]
                theta_r = VG_params["Sr"].iloc[ind]
                theta_s = VG_params["Ss"].iloc[ind]
                k = VG_params["Ks"].iloc[ind]

                new_namelist = infnamelist
                exnew_namelist = exnamelist
                make_dir(os.path.join(dir_path, f"test_case{i}"))
                
                # infiltration        
                new_namelist[64] = k_namelist.replace("0.01",f"{k}")
                new_namelist[214] = q_namelist.replace("0.001", f"{q}")
                new_namelist[142] = alfa_namelist1.replace("1.0", f"{alfa}")
                new_namelist[151] = alfa_namelist2.replace("1.", f"{alfa}")
                new_namelist[143] = n_namelist1.replace("2.0", f"{n}")
                new_namelist[152] = n_namelist2.replace("2.0", f"{n}")
                new_namelist[129] = porosity_namelist.replace("0.5", f"{theta_s}") # porosity = theta_s
                new_namelist[153] = Sres_namelist.replace("0.2", f"{theta_r/theta_s}") # theta_r/porosity
                #d = random.choice(d_range)
                new_namelist[54] = d_namelist.replace("4.0", f"{d}")
                new_namelist[31] = Nz_namelist.replace("40", f"{int(d/0.1)}")

                # exfiltration
                exnew_namelist[64] = k_namelist.replace("0.01",f"{k}")
                #exnew_namelist[214] = q_namelist.replace("0.001", f"{q}")
                exnew_namelist[142] = alfa_namelist1.replace("1.0", f"{alfa}")
                exnew_namelist[151] = alfa_namelist2.replace("1.", f"{alfa}")
                exnew_namelist[143] = n_namelist1.replace("2.0", f"{n}")
                exnew_namelist[152] = n_namelist2.replace("2.0", f"{n}")
                exnew_namelist[129] = porosity_namelist.replace("0.5", f"{theta_s}") # porosity = theta_s
                exnew_namelist[153] = Sres_namelist.replace("0.2", f"{theta_r/theta_s}") # theta_r/porosity
                #d = random.choice(d_range)
                exnew_namelist[54] = d_namelist.replace("4.0", f"{d}")
                exnew_namelist[31] = Nz_namelist.replace("40", f"{int(d/0.1)}")


                settings = ["k,q,alfa,n,theta_r,theta_s,d",f"{k}",f"{q}",f"{alfa}",f"{n}",f"{theta_r}", f"{theta_s}",f"{d}"]
                print(settings)
                with open(os.path.join(dir_path, f"test_case{i}", "settings.txt"),'w') as settingsfile:
                    settingsfile.write('\n'.join(settings))
                settingsfile.close()

                with open(os.path.join(dir_path, f"test_case{i}", "infiltration.tcl"), 'w') as tclfile:
                    tclfile.write('\n'.join(new_namelist))
                tclfile.close()

                with open(os.path.join(dir_path, f"test_case{i}", "exfiltration.tcl"), 'w') as extclfile:
                    extclfile.write('\n'.join(exnew_namelist))
                extclfile.close()
                i += 1
    return

def plot_pressure_profile():
    plt.rcParams.update({'font.size': 22})
    # plt.figure(figsize=(16, 9))
    fig, ax = plt.subplots(figsize=(16, 9))
    name=os.path.join(DIRPATH, "pressureprofileexample", "infiltration")
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
    fig.savefig(os.path.join(DIRPATH, "pressureprofileexample", "pressure_profile.png"))

def infexf_lambda_dependence():
    exf_cases_path = os.path.join(DIRPATH, "inputs", "drainage_inf_testcases_toplayer.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(exf_cases_path)

    d_column = "d"
    inf_t_column = "inf_time"
    exf_t_column = "exf_time"
    toplayer_pressure_column = "toplayer_pressure"

    # Prepare figure
    plt.figure(figsize=(16, 9))
    
    # Plot each group with different color
    for d_val, group in df.groupby(d_column):
        y = group[exf_t_column] / group[inf_t_column]
        x = group[d_column] - 0.05 - abs(group[toplayer_pressure_column])
        plt.scatter(x, y, label=f"d = {d_val}")

    # Plot formatting
    plt.grid(True)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("λ (m)")  # α
    plt.ylabel(r"$SST_{d} /SST_{i}  (-)$")
    plt.legend(title="d values")
    
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, "exfinf_lambda.png"), dpi=300)
    plt.close()

def infexf_dependence():    
    exf_cases_path = os.path.join(DIRPATH, "inputs", "drainage_inf_testcases_toplayer.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(exf_cases_path)

    d_column = "d"
    inf_t_column = "inf_time"
    exf_t_column = "exf_time"
    toplayer_pressure_column = "toplayer_pressure"

    # Create figure and twin y-axis
    fig, ax1 = plt.subplots(figsize=(16, 9))
    ax2 = ax1.twinx()
    # Create a color map for different d values
    unique_d = sorted(df[d_column].unique())
    colors = plt.cm.viridis(np.linspace(0, 1, len(unique_d)))
    # For building a combined legend
    handles = []
    labels = []

    # Plot each group with a different color
    for color, d_val in zip(colors, unique_d):
        group = df[df[d_column] == d_val]
        x = group[d_column] - 0.05 - abs(group[toplayer_pressure_column])
        yd = group[exf_t_column]
        yi = group[inf_t_column]
        # Plot exf_time (primary y-axis)
        s1 = ax1.scatter(x, yd, color=color, alpha=0.6, label=fr"$SST_d$ (d={d_val} m)")
        # Plot inf_time (secondary y-axis)
        s2 = ax2.scatter(x, yi, color=color, alpha=0.4, marker='x', label=fr"$SST_i$ (d={d_val} m)")
        handles.append((s1, s2))
        labels.append(fr"d={d_val}")

    ax1.set_xlabel("λ (m)")
    ax1.set_ylabel(r"$SST_{d}$ (hr)")#, color='blue')
    ax1.set_yscale('log')
    ax1.grid(True, linestyle='--', alpha=0.3)
    ax2.set_ylabel(r"$SST_{i}$ (hr)")#, color='green')
    ax2.set_yscale('log')
    # Legend
    ax1.legend(handles, labels, loc='lower right', fontsize=8, handler_map={tuple: HandlerTuple(ndivide=None)})
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, "exf_vs_inf.png"), dpi=300)
    plt.close()

def dexf_dependence():
    exf_cases_path = os.path.join(DIRPATH, "inputs", "drainage_inf_testcases_toplayer.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(exf_cases_path)

    d_column = "d"
    exf_t_column = "exf_time"

    # Prepare data
    x = df[d_column]
    y = df[exf_t_column]

    # Create plot
    plt.figure()
    plt.scatter(x, y, color='black')  # no coloring by group
    plt.grid(True)
    plt.yscale("log")
    plt.xlabel("d (m)")
    plt.ylabel(r"$SST_{d}$ (hr)")
    plt.tight_layout()

    # Save plot
    plt.savefig(os.path.join(output_path, r"tdr_vs_d.png"), dpi=300)
    plt.close()

def dinf_dependence():
    exf_cases_path = os.path.join(DIRPATH, "inputs", "drainage_inf_testcases_toplayer.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(exf_cases_path)

    d_column = "d"
    inf_t_column = "inf_time"

    # Prepare data
    x = df[d_column]
    y = df[inf_t_column]

    # Create plot
    plt.figure()
    plt.scatter(x, y, color='black')  # no coloring by group
    plt.grid(True)
    plt.yscale("log")
    plt.xlabel("d (m)")
    plt.ylabel(r"$SST_{i}$ (hr)")
    plt.tight_layout()

    # Save plot
    plt.savefig(os.path.join(output_path, r"tinf_vs_d.png"), dpi=300)
    plt.close()

def fitting_func(xlist, ylist):
    # linearize
    y_lin = np.log(ylist)
    x_lin = np.log(xlist)
    # Fit the function
    params, covariance = curve_fit(linear_law, x_lin, y_lin)
    a_fit, b_fit = params

    # fitting accuracy
    y_fit = [linear_law(x, a_fit, b_fit) for x in x_lin]
    R_square = r2_score(y_lin, y_fit)
    r = np.corrcoef(y_lin, y_fit)
    r2 = r[0][1]**2

    # back transform to power law
    a_fit = np.exp(a_fit)

    x_fit = [min(xlist), max(xlist)]
    y_fit = [powerlaw_func(x, a_fit, b_fit) for x in x_fit]
    return x_fit, y_fit, a_fit, b_fit, r2

def plot_exfsolution_1(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    x = alfa*d
    y = (exf_t*k*alfa)/(theta_s-theta_r)
    xlabel = r'$\alpha\cdot d (-)$'
    ylabel = r'$\frac{SST_{d}\cdot K_s\cdot \alpha}{\theta_s-\theta_r} (-)$'
    return x, y, xlabel, ylabel

def plot_exfsolution_2(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    x = alfa*d
    y = ((d-0.05)-abs(toplayer_pressure))/(exf_t*k)
    xlabel = r'$\alpha\cdot d (-)$'
    ylabel = r'$\frac{λ}{SST_{d}\cdot K_s} (-)$'
    return x, y, xlabel, ylabel

def exf_solution_plots():
    exf_cases_path = os.path.join(DIRPATH, "inputs", "drainage_inf_testcases_toplayer.csv")
    output_path = os.path.join(DIRPATH, "outputs")

    q_column = "q"
    k_column = "k"
    d_column = "d"
    n_column = "n"
    alfa_column = "alfa"
    theta_r_column = "theta_r"
    theta_s_column = "theta_s"
    inf_t_column = "inf_time"
    exf_t_column = "exf_time"
    toplayer_pressure_column = "toplayer_pressure"

    df = pd.read_csv(exf_cases_path)

    soil_types = list(df[k_column].unique())
    soil_types.sort()
    soil_types_colors = plt.cm.jet(np.linspace(0,1,len(soil_types)))
    colors_dic = {soil_types[i]:soil_types_colors[i] for i in range(len(soil_types))}

    list_markers = ["o", "^", "s", "P", "*", "X", "d", "p", "2", r"$\clubsuit$", (5,2), "x"]
    print(len(list_markers))
    print(len(soil_types))
    print(soil_types)
    markers = {soil_types[i]:list_markers[i] for i in range(len(soil_types))}
    print(markers)

    ax = plt
    ax.figure(figsize=(22,9))

    x_all = []
    y_all = []
    colors = []
    for soil in soil_types:
        x_axis = []
        y_axis = []
        df_soil = df[df[k_column]==soil]
        for i in range(len(df_soil)):
            q = df_soil[q_column].iloc[i]
            k = df_soil[k_column].iloc[i]
            d = df_soil[d_column].iloc[i]
            n = df_soil[n_column].iloc[i]
            alfa = df_soil[alfa_column].iloc[i]
            theta_r = df_soil[theta_r_column].iloc[i]
            theta_s = df_soil[theta_s_column].iloc[i]
            inf_t = df_soil[inf_t_column].iloc[i]
            exf_t = df_soil[exf_t_column].iloc[i]
            toplayer_pressure = df_soil[toplayer_pressure_column].iloc[i]
            if k>=q and exf_t!=0:# and k!=0.0025 and k!=0.0026 and k!=0.0012 and k!=0.002:
                x, y, xlabel, ylabel = plot_exfsolution_2(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure)
                y_axis.append(y)
                x_axis.append(x)
                
                ax.scatter(x, y, c=colors_dic[soil], cmap='jet',s=80, marker=markers[k])
                # Add text labels for each point
                #labels = f"q={q:.2e}"
                #plt.text(x, y, labels, fontsize=9, ha='right', va='bottom')  # Adjust alignment as needed


        x_all.extend(x_axis)
        y_all.extend(y_axis)
        colors.extend([soil]*len(x_axis))
    for k in markers.keys():
        ax.scatter([],[], c=colors_dic[k],s=80, marker=markers[k], label=f"{np.around(k, 4)}")

    ax.legend(title="Ks (m/hr)", loc='upper right', bbox_to_anchor=(1.01, 0.8))
    # ax.legend(title="Ks (m/hr)", loc='upper right', bbox_to_anchor=(1.008, 1.018))

    # linearize
    y_lin = np.log(y_all)
    x_lin = np.log(x_all)
    # Fit the function
    params, covariance = curve_fit(linear_law, x_lin, y_lin)
    a_fit, b_fit = params

    # fitting accuracy
    y_fit = [linear_law(x, a_fit, b_fit) for x in x_lin]
    R_square = r2_score(y_lin, y_fit)
    print(f"R2 = {R_square}")
    r = np.corrcoef(y_lin, y_fit)
    r2 = r[0][1]**2

    # back transform to power law
    a_fit = np.exp(a_fit)
    print(f"Fitted a: {a_fit} and b:{b_fit}")

    x_fit = [min(x_all), max(x_all)]
    y_fit = [powerlaw_func(x, a_fit, b_fit) for x in x_fit]
    plt.plot(x_fit, y_fit, color="k", label="fitted line", linewidth=5)
    #plt.annotate(f"R²={round(r2,2)}\nf(x)={round(a_fit, 2)}x^({round(b_fit, 2)})", xy=(min(x_all), max(y_all)/10), color="black")
    print("annotations:")
    print(f"r2:{round(r2,2)}")
    print(f"f(x)={round(a_fit, 2)}x^({round(b_fit, 2)}")
    ax.grid(True)
    ax.xscale("log")
    ax.yscale("log")
    ax.xlabel(f"{xlabel}", fontsize=32)
    ax.ylabel(f"{ylabel}", fontsize=32)
    #ax.savefig(os.path.join(output_path, f"exf_ad_vs_tka-thetasr_symbolstest.png"))
    ax.savefig(os.path.join(output_path, f"exf_ad_vs_v-k_symbolstest.png"))

def plot_infsolution(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    x = alfa*d*q/k
    y = ((d-0.05)-abs(toplayer_pressure))/(inf_t)
    xlabel = r'$\frac{\alpha\cdot d\cdot q}{K_s} (-)$'
    ylabel = r'$\frac{λ}{SST_{i}} (m/hr)$'
    return x, y, xlabel, ylabel

def inf_solution_plots():
    inf_cases_path = os.path.join(DIRPATH, "inputs", "infiltration_pressure_index.csv")
    output_path = os.path.join(DIRPATH, "outputs")

    q_column = "q"
    k_column = "k"
    d_column = "d"
    n_column = "n"
    alfa_column = "alfa"
    theta_r_column = "theta_r"
    theta_s_column = "theta_s"
    inf_t_column = "time"
    toplayer_pressure_column = "toplayer_pressure"

    df = pd.read_csv(inf_cases_path)

    soil_types = list(df[k_column].unique())
    soil_types.sort()
    soil_types_colors = plt.cm.jet(np.linspace(0,1,len(soil_types)))
    colors_dic = {soil_types[i]:soil_types_colors[i] for i in range(len(soil_types))}
    list_markers = ["o", "^", "s", "P", "*", "X", "d", "p", "2", r"$\clubsuit$", (5,2), "x"]

    print(len(soil_types))
    markers = {soil_types[i]:list_markers[i] for i in range(len(soil_types))}

    fig, ax = plt.subplots(figsize=(16, 9))
    #ax.figure(figsize=(24,11))
    print("plotting points...")
    x_all = []
    y_all = []
    colors = []
    for soil in soil_types:
        x_axis = []
        y_axis = []
        df_soil = df[df[k_column]==soil]
        for i in range(len(df_soil)):
            q = df_soil[q_column].iloc[i]
            k = df_soil[k_column].iloc[i]
            d = df_soil[d_column].iloc[i]
            n = df_soil[n_column].iloc[i]
            alfa = df_soil[alfa_column].iloc[i]
            theta_r = df_soil[theta_r_column].iloc[i]
            theta_s = df_soil[theta_s_column].iloc[i]
            inf_t = df_soil[inf_t_column].iloc[i]
            toplayer_pressure = df_soil[toplayer_pressure_column].iloc[i]
            if k>=q:# and k!=0.0025 and k!=0.0026 and k!=0.0012 and k!=0.002:# and (q/k)>=0.001:
                x, y, xlabel, ylabel = plot_infsolution(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)
                y_axis.append(y)
                x_axis.append(x)

                ax.scatter(x, y, c=colors_dic[soil], cmap='jet',s=80, marker=markers[k])#, label=k)

        x_all.extend(x_axis)
        y_all.extend(y_axis)
        colors.extend([soil]*len(x_axis))
    for k in markers.keys():
        # change logarithmic colors for every marker
        ax.scatter([],[], c=colors_dic[k],s=80, marker=markers[k], label=f"{np.around(k, 4)}")

    #ax.legend(title="Ks (m/hr)", loc='center right', bbox_to_anchor=(1.135, 0.5))
    ax.legend(title="Ks (m/hr)", loc='center right', bbox_to_anchor=(1.01, 0.41))
    #ax.tight_layout()

    # linearize
    y_lin = np.log(y_all)
    x_lin = np.log(x_all)
    # Fit the function
    params, covariance = curve_fit(linear_law, x_lin, y_lin)
    a_fit, b_fit = params

    # fitting accuracy
    y_fit = [linear_law(x, a_fit, b_fit) for x in x_lin]
    R_square = r2_score(y_lin, y_fit)
    print(f"R2 = {R_square}")
    r = np.corrcoef(y_lin, y_fit)
    r2 = r[0][1]**2
    print(f"pearson r2: {r2}")

    # back transform to power law
    a_fit = np.exp(a_fit)
    print(f"Fitted a: {a_fit} and b:{b_fit}")

    # check kolmogorov-smirnov
    #anderson_darling(y_all)
    #KS_fit(y_all, y_fit)

    x_fit = [min(x_all), max(x_all)]
    y_fit = [powerlaw_func(x, a_fit, b_fit) for x in x_fit]
    plt.plot(x_fit, y_fit, color="k", label="fitted line", linewidth=5)
    #plt.annotate(f"R²={round(r2,2)}\nf(x)={round(a_fit, 2)}x^({round(b_fit, 2)})", xy=(min(x_all), max(y_all)/10), color="black")
    print("annotations:")
    print(f"r2:{round(r2,2)}")
    print(f"f(x)={round(a_fit, 2)}x^({round(b_fit, 2)}")

    #ax.scatter(x_all, y_all,c=colors, cmap='jet', s=30, norm=matplotlib.colors.LogNorm(vmin=kmin, vmax=kmax))
    ax.grid(True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(f"{xlabel}", fontsize=32)
    ax.set_ylabel(f"{ylabel}", fontsize=32)
    #plt.gcf().subplots_adjust(bottom=0.4)
    plt.tight_layout()

    fig.savefig(os.path.join(output_path, "inf_adq-k_vs_v_symbolstest.png"))


def qfit():
    exf_cases_path = os.path.join(DIRPATH, "inputs", "exf_intercept.csv")
    output_path = os.path.join(DIRPATH, "outputs")

    q_column = "q"
    k_column = "k"
    d_column = "d"
    n_column = "n"
    alfa_column = "alfa"
    theta_r_column = "theta_r"
    theta_s_column = "theta_s"
    inf_t_column = "inf_time"
    exf_t_column = "exf_time"
    toplayer_pressure_column = "toplayer_pressure"

    df = pd.read_csv(exf_cases_path)

    soil_types = list(df[k_column].unique())
    soil_types.sort()

    markers = {0.002:"P", 0.0045:"d"}
    markers = {0.0045:"d"}

    ax = plt
    ax.figure(figsize=(16,9))

    df["x"] = df[d_column]*df[alfa_column]
    df["y"] = ((df[d_column]-0.05)-abs(df[toplayer_pressure_column]))/(df[exf_t_column]*df[k_column])
    # df["y"] = ((df[d_column]-0.05)-abs(df[toplayer_pressure_column]))#df[exf_t_column]
    alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"]
    ann_ab = {"0.00010.002": [r'$0.225x^{1.356}$', (0.73, 0.234386417)],
            "0.000150.002": [r'$0.234x^{1.385}$', (0.73, 0.338135512)],
            "0.0010.0045": [r'$0.304x^{1.717}$', (2, 0.09250808557112719)],
            "0.00010.0045": [r'$0.125x^{1.306}$', (2, 0.05052146488857635)],
            "0.000150.0045": [r'$0.156x^{1.409}$', (2, 0.05871341875690252)],
            "0.00050.0045": [r'$0.253x^{1.634}$', (2, 0.08165277976910569)]
            }
    ind=0
    colors = ["r", "g", "b", "darkorange"]
    for soil in soil_types:
        if soil == 0.002:
            continue
        x_axis = []
        y_axis = []
        df_soil = df[df[k_column]==soil]
        ax.scatter(df_soil["x"], df_soil["y"], c="k",s=80, marker=markers[soil])

        for q in df_soil[q_column].unique():
            df_q = df_soil[df_soil[q_column]==q]
            x_fit, y_fit, a_fit, b_fit, r2 = fitting_func(df_q["x"], df_q["y"])
            print(df_q["x"])
            print(type(df_q["x"]))
            #ax.plot(x_fit, y_fit, color="k", linewidth=1)
            ax.plot(x_fit, y_fit, color=colors[ind] ,linewidth=1)
            print(f"{q}{soil} this is {alphabet[ind]} and minx={min(x_fit)} and maxy={max(y_fit)}")
            #ax.annotate(f"R²={round(r2,2)}\nf(x)={round(a_fit, 2)}x^({round(b_fit, 2)})", xy=(min(x_fit), max(y_fit)), color="black")
            anx = min(x_fit)-0.1 if soil==0.0045 else min(x_fit)-0.03
            anx = min(x_fit) if soil==0.0045 else min(x_fit)-0.03
            any = max(y_fit)-0.05 if soil==0.002 and q==0.0001 else max(y_fit)
            #ax.annotate(ann_ab[f"{q}{soil}"][0], xy=ann_ab[f"{q}{soil}"][1], color="black", fontsize=15)
            ax.annotate(f"y={ann_ab[f'{q}{soil}'][0]}", xy=ann_ab[f"{q}{soil}"][1],color=colors[ind], fontsize=15)
            ax.annotate(f"q={q}", xy=(df_q["x"].to_list()[1]+0.1, df_q["y"].to_list()[1]-0.002),color=colors[ind], fontsize=15)
            ind += 1

    for k in markers.keys():
        ax.scatter([],[], c="k",s=80, marker=markers[k], label=f"{np.around(k, 4)}")

    ax.legend(title="Ks (m/hr)", loc='upper right')#, bbox_to_anchor=(1.135, 0.5))


    ax.grid(True)
    ax.xscale("log")
    ax.yscale("log")
    xlabel = r'$\alpha\cdot d (-)$'
    # xlabel = r'$d(m)$'
    ylabel = r'$\frac{λ}{SST_{d}\cdot K_s} (-)$'
    # ylabel = r'$d-|\psi_{toplayer}| (m)$'
    ax.xlabel(f"{xlabel}")
    ax.ylabel(f"{ylabel}")
    ax.savefig(os.path.join(output_path, "exf_intercept_qfit.png"))


def q_vs_slope():
    input_path = os.path.join(DIRPATH, "inputs", "q_vs_fittinglinesslope2.csv")
    output_path = os.path.join(DIRPATH, "outputs")


    # Extracting the data
    df = pd.read_csv(input_path)
    q = df['q (m/hr)'].values
    slope = df['fitting line slope'].values
    print(q)
    print(type(q))
    # Define a logarithmic model for fitting
    def log_model(x, a, b):
        return a * np.log(x) + b

    # Perform curve fitting
    popt, pcov = curve_fit(log_model, q, slope)

    # Get the fitted values
    slope_fitted = log_model(q, *popt)

    # Calculate R-squared for the fit
    r_squared = r2_score(slope, slope_fitted)

    # Plot the data points and the fitting curve
    plt.figure(figsize=(16,9))
    plt.scatter(q, slope, color='k', label='Data Points')
    plt.plot(q, slope_fitted, color='k', label=f'y = {popt[0]:.2f}ln(Δq) {popt[1]:.2f}\n$R^2$ = {r_squared:.2f}')
    plt.xscale('log')  # since we're dealing with logarithmic fitting
    plt.xlabel('Δq (m/hr)')
    plt.ylabel('Power laws slope')
    #plt.title('Logarithmic Fit of Fitting Line Slope vs. q')
    plt.legend()
    plt.grid(True)

    # Show plot
    plt.savefig(os.path.join(output_path, "q_vs_slope2.png"))
