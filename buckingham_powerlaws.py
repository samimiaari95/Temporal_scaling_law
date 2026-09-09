import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cmocean
from analysis import fitting_func
from utils import powerlaw_func
from file_io import DIRPATH

filepath = os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv")

# plt.rcParams.update({'font.size': 22})

filepath = os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv")

def infsolution(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    x = q/(k*alfa*d*n*theta_s)
    y = k * inf_t / d
    xlabel = r'$\frac{q}{\alpha\cdot d\cdot K_s\cdot \theta_s\cdot n} (-)$'
    ylabel = r'$\frac{K_s\cdot t_{inf}}{d} (-)$'
    return x, y, xlabel, ylabel


def exfsolution(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    x = k/(alfa * d*n*theta_s*q)
    y = q*exf_t/d
    xlabel = r'$\frac{K_s}{q\cdot \alpha\cdot d\cdot \phi\cdot n} (-)$'
    ylabel = r'$\frac{q\cdot t_{exf}}{d} (-)$'
    return x, y, xlabel, ylabel

############################################## infiltration solutions here ###########################################

def infsolution_1(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    x = q/(k)
    y = k * inf_t / d
    xlabel = r'$\frac{q}{K_s} (-)$'
    ylabel = r'$\frac{K_s\cdot t_{inf}}{d} (-)$'
    name = "1g1p1_p2"
    return x, y, xlabel, ylabel, name

def infsolution_2(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    x = q*alfa*d/(k)
    y = k * inf_t / d
    xlabel = r'$\frac{q\cdot \alpha\cdot d}{K_s} (-)$'
    ylabel = r'$\frac{K_s\cdot t_{inf}}{d} (-)$'
    name = "2g1p1_p2p3"
    return x, y, xlabel, ylabel, name

def infsolution_3(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    x = q/(k*alfa*d)
    y = k * inf_t / d
    xlabel = r'$\frac{q}{K_s\cdot \alpha\cdot d} (-)$'
    ylabel = r'$\frac{K_s\cdot t_{inf}}{d} (-)$'
    name = "3g1p1_p2-p3"
    return x, y, xlabel, ylabel, name

def infsolution_4(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    x = q/(k*alfa*d*n)
    y = k * inf_t / d
    xlabel = r'$\frac{q}{K_s\cdot \alpha\cdot d\cdot n} (-)$'
    ylabel = r'$\frac{K_s\cdot t_{inf}}{d} (-)$'
    name = "4g1p1_p2-p3p4"
    return x, y, xlabel, ylabel, name

def infsolution_5(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    x = q*alfa*d*n*theta_s/(k)
    y = k * inf_t / d
    xlabel = r'$\frac{q\cdot \alpha\cdot d\cdot n\cdot \phi}{K_s} (-)$'
    ylabel = r'$\frac{K_s\cdot t_{inf}}{d} (-)$'
    name = "5g1p1_p2p3p4p5"
    return x, y, xlabel, ylabel, name

def infsolution_6(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    x = q*n*theta_s/(k*alfa*d)
    y = k * inf_t / d
    xlabel = r'$\frac{q\cdot n\cdot \phi}{K_s\cdot \alpha\cdot d} (-)$'
    ylabel = r'$\frac{K_s\cdot t_{inf}}{d} (-)$'
    name = "6g1p1_p2p4p5-p3"
    return x, y, xlabel, ylabel, name

def infsolution_7(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    x = k/(q)
    y = q * inf_t / d
    xlabel = r'$\frac{K_s}{q} (-)$'
    ylabel = r'$\frac{q\cdot t_{inf}}{d} (-)$'
    name = "7g2p1_p2"
    return x, y, xlabel, ylabel, name

def infsolution_8(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    x = k/(q*alfa*d*theta_s*n)
    y = q * inf_t / d
    xlabel = r'$\frac{K_s}{q\cdot \alpha\cdot d\cdot n\cdot \phi} (-)$'
    ylabel = r'$\frac{q\cdot t_{inf}}{d} (-)$'
    name = "8g2p1_p2-p3p4p5"
    return x, y, xlabel, ylabel, name

############################################# exfiltration solutions here ###########################################

def exfsolution_1(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    x = q/(k)
    y = k * exf_t / d
    xlabel = r'$\frac{q}{K_s} (-)$'
    ylabel = r'$\frac{K_s\cdot t_{exf}}{d} (-)$'
    name = "1g1p1_p2"
    return x, y, xlabel, ylabel, name

def exfsolution_2(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    x = q*alfa*d/(k)
    y = k * exf_t / d
    xlabel = r'$\frac{q\cdot \alpha\cdot d}{K_s} (-)$'
    ylabel = r'$\frac{K_s\cdot t_{exf}}{d} (-)$'
    name = "2g1p1_p2p3"
    return x, y, xlabel, ylabel, name

def exfsolution_3(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    x = q/(k*alfa*d)
    y = k * exf_t / d
    xlabel = r'$\frac{q}{K_s\cdot \alpha\cdot d} (-)$'
    ylabel = r'$\frac{K_s\cdot t_{exf}}{d} (-)$'
    name = "3g1p1_p2-p3"
    return x, y, xlabel, ylabel, name

def exfsolution_4(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    x = q/(k*alfa*d*n)
    y = k * exf_t / d
    xlabel = r'$\frac{q}{K_s\cdot \alpha\cdot d\cdot n} (-)$'
    ylabel = r'$\frac{K_s\cdot t_{exf}}{d} (-)$'
    name = "4g1p1_p2-p3p4"
    return x, y, xlabel, ylabel, name

def exfsolution_5(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    x = q/(k*alfa*d*n*theta_s)
    y = k * exf_t / d
    xlabel = r'$\frac{q}{K_s\cdot \alpha\cdot d\cdot n\cdot \phi} (-)$'
    ylabel = r'$\frac{K_s\cdot t_{exf}}{d} (-)$'
    name = "5g1p1_p2-p3p4p5"
    return x, y, xlabel, ylabel, name

def exfsolution_6(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    x = q*n*theta_s/(k*alfa*d)
    y = k * exf_t / d
    xlabel = r'$\frac{q\cdot n\cdot \phi}{K_s\cdot \alpha\cdot d} (-)$'
    ylabel = r'$\frac{K_s\cdot t_{exf}}{d} (-)$'
    name = "6g1p1_p2p4p5-p3"
    return x, y, xlabel, ylabel, name

def exfsolution_7(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    x = k/(q)
    y = q * exf_t / d
    xlabel = r'$\frac{K_s}{q} (-)$'
    ylabel = r'$\frac{q\cdot t_{exf}}{d} (-)$'
    name = "7g2p1_p2"
    return x, y, xlabel, ylabel, name

def exfsolution_8(q, k, d, n, alfa, theta_r, theta_s, inf_t, exf_t, toplayer_pressure):
    x = k*alfa*d*theta_s*n/(q)
    y = q * exf_t / d
    xlabel = r'$\frac{K_s\cdot \alpha\cdot d\cdot n\cdot \phi}{q} (-)$'
    ylabel = r'$\frac{q\cdot t_{exf}}{d} (-)$'
    name = "8g2p1_p2p3p4p5"
    return x, y, xlabel, ylabel, name

############################################## inf combinations by chatgbt ###########################################

def inf_dKs_t_x1(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    x = q / (k*alfa*d*n*(theta_s-theta_r))
    y = k * inf_t / d
    xlabel = r"$\frac{q'}{d'\cdot n'} (-)$"
    ylabel = r"$K_s\cdot t_{i}' (-)$"
    name = "1_kt_qkalfadntheta"
    return x, y, xlabel, ylabel, name

def inf_dKs_t_x2(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    x = alfa * d
    y = k * inf_t / d
    xlabel = r"$d'$"
    ylabel = r"$K_s\cdot t_{i}'$"
    name = "2_dKs_alpha_d_vs_Ks_t_d"
    return x, y, xlabel, ylabel, name

def inf_dKs_t_x3(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = lam / d
    y = k * inf_t / d
    xlabel = r'$\lambda/d$'
    ylabel = r'$\frac{K_s t_{inf}}{d}$'
    name = "3_dKs_lambda_d_vs_Ks_t_d"
    return x, y, xlabel, ylabel, name

def inf_dKs_t_x4(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    x = (q / k) * (alfa * d)
    y = k * inf_t / d
    xlabel = r'$(q/K_s)\cdot(\alpha d)$'
    ylabel = r'$\frac{K_s t_{inf}}{d}$'
    name = "4_dKs_qKs_alpha_d"
    return x, y, xlabel, ylabel, name

def inf_dKs_t_x5(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (q / k) / (lam / d)
    y = k * inf_t / d
    xlabel = r'$(q/K_s)/(\lambda/d)$'
    ylabel = r'$\frac{K_s t_{inf}}{d}$'
    name = "5_dKs_qKs_over_lambda_d"
    return x, y, xlabel, ylabel, name

def inf_dKs_t_x6(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (alfa * d) * (lam / d)
    y = k * inf_t / d
    xlabel = r'$(\alpha d)(\lambda/d)$'
    ylabel = r'$\frac{K_s t_{inf}}{d}$'
    name = "6_dKs_alpha_d_lambda_d"
    return x, y, xlabel, ylabel, name

def inf_dKs_t_x7(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (q / k) * (alfa * d) / (lam / d)
    y = k * inf_t / d
    xlabel = r'$(q/K_s)(\alpha d)/(\lambda/d)$'
    ylabel = r'$\frac{K_s t_{inf}}{d}$'
    name = "7_dKs_qKs_alpha_d_lambda_d"
    return x, y, xlabel, ylabel, name

def inf_KsLam_t_x1(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = q / k
    y = k * inf_t / lam
    xlabel = r'$\frac{q}{K_s}$'
    ylabel = r'$\frac{K_s t_{inf}}{\lambda}$'
    name = "8_KsLam_qKs"
    return x, y, xlabel, ylabel, name

def inf_KsLam_t_x2(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = alfa * lam
    y = k * inf_t / lam
    xlabel = r'$\alpha \lambda$'
    ylabel = r'$\frac{K_s t_{inf}}{\lambda}$'
    name = "9_KsLam_alpha_lam"
    return x, y, xlabel, ylabel, name

def inf_KsLam_t_x3(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = d / lam
    y = k * inf_t / lam
    xlabel = r'$\frac{d}{\lambda}$'
    ylabel = r'$\frac{K_s t_{inf}}{\lambda}$'
    name = "10_KsLam_d_lam"
    return x, y, xlabel, ylabel, name

def inf_KsLam_t_x4(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (q / k) * (d / lam)
    y = k * inf_t / lam
    xlabel = r'$\frac{q\cdot d}{K_s\cdot \lambda}$'
    ylabel = r'$\frac{K_s t_{inf}}{\lambda}$'
    name = "11_KsLam_qKs_d_lam"
    return x, y, xlabel, ylabel, name

def inf_KsLam_t_x5(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (alfa * lam) / (q / k)
    y = k * inf_t / lam
    xlabel = r'$\frac{\alpha \lambda\cdot q}{K_s}$'
    ylabel = r'$\frac{K_s t_{inf}}{\lambda}$'
    name = "12_KsLam_alpha_over_qKs"
    return x, y, xlabel, ylabel, name

def inf_KsLam_t_x6(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (q / k) * (alfa * lam) / (d / lam)
    y = k * inf_t / lam
    xlabel = r'$\frac{q\cdot \alpha \lambda^2}{K_s\cdot d}$'
    ylabel = r'$\frac{K_s t_{inf}}{\lambda}$'
    name = "13_KsLam_full_combo"
    return x, y, xlabel, ylabel, name

def inf_qLam_t_x1(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = q / k
    y = lam / (q * inf_t)
    xlabel = r'$\frac{q}{K_s}$'
    ylabel = r'$\frac{q t_{inf}}{\lambda}$'
    name = "14_qLam_qKs"
    return x, y, xlabel, ylabel, name

def inf_qLam_t_x2(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = alfa * lam
    y = (q * inf_t) / lam
    xlabel = r'$\alpha \lambda$'
    ylabel = r'$\frac{q t_{inf}}{\lambda}$'
    name = "15_qLam_alpha_lam"
    return x, y, xlabel, ylabel, name

def inf_qLam_t_x3(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = d / lam
    y = (q * inf_t) / lam
    xlabel = r'$\frac{d}{\lambda}$'
    ylabel = r'$\frac{q t_{inf}}{\lambda}$'
    name = "16_qLam_d_lam"
    return x, y, xlabel, ylabel, name

def inf_qLam_t_x4(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (q / k) * (d / lam)
    y = (q * inf_t) / lam
    xlabel = r'$\frac{q\cdot d}{K_s\cdot \lambda}$'
    ylabel = r'$\frac{q t_{inf}}{\lambda}$'
    name = "17_qLam_qKs_d_lam"
    return x, y, xlabel, ylabel, name

def inf_qLam_t_x5(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (alfa * lam) * (q / k)
    y = (q * inf_t) / lam
    xlabel = r'$\frac{\alpha \lambda\cdot q}{K_s}$'
    ylabel = r'$\frac{q t_{inf}}{\lambda}$'
    name = "18_qLam_alpha_qKs"
    return x, y, xlabel, ylabel, name

def inf_qLam_t_x6(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (q / k) * (alfa * lam) / (d / lam)
    y = (q * inf_t) / lam
    xlabel = r'$\frac{(q\cdot \alpha \lambda^2)}{K_s\cdot d}$'
    ylabel = r'$\frac{q t_{inf}}{\lambda}$'
    name = "19_qLam_full_combo"
    return x, y, xlabel, ylabel, name

def inf_v_x1(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = alfa*d
    y = lam / (q * inf_t)
    xlabel = r'$\alpha\cdot d$'
    ylabel = r'$\frac{v_{inf}}{q}$'
    name = "20_v_alpha_d"
    return x, y, xlabel, ylabel, name

def inf_v_x2(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = alfa*d
    y = lam / (k * inf_t)
    xlabel = r'$\alpha\cdot d$'
    ylabel = r'$\frac{v_{inf}}{K_s}$'
    name = "21_v_alpha_d"
    return x, y, xlabel, ylabel, name

def run_infsolution_case(case_id,
                         q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure):
    if case_id == 1:
        return inf_dKs_t_x1(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 2:
        return inf_dKs_t_x2(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 3:
        return inf_dKs_t_x3(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 4:
        return inf_dKs_t_x4(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 5:
        return inf_dKs_t_x5(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 6:
        return inf_dKs_t_x6(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 7:
        return inf_dKs_t_x7(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    # -----------------------------
    # Ks, lambda group
    # -----------------------------

    if case_id == 8:
        return inf_KsLam_t_x1(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 9:
        return inf_KsLam_t_x2(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 10:
        return inf_KsLam_t_x3(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 11:
        return inf_KsLam_t_x4(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 12:
        return inf_KsLam_t_x5(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 13:
        return inf_KsLam_t_x6(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    # -----------------------------
    # q, lambda group
    # -----------------------------

    elif case_id == 14:
        return inf_qLam_t_x1(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 15:
        return inf_qLam_t_x2(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 16:
        return inf_qLam_t_x3(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 17:
        return inf_qLam_t_x4(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 18:
        return inf_qLam_t_x5(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 19:
        return inf_qLam_t_x6(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)
    
    elif case_id == 20:
        return inf_v_x1(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    elif case_id == 21:
        return inf_v_x2(q, k, d, n, alfa, theta_r, theta_s, inf_t, toplayer_pressure)

    else:
        raise ValueError("case_id must be between 1 and 21")


############################################## exf combinations by chatgbt ###########################################

def exf_dKs_t_x1(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    x = alfa*d
    y = (d-abs(toplayer_pressure)) / (k * exf_t)
    xlabel = r"$d' (-)$"
    ylabel = r'$\frac{v_{d}}{K_s} (-)$'
    name = "1_vKs_alpha_d"
    return x, y, xlabel, ylabel, name

def exf_dKs_t_x2(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    x = (q*alfa*d*n*(theta_s-theta_r))/k
    y = q * exf_t / d
    xlabel = r"$q'\cdot d'\cdot n'$ (-)"
    ylabel = r"$q\cdot t_{d}' (-)$"
    name = "2_qt_qkalfadntheta"
    return x, y, xlabel, ylabel, name

def exf_dKs_t_x3(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    x = alfa*d
    y = exf_t*k*alfa / (theta_s-theta_r)
    xlabel = r"$d' (-)$"
    ylabel = r'$\frac{t_{d}\cdot K_s\cdot \alpha}{\phi} (-)$'
    name = "3_tkalfa_d"
    return x, y, xlabel, ylabel, name

def exf_dKs_t_x4(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    x = alfa * d
    y = exf_t
    xlabel = r"$d' (-)$"
    ylabel = r"$t_{d} (m/hr)$"
    name = "4_t_d"
    return x, y, xlabel, ylabel, name

def exf_dKs_t_x5(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (q / k) / (lam / d)
    y = k * exf_t / d
    xlabel = r'$(q/K_s)/(\lambda/d)$'
    ylabel = r'$\frac{K_s t_{dr}}{d}$'
    name = "5_dKs_qKs_over_lambda_d"
    return x, y, xlabel, ylabel, name

def exf_dKs_t_x6(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (alfa * d) * (lam / d)
    y = k * exf_t / d
    xlabel = r'$(\alpha d)(\lambda/d)$'
    ylabel = r'$\frac{K_s t_{dr}}{d}$'
    name = "6_dKs_alpha_d_lambda_d"
    return x, y, xlabel, ylabel, name

def exf_dKs_t_x7(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (q / k) * (alfa * d) / (lam / d)
    y = k * exf_t / d
    xlabel = r'$(q/K_s)(\alpha d)/(\lambda/d)$'
    ylabel = r'$\frac{K_s t_{dr}}{d}$'
    name = "7_dKs_full_combo_1"
    return x, y, xlabel, ylabel, name

def exf_KsLam_t_x1(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = q / k
    y = k * exf_t / lam
    xlabel = r'$\frac{q}{K_s}$'
    ylabel = r'$\frac{K_s t_{dr}}{\lambda}$'
    name = "8_KsLam_qKs"
    return x, y, xlabel, ylabel, name

def exf_KsLam_t_x2(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = alfa * lam
    y = k * exf_t / lam
    xlabel = r'$\alpha \lambda$'
    ylabel = r'$\frac{K_s t_{dr}}{\lambda}$'
    name = "9_KsLam_alpha_lam"
    return x, y, xlabel, ylabel, name

def exf_KsLam_t_x3(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = d / lam
    y = k * exf_t / lam
    xlabel = r'$\frac{d}{\lambda}$'
    ylabel = r'$\frac{K_s t_{dr}}{\lambda}$'
    name = "10_KsLam_d_lam"
    return x, y, xlabel, ylabel, name

def exf_KsLam_t_x4(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (q / k) * (d / lam)
    y = k * exf_t / lam
    xlabel = r'$\frac{q\cdot d}{K_s\cdot \lambda}$'
    ylabel = r'$\frac{K_s t_{dr}}{\lambda}$'
    name = "11_KsLam_qKs_d_lam"
    return x, y, xlabel, ylabel, name

def exf_KsLam_t_x5(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (alfa * lam) / (q / k)
    y = k * exf_t / lam
    xlabel = r'$\frac{\alpha \lambda\cdot q}{K_s}$'
    ylabel = r'$\frac{K_s t_{dr}}{\lambda}$'
    name = "12_KsLam_alpha_over_qKs"
    return x, y, xlabel, ylabel, name

def exf_KsLam_t_x6(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (q / k) * (alfa * lam) / (d / lam)
    y = k * exf_t / lam
    xlabel = r'$\frac{q\cdot \alpha \lambda^2}{K_s\cdot d}$'
    ylabel = r'$\frac{K_s t_{dr}}{\lambda}$'
    name = "13_KsLam_full_combo"
    return x, y, xlabel, ylabel, name

def exf_qLam_t_x1(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = q / k
    y = lam / (q * exf_t)
    xlabel = r'$\frac{q}{K_s}$'
    ylabel = r'$\frac{q t_{dr}}{\lambda}$'
    name = "14_qLam_qKs"
    return x, y, xlabel, ylabel, name

def exf_qLam_t_x2(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = alfa * lam
    y = (q * exf_t) / lam
    xlabel = r'$\alpha \lambda$'
    ylabel = r'$\frac{q t_{dr}}{\lambda}$'
    name = "15_qLam_alpha_lam"
    return x, y, xlabel, ylabel, name

def exf_qLam_t_x3(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = d / lam
    y = (q * exf_t) / lam
    xlabel = r'$\frac{d}{\lambda}$'
    ylabel = r'$\frac{q t_{dr}}{\lambda}$'
    name = "16_qLam_d_lam"
    return x, y, xlabel, ylabel, name

def exf_qLam_t_x4(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (q / k) * (d / lam)
    y = (q * exf_t) / lam
    xlabel = r'$\frac{q\cdot d}{K_s\cdot \lambda}$'
    ylabel = r'$\frac{q t_{dr}}{\lambda}$'
    name = "17_qLam_qKs_d_lam"
    return x, y, xlabel, ylabel, name

def exf_qLam_t_x5(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (alfa * lam) * (q / k)
    y = (q * exf_t) / lam
    xlabel = r'$\frac{\alpha \lambda\cdot q}{K_s}$'
    ylabel = r'$\frac{q t_{dr}}{\lambda}$'
    name = "18_qLam_alpha_qKs"
    return x, y, xlabel, ylabel, name

def exf_qLam_t_x6(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = (q / k) * (alfa * lam) / (d / lam)
    y = (q * exf_t) / lam
    xlabel = r'$\frac{(q\cdot \alpha \lambda^2)}{K_s\cdot d}$'
    ylabel = r'$\frac{q t_{dr}}{\lambda}$'
    name = "19_qLam_full_combo"
    return x, y, xlabel, ylabel, name

def exf_v_x1(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = alfa*d
    y = lam / (q * exf_t)
    xlabel = r'$\alpha\cdot d$'
    ylabel = r'$\frac{v_{dr}}{q}$'
    name = "20_v_alpha_d"
    return x, y, xlabel, ylabel, name

def exf_v_x2(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):
    lam = d - abs(toplayer_pressure)
    x = alfa*d
    y = lam / (k * exf_t)
    xlabel = r'$\alpha\cdot d$'
    ylabel = r'$\frac{v_{dr}}{K_s}$'
    name = "21_v_alpha_d"
    return x, y, xlabel, ylabel, name

def run_exfsolution_case(case_id,
                         q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure):

    if case_id == 1:
        return exf_dKs_t_x1(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 2:
        return exf_dKs_t_x2(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 3:
        return exf_dKs_t_x3(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 4:
        return exf_dKs_t_x4(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 5:
        return exf_dKs_t_x5(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 6:
        return exf_dKs_t_x6(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 7:
        return exf_dKs_t_x7(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    # -----------------------------
    # Ks, lambda group
    # -----------------------------

    if case_id == 8:
        return exf_KsLam_t_x1(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 9:
        return exf_KsLam_t_x2(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 10:
        return exf_KsLam_t_x3(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 11:
        return exf_KsLam_t_x4(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 12:
        return exf_KsLam_t_x5(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 13:
        return exf_KsLam_t_x6(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    # -----------------------------
    # q, lambda group
    # -----------------------------

    elif case_id == 14:
        return exf_qLam_t_x1(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 15:
        return exf_qLam_t_x2(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 16:
        return exf_qLam_t_x3(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 17:
        return exf_qLam_t_x4(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 18:
        return exf_qLam_t_x5(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 19:
        return exf_qLam_t_x6(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)
    
    elif case_id == 20:
        return exf_v_x1(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)

    elif case_id == 21:
        return exf_v_x2(q, k, d, n, alfa, theta_r, theta_s, exf_t, toplayer_pressure)
    
    else:
        raise ValueError("case_id must be between 1 and 21")

def inf_solution_plots(solind):
    """
    Plots the non-dimensional infiltration velocity (lambda / SST_i)
    against the non-dimensional infiltration parameter (alpha * d * q / Ks).
    Includes power-law fit to the data.
    """
    inf_cases_path = os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(inf_cases_path)
    
    # Get unique soil types
    soil_types = sorted(df["k"].unique())
    
    # Create color map and markers
    soil_types_colors = cmocean.cm.balance(np.linspace(0, 1, len(soil_types)))
    colors_dic = {soil_types[i]: soil_types_colors[i] for i in range(len(soil_types))}
    
    list_markers = ["o", "^", "s", "P", "*", "X", "d", "p", "2", r"$\clubsuit$", (5,2), "x"]
    markers = {soil_types[i]: list_markers[i] for i in range(len(soil_types))}
    
    print(f"Number of soil types: {len(soil_types)}")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(7.09, 3.54))  # A4 size in inches (width, height)
    
    print("Plotting points...")
    x_all = []
    y_all = []
    print(f"Running infiltration solution case {solind}...")

    # Plot each soil type
    for soil in soil_types:
        df_soil = df[df["k"] == soil]
        
        for i in range(len(df_soil)):
            q = df_soil["q"].iloc[i]
            k = df_soil["k"].iloc[i]
            d = df_soil["d"].iloc[i]
            alfa = df_soil["alfa"].iloc[i]
            inf_t = df_soil["inf_time"].iloc[i]
            toplayer_pressure = df_soil["toplayer_pressure"].iloc[i]
            
            # Filter: only cases where k >= q
            if k >= q:
                # Compute non-dimensional variables
                x, y, xlabel, ylabel, name = run_infsolution_case(
                    solind,  # case_id
                    q, k, d, df_soil["n"].iloc[i], alfa,
                    df_soil["theta_r"].iloc[i], df_soil["theta_s"].iloc[i],
                    inf_t, toplayer_pressure
                )
                
                x_all.append(x)
                y_all.append(y)
                
                ax.scatter(x, y, c=colors_dic[soil], cmap='jet', s=10, marker=markers[k])
    
    # Create legend
    for k in markers.keys():
        ax.scatter([], [], c=colors_dic[k], s=10, marker=markers[k],
                  label=f"{np.around(k, 4)}")
    
    # ax.legend(title="Ks (m/hr)", loc='best', framealpha=0.6)#, bbox_to_anchor=(1.01, 0.41))
    ax.legend(
    title="Ks (m/hr)", 
    loc='upper right',      # Places it in the empty bottom-left area
    ncols=3,               # Splits the legend into 3 columns
    fontsize='small',      # Shrinks the text slightly
    labelspacing=0.3,      # Reduces vertical space between rows
    handletextpad=0.5,      # Reduces space between marker and text
    framealpha=0.6         # Makes the legend background slightly transparent
    )
    
    # Perform power-law fit
    x_fit, y_fit, a_fit, b_fit, r2 = fitting_func(x_all, y_all)
    
    plt.plot(x_fit, y_fit, color="k", label="fitted line", linewidth=2)
    
    print(f"Power-law fit: y = {round(a_fit, 2)} * x^{round(b_fit, 2)}")
    print(f"R² = {round(r2, 2)}")
    
    # Format plot
    ax.grid(True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(f"{xlabel}", fontsize=12)
    ax.set_ylabel(f"{ylabel}", fontsize=12)

    # ax.set_title(rf"$y = {round(a_fit, 2)} x^{{{round(b_fit, 2)}}}$  $\mathit{{R}}^2 = {round(r2, 2)}$", fontsize=12)

    plt.tight_layout()
    
    # Save figure
    fig.savefig(os.path.join(output_path, f"inf_{name}.pdf"), dpi=500)
    plt.close()

def exf_solution_plots(solind):
    """
    Create comprehensive plot showing drainage time scaling relationships.
    
    Plots the non-dimensional drainage velocity (lambda / (SST_d * Ks))
    against the non-dimensional depth (alpha * d) for different soil types.
    Includes power-law fit to the data.
    """
    exf_cases_path = os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(exf_cases_path)
    
    # Get unique soil types (by hydraulic conductivity)
    soil_types = sorted(df["k"].unique())
    
    # Create color map and marker styles
    soil_types_colors = cmocean.cm.balance(np.linspace(0, 1, len(soil_types)))
    colors_dic = {soil_types[i]: soil_types_colors[i] for i in range(len(soil_types))}
    
    list_markers = ["o", "^", "s", "P", "*", "X", "d", "p", "2", r"$\clubsuit$", (5,2), "x"]
    markers = {soil_types[i]: list_markers[i] for i in range(len(soil_types))}
    
    print(f"Number of soil types: {len(soil_types)}")
    print(f"Soil types (Ks values): {soil_types}")
    
    # Create figure
    ax = plt
    ax.figure(figsize=(7.09, 3.54))
    
    # Collect all data for fitting
    x_all = []
    y_all = []
    
    # Plot each soil type
    for soil in soil_types:
        df_soil = df[df["k"] == soil]
        
        for i in range(len(df_soil)):
            q = df_soil["q"].iloc[i]
            k = df_soil["k"].iloc[i]
            d = df_soil["d"].iloc[i]
            alfa = df_soil["alfa"].iloc[i]
            exf_t = df_soil["exf_time"].iloc[i]
            toplayer_pressure = df_soil["toplayer_pressure"].iloc[i]
            
            # Filter: only cases where k >= q and exf_t is valid
            if k >= q and exf_t != 0:
                # Compute non-dimensional variables
                x, y, xlabel, ylabel, name = run_exfsolution_case(solind,
                    q, k, d, df_soil["n"].iloc[i], alfa,
                    df_soil["theta_r"].iloc[i], df_soil["theta_s"].iloc[i],
                    exf_t, toplayer_pressure
                )
                
                x_all.append(x)
                y_all.append(y)
                
                ax.scatter(x, y, c=colors_dic[soil], cmap='jet', s=10, marker=markers[k])
    
    # Create legend for soil types
    for k in markers.keys():
        ax.scatter([], [], c=colors_dic[k], s=10, marker=markers[k], 
                  label=f"{np.around(k, 4)}")
    
    # ax.legend(title="Ks (m/hr)", loc='upper right', bbox_to_anchor=(1.01, 0.8))
    ax.legend(
    title="Ks (m/hr)", 
    loc='best',      # Places it in the empty bottom-left area
    ncols=3,               # Splits the legend into 3 columns
    fontsize='small',      # Shrinks the text slightly
    labelspacing=0.3,      # Reduces vertical space between rows
    handletextpad=0.5,      # Reduces space between marker and text
    framealpha=0.6         # Makes the legend background slightly transparent
    )
    # Perform power-law fit
    x_fit, y_fit, a_fit, b_fit, r2 = fitting_func(x_all, y_all)
    
    plt.plot(x_fit, y_fit, color="k", label="fitted line", linewidth=2) if solind!=4 else None
    
    print(f"Power-law fit: y = {round(a_fit, 2)} * x^{round(b_fit, 2)}")
    print(f"R² = {round(r2, 2)}")
    
    # Format plot
    ax.grid(True)
    ax.xscale("log")
    ax.yscale("log")
    ax.xlabel(f"{xlabel}", fontsize=12)
    ax.ylabel(f"{ylabel}", fontsize=12)
    # ax.title(rf"$y = {round(a_fit, 2)} x^{{{round(b_fit, 2)}}}$  $\mathit{{R}}^2 = {round(r2, 2)}$", fontsize=12)
    plt.tight_layout()
    
    # Save figure
    figname = f"exf_{name}.pdf"
    ax.savefig(os.path.join(output_path, figname), dpi=500)
    plt.close()

def _powerlaw_r2_from_curve(x_values, y_values, a_fit, b_fit):
    """
    Evaluate a power-law curve against observed data in log space.

    This mirrors the log-log comparison used by fitting_func so the returned
    score is directly comparable to the original fit.
    """
    x_values = np.asarray(x_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)

    valid_mask = np.isfinite(x_values) & np.isfinite(y_values) & (x_values > 0) & (y_values > 0)
    x_values = x_values[valid_mask]
    y_values = y_values[valid_mask]

    if x_values.size < 2:
        raise ValueError("At least two positive data points are required to compute R².")

    y_pred = powerlaw_func(x_values, a_fit, b_fit)
    corr_matrix = np.corrcoef(np.log(y_values), np.log(y_pred))

    if np.isnan(corr_matrix).any():
        return np.nan

    return corr_matrix[0, 1] ** 2

def plot_solution_transferability_across_datasets(output_path=None):
    """
    Create one vertically stacked figure for each of 4 solutions across 2 CSV inputs.

    For each of inf solution 1, exf solution 1, exf solution 2, and exf
    solution 3, this function creates one figure with two subplots: the top
    subplot uses inf_exf_times_config.csv and the bottom subplot uses
    tolerance_05_inf_exf_times.csv with the same fitted curve from the first
    CSV.
    """
    if output_path is None:
        output_path = os.path.join(DIRPATH, "outputs")

    config_df = pd.read_csv(os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv"))
    tolerance_df = pd.read_csv(os.path.join(DIRPATH, "inputs", "tolerance_05_inf_exf_times.csv"))

    solution_specs = [
        {"display_name": "inf solution 1", "slug": "inf_solution_1", "kind": "inf", "case_id": 1},
        {"display_name": "exf solution 1", "slug": "exf_solution_1", "kind": "exf", "case_id": 1},
        {"display_name": "exf solution 2", "slug": "exf_solution_2", "kind": "exf", "case_id": 2},
        {"display_name": "exf solution 3", "slug": "exf_solution_3", "kind": "exf", "case_id": 3},
    ]

    def collect_rows(df, solution_spec):
        soil_types = sorted(df["k"].unique())
        soil_colors = cmocean.cm.balance(np.linspace(0, 1, len(soil_types)))
        colors_dic = {soil_types[i]: soil_colors[i] for i in range(len(soil_types))}
        list_markers = ["o", "^", "s", "P", "*", "X", "d", "p", "2", r"$\clubsuit$", (5, 2), "x"]
        markers = {soil_types[i]: list_markers[i % len(list_markers)] for i in range(len(soil_types))}

        rows = []
        xlabel = None
        ylabel = None

        for soil in soil_types:
            df_soil = df[df["k"] == soil]

            for i in range(len(df_soil)):
                q = df_soil["q"].iloc[i]
                k = df_soil["k"].iloc[i]

                if k < q:
                    continue

                d = df_soil["d"].iloc[i]
                n = df_soil["n"].iloc[i]
                alfa = df_soil["alfa"].iloc[i]
                theta_r = df_soil["theta_r"].iloc[i]
                theta_s = df_soil["theta_s"].iloc[i]
                top_layer_pressure = df_soil["toplayer_pressure"].iloc[i]

                if solution_spec["kind"] == "inf":
                    x, y, xlabel, ylabel, _ = run_infsolution_case(
                        solution_spec["case_id"],
                        q, k, d, n, alfa, theta_r, theta_s,
                        df_soil["inf_time"].iloc[i],
                        top_layer_pressure,
                    )
                else:
                    x, y, xlabel, ylabel, _ = run_exfsolution_case(
                        solution_spec["case_id"],
                        q, k, d, n, alfa, theta_r, theta_s,
                        df_soil["exf_time"].iloc[i],
                        top_layer_pressure,
                    )

                rows.append((x, y, soil))

        return rows, xlabel, ylabel, colors_dic, markers

    def draw_subplot(ax, rows, colors_dic, markers, xlabel, ylabel, curve_x, curve_y, annotation_text,
                     x_limits=None, y_limits=None, annotation_loc='upper left', panel_label=None):
        for x, y, soil in rows:
            ax.scatter(x, y, c=[colors_dic[soil]], s=10, marker=markers[soil])

        for soil in sorted(colors_dic.keys()):
            ax.scatter([], [], c=[colors_dic[soil]], s=10, marker=markers[soil], label=f"{np.around(soil, 4)}")

        ax.plot(curve_x, curve_y, color="k", linewidth=2)
        ax.legend(
            title="Ks (m/hr)",
            loc="best",
            ncols=3,
            fontsize='small',
            labelspacing=0.3,
            handletextpad=0.5,
            framealpha=0.6,
        )
        ax.grid(True)
        ax.set_xscale("log")
        ax.set_yscale("log")
        if x_limits is not None:
            ax.set_xlim(x_limits)
        if y_limits is not None:
            ax.set_ylim(y_limits)
        ax.set_xlabel(f"{xlabel}", fontsize=12)
        ax.set_ylabel(f"{ylabel}", fontsize=12)
        if panel_label is not None:
            ax.text(
                0.02,
                0.98,
                panel_label,
                transform=ax.transAxes,
                va='top',
                ha='left',
                fontsize=12,
                fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='none', alpha=0.65),
            )
        ax.text(
            0.97 if annotation_loc.endswith('right') else 0.03,
            0.97 if annotation_loc.startswith('upper') else 0.03,
            annotation_text,
            transform=ax.transAxes,
            va='top' if annotation_loc.startswith('upper') else 'bottom',
            ha='right' if annotation_loc.endswith('right') else 'left',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.6),
        )

    results = {}

    for solution_spec in solution_specs:
        print(f"Processing {solution_spec['display_name']}...")
        config_rows, xlabel, ylabel, config_colors, config_markers = collect_rows(config_df, solution_spec)
        tolerance_rows, _, _, _, _ = collect_rows(tolerance_df, solution_spec)

        if len(config_rows) < 2:
            raise ValueError(f"Not enough config data points to fit {solution_spec['display_name']}.")
        if len(tolerance_rows) < 2:
            raise ValueError(f"Not enough tolerance data points to evaluate {solution_spec['display_name']}.")

        x_config = np.asarray([row[0] for row in config_rows], dtype=float)
        y_config = np.asarray([row[1] for row in config_rows], dtype=float)
        x_tolerance = np.asarray([row[0] for row in tolerance_rows], dtype=float)
        y_tolerance = np.asarray([row[1] for row in tolerance_rows], dtype=float)

        x_all = np.concatenate([x_config, x_tolerance])
        y_all = np.concatenate([y_config, y_tolerance])

        x_min = np.min(x_all)
        x_max = np.max(x_all)
        y_min = np.min(y_all)
        y_max = np.max(y_all)

        x_limits = (x_min, x_max)
        y_limits = (y_min, y_max)

        _, _, a_fit, b_fit, original_r2 = fitting_func(x_config, y_config)
        tolerance_r2 = _powerlaw_r2_from_curve(x_tolerance, y_tolerance, a_fit, b_fit)
        print(f"{solution_spec['display_name']} | config R² = {original_r2:.3f}, tolerance R² = {tolerance_r2:.3f}")
        print(f"Fitted curve: y = {a_fit:.3g} * x^{b_fit:.3g}")

        x_config_curve = np.logspace(np.log10(np.min(x_config)), np.log10(np.max(x_config)), 200)
        y_config_curve = powerlaw_func(x_config_curve, a_fit, b_fit)
        x_tolerance_curve = np.logspace(np.log10(np.min(x_tolerance)), np.log10(np.max(x_tolerance)), 200)
        y_tolerance_curve = powerlaw_func(x_tolerance_curve, a_fit, b_fit)

        fig, axes = plt.subplots(2, 1, figsize=(7.09, 7.08), sharex=False)

        draw_subplot(
            axes[0],
            config_rows,
            config_colors,
            config_markers,
            xlabel,
            ylabel,
            x_config_curve,
            y_config_curve,
            rf"$R^2$ (Solver residual tolerance $10^{{-7}}$) = {original_r2:.3f}",
            x_limits=x_limits,
            y_limits=y_limits,
            annotation_loc='lower left' if solution_spec['slug'] == 'exf_solution_1' else 'upper right',
            panel_label='(a)',
        )

        draw_subplot(
            axes[1],
            tolerance_rows,
            config_colors,
            config_markers,
            xlabel,
            ylabel,
            x_tolerance_curve,
            y_tolerance_curve,
            # rf"$R^2$ (Solver residual tolerance $10^{{-7}}$) = {original_r2:.3f}"
            # "\n"
            rf"$R^2$ (Solver residual tolerance $10^{{-5}}$) = {tolerance_r2:.3f}",
            x_limits=x_limits,
            y_limits=y_limits,
            annotation_loc='lower left' if solution_spec['slug'] == 'exf_solution_1' else 'upper right',
            panel_label='(b)',
        )

        # fig.suptitle(solution_spec['display_name'], fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        fig.savefig(os.path.join(output_path, f"{solution_spec['slug']}_transferability.png"), dpi=500)
        plt.close(fig)

        results[solution_spec["slug"]] = {
            "a_fit": a_fit,
            "b_fit": b_fit,
            "config_r2": original_r2,
            "tolerance_r2": tolerance_r2,
        }

    return results

def evaluate_residuals_and_systematic_errors(x_values=None, y_values=None, a_fit=None, b_fit=None, xlabel="x",
                                             output_path=None, figname="residuals_vs_x.png"):
    """
    Evaluate residuals and systematic errors for a fitted power-law relationship.

    Residuals are defined as observed y minus fitted y. The fitted relationship is
    y = a_fit * x ** b_fit.

    Parameters
    ----------
    x_values : array-like
        Observed x data.
    y_values : array-like
        Observed y data.
    a_fit : float or None
        Power-law coefficient. Must be provided.
    b_fit : float or None
        Power-law exponent. Must be provided.
    output_path : str or None
        Directory where the plot will be saved. Defaults to the module output folder.
    figname : str
        Name of the saved residual plot.

    Returns
    -------
    dict
        Dictionary containing x values, observed y, fitted y, and residuals.
    """
    if x_values is None or y_values is None:
        raise ValueError("x_values and y_values must be provided to compute residuals.")
    if a_fit is None or b_fit is None:
        raise ValueError("a_fit and b_fit must be provided to evaluate the fitted power-law.")

    x_values = np.asarray(x_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)

    if x_values.shape != y_values.shape:
        raise ValueError("x_values and y_values must have the same shape.")

    if np.any(x_values <= 0):
        raise ValueError("x_values must be positive to evaluate and plot the power-law on log scales.")

    y_fitted = powerlaw_func(x_values, a_fit, b_fit)
    residuals = np.log(y_values) - np.log(y_fitted)

    if output_path is None:
        output_path = os.path.join(DIRPATH, "outputs")

    # os.makedirs(output_path, exist_ok=True)

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.scatter(y_fitted, residuals, s=70, alpha=0.75)
    # ax.scatter(y_fitted, y_values, s=70, alpha=0.75)
    ax.axhline(0.0, color="k", linestyle="--", linewidth=1.5)
    ax.set_xscale("log")
    # ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Residuals (observed - fitted)")
    ax.set_title("Residuals vs x")
    ax.grid(True, which="both", alpha=0.3)
    plt.tight_layout()

    fig.savefig(os.path.join(output_path, figname), dpi=300)
    plt.close(fig)

    return {
        "x_values": x_values,
        "y_observed": y_values,
        "y_fitted": y_fitted,
        "residuals": residuals,
    }

def run_residual_analysis_inf():
    """
    Run residual analysis for the fitted power-law relationship.
    """
    cases_path = os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(cases_path)
    df["phi"] = df["theta_s"] - df["theta_r"]

    # inf solution
    df["inf_xd"] = df[["k", "alfa", "d", "n", "phi"]].apply(
        lambda row: row["k"] * row["alfa"] * row["d"] * row["n"] * row["phi"], axis=1
    )
    df["inf_y"] = df.apply(lambda row: row["k"] * row["inf_time"] / row["d"], axis=1)
    df["inf_x"] = df["q"] / df["inf_xd"]
    x_values = df["inf_x"].values
    y_values = df["inf_y"].values
    print(df)

    x_fit, y_fit, a_fit, b_fit, r2 = fitting_func(x_values, y_values)
    print(a_fit, b_fit, r2)
    # a_fit = 0.07
    # b_fit = -0.65
    y_fitted = powerlaw_func(x_values, a_fit, b_fit)

    results = evaluate_residuals_and_systematic_errors(
        x_values=x_values,
        y_values=y_values,
        a_fit=a_fit,
        b_fit=b_fit,
        xlabel=r"$K_s\cdot t'$",
        output_path=output_path,
        figname="residuals_infiltration.png"
    )    
    return results

def run_residual_analysis_dr_Bpi():
    """
    Run residual analysis for the fitted power-law relationship.
    """
    cases_path = os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(cases_path)
    df["phi"] = df["theta_s"] - df["theta_r"]

    # exf solution
    df["exf_xd"] = df[["q", "alfa", "d", "n", "phi"]].apply(
        lambda row: row["q"] * row["alfa"] * row["d"] * row["n"] * row["phi"], axis=1
    )
    df["exf_y"] = df.apply(lambda row: row["q"] * row["exf_time"] / row["d"], axis=1)
    df["exf_x"] = df["k"] / df["exf_xd"]
    x_values = df["exf_x"].values
    y_values = df["exf_y"].values
    print(df)

    x_fit, y_fit, a_fit, b_fit, r2 = fitting_func(x_values, y_values)
    print(a_fit, b_fit, r2)
    # a_fit = 0.07
    # b_fit = -0.65
    y_fitted = powerlaw_func(x_values, a_fit, b_fit)

    results = evaluate_residuals_and_systematic_errors(
        x_values=x_values,
        y_values=y_values,
        a_fit=a_fit,
        b_fit=b_fit,
        xlabel=r"$q\cdot t'$",
        output_path=output_path,
        figname="residuals_drainage_Bpi.png"
    )    
    return results

def run_residual_analysis_dr_v():
    """
    Run residual analysis for the fitted power-law relationship.
    """
    cases_path = os.path.join(DIRPATH, "inputs", "inf_exf_times_config.csv")
    output_path = os.path.join(DIRPATH, "outputs")
    df = pd.read_csv(cases_path)
    df["phi"] = df["theta_s"] - df["theta_r"]

    # exf solution
    # v = (d-abs(toplayer_pressure)) / exf_time
    df["v"] = (df["d"] - df["toplayer_pressure"].abs()) / df["exf_time"]
    df["exf_x"] = df["alfa"] * df["d"]
    df["exf_y"] = df["v"] / df["k"]
    x_values = df["exf_x"].values
    y_values = df["exf_y"].values
    print(df)

    x_fit, y_fit, a_fit, b_fit, r2 = fitting_func(x_values, y_values)
    print(a_fit, b_fit, r2)
    # a_fit = 0.07
    # b_fit = -0.65
    y_fitted = powerlaw_func(x_values, a_fit, b_fit)

    results = evaluate_residuals_and_systematic_errors(
        x_values=x_values,
        y_values=y_values,
        a_fit=a_fit,
        b_fit=b_fit,
        xlabel=r"$\frac{v_d}{k_s}$",
        output_path=output_path,
        figname="residuals_drainage_v.png"
    )    
    return results

if __name__ == "__main__":
    # inf_solution_plots(1)

    # exf_solution_plots(1)
    # exf_solution_plots(2)
    # exf_solution_plots(3)

    plot_solution_transferability_across_datasets()
    # run_residual_analysis_inf()
    # run_residual_analysis_dr_Bpi()
    # run_residual_analysis_dr_v()