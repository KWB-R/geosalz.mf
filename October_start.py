# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 16:47:31 2022

@author: dbaldw
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import flopy
import flopy.utils.binaryfile as bf
import pandas as pd

sys.path.append('modules')

sys.path.append(os.path.join("..", "common"))

plt.rcParams.update({'font.size': 3})

figure_size = (12, 12)

parameter_units = {"recharge": "$m^{3} days^{-1} m^{-2}$", "slt_cnc": "$kg m^{-3}$", "slt_strt": "$kg m^{-3}$",
                   "perlen": "$days$",
                   "perlent": "$days$"}


lay_disc = 1  # layer discretization option 
           

perlen = 100*365  # Zeit für Flow-Modell
perlent = 5000*365  # Zeit für Transportmodell
steady = [False]  # instationär // transient

nstpf = 1
nstpt = 1

nper = 1  # Number of periods

riv_bed = -3.0  # Flusstiefe
flt1_top = -10.0  # Tiefe Oberkante 1. Filter
flt2_top = -21.0  # Tiefe Oberkante 2. Filter
top_hol = -35.0  # Entspricht Top Holstein

ncol = 50  # Number of columns
nrow = 50  # Number of rows
delr = 10.0  # Column width ($m$)
delc = 10.0  # Row width ($m$)
top = 0.0  # Top of the model ($m$)
slt_cnc = 4.0  # Konzentration an Holstein-Fehlstelle
porosity = 0.30
al = 20.0  # longitudinale Dispersivität
rch = 0.1/365  # GW Neubildungsrate
strt_head_u = 0.065
rch_slt = 0.025
crch=rch_slt

# Model layers

if lay_disc == 2:
    nlay = 24
    lay1 = -3
    lay2 = -5
    lay3 = -7
    lay4 = -9
    lay5 = -11
    lay6 = -13
    lay7 = -15
    lay8 = -17
    lay9 = -19
    lay10 = -21
    lay11 = -23
    lay12 = -25
    lay13 = -27
    lay14 = -29
    lay15 = -31
    lay16 = -33
    lay17 = -35
    lay18 = -37
    lay19 = -39
    lay20 = -41
    lay21 = -50
    lay22 = -70
    lay23 = -90
    lay24 = -100
    botm = [lay1, lay2, lay3, lay4, lay5, lay6, lay7, lay8, lay9, lay10,
            lay11, lay12, lay13, lay14, lay15, lay16, lay17, lay18, lay19,
            lay20, lay21, lay22, lay23, lay24]  # Layer bottom elevations ($m$)
    topm = [top, lay1, lay2, lay3, lay4, lay5, lay6, lay7, lay8, lay9, lay10,
            lay11, lay12, lay13, lay14, lay15, lay16, lay17, lay18, lay19,
            lay20, lay21, lay22, lay23]  # Layer top elevations
    botm_calc = botm
    topm_calc = topm
    lay_thick = list()
    for topm_calc, botm_calc in zip(topm_calc, botm_calc):
        lay_thick.append(topm_calc-botm_calc)

if lay_disc == 1:
    nlay = 30
    lay1 = -3
    lay2 = -5
    lay3 = -7
    lay4 = -9
    lay5 = -11
    lay6 = -13
    lay7 = -15
    lay8 = -17
    lay9 = -19
    lay10 = -21
    lay11 = -23
    lay12 = -25
    lay13 = -27
    lay14 = -29
    lay15 = -31
    lay16 = -33
    lay17 = -35
    lay18 = -37
    lay19 = -39
    lay20 = -41
    lay21 = -50
    lay22 = -55
    lay23 = -60
    lay24 = -65
    lay25 = -70
    lay26 = -75
    lay27 = -80
    lay28 = -85
    lay29 = -90
    lay30 = -100
    botm = [lay1, lay2, lay3, lay4, lay5, lay6, lay7, lay8, lay9, lay10,
            lay11, lay12, lay13, lay14, lay15, lay16, lay17, lay18, lay19,
            lay20, lay21, lay22, lay23, lay24, lay25, lay26, lay27, lay28, lay29, lay30]  # Layer bottom elevations ($m$)
    topm = [top, lay1, lay2, lay3, lay4, lay5, lay6, lay7, lay8, lay9, lay10,
            lay11, lay12, lay13, lay14, lay15, lay16, lay17, lay18, lay19,
            lay20, lay21, lay22, lay23, lay24, lay25, lay26, lay27, lay28, lay29]  # Layer top elevations
    botm_calc = botm
    topm_calc = topm
    lay_thick = list()
    for topm_calc, botm_calc in zip(topm_calc, botm_calc):
        lay_thick.append(topm_calc-botm_calc)





else:
    raise Exception('layer discretization incorrectly defined')
    sys.exit()


# Parameters defined after grid definition
k1 = 25.92  # Horizontale hydraulische Leitfähigkeit ($m/day$)
k3 = 2.592  # Vertikale hydraulische Leitfähigkeit ($m/day$)
k11 = k1 * np.ones((nlay, nrow, ncol), dtype=float)
# Horizontale hydraulische Leitfähigkeit für qhol in unterster Layer ($m/day$)
k2 = 0.002592
# Vertikale hydraulische Leitfähigkeit für qhol in unterster Layer ($m/day$)
k4 = 0.0002592
riv_bank_k = 86400 * 5 * 10** -5
riv_bed_k = 86400 * 1 * 10 ** -9

k33 = k3 * np.ones((nlay, nrow, ncol), dtype=float)

if lay_disc == 1:
    k11[17:20, 0:20, 0:50] = k2
    k11[17:20, 30:50, 0:50] = k2
    k33[17:20, 0:20, 0:50] = k4
    k33[17:20, 30:50, 0:50] = k4
    ghb_spd = []
    for k in range(49, 50):
        for l in range(0, 12):
            ghb_spd += [[l, i, k, 0, k1*lay_thick[l]
                         * (delr/500)] for i in range(nrow)]

elif lay_disc == 2:
    k11[17:20, 0:20, 0:50] = k2
    k11[17:20, 30:50, 0:50] = k2
    k33[17:20, 0:20, 0:50] = k4
    k33[17:20, 30:50, 0:50] = k4
    ghb_spd = []
    for k in range(49, 50):
        for l in range(0, 12):
            ghb_spd += [[l, i, k, 0, k1*lay_thick[l]
                         * (delr/500)] for i in range(nrow)]


# Well network


flt1_len = 5.0  # Länge Filter flache Brunnen
flt2_len = 6.0  # Länge Filter tiefer Brunnen

frate_fl = 0##-720 # Förderrate flache Brunnen
frate_tf = 0##-600 # Förderrate tiefer Brunnen

# Well arrays

if lay_disc == 1:
 wel_spd1 = [
     (5, 0, 24, frate_fl / 6),
     (6, 0, 24, frate_fl / 6),
     (7, 0, 24, frate_fl / 6),
     (5, 12, 24, frate_fl / 3),
     (6, 12, 24, frate_fl / 3),
     (7, 12, 24, frate_fl / 3),
     (12, 24, 24, frate_tf / 3),
     (13, 24, 24, frate_tf / 3),
     (14, 24, 24, frate_tf / 3),
     (5, 36, 24, frate_fl / 3),
     (6, 36, 24, frate_fl / 3),
     (7, 36, 24, frate_fl / 3),
     (5, 49, 24, frate_fl / 6),
     (6, 49, 24, frate_fl / 6),
     (7, 49, 24, frate_fl / 6),]
elif lay_disc == 2:
    wel_spd1 = [
        (5, 0, 24, frate_fl / 6),
        (6, 0, 24, frate_fl / 6),
        (7, 0, 24, frate_fl / 6),
        (5, 12, 24, frate_fl / 3),
        (6, 12, 24, frate_fl / 3),
        (7, 12, 24, frate_fl / 3),
        (12, 24, 24, frate_tf / 3),
        (13, 24, 24, frate_tf / 3),
        (14, 24, 24, frate_tf / 3),
        (5, 36, 24, frate_fl / 3),
        (6, 36, 24, frate_fl / 3),
        (7, 36, 24, frate_fl / 3),
        (5, 49, 24, frate_fl / 6),
        (6, 49, 24, frate_fl / 6),
        (7, 49, 24, frate_fl / 6),]
    

wells = wel_spd1  # + wel_spd2
wel_spd12 = {0: wells}

# Flow model section
modelname_mf = "flow"

# Generates the model folder containing all model files
ws = './October_start'+str(nlay)+'test_3'
orig = './October_start'+str(nlay)+'test_3'

mf = flopy.modflow.Modflow(
    modelname=modelname_mf, model_ws=ws, exe_name="C:/WRDD/MF2005.1_12/bin/mf2005.exe")

flopy.modflow.ModflowDis(
    mf,
    nlay=nlay,
    nrow=nrow,
    ncol=ncol,
    delr=delr,
    delc=delc,
    top=top,
    botm=botm,
    perlen=perlen,
    steady=steady,
    nstp=nstpf
)

# Flow boundary conditions

if lay_disc == 2:
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    ibound[23:24, 0:49, 0] = -1
    ibound[23:24, 0:49, 49] = -1
    strt = np.ones((nlay, nrow, ncol), dtype=np.float32)
    strt[23:24, 0:50, 0] = strt_head_u
    strt[23:24, 0:50, 49] = strt_head_u  # Starting head ($m$)
    
elif lay_disc == 1:
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
  #  ibound[29:30, 0:50, 0] = -1
    ibound[29:30, 49, 0:50] = -1
    strt = np.ones((nlay, nrow, ncol), dtype=np.float32)
  #  strt[29:30, 0:50, 0] = strt_head_u
    strt[29:30, 49, 0:50] = strt_head_u  # Starting head ($m$)

# River boundary conditions array
riv_spd1 = []
for k in range(0, 6):
    riv_spd1 += [[0, i, k, 0, riv_bed_k, riv_bed] for i in range(nrow)]

riv_spd2 = []
for k in range(6, 10):
    riv_spd2 += [[0, i, k, 0, riv_bank_k, riv_bed] for i in range(nrow)]

riv_spd = riv_spd1 + riv_spd2

riv_spd = {0: riv_spd}


# Instantiate layer property bas package
flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)

# Instantiate layer property flow package
flopy.modflow.ModflowLpf(mf, hk=k11, vka=k33, ipakcb=53, ss=5 * 10e-4)

# Instantiate well package
flopy.modflow.ModflowWel(mf, stress_period_data=wel_spd12)

# Instantiate recharge package
flopy.modflow.ModflowRch(mf, ipakcb=1, rech=rch)

# Instantiate river package
flopy.modflow.ModflowRiv(mf, stress_period_data=riv_spd)

# Instantiate solver package
flopy.modflow.ModflowSip(mf, mxiter=90000000)

# Instantiate link mass transport package (for writing linker file)
flopy.modflow.ModflowLmt(mf)

# Instatiate GHB package
flopy.modflow.ModflowGhb(mf, stress_period_data=ghb_spd)

spd = {(0, 0): ["print head", "print budget", "save head", "save budget"]}
oc = flopy.modflow.ModflowOc(
    mf, stress_period_data=spd, compact=True, save_specific_discharge=True)


# Write model inputs
mf.write_input()

# Run flow model
mf.run_model()


# Extract the heads for hydraulic plots
hds = bf.HeadFile(os.path.join(ws, modelname_mf + '.hds'))
times = hds.get_times()
head = hds.get_data(totim=times[-1])

cbb = bf.CellBudgetFile(os.path.join(ws, modelname_mf + '.cbc'))
cbb.get_data(kstpkper=(0, 0))
kstpkper_list = cbb.get_kstpkper()
frf = cbb.get_data(text="FLOW RIGHT FACE", totim=times[-1])[0]
fff = cbb.get_data(text="FLOW FRONT FACE", totim=times[-1])[0]
flf = cbb.get_data(text="FLOW LOWER FACE", totim=times[-1])[0]

# Transport model section
if lay_disc == 2:
    obs_spd = [
        (7, 0, 24),
        (7, 12, 24),
        (14, 24, 24),
        (7, 36, 24),
        (7, 49, 24),]

if lay_disc == 2:
        slt_strt = []
        for l in range (0,19): 
            slt_strt += [(0.05)]
        for l in range(19,23):
            slt_strt += [(slt_cnc/2)]
        for l in range(23,24):
            slt_strt += [(slt_cnc)]

        slt_strt = slt_strt  # Hintergrundkonzentration im Modell

if lay_disc == 1:
    obs_spd = [
        (7, 0, 24),
        (7, 12, 24),
        (14, 24, 24),
        (7, 36, 24),
        (7, 49, 24),]

if lay_disc == 1:
        slt_strt = []
        for l in range (0,19): 
            slt_strt += [(0.05)]
        for l in range(19,29):
            slt_strt += [(slt_cnc/2.1)]
        for l in range(29,30):
            slt_strt += [(slt_cnc-1)]

        slt_strt = slt_strt  # Hintergrundkonzentration im Modell


# Extract concentrations of all wells from .UCN files and create model cross section with concentration overlay
# fname = os.path.join(ws, "MT3D001.UCN")
# ucnobj = flopy.utils.UcnFile(fname)
# times = ucnobj.get_times()
# concentration2 = ucnobj.get_data(totim=times[-1])
# slt_strt = concentration2

modelname_mt = "transport"
model_ws = os.path.join(ws, "mfgwt")
mt = flopy.mt3d.Mt3dms(
    modelname=modelname_mt,
    model_ws=ws,
    exe_name="C:/WRDD/mt3dusgs1.1.0/bin/mt3d-usgs_1.1.0_64.exe",
    modflowmodel=mf,
)

btn = flopy.mt3d.Mt3dBtn(
    mt,
    icbund=1,
    prsity=porosity,
    sconc=slt_strt,
    perlen=perlent,
    nper=1,
    nstp=nstpt,
    obs=obs_spd
)

dceps = 1.0e-5
nplane = 1
npl = 0
nph = 4
npmin = 0
npmax = 8
nlsink = nplane
npsink = nph

mixelm = -1  # -1 bedeutet TVD Solver
adv = flopy.mt3d.Mt3dAdv(mt, mixelm=mixelm)
# dceps=dceps,
# nplane=nplane,
# npl=npl,
# nph=nph,
# npmin=npmin,
# npmax=npmax,
# nlsink=nlsink,
# npsink=npsink,
# percel=0.5,
# )

dsp = flopy.mt3d.Mt3dDsp(mt, al=al)




if lay_disc == 2:
    slt_spd = []
    for k in range(ncol):
       for l in range(23,24):
           slt_spd += [(l, i, 0, slt_cnc, -1) for i in range(nrow)]
           slt_spd += [(l, i, 50, slt_cnc, -1) for i in range(nrow)]
                
    criv_spd = []
    for k in range(10):
        criv_spd += [[0, i, k, 0.05, -1] for i in range(nrow)]
elif lay_disc ==1:         
    slt_spd = []
    for k in range(ncol):
       for l in range(29,30):
        #   slt_spd += [(l, i, 0, slt_cnc, -1) for i in range(nrow)]
           slt_spd += [(l, 49, k, slt_cnc, -1) for i in range(nrow)]
                
    criv_spd = []
    for k in range(10):
        criv_spd += [[0, i, k, 0.05, -1] for i in range(nrow)]



all_spd = slt_spd + criv_spd 
all_spd = {0: all_spd}

ssm = flopy.mt3d.Mt3dSsm(mt, crch=crch, stress_period_data=all_spd)
gcg = flopy.mt3d.Mt3dGcg(mt)


# Section for seawat module, currently not used in results
swt = flopy.seawat.Seawat(
    modflowmodel=mf,
    mt3dmodel=mt,
    modelname=modelname_mt,
    namefile_ext="nam_swt",
    model_ws=ws,
    exe_name="C:WRDD/swt_v4_00_05/exe/swt_v4.exe"
)
vdf = flopy.seawat.SeawatVdf(
    swt,
    mtdnconc=0,
    iwtable=0,
    indense=-1,
    densemin=0,
    densemax=0,
    denseref=1000.0,
    denseslp=0.7143,
    firstdt=1e-3
)

mt.write_input()
# Run transport model
mt.run_model()


  ##Seawat specific code section for seawat files
# fname = modelname_mt + ".vdf"
# f = open(os.path.join(ws, fname), "r")
# lines = f.readlines()
# f.close()
# f = open(os.path.join(ws, fname), "w")
# for line in lines:
#     f.write(line)
# for kper in range(nper):
#     f.write("-1\n")
#     f.close()

def plot_cross():
    fig = plt.figure(figsize=(15, 5))
    ax = fig.add_subplot(1, 1, 1)
    # Next we create an instance of the PlotCrossSection class
    xsect = flopy.plot.PlotCrossSection(model=mf, line={"ROW": 24})
    # Then we can use the plot_grid() method to draw the grid
    # The return value for this function is a matplotlib LineCollection object,
    # which could be manipulated (or used) later if necessary.
    xsect.plot_grid(lw=0.5)
    hl = xsect.plot_array(head, head=head, alpha=0.4)
    # xsect.plot_array(hk, masked_values=[hk[0, 0, 0]], alpha=0.2)
    # xsect.plot_array(vka, masked_values=[hk[0, 0, 0]], alpha=0.2)
    #xsect.plot_array(ibound, masked_values=[ibound[0, 0, 0]], alpha=0.2)
    #xsect.plot_array(strt, masked_values=[strt[0, 0, 0]], alpha=0.2)
    #xsect.plot_array(head, head=head)
    hl
    ha = xsect.contour_array(head, head=head, linewidths=4,
                             cmap="viridis", ) ##levels=[ -2.5, -2, -1.5, -1., -0.5, 0, 0.5])
    ha
    xsect.plot_bc("WEL", color="red")
    xsect.plot_bc("GHB", color="orange")
    xsect.plot_bc("RIV", color="blue")
    # quiver = xsect.plot_discharge(frf, fff, flf, head=head, normalize=True, color="0.75", scale=50, headwidth=3,
    # linewidths=0.1)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    ax.set_title("West-Ost Profilschnitt \n Qmodel = {} [m\u00b3d-\u00b9]".format(-1*(frate_tf+frate_fl*3))+ "  Qtot = {} [m\u00b3d-\u00b9]".format(-1*(frate_tf+frate_fl*4))+"  CHD = {} [m]".format(strt_head_u), fontsize=15)
    plt.xlabel("Länge des Modells [m]", fontsize=15)
    plt.ylabel("Tiefe unter Geländeoberkante [m]", fontsize=15)
    plt.rcParams.update({'font.size': 15})
    cbar = plt.colorbar(hl, orientation='vertical', cmap='viridis')
    cbar.add_lines(ha)
    #cbar.ax.tick_params(size=0)
    cbar.ax.set_ylabel("Head [m]")
    plt.show()


plot_cross()

def plot_cross_NS():
    fig = plt.figure(figsize=(15, 5))
    ax = fig.add_subplot(1, 1, 1)
    # Next we create an instance of the PlotCrossSection class
    xsect = flopy.plot.PlotCrossSection(model=mf, line={"COLUMN": 24})
    # Then we can use the plot_grid() method to draw the grid
    # The return value for this function is a matplotlib LineCollection object,
    # which could be manipulated (or used) later if necessary.
    xsect.plot_grid(lw=0.5)
    hl = xsect.plot_array(head, head=head, alpha=0.4)
    # xsect.plot_array(hk, masked_values=[hk[0, 0, 0]], alpha=0.2)
    # xsect.plot_array(vka, masked_values=[hk[0, 0, 0]], alpha=0.2)
    #xsect.plot_array(ibound, masked_values=[ibound[0, 0, 0]], alpha=0.2)
    #xsect.plot_array(strt, masked_values=[strt[0, 0, 0]], alpha=0.2)
    #xsect.plot_array(head, head=head)
    hl
    ha = xsect.contour_array(head, head=head, linewidths=4,
                             cmap="viridis", ) ##levels=[ -2.5, -2, -1.5, -1., -0.5, 0, 0.5])
    ha
    xsect.plot_bc("WEL", color="red")
    xsect.plot_bc("GHB", color="orange")
    xsect.plot_bc("RIV", color="blue")
    # quiver = xsect.plot_discharge(frf, fff, flf, head=head, normalize=True, color="0.75", scale=50, headwidth=3,
    # linewidths=0.1)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    ax.set_title("Nord-Sud Profilschnitt \n Qmodel = {} [m\u00b3d-\u00b9]".format(-1*(frate_tf+frate_fl*3))+ "  Qtot = {} [m\u00b3d-\u00b9]".format(-1*(frate_tf+frate_fl*4))+"  CHD = {} [m]".format(strt_head_u), fontsize=15)
    plt.xlabel("Länge des Modells [m]", fontsize=15)
    plt.ylabel("Tiefe unter Geländeoberkante [m]", fontsize=15)
    plt.rcParams.update({'font.size': 15})
    cbar = plt.colorbar(hl, orientation='vertical', cmap='viridis')
    cbar.add_lines(ha)
    #cbar.ax.tick_params(size=0)
    cbar.ax.set_ylabel("Head [m]")
    plt.show()

plot_cross_NS()

def plot_cross_boundaries():
    fig = plt.figure(figsize=(15, 5))
    ax = fig.add_subplot(1, 1, 1)
    xsect = flopy.plot.PlotCrossSection(model=mf, line={"ROW": 24})
    xsect.plot_grid(lw=0.5)
    xsect.plot_array(ibound, masked_values=[ibound[0, 0, 0]], alpha=0.8)
    xsect.plot_bc("WEL", color="red")
    xsect.plot_bc("GHB", color="orange")
    xsect.plot_bc("RIV", color="blue")
    ##xsect.plot_array(k11, alpha=1)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    ax.set_title("Boundary Conditions N-S")
    plt.xlabel("Länge des Modells [m]", fontsize=15)
    plt.ylabel("Tiefe unter Geländeoberkante [m]", fontsize=15)
    plt.rcParams.update({'font.size': 15})
    plt.show()

plot_cross_boundaries()


# Extract concentrations of all wells from .UCN files and create concentration plot
fname = os.path.join(ws, "MT3D001.UCN")
ucnobj = flopy.utils.UcnFile(fname)
times = ucnobj.get_times()
conc = ucnobj.get_alldata()

fname = os.path.join(ws, "MT3D001.OBS")
if os.path.isfile(fname):
    cvt = mt.load_obs(fname)
else:
    cvt = None

fname = os.path.join(ws, "MT3D001.MAS")
mvt = mt.load_mas(fname)

plt.figure(figsize=(8,6))

# x = cvt["time"]
# y = cvt["(8, 1, 25)"] * (-1*(frate_fl) / 2) \
#     + cvt["(8, 13, 25)"] * (-1*(frate_fl)) \
#         + cvt["(8, 37, 25)"] * (-1*(frate_fl)) \
#             + cvt["(8, 50, 25)"] * (-1*(frate_fl) / 2) \
#                 +cvt["(15, 25, 25)"] *  (-1* frate_tf)
# plt.plot(x, y, label="Total Salt", color = 'black', linestyle = ':', linewidth = 4)

# x = cvt["time"]
# y = cvt["(8, 1, 25)"] * (-1*(frate_fl) / 2) \
#     + cvt["(8, 13, 25)"] *(-1*(frate_fl)) \
#         + cvt["(8, 37, 25)"] *(-1*(frate_fl)) \
#             + cvt["(8, 50, 25)"] * (-1*(frate_fl) / 2)
# plt.plot(x, y, label="Shallow well Salt", color = 'blue', linestyle = '-.')

# x = cvt["time"]
# y =  cvt["(15, 25, 25)"] * (-1* frate_tf) 
# plt.plot(x, y, label="Deep well salt", color = 'red', linestyle = '-.')

x = cvt["time"]
y = cvt["(8, 1, 25)"] * 1000 
plt.plot(x, y, label="Brunnen 1", color = "orange", linestyle = '--')

x = cvt["time"]
y = cvt["(8, 13, 25)"] * 1000
plt.plot(x, y, label="Brunnen 2", color = "violet", linestyle = '--')

x = cvt["time"]
y = cvt["(15, 25, 25)"] * 1000
plt.plot(x, y, label="Brunnen 3", color = 'blue')

x = cvt["time"]
y = cvt["(8, 37, 25)"] * 1000 
plt.plot(x, y, label="Brunnen 4", color = 'purple', linestyle = '-.')

x = cvt["time"]
y = cvt["(8, 50, 25)"] * 1000  
plt.plot(x, y, label="Brunnen 5", color = 'brown', linestyle = ":")

plt.xticks(fontsize=15)
plt.yticks(fontsize=12)
# plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.])
# plt.yticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400,
#             1500, 1600, 1700, 1800, 1900, 2000])
plt.xlim(0, perlent)
#plt.ylim(0, 24000)
plt.xlabel("Zeit (day)", fontsize=15)
plt.ylabel(("Chlorid (mg/L)"), fontsize=15)
plt.rcParams.update({'font.size': 15})
plt.grid(visible=True)
# title1 = "Chlorid-Konzentrationen\n in der GHBvariante des KHM Friedrichshagen"
# title2 = "Filterlänge: Brunnen A, B, D, E = {} m, Brunnen C = {} m und " \
#         "Tiefe Filteroberkante: Brunnen A, B, D, E = {} m, Brunnen C = {} m".format(flt1_len, flt2_len, flt1_top,
#                                                                                      flt2_top)
#plt.title(title1, pad=10, fontsize=16)
#plt.title(title2, pad=6, fontsize=12)
plt.legend(fontsize=15)

plt.show()


plt.figure(figsize=(8,6))

# x = cvt["time"]
# y = cvt["(8, 1, 25)"] * (-1*(frate_fl) / 2) \
#     + cvt["(8, 13, 25)"] * (-1*(frate_fl)) \
#         + cvt["(8, 37, 25)"] * (-1*(frate_fl)) \
#             + cvt["(8, 50, 25)"] * (-1*(frate_fl) / 2) \
#                 +cvt["(15, 25, 25)"] *  (-1* frate_tf)
# plt.plot(x, y, label="Total Salt", color = 'black', linestyle = ':', linewidth = 4)

# x = cvt["time"]
# y = cvt["(8, 1, 25)"] * (-1*(frate_fl) / 2) \
#     + cvt["(8, 13, 25)"] *(-1*(frate_fl)) \
#         + cvt["(8, 37, 25)"] *(-1*(frate_fl)) \
#             + cvt["(8, 50, 25)"] * (-1*(frate_fl) / 2)
# plt.plot(x, y, label="Shallow well Salt", color = 'blue', linestyle = '-.')

# x = cvt["time"]
# y =  cvt["(15, 25, 25)"] * (-1* frate_tf) 
# plt.plot(x, y, label="Deep well salt", color = 'red', linestyle = '-.')

x = cvt["time"] / 365
y = cvt["(8, 1, 25)"] * 1000 
plt.plot(x, y, label="Brunnen 1", color = "orange", linestyle = '--')

x = cvt["time"] / 365
y = cvt["(8, 13, 25)"] * 1000
plt.plot(x, y, label="Brunnen 2", color = "violet", linestyle = '--')

x = cvt["time"] / 365
y = cvt["(15, 25, 25)"] * 1000
plt.plot(x, y, label="Brunnen 3", color = 'blue')

x = cvt["time"] / 365
y = cvt["(8, 37, 25)"] * 1000 
plt.plot(x, y, label="Brunnen 4", color = 'purple', linestyle = '-.')

x = cvt["time"] / 365
y = cvt["(8, 50, 25)"] * 1000  
plt.plot(x, y, label="Brunnen 5", color = 'brown', linestyle = ":")

plt.xticks(fontsize=15)
plt.yticks(fontsize=12)
# plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.])
# plt.yticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400,
#             1500, 1600, 1700, 1800, 1900, 2000])
plt.xlim(0, perlent / 365)
#plt.ylim(0, 24000)
plt.xlabel("Zeit (Jahr)", fontsize=15)
plt.ylabel(("Chlorid (mg/l)"), fontsize=15)
plt.rcParams.update({'font.size': 15})
plt.grid(visible=True)
# title1 = "Chlorid-Konzentrationen\n in der GHBvariante des KHM Friedrichshagen"
# title2 = "Filterlänge: Brunnen A, B, D, E = {} m, Brunnen C = {} m und " \
#         "Tiefe Filteroberkante: Brunnen A, B, D, E = {} m, Brunnen C = {} m".format(flt1_len, flt2_len, flt1_top,
#                                                                                      flt2_top)
#plt.title(title1, pad=10, fontsize=16)
#plt.title(title2, pad=6, fontsize=12)
plt.legend(fontsize=15)

plt.show()


# Extract concentrations of all wells from .UCN files and create model cross section with concentration overlay
fname = os.path.join(ws, "MT3D001.UCN")
ucnobj = flopy.utils.UcnFile(fname)
times = ucnobj.get_times()
# conc = ucnobj.get_alldata()
concentration = ucnobj.get_data(totim=times[-1])

fig = plt.figure(figsize=(15, 5))
ax = fig.add_subplot(1, 1, 1)
xsect = flopy.plot.PlotCrossSection(model=mf, ax=ax, line={"COLUMN": 24})
arr = xsect.plot_array(concentration*1000)
xsect.plot_grid(lw=0.2)
xsect.plot_bc("WEL", color="red")
#xsect.plot_array(hk, masked_values=[hk[0, 0, 0]], alpha=0.2)
cv = xsect.contour_array(concentration*1000, levels=[50, 250, 2000, 3000, 4000], linewidths=2,
                         colors="white", )
def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s}" if plt.rcParams["text.usetex"] else f"{s}"

ax.clabel(cv, cv.levels, inline=True, fontsize =12, fmt = fmt)
ax.set_title("Simulierte Konzentration T = {}".format(perlent/365)+" a \n Qmodel = {} [m\u00b3d-\u00b9]".format(-1*(frate_tf+frate_fl*3))+ "  Qtot = {} [m\u00b3d-\u00b9]".format(-1*(frate_tf+frate_fl*4))+"  CHD = {} [m]".format(strt_head_u), pad=15)
#quiver = xsect.plot_discharge(frf, fff, flf, head=head, normalize=True, color="0.75", scale=50, headwidth=3,
 #                             linewidths=0.1)
plt.rcParams.update({'font.size': 15})
cbar = plt.colorbar(arr, orientation='vertical', cmap='viridis')
cbar.add_lines(cv)
cbar.ax.set_ylabel(r"Konzentration in mg/L")
plt.show()

