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

nlay = 9  # Number of layers
ncol = 50  # Number of columns
nrow = 50  # Number of rows
delr = 10.0  # Column width ($m$)
delc = 10.0  # Row width ($m$)
top = 0.0  # Top of the model ($m$)
slt_cnc = 2.0  # Konzentration an Holstein-Fehlstelle
slt_strt = 0.00  # Hintergrundkonzentration im Modell
porosity = 0.30
al = 20.0  # longitudinale Dispersivität
rch = 2.7e-4  # GW Neubildungsrate
k1 = 25.92  # Horizontale hydraulische Leitfähigkeit ($m/day$)
k3 = 2.592  # Vertikale hydraulische Leitfähigkeit ($m/day$)
k11 = k1 * np.ones((nlay, nrow, ncol), dtype=float)
k2 = 0.002592  # Horizontale hydraulische Leitfähigkeit für qhol in unterster Layer ($m/day$)
k4 = 0.0002592  # Vertikale hydraulische Leitfähigkeit für qhol in unterster Layer ($m/day$)
#k5 = 50 # Horizontale hydraulische LF in Fehlstelle
#k6 = 50 # Horizontale hydraulische LF in unterster Layer
#k7 = 5 # Vertikale hydraulische LF in Fehlstelle
#k8 = 5 # Vertikale hydraulische LF in unterster Layer


k33 = k3 * np.ones((nlay, nrow, ncol), dtype=float)
k11[7, 0:20, 0:50] = k2
k11[7, 30:50, 0:50] = k2
k33[7, 0:20, 0:50] = k4
k33[7, 30:50, 0:50] = k4
#k11[7, 20:30, 0:50] = k5
#k11[8, 0:50, 0:50] = k6
#k33[7, 20:30, 0:50] = k7
#k33[8, 0:50, 0:50] = k8

perlen = 365  # Zeit für Flow-Modell
perlent = 3500  # Zeit für Transportmodell
steady = [False]  # instationär

nstpf = 1
nstpt = 1

nper = 1  # Number of periods

riv_bed = -3.0  # Flusstiefe
flt1_top = -10.0  # Tiefe Oberkante 1. Filter
flt2_top = -21.0  # Tiefe Oberkante 2. Filter
top_hol = -35.0  # Entspricht Top Holstein

flt1_len = 5.0  # Länge Filter flache Brunnen
flt2_len = 6.0  # Länge Filter tiefer Brunnen

frate_fl = -720.0  # Förderrate flache Brunnen
frate_tf = -600.0  # Förderrate tiefer Brunnen

#  pot_dif = (3.0e-5 * 500 * 100 * (0.05/40))/500 * 86400

#  print(pot_dif)

#  ws = './Szenario FR_flach {} FR_tief {}'.format(frate_fl, frate_tf)
ws = './test'  # Generates the model folder containing all model files

# Model layer geometry
if flt1_top > flt2_top:
    lay1 = top + riv_bed
    lay2 = flt1_top
    lay3 = lay2 - flt1_len
    lay4 = flt2_top
    lay5 = lay4 - flt2_len
    lay6 = lay5 + ((top_hol - lay5) / 2)
    lay7 = top_hol
    lay8 = -40
    lay9 = -50
elif flt1_top == flt2_top:
    lay1 = top + riv_bed
    lay2 = flt1_top
    lay3 = lay2 - flt1_len
    lay4 = lay3 - 5
    lay5 = lay4 - 5
    lay6 = lay5 + ((top_hol - lay5) / 2)
    lay7 = top_hol
    lay8 = -40
    lay9 = -50

if lay1 < lay2 or lay2 < lay3 or lay3 < lay4 or lay1 == lay2 or lay2 == lay3 or lay3 == lay4:
    raise Exception(
        "Layer Überschneidung! Werte für River Bed, Filterlängen oder Filter-Tops verursachen Unterschreitung oder Überlagerung der Basis von Layer 2, 3 oder 4 durch die Basis der Layer darüber.")
    sys.exit(1)

botm = [lay1, lay2, lay3, lay4, lay5, lay6, lay7, lay8, lay9]  # Layer bottom elevations ($m$)

icelltype_str = "1, 0, 0, 0, 0, 0, 0"  # Cell conversion type

# Well arrays
if flt1_top > flt2_top:
    wel_spd1 = [
        (2, 0, 24, frate_fl / 2),
        (2, 12, 24, frate_fl),
        (4, 24, 24, frate_tf),
        (2, 36, 24, frate_fl),
        (2, 49, 24, frate_fl / 2), ]

    wel_spd2 = []
    for k in range(49,50):
        for l in range(0, 7):
            wel_spd2 += [(l, i, k, 3.2) for i in range(nrow)]


    wells = wel_spd1 + wel_spd2
    wel_spd12 = {0: wells}

elif flt1_top == flt2_top:
    wel_spd1 = [
        (2, 0, 24, frate_fl / 2),
        (2, 12, 24, frate_fl),
        (2, 24, 24, frate_tf),
        (2, 36, 24, frate_fl),
        (2, 49, 24, frate_fl / 2), ]

wel_spd2 = []
for k in range(49, 50):
    for l in range(0, 7):
        wel_spd2 += [(l, i, k, 3.2) for i in range(nrow)]

wells = wel_spd1 + wel_spd2
wel_spd12 = {0: wells}

# River array
riv_spd = []
for k in range(10):
    riv_spd += [[0, i, k, top, 43.2, riv_bed] for i in range(nrow)]
riv_spd = {0: riv_spd}

if flt1_top > flt2_top:
    obs_spd = [
        (2, 0, 24),
        (2, 12, 24),
        (4, 24, 24),
        (2, 36, 24),
        (2, 49, 24)]

elif flt1_top == flt2_top:
    obs_spd = [
        (2, 0, 24),
        (2, 12, 24),
        (2, 24, 24),
        (2, 36, 24),
        (2, 49, 24)]

# Flow model section
modelname_mf = "flow"

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
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
ibound[8:, 0:49, 0] = -1
ibound[8:, 0:49, 49] = -1
#ibound[0:8, 0:49, 49] = 1
strt = np.ones((nlay, nrow, ncol), dtype=np.float32)
strt[8:, 0:49, 0] = 0.0
strt[8:, 0:49, 49] = 0.0  # Starting head ($m$)

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
flopy.modflow.ModflowSip(mf)

# Instantiate link mass transport package (for writing linker file)
flopy.modflow.ModflowLmt(mf)

spd = {(0, 0): ["print head", "print budget", "save head", "save budget"]}
oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True, save_specific_discharge=True)

# Transport model section
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

slt_spd = []
for k in range(ncol):
    for l in range(8,9):
        slt_spd += [(l, i, k, slt_cnc, -1) for i in range(nrow)]


criv_spd = []
for k in range(10):
    criv_spd += [[0, i, k, 0.05, -1] for i in range(nrow)]


const_spd = []
for i in range (nrow):
    for k in range(49,50):
        const_spd += [[l, i, k, 0.05, -1] for l in range (0,7)]

all_spd = slt_spd+criv_spd+const_spd
all_spd = {0: all_spd}

print(riv_spd)
print(slt_spd)

ssm = flopy.mt3d.Mt3dSsm(mt, stress_period_data=all_spd)
gcg = flopy.mt3d.Mt3dGcg(mt)

# Section for seawat module, currently not used in results
swt = flopy.seawat.Seawat(
    modflowmodel=mf,
    mt3dmodel=mt,
    modelname=modelname_mt,
    namefile_ext="nam_swt",
    model_ws=ws,
    exe_name="C:/development/bin/bin/swt_v4.exe"
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

# Write model inputs
mf.write_input()
mt.write_input()

#  Seawat specific code section for seawat files
#  fname = modelname_mt + ".vdf"
#  f = open(os.path.join(ws, fname), "r")
#  lines = f.readlines()
#  f.close()
#  f = open(os.path.join(ws, fname), "w")
#  for line in lines:
#  f.write(line)
#  for kper in range(nper):
#  f.write("-1\n")
#  f.close()

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

# Function for hydraulic plot
def plot_mf():
    fig, axes = plt.subplots(3, 3, figsize=figure_size, dpi=300, constrained_layout=True, )
    extents = (0, ncol * delc, 0, nrow * delr)
    vmin, vmax = -1, 1

    for ax in axes.flatten():
        ax.set_aspect("equal")
        ax.set_xlim(extents[:2])
        ax.set_ylim(extents[:2])
    for idx, ax in enumerate(axes.flatten()[:nlay]):
        fmp = flopy.plot.PlotMapView(model=mf, ax=ax, layer=idx, extent=extents)
        fmp.plot_grid(lw=0.2)
        plot_obj = fmp.plot_array(head, vmin=vmin, vmax=vmax)
        # fmp.plot_bc("RIV", color="red", plotAll=True)
        # fmp.plot_bc("GHB", color="red", plotAll=True)
        fmp.plot_bc("WEL", color="red")
        #cv = fmp.contour_array(head, levels=[-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3], linewidths=0.2,
                               #colors="black", )
        #plt.clabel(cv, fmt="%1.0f", fontsize=4)
        # fmp.plot_specific_discharge(spdis, normalize=True, color="0.75", scale=15, istep=2, jstep=4, headwidth=6)
        title = "Model Layer {}".format(idx + 1)
        letter = chr(ord("@") + idx + 1)
        ax.set_title(title, fontsize=4, pad=2)
        ax.tick_params(width=0.2, length=0.5)
        plt.setp(ax.spines.values(), linewidth=0.2)
        ax.set_ylabel(r"$Meter$", fontsize=4)
        quiver = fmp.plot_discharge(frf, fff, head=head, normalize=True, color="0.75", scale=15, istep=2, jstep=4,
                                    headwidth=6)
        # fmp.plot_bc(slt_spd, color="red", plotAll=True)
    cbar = plt.colorbar(plot_obj, shrink=0.3, orientation="horizontal", )
    cbar.outline.set_linewidth(0.2)
    cbar.ax.tick_params(size=0)
    cbar.ax.set_xlabel(r"Hydraulisches Potential, $m$", fontsize=5)

    plt.show()

# Call function for hydraulic plot
plot_mf()

print(wel_spd1)

# Array of different hydraulic conductivities in the model
hk = mf.lpf.hk.array
vka = mf.lpf.vka.array

# Function to plot cross section with hydraulic conductivities and wells etc.
def plot_cross():
    fig = plt.figure(figsize=(15, 5))
    ax = fig.add_subplot(1, 1, 1)
    vmin, vmax = -1, 1
    # Next we create an instance of the PlotCrossSection class
    xsect = flopy.plot.PlotCrossSection(model=mf, line={"ROW": 24})

    # Then we can use the plot_grid() method to draw the grid
    # The return value for this function is a matplotlib LineCollection object,
    # which could be manipulated (or used) later if necessary.
    xsect.plot_grid(lw=0.2)
    xsect.plot_array(hk, masked_values=[hk[0, 0, 0]], alpha=0.2)
    # xsect.plot_array(vka, masked_values=[hk[0, 0, 0]], alpha=0.2)
    xsect.plot_array(ibound, masked_values=[ibound[0, 0, 0]], alpha=0.2)
    #xsect.plot_array(strt, masked_values=[strt[0, 0, 0]], alpha=0.2)
    # plot_obj = xsect.plot_array(head, vmin=vmin, vmax=vmax)
    #cv = xsect.contour_array(head, levels=[-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3], linewidths=0.2,
                             #colors="black", )
    xsect.plot_bc("WEL", color="red")
    #quiver = xsect.plot_discharge(frf, fff, flf, head=head, normalize=True, color="0.75", scale=50, headwidth=3,
                                  #linewidths=0.1)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    t = ax.set_title("Ost-West Querschnitt in der Mitte des KHM", fontsize=15)
    plt.xlabel("Länge des Modells in Meter", fontsize=15)
    plt.ylabel("Tiefe unter Geländeoberkante in Meter", fontsize=15)
    plt.rcParams.update({'font.size': 15})
    #  cbar = plt.colorbar(plot_obj, shrink=0.3, orientation="horizontal")
    #  cbar.outline.set_linewidth(0.2)
    #  cbar.ax.tick_params(size=0)
    #  cbar.ax.set_xlabel(r"Hydraulisches Potential, $m$")
    plt.show()

# Call cross section function
plot_cross()


# Layer budget section of flow model, not really needed but maybe useful later
def get_layerbudget(modelname,
                    nper,
                    perlen,
                    nlay,
                    model_ws='.',
                    debug=True
                    ):
    perlen = np.array(perlen, ndmin=1, copy=False)

    bud_agg = pd.DataFrame(columns=['stress_period',
                                    'time_step',
                                    'layer',
                                    'STORAGE_IN',
                                    'STORAGE_OUT',
                                    'CONSTANT_HEAD_IN',
                                    'CONSTANT_HEAD_OUT',
                                    'FLOW_RIGHT_FACE',
                                    'FLOW_FRONT_FACE',
                                    'FLOW_LOWER_FACE', ])
    cbb = bf.CellBudgetFile(os.path.join(model_ws, modelname + '.cbc'))
    for stress_period in np.arange(0, nper).astype('int'):
        for time_step in np.arange(0, perlen[stress_period]).astype('int'):
            stressperiod_timestep = cbb.get_kstpkper()

    for stress_period in [item[0] for item in stressperiod_timestep]:
        for time_step in [item[1] for item in stressperiod_timestep]:
            bud = cbb.get_data(kstpkper=(time_step, stress_period),
                               full3D=True)
            for layer in np.arange(0, nlay).astype('int'):
                if debug:
                    print("Stress period: " + str(stress_period + 1) + ", Time step: " + str(
                        time_step + 1) + ", Layer: " + str(layer + 1))
                tmp = pd.DataFrame([[stress_period,
                                     time_step + 1,
                                     layer,
                                     bud[0][layer][bud[0][layer] > 0].sum(),
                                     bud[0][layer][bud[0][layer] < 0].sum(),
                                     bud[1][layer][bud[1][layer] > 0].sum(),
                                     bud[1][layer][bud[1][layer] < 0].sum(),
                                     bud[2][layer].sum(),
                                     bud[3][layer].sum(),
                                     bud[4][layer].sum()
                                     ]],
                                   columns=['stress_period',
                                            'time_step',
                                            'layer',
                                            'STORAGE_IN',
                                            'STORAGE_OUT',
                                            'CONSTANT_HEAD_IN',
                                            'CONSTANT_HEAD_OUT',
                                            'FLOW_RIGHT_FACE',
                                            'FLOW_FRONT_FACE',
                                            'FLOW_LOWER_FACE'])
                bud_agg = bud_agg.append(tmp, ignore_index=True)
    # bud_agg.loc[:,['CONSTANT_HEAD_IN']] = bud_agg['CONSTANT_HEAD_IN'].as_matrix().astype("float32")
    # bud_agg.loc[:,['CONSTANT_HEAD_OUT']] = bud_agg['CONSTANT_HEAD_OUT'].as_matrix().astype("float32")
    # bud_agg.loc[np.isnan(bud_agg['CONSTANT_HEAD_IN']),['CONSTANT_HEAD_IN']] = 0
    # bud_agg.loc[np.isnan(bud_agg['CONSTANT_HEAD_OUT']),['CONSTANT_HEAD_OUT']] = 0
    # bud_agg.loc[:,['CONSTANT_HEAD_IN']] = bud_agg['CONSTANT_HEAD_IN'].as_matrix().astype("float32")
    # bud_agg.loc[:,['CONSTANT_HEAD_OUT']] = bud_agg['CONSTANT_HEAD_OUT'].as_matrix().astype("float32")
    # bud_agg.loc[np.isnan(bud_agg['CONSTANT_HEAD_IN']),['CONSTANT_HEAD_IN']] = 0
    # bud_agg.loc[np.isnan(bud_agg['CONSTANT_HEAD_OUT']),['CONSTANT_HEAD_OUT']] = 0
    return bud_agg


### Get layer based budget for each time step
layer_budget = get_layerbudget(modelname=modelname_mf,
                               nlay=mf.dis.nlay,
                               model_ws=ws, nper=nper, perlen=perlen)

### Aggregate budget for layer for whole simulation
#  layer_budget_perLayer = layer_budget.groupby(['layer']).sum().reset_index()

### Aggregate budget for stress period & layer
#  layer_budget_perStressPeriod = layer_budget.groupby(['layer', 'stress_period']).sum().reset_index()

### Filter only lowest layer
#  layer3_budget_perStressPeriod = layer_budget_perStressPeriod[layer_budget_perStressPeriod['layer'] == 6]

#  mf_list = flopy.utils.MfListBudget(os.path.join(ws, modelname_mf + ".list"))
#  budget_incremental, budget_cumulative = mf_list.get_dataframes(start_datetime='31-12-2006')

# Show layer budget plot
layer_budget[[
    'CONSTANT_HEAD_IN',
    'CONSTANT_HEAD_OUT',
    'FLOW_RIGHT_FACE',
    'FLOW_FRONT_FACE',
    'FLOW_LOWER_FACE']].plot(kind="bar")

plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.legend(fontsize=10)
plt.show()

# layer3_budget_perStressPeriod['MNW2_IN'] = np.append(budget_cumulative['MNW2_IN'][0],
# budget_cumulative['MNW2_IN'].diff().as_matrix()[1:])

# layer3_budget_perStressPeriod['MNW2_OUT'] = np.append(budget_cumulative['MNW2_OUT'][0],
# budget_cumulative['MNW2_OUT'].diff().as_matrix()[1:])

# Run transport model
mt.run_model()

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

x = cvt["time"]
y = cvt["(3, 1, 25)"] * 1000 + 50.0
plt.plot(x, y, label="Brunnen A")

x = cvt["time"]
y = cvt["(3, 13, 25)"] * 1000 + 50.0
plt.plot(x, y, label="Brunnen B")

x = cvt["time"]
y = cvt["(5, 25, 25)"] * 1000 + 50.0
plt.plot(x, y, label="Brunnen C")

x = cvt["time"]
y = cvt["(3, 37, 25)"] * 1000 + 50.0
plt.plot(x, y, label="Brunnen D")

x = cvt["time"]
y = cvt["(3, 50, 25)"] * 1000 + 50.0
plt.plot(x, y, label="Brunnen E")

plt.xticks(fontsize=15)
plt.yticks(fontsize=12)
#  plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.])
plt.yticks([50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400,
            1500, 1600, 1700, 1800, 1900, 2000])
plt.xlim(0, perlent)
plt.ylim(0, 2000)
plt.xlabel("ZEIT IN TAGEN", fontsize=15)
plt.ylabel(("Chlorid-Konzentration in mg/L").format(slt_cnc), fontsize=15)
plt.rcParams.update({'font.size': 15})
plt.grid(visible=True)
title1 = "Chlorid-Konzentrationen in der Basisvariante des KHM Friedrichshagen"
title2 = "Filterlänge: Brunnen A, B, D, E = {} m, Brunnen C = {} m und " \
         "Tiefe Filteroberkante: Brunnen A, B, D, E = {} m, Brunnen C = {} m".format(flt1_len, flt2_len, flt1_top,
                                                                                     flt2_top)
plt.suptitle(title1, fontsize=18)
plt.title(title2, pad=15, fontsize=12)
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
arr = xsect.plot_array(concentration)
xsect.plot_grid(lw=0.2)
xsect.plot_bc("WEL", color="red")
xsect.plot_array(hk, masked_values=[hk[0, 0, 0]], alpha=0.2)
cv = xsect.contour_array(head, levels=[-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3], linewidths=0.2,
                         colors="black", )
ax.set_title("Simulierte Konzentration bei Pumprate tiefer Brunnen = {}".format(frate_tf), pad=15)
quiver = xsect.plot_discharge(frf, fff, flf, head=head, normalize=True, color="0.75", scale=50, headwidth=3,
                              linewidths=0.1)
plt.rcParams.update({'font.size': 15})
cbar = plt.colorbar(arr, shrink=0.3, orientation="horizontal")
cbar.outline.set_linewidth(0.2)
cbar.ax.tick_params(size=0)
cbar.ax.set_xlabel(r"Konzentration in mg/L")
plt.show()



