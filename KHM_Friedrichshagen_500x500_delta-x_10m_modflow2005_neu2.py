import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import flopy
import flopy.utils.binaryfile as bf

sys.path.append(os.path.join("..", "common"))

plt.rcParams.update({'font.size': 3})

# plt.rcParams.update({'lines.linewidth': 0.2})
# from figspecs import USGSFigure
# sys.path.append(os.path.join("..", "common"))

figure_size = (12, 12)
ws = './flow_model'
parameter_units = {"recharge": "$m/s$", "slt_cnc": "$kg m^{-3}$", "slt_strt": "$kg m^{-3}$", "perlen": "$s$"}
slt_cnc = 3.0
slt_strt = 0.2
porosity = 0.30
#length_units = "meter"
#time_units = "years"
perlen = 315360
perlent = 5*864000
steady = [False]
nstp = 10

nper = 1  # Number of periods
nlay = 6  # Number of layers
ncol = 50  # Number of columns
nrow = 50  # Number of rows
delr = 10.0  # Column width ($m$)
delc = 10.0  # Row width ($m$)
top = 0.0  # Top of the model ($m$)

riv_bed = -3.0  # Flusstiefe
flt1_top = -20.0  # Tiefe Oberkante 1. Filter
flt2_top = -25.0  # Tiefe Oberkante 2. Filter
top_hol = -35.0

flt1_len = 2.0  # Länge 1. Filter
flt2_len = 5.0  # Länge 2. Filter

if flt1_top > flt2_top:
    lay1 = top + riv_bed
    lay2 = flt1_top
    lay3 = lay2 - flt1_len
    lay4 = flt2_top
    lay5 = lay4 - flt2_len
    lay6 = top_hol
elif flt1_top == flt2_top:
    lay1 = top + riv_bed
    lay2 = flt1_top
    lay3 = lay2 - flt1_len
    lay4 = lay3 - 10
    lay5 = lay4 - 10
    lay6 = top_hol

if lay1 < lay2 or lay2 < lay3 or lay3 < lay4 or lay1 == lay2 or lay2 == lay3 or lay3 == lay4:
    raise Exception(
        "Layer Überschneidung! Werte für River Bed, Filterlängen oder Filter-Tops verursachen Unterschreitung oder Überlagerung der Basis von Layer 2, 3 oder 4 durch die Basis der Layer darüber.")
    sys.exit(1)

botm = [lay1, lay2, lay3, lay4, lay5, lay6]  # Layer bottom elevations ($m$)

icelltype_str = "1, 0, 0, 0, 0, 0"  # Cell conversion type
k11 = 3.0e-4  # Horizontal hydraulic conductivity ($m/s$)
k33 = 3.0e-5  # Vertical hydraulic conductivity ($m/s$)
recharge = 3e-8  # Recharge rate ($m/s$)

tdis_ds = ((31536000, 10,
            1),)  # (length of a stress period, number of time steps in a stress period, multiplier for the length of successive time steps)

icelltype = tuple([int(value) for value in icelltype_str.split(",")])

chd_spd = []
for k in range(6):
    chd_spd += [[(k, i, 49), 0.0] for i in range(nrow)]
chd_spd = {0: chd_spd}

if flt1_top > flt2_top:
    wel_spd = {
        0: [
            [2, 0, 24, -0.0058],
            [2, 12, 24, -0.0116],
            [4, 24, 24, -0.01116],
            [2, 36, 24, -0.0116],
            [2, 49, 24, -0.0058],
        ]
    }
elif flt1_top == flt2_top:
    wel_spd = {
        0: [
            [2, 0, 24, -0.0058],
            [2, 12, 24, -0.0116],
            [2, 24, 24, -0.01116],
            [2, 36, 24, -0.0116],
            [2, 49, 24, -0.0058],
        ]
    }

riv_spd = []
for k in range(10):
    riv_spd += [[0, i, k, top, 10000.0, top - 3] for i in range(nrow)]
riv_spd = {0: riv_spd}

# Flow model

modelname_mf = "flow"
# model_ws = os.path.join(ws, "mfgwf")

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
    nstp=nstp
)

ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
# ibound[:, :, 0] = -1
ibound[:, :, -1] = -1
strt = np.ones((nlay, nrow, ncol), dtype=np.float32)
# strt[:, :, 0] = 0.0
strt[:, :, -1] = 0.0  # Starting head ($m$)

bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)

# Instantiate layer property flow package
flopy.modflow.ModflowLpf(mf, hk=k11, vka=k33, ipakcb=53)

# Instantiate well package
flopy.modflow.ModflowWel(mf, stress_period_data=wel_spd)

# Instantiate river package
flopy.modflow.ModflowRiv(mf, stress_period_data=riv_spd)

# Instantiate solver package
flopy.modflow.ModflowSip(mf)

# Instantiate link mass transport package (for writing linker file)
flopy.modflow.ModflowLmt(mf)

spd = {(0, 0): ["print head", "print budget", "save head", "save budget"]}
oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True, save_specific_discharge=True)

mf.write_input()
mf.run_model()

# Extract the heads

hds = bf.HeadFile("flow_model/flow.hds")
times = hds.get_times()
head = hds.get_data(totim=times[-1])

cbb = bf.CellBudgetFile("flow_model/flow.cbc")
kstpkper_list = cbb.get_kstpkper()
frf = cbb.get_data(text="FLOW RIGHT FACE", totim=times[-1])[0]
fff = cbb.get_data(text="FLOW FRONT FACE", totim=times[-1])[0]

def plot_mf():
    fig, axes = plt.subplots(2, 3, figsize=figure_size, dpi=300, constrained_layout=True, )
    extents = (0, ncol * delc, 0, nrow * delr)
    vmin, vmax = -3, 1

    for ax in axes.flatten():
        ax.set_aspect("equal")
        ax.set_xlim(extents[:2])
        ax.set_ylim(extents[:2])
    for idx, ax in enumerate(axes.flatten()[:nlay]):
        fmp = flopy.plot.PlotMapView(model=mf, ax=ax, layer=idx, extent=extents)
        fmp.plot_grid(lw=0.2)
        plot_obj = fmp.plot_array(head, vmin=vmin, vmax=vmax)
        fmp.plot_bc("RIV", color="red", plotAll=True)
        #fmp.plot_bc("GHB", color="red", plotAll=True)
        fmp.plot_bc("WEL", color="red")
        cv = fmp.contour_array(head, levels=[-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3], linewidths=0.2,
                           colors="black", )
        plt.clabel(cv, fmt="%1.0f", fontsize=4)
        #fmp.plot_specific_discharge(spdis, normalize=True, color="0.75", scale=15, istep=2, jstep=4, headwidth=6)
        title = "Model Layer {}".format(idx + 1)
        letter = chr(ord("@") + idx + 1)
        ax.set_title(title, fontsize=4)
        ax.tick_params(width=0.2, length=0.5)
        plt.setp(ax.spines.values(), linewidth=0.2)
        ax.set_ylabel(r"$m$", fontsize=4)
        quiver = fmp.plot_discharge(frf, fff, head=head, normalize=True, color="0.75", scale=15, istep=2, jstep=4,
                                 headwidth=6)

    plt.show()

#plot_mf()

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
    nper=1,
    nstp=nstp,
    perlen=perlent,
    obs=[(0, 24, 11)]
)

     
dceps = 1.0e-5
nplane = 1
npl = 0
nph = 4
npmin = 0
npmax = 8
nlsink = nplane
npsink = nph
adv = flopy.mt3d.Mt3dAdv(mt)
    #dceps=dceps,
    #nplane=nplane,
    #npl=npl,
    #nph=nph,
    #npmin=npmin,
    #npmax=npmax,
    #nlsink=nlsink,
    #npsink=npsink,
    #percel=0.5,
#)

al = 0.5
dsp = flopy.mt3d.Mt3dDsp(mt, al=al)

mt.write_input()

slt_spd = []
for k in range(30, 40):
    slt_spd += [[5, i, k, slt_cnc, -1] for i in range(30, 50)]
slt_spd = {0: slt_spd}

criv_spd = []
for k in range(10):
    criv_spd += [[0, i, k, slt_cnc, -1] for i in range(nrow)]
criv_spd = {0: criv_spd}


print(riv_spd)
#print(slt_spd)

ssm = flopy.mt3d.Mt3dSsm(mt, stress_period_data=criv_spd)
gcg = flopy.mt3d.Mt3dGcg(mt)
mt.write_input()

fname = os.path.join(model_ws, "MT3D001.UCN")
if os.path.isfile(fname):
    os.remove(fname)
mt.run_model()

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

# mf, mt, conc, cvt, mvt = p02("freundlich", 2, 0.3, 0.7, -1)
x = cvt["time"]/86400
y = cvt["(1, 25, 12)"]
plt.rcParams.update({'font.size': 10})
plt.plot(x, y, label="Salt")
plt.xlim(0, perlent)
plt.ylim(0, 4.0)
plt.xlabel("TIME IN SECONDS", fontsize=10)
plt.ylabel("CONCENTRATION", fontsize=10)
plt.legend()
plt.show()

