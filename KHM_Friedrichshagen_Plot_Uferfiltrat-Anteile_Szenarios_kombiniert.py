import os
import matplotlib.pyplot as plt
import flopy

ws1 = './Szenario 1 FOK_flach -15.0 FOK_tief -21.0 Fluss'  # Pfad zum Modellordner des Szenarios
ws2 = './Szenario 0 Basismodell Fluss'  # Pfad zum Modellordner des
szn_nr = 1  # Nummerierung des Szenarios für die Legende des Plots

f1 = os.path.join(ws1, 'transport.nam')
mt1 = flopy.mt3d.mt.Mt3dms.load(f1)

f2 = os.path.join(ws2, 'transport.nam')
mt2 = flopy.mt3d.mt.Mt3dms.load(f1)

fname = os.path.join(ws1, "MT3D001.UCN")
ucnobj = flopy.utils.UcnFile(fname)
times = ucnobj.get_times()
conc = ucnobj.get_alldata()

fname = os.path.join(ws1, "MT3D001.OBS")
if os.path.isfile(fname):
    cvt = mt1.load_obs(fname)
else:
    cvt = None

fname = os.path.join(ws1, "MT3D001.MAS")
mvt = mt1.load_mas(fname)

fname2 = os.path.join(ws2, "MT3D001.UCN")
ucnobj2 = flopy.utils.UcnFile(fname2)
times2 = ucnobj2.get_times()
conc2 = ucnobj2.get_alldata()

fname2 = os.path.join(ws2, "MT3D001.OBS")
if os.path.isfile(fname2):
    cvt2 = mt2.load_obs(fname2)
else:
    cvt2 = None

fname2 = os.path.join(ws2, "MT3D001.MAS")
mvt2 = mt2.load_mas(fname2)

plt.rcParams.update({'font.size': 20})
x = cvt["time"]
y = cvt["(3, 1, 25)"]
plt.plot(x, y, label="Flache Brunnen Szenario {}".format(szn_nr),  color='blue')

x = cvt2["time"]
y = cvt2["(3, 1, 25)"]
plt.plot(x, y, label="Flache Brunnen Basismodell", color='blue', linestyle='dashed')

x = cvt["time"]
y = cvt["(5, 25, 25)"]
plt.plot(x, y, label="Tiefer Brunnen Szenario {}".format(szn_nr), color='red')


x = cvt2["time"]
y = cvt2["(5, 25, 25)"]
plt.plot(x, y, label="Tiefer Brunnen Basismodell", color='red', linestyle='dashed')

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#  plt.yticks([50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400,
           # 1500, 1600, 1700, 1800, 1900, 2000])
plt.xlim(0, 3500)
plt.ylim(0, 1)
plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.])
plt.xlabel("ZEIT IN TAGEN", fontsize=18)
plt.ylabel(("Konzentration in mg/L"), fontsize=18)
#  plt.rcParams.update({'font.size': 18})
plt.grid(visible=True)
title1 = "Uferfiltrat-Anteile"

plt.legend(loc=2)

plt.show()
