
import os
from root_pandas import *

blue = (0.21568627450980393, 0.47058823529411764, 0.7490196078431373)

df = read_root('../store/tmp/DATA_B_Dmumu_ALL.Reduce.Cut.AddHypos.root', columns=['B_M', 'Jpsi_M', 'Dstar_M'])
mc = read_root('../store/tmp/SIM_B_Dmumu_ALL.Reduce.Cut.AddHypos.root', columns=['B_M', 'Jpsi_M', 'Dstar_M'])

from scale import mksize
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

if not os.path.exists('../store/tmp/vetoes'):
    os.makedirs('../store/tmp/vetoes')

print('Documenting D* veto')
plt.figure(figsize=mksize(0.8))
plt.hist(df.Dstar_M, range=(1900, 2100), bins=100, histtype='stepfilled', color=blue)
plt.axvspan(1990, 2030, alpha=0.2, color='red')
plt.axvline(1990, alpha=0.8, color='red')
plt.axvline(2030, alpha=0.8, color='red')
plt.xlabel('$m(K^+\\pi^-(\\mu^- \\to \\pi^-))\\ /\\ \\mathrm{MeV}$')
plt.tight_layout()
plt.savefig('../store/tmp/vetoes/Dstar.pdf')
plt.savefig('../store/tmp/vetoes/Dstar.pgf')
plt.clf()

print('Documenting J/psi vetoes')
q2 = (df.Jpsi_M / 1000)**2
plt.figure(figsize=mksize(0.8, 1.0))
col = plt.hexbin(df.B_M, q2, cmap=plt.cm.Blues, gridsize=100, linewidths=(0.15,), norm=LogNorm())
ax = plt.gca()
plt.xlabel('$m(K^+\\pi^-\\mu^+\\mu^-)\ /\ \mathrm{MeV}$')
plt.ylabel('$q^2(\\mu^+\\mu^-)\ /\ \mathrm{GeV}^2$')
lower = 2.9**2
upper = 3.2**2
edge = 3.5**2
plt.axhspan(lower, upper, color='red', alpha=0.1)
plt.axhline(lower, color='red')
plt.axhline(upper, color='red')
plt.axhspan(edge, q2.max(), color='red', alpha=0.1)
plt.axhline(edge, color='red')
#plt.ylim(q2.min(), q2.max())
ax = plt.gca()
#ax.set_aspect('equal')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
bar = plt.colorbar(cax=cax, norm=col)
bar.set_label('Candidates')
plt.savefig('../store/tmp/vetoes/B_vs_Jpsi.pdf')
plt.savefig('../store/tmp/vetoes/B_vs_Jpsi.pgf')
plt.clf()

plt.figure(figsize=mksize(0.8))
plt.hist((mc.Jpsi_M/1000)**2, bins=100, histtype='stepfilled', color=blue)
plt.axvspan(lower, upper, color='red', alpha=0.1)
plt.axvspan(edge, 14, color='red', alpha=0.1)
plt.axvline(lower, color='red')
plt.axvline(upper, color='red')
plt.axvline(edge, color='red')
plt.ylabel('Simulated candidates')
plt.xlabel('$q^2(\\mu^+\\mu^-)\\ /\\ \\mathrm{GeV}^2$')
plt.savefig('../store/tmp/vetoes/MC_Jpsi.pdf')
plt.savefig('../store/tmp/vetoes/MC_Jpsi.pgf')
plt.clf()

print('Documenting swapped J/psiK* veto')

import numpy as np

df = read_root('../store/tmp/DATA_B_JpsiKst.Reduce.Cut.AddHypos.Cut.CalcSWeights.root', columns=['*_P{X,Y,Z}', '*_M', 'sigYield_sw'])

kaon_m = 493.667
pion_m = 139.570
muon_m = 105.659
dzero_m = 1864.84

PX = df.Kplus_PX + df.piminus_PX + df.muplus_PX + df.muminus_PX
PY = df.Kplus_PY + df.piminus_PY + df.muplus_PY + df.muminus_PY
PZ = df.Kplus_PZ + df.piminus_PZ + df.muplus_PZ + df.muminus_PZ

PE  = np.sqrt(kaon_m**2 + df.Kplus_PX**2 + df.Kplus_PY**2 + df.Kplus_PZ**2)
PE += np.sqrt(muon_m**2 + df.piminus_PX**2 + df.piminus_PY**2 + df.piminus_PZ**2)
PE += np.sqrt(muon_m**2 + df.muplus_PX**2 + df.muplus_PY**2 + df.muplus_PZ**2)
PE += np.sqrt(pion_m**2 + df.muminus_PX**2 + df.muminus_PY**2 + df.muminus_PZ**2)

m = np.sqrt(PE**2 - PX**2 - PY**2 - PZ**2)

Meson_PX = df.Kplus_PX + df.muminus_PX
Meson_PY = df.Kplus_PY + df.muminus_PY
Meson_PZ = df.Kplus_PZ + df.muminus_PZ

Meson_PE = np.sqrt(kaon_m**2 + df.Kplus_PX**2 + df.Kplus_PY**2 + df.Kplus_PZ**2)
Meson_PE += np.sqrt(pion_m**2 + df.muminus_PX**2 + df.muminus_PY**2 + df.muminus_PZ**2)

Meson_m = np.sqrt(Meson_PE**2 - Meson_PX**2 - Meson_PY**2 - Meson_PZ**2)

plt.figure(figsize=mksize(0.8))
plt.hist(Meson_m.ravel(), bins=100, histtype='stepfilled', color=blue, weights=df.sigYield_sw.ravel())
plt.axvspan(500, dzero_m - 100, color='red', alpha=0.1)
plt.axvspan(dzero_m + 100, 4500, color='red', alpha=0.1)
plt.axvline(dzero_m - 100, color='red')
plt.axvline(dzero_m + 100, color='red')
plt.xlim(500, 4500)
plt.ylim(bottom=0)
plt.ylabel('Candidates')
plt.xlabel('$m(K^+\\mu^-_\\pi)\\ /\\ \\mathrm{MeV}$')
plt.savefig('../store/tmp/vetoes/Meson_swapped.pdf')
plt.savefig('../store/tmp/vetoes/Meson_swapped.pgf')
plt.clf()

select = (Meson_m > dzero_m - 100) & (Meson_m < dzero_m + 100)

plt.figure(figsize=mksize(0.8))
plt.hist(m[select].ravel(), bins=100, histtype='stepfilled', color=blue, weights=df.sigYield_sw[select].ravel())
plt.ylim(bottom=0)
plt.ylabel('Candidates')
plt.xlabel('$m(K^+\\mu^-_\\pi)\\ /\\ \\mathrm{MeV}$')
plt.savefig('../store/tmp/vetoes/B_swapped.pdf')
plt.savefig('../store/tmp/vetoes/B_swapped.pgf')
plt.clf()


