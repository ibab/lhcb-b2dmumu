from root_pandas import *

simple_mc = read_root('../store/tmp/SIM_B_JpsiKst_ALL.Reduce.Cut.AddHypos.Cut.KFoldTrainAndApplyExtra.root', columns=['clf'])
#signal_mc = read_root('../store/tmp/SIM_B_Dmumu_ALL.Reduce.Cut.AddHypos.Cut.KFoldTrainAndApplyExtra.root', columns=['clf'])
resampled_mc = read_root('../store/tmp/SIM_B_JpsiKst_ALL.Reduce.Cut.AddHypos.Cut.ResamplePID.CalcSimpleWeights.KFoldTrainAndApplyExtra.root', columns=['clf', 'SimpleWeight'])
reweighted_mc = read_root('../store/tmp/SIM_B_JpsiKst_ALL.Reduce.Cut.AddHypos.Cut.CalcSuperWeights.KFoldTrainAndApplyExtra.root', columns=['SuperWeight', 'clf'])
#reweighted_signal_mc = read_root('../store/tmp/SIM_B_Dmumu_ALL.Reduce.Cut.AddHypos.Cut.ApplySuperWeights.KFoldTrainAndApplyExtra.root', columns=['SuperWeight', 'clf'])
norm_data = read_root('../store/tmp/DATA_B_JpsiKst.Reduce.Cut.AddHypos.Cut.CalcSWeights.KFoldTrainAndApplyExtra.root', columns=['sigYield_sw', 'clf'])

from scale import mksize
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from scipy.interpolate import interp1d
import numpy as np

datasets = [
        (simple_mc, None, 'b'),
#        (signal_mc, None),
        #(resampled_mc, resampled_mc.SimpleWeight.ravel(), 'r'),
        (resampled_mc, resampled_mc.SimpleWeight.ravel(), 'r'),
        (reweighted_mc, reweighted_mc.SuperWeight.ravel(), 'g'),
        #(reweighted_signal_mc, reweighted_signal_mc.SuperWeight.ravel(), 'black'),
]

print(norm_data)
_, y_, x_ = roc_curve(np.ones(len(norm_data)), norm_data.clf, sample_weight=norm_data.sigYield_sw)
curve_ = interp1d(x_, y_)

curves = []

print(curve_(3.905))

for d in datasets:
    _, y, x = roc_curve(np.ones(len(d[0])), d[0].clf, sample_weight=d[1])
    curves.append((interp1d(x, y), d[2]))
    print(interp1d(x, y)(3.905))

left = -2
right = 4.5

plt.figure(figsize=mksize(0.8))
x = np.linspace(left, right, 300)
plt.axhline(1, color='black')
for c, color in curves:
    plt.plot(x, c(x) / curve_(x), color=color)
plt.axvline(3.905, color='black', linestyle='dashed')

plt.xlim(left, right)
plt.ylim(0.95, 1.4)
plt.xlabel('Classifier threshold')
plt.ylabel('Signal efficiency relative to data')

plt.savefig('../store/tmp/data_mc_reweighted.pdf')
plt.savefig('../store/tmp/data_mc_reweighted.pgf')

