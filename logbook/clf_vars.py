
import os
from root_pandas import *

blue = (0.21568627450980393, 0.47058823529411764, 0.7490196078431373)

norm_data = read_root('../store/tmp/DATA_B_JpsiKst.Reduce.Cut.AddHypos.Cut.CalcSWeights.KFoldTrainAndApplyExtra.root')
norm_mc = read_root('../store/tmp/SIM_B_JpsiKst_ALL.Reduce.Cut.AddHypos.Cut.KFoldTrainAndApplyExtra.root')
norm_resampled_mc = read_root('../store/tmp/SIM_B_JpsiKst_ALL.Reduce.Cut.AddHypos.Cut.ResamplePID.CalcSimpleWeights.KFoldTrainAndApplyExtra.root')
norm_reweighted_mc = read_root('../store/tmp/SIM_B_JpsiKst_ALL.Reduce.Cut.AddHypos.Cut.CalcSuperWeights.KFoldTrainAndApplyExtra.root')
sig_mc = read_root('../store/tmp/SIM_B_Dmumu_ALL.Reduce.Cut.AddHypos.Cut.ResamplePID.ApplySimpleWeights.KFoldTrainAndApplyExtra.root')
bkg_data = read_root('../store/tmp/DATA_B_Dmumu_ALL.Reduce.Cut.AddHypos.Cut.KFoldTrainAndApplyExtra.root')

from scale import mksize
import matplotlib as mpl
import matplotlib.pyplot as plt

# Plot all the variables!
variables = {
        'B_DiraAngle': (0, 0.014),
        'B_TAU': None,
        'B_ENDVERTEX_CHI2_NDOF': None,
        'B_ISOLATION_BDT_Soft': None,
        'Kplus_PIDK': (-15, 130),
        'Kplus_PIDmu': (-15, 15),
        'piminus_PIDK': (-130, 30),
        'piminus_PIDmu': (-15, 15),
        'muplus_PIDmu': None,
        'muminus_PIDmu': None,
        'muplus_PIDK': None,
        'muminus_PIDK': None,
        'B_TAU': (0, 10),
        'clf': None,
        'nSPDHits': None,
        'nTracks': None,
        'B_PT': None,
        'Kplus_P': None,
        'piminus_P': None,
        'muplus_P': None,
        'muminus_P': None,
        'Kplus_TRACK_CHI2NDOF': None,
        'piminus_TRACK_CHI2NDOF': None,
        'muplus_TRACK_CHI2NDOF': None,
        'muminus_TRACK_CHI2NDOF': None,
        #'Meson_TAU': None,
}

labels = {
        'clf': 'classifier response',
        'B_DiraAngle': r'$\mathrm{cos}(\mathrm{DIRA\ angle})$',
        'B_TAU': r'$t_{B^0}$',
        'B_ENDVERTEX_CHI2_NDOF': r'$B^0\ \mathrm{vertex}\ \chi^2 / \mathrm{ndf}$',
        'B_ISOLATION_BDT_Soft': 'muon isolation BDT response',
        'Kplus_PIDK': r'$\mathrm{DLL}_{K/\pi}(K^+)$',
        'Kplus_PIDmu': r'$\mathrm{DLL}_{\mu/\pi}(K^+)$',
        'piminus_PIDK': r'$\mathrm{DLL}_{K/\pi}(\pi^-)$',
        'piminus_PIDmu': r'$\mathrm{DLL}_{\mu/\pi}(\pi^-)$',
        'muplus_PIDmu': r'$\mathrm{DLL}_{\mu/\pi}(\mu^+)$',
        'muminus_PIDmu': r'$\mathrm{DLL}_{\mu/\pi}(\mu^-)$',
}

scale = {
        'clf': 0.8,
        }

if not os.path.exists('../store/tmp/variables'):
    os.makedirs('../store/tmp/variables')

mpl.rcParams['figure.figsize'] = (8, 5)
for v, rng in variables.items():
    print('Plotting {}'.format(v))
    kwargs = dict()
    if rng:
        kwargs['range'] = rng

    if v in labels:
        xlabel = labels[v]
    else:
        xlabel = None

    if v in scale:
        scl = scale[v]
    else:
        scl = 0.55
        plt.locator_params(axis='y', nbins=5)
        plt.locator_params(axis='x', nbins=6)
    
    plt.figure(figsize=mksize(scl, aspect=0.70))
    # Plot data-resampling comparison
    N, bins, _ = plt.hist(norm_data[v].ravel(), bins=50, histtype='stepfilled', color=blue, lw=0, weights=norm_data['sigYield_sw'].ravel(), normed=True, **kwargs)
    plt.hist(norm_mc[v].ravel(), bins=bins, histtype='step', color='red', lw=0.5, normed=True)
    plt.hist(norm_resampled_mc[v].ravel(), bins=bins, histtype='step', color='orange', lw=0.5, weights=norm_resampled_mc['SimpleWeight'].ravel(), normed=True)

    if xlabel:
        plt.xlabel(xlabel)
    else:
        plt.xlabel(v)

    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.savefig('../store/tmp/variables/DATA_MC_{}.pdf'.format(v))
    plt.savefig('../store/tmp/variables/DATA_MC_{}.pgf'.format(v))
    plt.clf()

    # Plot data-reweighting comparison
    N, bins, _ = plt.hist(norm_data[v].ravel(), bins=50, lw=0, histtype='stepfilled', color=blue, weights=norm_data['sigYield_sw'].ravel(), normed=True, **kwargs)
    plt.hist(norm_mc[v].ravel(), bins=bins, histtype='step', lw=0.5, color='red', normed=True)
    plt.hist(norm_reweighted_mc[v].ravel(), bins=bins, histtype='step', color='orange', weights=norm_reweighted_mc['SuperWeight'].ravel(), normed=True)

    if xlabel:
        plt.xlabel(xlabel)
    else:
        plt.xlabel(v)

    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.savefig('../store/tmp/variables/DATA_MC_REWEIGHT_{}.pdf'.format(v))
    plt.savefig('../store/tmp/variables/DATA_MC_REWEIGHT_{}.pgf'.format(v))
    plt.clf()

    # Plot bkg-sig comparison
    N, bins, _ = plt.hist(sig_mc[v].ravel(), bins=50, histtype='stepfilled', lw=0, color=blue, weights=sig_mc['SimpleWeight'].ravel(), normed=True, **kwargs)
    plt.hist(bkg_data[v].ravel(), bins=bins, histtype='step', color='red', lw=0.5, normed=True)

    if xlabel:
        plt.xlabel(xlabel)
    else:
        plt.xlabel(v)

    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.savefig('../store/tmp/variables/SIG_BKG_{}.pdf'.format(v))
    plt.savefig('../store/tmp/variables/SIG_BKG_{}.pgf'.format(v))
    plt.clf()

