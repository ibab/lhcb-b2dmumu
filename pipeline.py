#!/usr/bin/env python2
"""
Analysis pipeline:

 *
 |
 Reduce dataset (Throw away unused variables)
 |
 Blind signal peak (Remove events with B_M close to B0 mass)
 |
 Apply cut-based selection
 |
 |--------> Plot variables
 |
 ?

"""


import logging
from log import setup_logging; setup_logging()

import matplotlib
matplotlib.use('Agg')
from matplotlib.mlab import rec_append_fields as append_fields

# Simulated data fetched from the grid
input_mc = [
    './storage/MC/2012/Stripping20/AllStreams/DVBd2MuMuD0_MC/combined/Bd2MuMuD0.root'
]

# Real data fetched from the grid
input_data = [
    './storage/Data/AllYears/Stripping20/Dimuon/DVBd2MuMuD0_data/combined/Bd2MuMuD0.root',
]

# Variables that are used in the analysis
variables = [
        'B_M',
        'B_DIRA_OWNPV',
        'B_FD_OWNPV',
        'B_ENDVERTEX_CHI2',
        'B_P',
        'B_PT',
        'B_ISOLATION_BDT_Hard',
        'B_ISOLATION_BDT_Soft',
        'B_OWNPV_CHI2',
        'Psi_M',
        'Psi_FD_ORIVX',
        'Psi_FDCHI2_ORIVX',
        'D~0_M',
        'D~0_FD_ORIVX',
        'D~0_CosTheta',
        'D~0_DIRA_OWNPV',
        'muplus_isMuon',
        'muplus_ProbNNpi',
        'muplus_ProbNNmu',
        'muminus_isMuon',
        'muminus_ProbNNpi',
        'muminus_ProbNNmu',
        'Kplus_hasRich',
        'Kplus_ProbNNpi',
        'Kplus_ProbNNk',
        'Kplus_ProbNNmu',
        'Kplus_ProbNNp',
        'Kplus_PX',
        'Kplus_PY',
        'Kplus_PZ',
        'Kplus_PT',
        'piminus_PX',
        'piminus_PY',
        'piminus_PZ',
        'piminus_PT',
        'Psi_PX',
        'Psi_PY',
        'Psi_PZ',
        'Psi_PE',
        'piminus_hasRich',
        'piminus_ProbNNpi',
        'piminus_ProbNNk',
        'piminus_ProbNNmu',
        'piminus_ProbNNp',
        'B_L0MuonDecision_TOS',
        'B_Hlt1TrackAllL0Decision_TOS',
        'B_Hlt1TrackMuonDecision_TOS',
        'B_Hlt2Topo2BodyBBDTDecision_TOS',
        'B_Hlt2Topo3BodyBBDTDecision_TOS',
        'B_Hlt2Topo4BodyBBDTDecision_TOS',
        'B_Hlt2TopoMu2BodyBBDTDecision_TOS',
        'B_Hlt2TopoMu3BodyBBDTDecision_TOS',
        'B_Hlt2TopoMu4BodyBBDTDecision_TOS',
        'Psi_Hlt2SingleMuonDecision_TOS',
        'B_Hlt2DiMuonDetachedDecision_TOS',
        'Kplus_PIDK',
        'piminus_PIDK',
        'Kplus_TRACK_GhostProb',
        'piminus_TRACK_GhostProb',
        'muplus_TRACK_GhostProb',
        'muminus_TRACK_GhostProb',
        'Kplus_TRACK_CHI2NDOF',
        'piminus_TRACK_CHI2NDOF',
        'muplus_TRACK_CHI2NDOF',
        'muminus_TRACK_CHI2NDOF',
        'Kplus_isMuonLoose',
        'piminus_isMuonLoose',
]

# Variables that are used in the multivariate classification
bdt_variables = [
        'B_TAU',
        'B_ISOLATION_BDT_Soft',
        'B_ENDVERTEX_CHI2',
        'B_DIRA_OWNPV',
        'D~0_DIRA_OWNPV',
        'muplus_ProbNNmu',
        'muminus_ProbNNmu',
        'Kplus_ProbNNpi',
        'Kplus_ProbNNk',
        'Kplus_ProbNNp',
        'piminus_ProbNNpi',
        'piminus_ProbNNp',
        'Psi_FD_ORIVX',
        'Kplus_TRACK_GhostProb',
        'piminus_TRACK_GhostProb',
        'muplus_TRACK_GhostProb',
        'muminus_TRACK_GhostProb',
]

mc_variables = [
        'B_BKGCAT',
]

mc_selection = [
        'B_BKGCAT == 10',
]

mass = {
        'Kaon': 493.67,
        'Pion': 139.57,
        'Muon': 105.66,
        'Proton': 938.27,
}

from ruffus import *

@transform(input_data, suffix('.root'), '.reduced.root')
def reduce(infile, outfile):
    from root_numpy import root2array, array2root
    arr = root2array(infile, 'B2XMuMu_Line_TupleDST/DecayTree', variables)

    arr = append_fields(arr, 'B_TAU', calc_tau(arr))
    
    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

@transform(input_mc, suffix('.root'), '.reduced.root')
def reduce_mc(infile, outfile):
    from root_numpy import root2array, array2root
    arr = root2array(infile, 'B2XMuMu_Line_TupleMC/DecayTree', variables + mc_variables, selection=prepare_sel(mc_selection))
    arr = append_fields(arr, 'B_TAU', calc_tau(arr))
    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

@transform(reduce, suffix('.root'), '.blinded.root')
def blind_signalpeak(infile, outfile):
    from root_numpy import root2array, array2root
    B_mass = 5279
    width = 50

    selstr = '(B_M < {}) || (B_M > {})'.format(B_mass - width, B_mass + width)
    arr = root2array(infile, 'B2dD0MuMu', selection=selstr)
    logging.info('Events in blinded dataset: ' + str(len(arr)))
    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

@transform(reduce_mc, suffix('.root'), '.mc_cut.root')
def select_mc(infile, outfile):
    select(infile, outfile, plots=None)

@transform(blind_signalpeak, suffix('.root'), '.cut.root')
def select(infile, outfile, plots='plots/select.pdf'):
    from root_numpy import root2array, array2root

    selection = [
        'Psi_M < 2860 || Psi_M > 3200',
        'Psi_M < 3500',
        'Kplus_isMuonLoose == 0',
        'piminus_isMuonLoose == 0',
    ]

    trigger_lines = [
        'B_L0MuonDecision_TOS == 1',
        'B_Hlt1TrackAllL0Decision_TOS == 1',
        'B_Hlt1TrackMuonDecision_TOS == 1',
        'B_Hlt2Topo2BodyBBDTDecision_TOS == 1',
        'B_Hlt2Topo3BodyBBDTDecision_TOS == 1',
        'B_Hlt2Topo3BodyBBDTDecision_TOS == 1',
        'B_Hlt2Topo4BodyBBDTDecision_TOS == 1',
        'B_Hlt2TopoMu2BodyBBDTDecision_TOS == 1',
        'B_Hlt2TopoMu3BodyBBDTDecision_TOS == 1',
        'B_Hlt2TopoMu4BodyBBDTDecision_TOS == 1',
        'Psi_Hlt2SingleMuonDecision_TOS == 1',
        'B_Hlt2DiMuonDetachedDecision_TOS == 1',
    ]

    pid_cuts = [
        'Kplus_PIDK > 0',
        'piminus_PIDK < 0',
    ]

    trigger_cut = '(' + ' || '.join(trigger_lines) + ')'

    arr = root2array(infile, 'B2dD0MuMu', selection=prepare_sel(selection) + ' && ' + trigger_cut)

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import seaborn as sns
    sns.set_palette("deep", desat=.6)

    logging.info('Plotting Kplus_PIDK vs. piminus_PIDK')

    mask = (arr['Kplus_PIDK'] > -999) & (arr['piminus_PIDK'] > -999)
    plt.hist2d(arr['Kplus_PIDK'][mask], arr['piminus_PIDK'][mask], bins=50, norm=LogNorm())
    #jp = sns.jointplot(arr['Kplus_PIDK'], arr['piminus_PIDK'], kind='hex', joint_kws={'norm': LogNorm()})
    #jp.set_axis_labels('$K^{+}_\\mathrm{DLLk}}$', '$\\pi^{-}_\\mathrm{DLLk}}$')
    plt.xlabel('$K^{+}_{\\mathrm{DLL}K\\pi}}$')
    plt.ylabel('$\\pi^{-}_{\\mathrm{DLL}K\\pi}$')
    plt.axhline(0, color='r', ls='-')
    plt.axvline(0, color='r', ls='-')
    # idea: diagonal cut
    #plt.plot([-16, 80], [-10, 50], 'r--', alpha=0.6)
    plt.xlim(-140, 140)
    plt.tight_layout()
    plt.text(50, 120, 'True $K$ & False $\\pi$', color='r', fontsize=18)
    plt.text(50, -140, 'True $K$ & True $\\pi$', color='r', fontsize=18)
    plt.text(-130, 120, 'False $K$ & False $\\pi$', color='r', fontsize=18)
    plt.text(-130, -140, 'False $K$ & True $\\pi$', color='r', fontsize=18)
    plt.savefig('plots/pid_plot.pdf')
    plt.clf()

    arr = root2array(infile, 'B2dD0MuMu', selection=prepare_sel(selection + pid_cuts) + ' && ' + trigger_cut)
    logging.info('Events after selection: ' + str(len(arr)))
    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

    if plots:
        plot(outfile, plots)

#@transform(select, suffix('.root'), '.misid.root')
def add_misid(infile, outfile):
    from root_numpy import root2array, array2root
    import numpy as np
    from itertools import product

    arr = root2array(infile, 'B2dD0MuMu')
    K_px = arr['Kplus_PX']
    K_py = arr['Kplus_PY']
    K_pz = arr['Kplus_PZ']

    pi_px = arr['piminus_PX']
    pi_py = arr['piminus_PY']
    pi_pz = arr['piminus_PZ']

    Psi_px = arr['Psi_PX']
    Psi_py = arr['Psi_PY']
    Psi_pz = arr['Psi_PZ']
    Psi_pe = arr['Psi_PE']

    hypotheses = ['Pion', 'Kaon', 'Muon', 'Proton']

    newnames = []
    newfields = []

    for m_K, m_pi in product(hypotheses, hypotheses):
        logging.info('Hypothesis: {}-{}'.format(m_K, m_pi))
        K_pe  = np.sqrt(K_px**2  + K_py**2  + K_pz**2  + mass[m_K]**2)
        pi_pe = np.sqrt(pi_px**2 + pi_py**2 + pi_pz**2 + mass[m_pi]**2)

        D_px = K_px + pi_px
        D_py = K_py + pi_py
        D_pz = K_pz + pi_pz
        D_pe = K_pe + pi_pe

        B_px = D_px + Psi_px
        B_py = D_py + Psi_py
        B_pz = D_pz + Psi_pz
        B_pe = D_pe + Psi_pe

        D_m = np.sqrt(D_pe**2 - D_px**2 - D_py**2 - D_pz**2)
        B_m = np.sqrt(B_pe**2 - B_px**2 - B_py**2 - B_pz**2)

        newfields.append(D_m)
        newnames.append('D_M_HYPO_K={},pi={}'.format(m_K, m_pi))
        newfields.append(B_m)
        newnames.append('B_M_HYPO_K={},pi={}'.format(m_K, m_pi))

    # Hypothesis: Missing Kaon
    B_px = K_px + Psi_px
    B_py = K_py + Psi_py
    B_pz = K_pz + Psi_pz
    B_pe = K_pe + Psi_pe
    B_m = np.sqrt(B_pe**2 - B_px**2 - B_py**2 - B_pz**2)
    newfields.append(B_m)
    newnames.append('B_M_HYPO_K=Missing,pi=Pion')

    # Hypothesis: Missing Pion
    B_px = pi_px + Psi_px
    B_py = pi_py + Psi_py
    B_pz = pi_pz + Psi_pz
    B_pe = pi_pe + Psi_pe
    B_m = np.sqrt(B_pe**2 - B_px**2 - B_py**2 - B_pz**2)
    newfields.append(B_m)
    newnames.append('B_M_HYPO_K=Kaon,pi=Missing')

    arr = append_fields(arr, newnames, newfields)
    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

@transform(select, suffix('.root'), add_inputs(select_mc), '.classified.root')
def classify(inputs, output):
    from root_numpy import root2array, array2root
    import numpy as np
    fname = inputs[0]
    mcname = inputs[1]

    select_sidebands = [
            'B_M > 5300'
    ]

    step = 50

    mcname_new = mcname.replace('.root', '.classified.root')
    sidebands = root2array(fname, 'B2dD0MuMu', bdt_variables, selection=prepare_sel(select_sidebands), step=step)
    mc = root2array(mcname, 'B2dD0MuMu', bdt_variables, step=step)

    from sklearn.tree import DecisionTreeClassifier
    from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
    from sklearn.cross_validation import StratifiedKFold
    from sklearn.metrics import roc_curve, auc
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    X = np.append(np.array(sidebands), np.array(mc))
    X = np.array(X.tolist())
    y = np.append(np.zeros(len(sidebands)), np.ones(len(mc)))

    clf = AdaBoostClassifier(
            DecisionTreeClassifier(max_depth=3),
            algorithm='SAMME.R',
            n_estimators=300,
            learning_rate=0.5
    )

    mean_fpr = np.linspace(0, 1, 200)

    with PdfPages('plots/classifier.pdf') as pdf:

        logging.info('Running x-validation...')
        skf = StratifiedKFold(y, 10)

        mean_fpr = np.linspace(0, 1, 200)
        mean_tpr = np.zeros(200)

        for i, (train, test) in enumerate(skf):
            logging.info('Running fold #{}'.format(i + 1))
            probs = clf.fit(X[train], y[train]).predict_proba(X[test])
            fpr, tpr, thresholds = roc_curve(y[test], probs[:,1])
            plt.plot(1 - fpr, tpr)
        pdf.savefig()
        plt.clf()

        logging.info("Variable importances:")
        imp = sorted(zip(sidebands.dtype.names, clf.feature_importances_), key=lambda x: -x[1])
        for n, i in imp:
            logging.info("{} - {}".format(n, i))
        plt.bar(np.arange(len(imp)), [entr[1] for entr in imp], 0.35)
        plt.gca().set_xticklabels([entr[0] for entr in imp])
        pdf.savefig()
        plt.clf()

        logging.info('Fit classifier model...')
        clf.fit(X, y)
        import pickle
        logging.info('Dump classifier to disk...')
        s = pickle.dumps(clf)
        with open('classifier.pkl', 'wb') as f:
            f.write(s)

    data_vars = root2array(fname, 'B2dD0MuMu', bdt_variables)
    data_vars = np.array(data_vars.tolist())
    data = root2array(fname, 'B2dD0MuMu')

    logging.info('Apply classifier to data...')
    pred = clf.decision_function(data_vars)
    data = append_fields(data, 'classifier', pred)
    array2root(data, output, 'B2dD0MuMu', 'recreate')

    logging.info('Apply classifier to MC...')
    mc_vars = root2array(mcname, 'B2dD0MuMu', bdt_variables)
    mc_vars = np.array(mc_vars.tolist())
    mc = root2array(mcname, 'B2dD0MuMu')
    pred = clf.decision_function(mc_vars)
    mc = append_fields(mc, 'classifier', pred)
    array2root(mc, mcname_new, 'B2dD0MuMu', 'recreate')

@transform(classify, formatter(), add_inputs(select_mc), 'plots/final.pdf')
def plot_final(infile, outfile):
    cuts = [
        'B_M > 5200 && B_M < 5450',
        'D~0_M > 1800 && D~0_M < 1940',
        'classifier > 0.03',
    ]
    plot(infile[0], outfile, mcfile=infile[1], cuts=cuts)

def plot(data, plotfile, mcfile=None, cuts=None):
    import numpy as np
    from root_numpy import root2array
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
    sns.set_palette("deep", desat=.6)

    if cuts is None:
        cuts = []

    arr = root2array(data, 'B2dD0MuMu', selection=prepare_sel(cuts))

    if mcfile:
        mc = mcfile.replace('.root', '.classified.root')
        arr_mc = root2array(mc, 'B2dD0MuMu', selection=prepare_sel(cuts))
        total_mc = len(root2array(mc, 'B2dD0MuMu'))
        factor = float(50) / total_mc

    logging.info('Saving plots to {}'.format(plotfile))
    with PdfPages(plotfile) as pdf:
        for vname in arr.dtype.names:
            logging.debug('Plotting ' + vname)
            #x = arr[vname][arr[vname] > -1000]
            #x_mc = arr_mc[vname][arr_mc[vname] > -1000]
            x = arr[vname]
            n, bins, _ = plt.hist(x, histtype='stepfilled', bins=40, alpha=0.7, normed=True)

            if mcfile:
                x_mc = arr_mc[vname]
                if vname in arr_mc.dtype.names:
                    n_mc, edges = np.histogram(arr_mc[vname], bins)
                    binned_hist(plt.gca(), factor * n_mc, edges, histtype='stepfilled', alpha=0.7, normed=True)

                    #plt.hist(x_mc, histtype='stepfilled', bins=bins, alpha=0.8, normed=True)
            #plt.yscale('log')
            if 'B_M' in vname:
                plt.axvline(5279, color='b', ls='--')
            plt.xlabel(vname)
            plt.ylim(0, max(n) * 1.05)
            pdf.savefig()
            plt.clf()

        logging.info('Plotting m_B vs. q^2')
        jp = sns.jointplot(arr['B_M'], arr['Psi_M'], kind='hex', joint_kws={'norm': LogNorm()})
        jp.set_axis_labels('$m_{B^0}$', '$q^2_{\\mu\\mu}$')
        plt.tight_layout()
        pdf.savefig()
        plt.clf()

        logging.info('Plotting m_D vs. q^2')
        jp = sns.jointplot(arr['D~0_M'], arr['Psi_M'], kind='hex', joint_kws={'norm': LogNorm()})
        jp.set_axis_labels('$m_{\\bar{D}^0}$', '$q^2_{\\mu\\mu}$')
        plt.tight_layout()
        pdf.savefig()
        plt.clf()

def prepare_sel(selections):
    return ' && '.join(map(lambda x: '(' + x + ')', selections))

def calc_tau(arr):
    from scipy.constants import c
    arr = arr['B_FD_OWNPV'] * arr['B_M'] / (arr['B_P'] * c * 10**3) * 10**12
    return arr

def binned_hist(ax, data, binedges, *args, **kwargs):
    #The dataset values are the bin centres
    x = (binedges[1:] + binedges[:-1]) / 2.0
    #The weights are the y-values of the input binned data
    weights = data
    return ax.hist(x, bins=binedges, weights=weights, *args, **kwargs)

if __name__ == '__main__':
    import sys
    pipeline_printout_graph("flow.pdf", forcedtorun_tasks = [reduce, reduce_mc], no_key_legend = True)
    pipeline_run(forcedtorun_tasks=sys.argv[1:])
