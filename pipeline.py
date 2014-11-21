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
        'Kplus_PIDk',
        'piminus_PIDk',
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
        # Ghost probability
        # no isMuonLoose on K,pi
]

selection = [
        'Psi_M < 2860 || Psi_M > 3200',
        'Psi_M < 3500',
        'Kplus_PIDk > 0',
        'piminus_PIDk < 0',
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
    try:
        arr = root2array(infile, 'B2XMuMu_Line_TupleDST/DecayTree', variables)
    except:
        arr = root2array(infile, 'DecayTree', variables)

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

@transform(blind_signalpeak, suffix('.root'), '.cut.root')
def select(infile, outfile):
    from root_numpy import root2array, array2root

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

    arr = root2array(infile, 'B2dD0MuMu', selection=prepare_sel(selection) + ' && (' + ' || '.join(trigger_lines) + ')')
    logging.info('Events after selection: ' + str(len(arr)))
    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

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

@transform(select, suffix('.root'), add_inputs(reduce_mc), '.classified.root')
def classify(inputs, output):
    from root_numpy import root2array, array2root
    import numpy as np
    fname = inputs[0]
    mcname = inputs[1]

    select_sidebands = [
            'B_M > 5300'
    ]

    mcname_new = mcname.replace('.root', '.classified.root')
    sidebands = root2array(fname, 'B2dD0MuMu', bdt_variables, selection=prepare_sel(select_sidebands), step=10)
    mc = root2array(mcname, 'B2dD0MuMu', bdt_variables, step=10)

    from sklearn.tree import DecisionTreeClassifier
    from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
    from sklearn.cross_validation import cross_val_score
    from sklearn.metrics import roc_curve, auc

    X = np.append(np.array(sidebands), np.array(mc))
    X = np.array(X.tolist())
    y = np.append(np.zeros(len(sidebands)), np.ones(len(mc)))

    clf = AdaBoostClassifier(
            DecisionTreeClassifier(max_depth=3),
            algorithm='SAMME',
            n_estimators=700,
            learning_rate=0.5
    )

    logging.info('Skip x-validation.')
    #logging.info('Running x-validation...')
    #scores = cross_val_score(clf, X, y, cv=10, n_jobs=4)
    #from scipy.stats import sem
    #logging.info('Scores: {} +/- {}'.format(np.mean(scores), sem(scores)))
    
    logging.info('Fit classifier model...')
    clf.fit(X, y)
    import pickle
    logging.info('Dump classifier to disk...')
    s = pickle.dumps(clf)
    with open('classifier.pkl', 'wb') as f:
        f.write(s)

    logging.info("Variable importances:")
    for n, i in sorted(zip(sidebands.dtype.names, clf.feature_importances_), key=lambda x: -x[1]):
        logging.info("{} - {}".format(n, i))

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

@transform(classify, formatter(), add_inputs(reduce_mc), 'plots.pdf')
def plot_vars(infile, outfile):
    import numpy as np
    from root_numpy import root2array
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
    sns.set_palette("deep", desat=.6)

    cuts = [
        #'D~0_M > 1800 && D~0_M < 1920',
        'B_M > 5200 && B_M < 5450',

        'D~0_M > 1800 && D~0_M < 1940',
        #'classifier > 0.13',
        #'classifier > 0.13',
        #'B_DIRA_OWNPV > 0.999996',

        #'B_ISOLATION_BDT_Soft < 0',
        #'B_ENDVERTEX_CHI2 < 10',
        ##'Kplus_PT + piminus_PT > 1000',

        #'D~0_DIRA_OWNPV > 0.9980',

        #'muplus_ProbNNmu > 0.6',
        #'muminus_ProbNNmu > 0.5',

        #'Kplus_ProbNNpi < 0.05',
        #'Kplus_ProbNNk > 0.8',
        #'Kplus_ProbNNp < 0.02',

        #'piminus_ProbNNpi > 0.9',
        #'piminus_ProbNNp < 0.002',
        #'piminus_ProbNNpi > 0.4',

        #'Psi_FD_ORIVX < 1.5',
    ]

    data = infile[0]
    mc = infile[1].replace('.root', '.classified.root')

    arr = root2array(data, 'B2dD0MuMu', selection=prepare_sel(cuts))
    total_mc = len(root2array(mc, 'B2dD0MuMu'))

    arr_mc = root2array(mc, 'B2dD0MuMu', selection=prepare_sel(cuts))

    #factor = float(30) / total_mc
    factor = float(50) / total_mc

    with PdfPages(outfile) as pdf:
        for vname in arr.dtype.names:
            logging.info('Plotting ' + vname)
            x = arr[vname][arr[vname] > -1000]
            x_mc = arr_mc[vname][arr_mc[vname] > -1000]
            n, bins, _ = plt.hist(x, histtype='stepfilled', bins=40, alpha=0.7, normed=True)
            if vname in arr_mc.dtype.names:
                n_mc, edges = np.histogram(arr_mc[vname], bins)
                #binned_hist(plt.gca(), factor * n_mc, edges, histtype='stepfilled', alpha=0.7, normed=True)

                plt.hist(x_mc, histtype='stepfilled', bins=bins, alpha=0.8, normed=True)
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

        # TODO plot K_ProbNN vs pi_ProbNN

def prepare_sel(selections):
    return ' && '.join(map(lambda x: '(' + x + ')', selections))

def calc_tau(arr):
    from scipy.constants import c
    arr = arr['B_FD_OWNPV'] * arr['B_M'] / (arr['B_P'] * c * 10**3) * 10**12
    print(arr)
    print(arr.dtype)
    return arr

def binned_hist(ax, data, binedges, *args, **kwargs):
    #The dataset values are the bin centres
    x = (binedges[1:] + binedges[:-1]) / 2.0
    #The weights are the y-values of the input binned data
    weights = data
    return ax.hist(x, bins=binedges, weights=weights, *args, **kwargs)

if __name__ == '__main__':
    import sys
    pipeline_run(forcedtorun_tasks=sys.argv[1:])
