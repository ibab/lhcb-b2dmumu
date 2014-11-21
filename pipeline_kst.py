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
    './storage/Data/AllYears/Stripping20/Dimuon/DVBd2Kstmumu_data/combined/Bd2Kstmumu.root',
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
        'Kstar_M',
        'Kstar_FD_ORIVX',
        'Kstar_CosTheta',
        'Kstar_DIRA_OWNPV',
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
        'B_Hlt2SingleMuonDecision_TOS',
        'B_Hlt2DiMuonDetachedDecision_TOS',
]

# Variables that are used in the multivariate classification
bdt_variables = [
        'B_TAU',
        'B_ISOLATION_BDT_Soft',
        'B_ENDVERTEX_CHI2',
        'B_DIRA_OWNPV',
        'Kstar_DIRA_OWNPV',
        'muplus_ProbNNmu',
        'muminus_ProbNNmu',
        'Kplus_ProbNNpi',
        'Kplus_ProbNNk',
        'Kplus_ProbNNp',
        'piminus_ProbNNpi',
        'piminus_ProbNNp',
        'Psi_FD_ORIVX',
]

selection = [
        'Psi_M < 2860 || Psi_M > 3200',
        'Psi_M < 3500',
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
    
    array2root(arr, outfile, 'Bd2Kstarmumu', 'recreate')

@transform(reduce, suffix('.root'), '.cut.root')
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
        'B_Hlt2SingleMuonDecision_TOS == 1',
        'B_Hlt2DiMuonDetachedDecision_TOS == 1',
    ]

    arr = root2array(infile, 'Bd2Kstarmumu', selection=prepare_sel(selection) + ' && (' + ' || '.join(trigger_lines) + ')')
    logging.info('Events after selection: ' + str(len(arr)))
    array2root(arr, outfile, 'Bd2Kstarmumu', 'recreate')

@transform(select, suffix('.root'), '.classified.root')
def classify(inputs, output):
    from root_numpy import root2array, array2root
    import numpy as np
    fname = inputs

    import pickle
    with open('classifier.pkl', 'r') as f:
        clf = pickle.loads(f.read())

    data_vars = root2array(fname, 'Bd2Kstarmumu', bdt_variables)
    data_vars = np.array(data_vars.tolist())
    data = root2array(fname, 'Bd2Kstarmumu')

    pred = clf.decision_function(data_vars)
    data = append_fields(data, 'classifier', pred)
    array2root(data, output, 'Bd2Kstarmumu', 'recreate')

@transform(classify, formatter(), 'plots_kst.pdf')
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
        #'Kstar_M > 1800 && Kstar_M < 1920',
        'B_M > 5150 && B_M < 5700',

        'Kstar_M > 792 && Kstar_M < 992',
        #'Kstar_M >  && Kstar_M < 1890',
        'classifier > 0.03',
        #'classifier > -0.15',
        #'B_DIRA_OWNPV > 0.999996',

        #'B_ISOLATION_BDT_Soft < 0',
        #'B_ENDVERTEX_CHI2 < 10',
        ##'Kplus_PT + piminus_PT > 1000',

        #'Kstar_DIRA_OWNPV > 0.9980',

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

    data = infile

    arr = root2array(data, 'Bd2Kstarmumu', selection=prepare_sel(cuts))

    with PdfPages(outfile) as pdf:
        for vname in arr.dtype.names:
            logging.info('Plotting ' + vname)
            n, bins, _ = plt.hist(arr[vname], histtype='stepfilled', bins=50, alpha=0.7)

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
        jp = sns.jointplot(arr['Kstar_M'], arr['Psi_M'], kind='hex', joint_kws={'norm': LogNorm()})
        jp.set_axis_labels('$m_{K^*}$', '$q^2_{\\mu\\mu}$')
        plt.tight_layout()
        pdf.savefig()
        plt.clf()

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
