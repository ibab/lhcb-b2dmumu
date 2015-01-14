#!/usr/bin/env python2

import logging
from log import setup_logging; setup_logging()
import plotting
from util import *

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
        'Psi_M',
        'Psi_FD_ORIVX',
        'Psi_FDCHI2_ORIVX',
        'D~0_M',
        'D~0_FD_ORIVX',
        'D~0_CosTheta',
        'D~0_DIRA_OWNPV',
        'Kplus_hasRich',
        'Kplus_ProbNNpi',
        'Kplus_ProbNNk',
        'Kplus_ProbNNmu',
        'Kplus_ProbNNp',
        'Kplus_PX',
        'Kplus_PY',
        'Kplus_PZ',
        'Kplus_PT',
        'Kplus_P',
        'Kplus_TRACK_GhostProb',
        'piminus_PX',
        'piminus_PY',
        'piminus_PZ',
        'piminus_PT',
        'piminus_P',
        'piminus_hasRich',
        'piminus_ProbNNpi',
        'piminus_ProbNNk',
        'piminus_ProbNNmu',
        'piminus_ProbNNp',
        'piminus_TRACK_GhostProb',
        'Psi_PX',
        'Psi_PY',
        'Psi_PZ',
        'Psi_PE',
        'muplus_TRACK_GhostProb',
        'muminus_TRACK_GhostProb',
        'muminus_ProbNNpi',
        'muminus_ProbNNk',
        'muminus_ProbNNmu',
        'muminus_ProbNNp',
        'muminus_P',
        'muminus_PT',
        'Kplus_TRACK_CHI2NDOF',
        'piminus_TRACK_CHI2NDOF',
        'muplus_TRACK_CHI2NDOF',
        'muminus_TRACK_CHI2NDOF',
        'muplus_ProbNNpi',
        'muplus_ProbNNk',
        'muplus_ProbNNmu',
        'muplus_ProbNNp',
        'muplus_P',
        'muplus_PT',
        'Kplus_isMuonLoose',
        'piminus_isMuonLoose',
        'muplus_isMuon',
        'muminus_isMuon',
        'nTracks',
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
        'Kplus_TRUEID',
        'piminus_TRUEID',
        'muplus_TRUEID',
        'muminus_TRUEID',
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

    plotting.plot(outfile, 'plots/blinded.pdf', bins=200, variables=['B_M'])

@transform(reduce_mc, suffix('.root'), '.pid_resampled.root')
def resample_pid(infile, outfile):
    import numpy as np
    import pid_resample
    resample = pid_resample.create_resampler()

    from root_numpy import root2array, array2root
    pids = root2array(infile, 'B2dD0MuMu')
    arr = root2array(infile, 'B2dD0MuMu')

    for part in ['Kplus', 'piminus', 'muplus', 'muminus']:
        for pid in ['K', 'mu']:
            new = []
            for id, p, pt, ntracks in zip(arr[part+'_TRUEID'], arr[part+'_P'], arr[part+'_PT'], arr['nTracks']):
                new.append(resample('DLL' + pid + 'Down', id, p, pt, ntracks))
            arr = append_fields(arr, part + '_ResampledProbNN' + pid, np.array(new))

    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

@transform(resample_pid, suffix('.root'), '.mc_cut.root')
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
        'B_Hlt2SingleMuonDecision_TOS == 1',
        'B_Hlt2DiMuonDetachedDecision_TOS == 1',
    ]

    pid_cuts = [
        'Kplus_PIDK > 0',
        'piminus_PIDK < 0',
    ]

    trigger_cut = '(' + ' || '.join(trigger_lines) + ')'

    arr = root2array(infile, 'B2dD0MuMu', selection=prepare_sel(selection) + ' && ' + trigger_cut)

    if plots:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        import seaborn as sns
        sns.set_palette("deep", desat=.6)
        sns.set_context('talk')

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
        plotting.plot(outfile, plots, bins=100)

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

    X = np.append(np.array(sidebands), np.array(mc))
    X = np.array(X.tolist())
    y = np.append(np.zeros(len(sidebands)), np.ones(len(mc)))

    clf = AdaBoostClassifier(
            DecisionTreeClassifier(max_depth=3),
            algorithm='SAMME.R',
            n_estimators=400,
            learning_rate=0.5
    )

    logging.info('Validating classifier...')
    from classification import validate_classifier
    validate_classifier(clf, X, y, 'plots')

    logging.info('Fit classifier to data...')
    clf.fit(X, y)

    logging.info('Dump classifier to disk...')
    import pickle
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
        'classifier > 0.015',
    ]
    plotting.plot(infile[0], outfile, mcfile=infile[1], cuts=cuts, variables=variables)

if __name__ == '__main__':
    import sys
    pipeline_printout_graph("flow.pdf", forcedtorun_tasks = [reduce, reduce_mc], no_key_legend = True)
    pipeline_run(forcedtorun_tasks=sys.argv[1:])

