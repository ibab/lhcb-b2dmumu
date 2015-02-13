#!/usr/bin/env python2

import os
import logging
import numpy as np
from pandas import DataFrame, Series
from root_pandas import read_root
import matplotlib
matplotlib.use('Agg')
from matplotlib.mlab import rec_append_fields as append_fields

from analysis.log import setup_logging, setup_roofit
setup_logging();
setup_roofit()
from analysis import plotting
from analysis.util import *

DATASTORE='/fhgfs/users/ibabuschkin/DataStore'
BLINDED = True

# Simulated data fetched from the grid
input_mc = [
    os.path.join(DATASTORE, 'MC/2012/Stripping20/AllStreams/DVBd2MuMuD0_MC/combined/Bd2MuMuD0.root'),
]

# Real data fetched from the grid
input_data = [
    os.path.join(DATASTORE, 'Data/AllYears/Stripping20/Dimuon/DVBd2MuMuD0_data/combined/Bd2MuMuD0.root'),
]

# Variables that are used in the analysis
variables = [
        '{B,D~0,Psi}_M',
        '{B,D~0,Psi}_P',
        '{B,D~0,Psi}_PT',
        'B_{DIRA,FD}_OWNPV',
        'B_{OWNPV,ENDVERTEX}_CHI2',
        'B_ISOLATION_BDT_{Hard,Soft}',

        'B_L0MuonDecision_TOS',
        'B_Hlt1TrackAllL0Decision_TOS',
        'B_Hlt1TrackMuonDecision_TOS',
        'B_Hlt2Topo{2,3,4}BodyBBDTDecision_TOS',
        'B_Hlt2TopoMu{2,3,4}BodyBBDTDecision_TOS',
        'B_Hlt2SingleMuonDecision_TOS',
        'B_Hlt2DiMuonDetachedDecision_TOS',

        'Psi_FD_ORIVX',
        'Psi_FDCHI2_ORIVX',

        'D~0_FD_ORIVX',
        'D~0_CosTheta',
        'D~0_DIRA_OWNPV',

        '{Kplus,piminus}_ProbNN*',
        '{Kplus,piminus}_PID*',
        '{Kplus,piminus}_hasRich',
        '{Kplus,piminus}_TRACK_GhostProb',
        '{Kplus,piminus}_TRACK_CHI2NDOF',
        '{Kplus,piminus}_isMuonLoose',

        'mu{plus,minus}_ProbNN*',
        'mu{plus,minus}_TRACK_GhostProb',
        'mu{plus,minus}_TRACK_CHI2NDOF',
        'mu{plus,minus}_isMuon',

        'nTracks',
]

# Variables used in the multivariate classification
bdt_variables = [
        'B_TAU',
        'B_ISOLATION_BDT_Soft',
        'B_ENDVERTEX_CHI2',
        'B_DIRA_OWNPV',
        'D~0_DIRA_OWNPV',
        'Psi_FD_ORIVX',
        'Kplus_TRACK_GhostProb',
        'piminus_TRACK_GhostProb',
        'muplus_TRACK_GhostProb',
        'muminus_TRACK_GhostProb',
        #'Kplus_ProbNNk',
        #'piminus_ProbNNk',
        #'muminus_ProbNNk',
        #'muplus_ProbNNk',
        #'Kplus_ProbNNmu',
        #'piminus_ProbNNmu',
        #'muminus_ProbNNmu',
        #'muplus_ProbNNmu',
]

mc_variables = [
        'B_BKGCAT',
        '*_TRUEID',
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

    if BLINDED:
        B_mass = 5279
        width = 50
        blindcut = '(B_M < {}) || (B_M > {})'.format(B_mass - width, B_mass + width)
    else:
        blindcut = None

    df = read_root(infile, 'B2XMuMu_Line_TupleDST/DecayTree', columns=variables, where=blindcut)

    # Add missing decay time
    df['B_TAU'] = Series(calc_tau(df), index=df.index)

    df.to_root(outfile, mode='recreate')

@transform(input_mc, suffix('.root'), '.reduced_mc.root')
def reduce_mc(infile, outfile):
    df = read_root(infile, 'B2XMuMu_Line_TupleMC/DecayTree', columns=variables + mc_variables, where=prepare_sel(mc_selection))

    # Add missing decay time
    df['B_TAU'] = Series(calc_tau(df), index=df.index)
    df.to_root(outfile, mode='recreate')

@transform(reduce, suffix('.root'), '.imputed.root')
def impute_data(infile, outfile):
    df = read_root(infile)
    for var in df.columns:
        if var == 'B_ISOLATION_BDT_Hard' or var == 'B_ISOLATION_BDT_Soft':
            df[df[var] == -2] = np.nan
        if 'PID' in var:
            df[df[var] == -1000] = np.nan
    df.to_root(outfile, mode='recreate')

@transform(reduce, suffix('.root'), '.imputed_mc.root')
def impute_mc(infile, outfile):
    impute_data(infile, outfile)

@transform(reduce_mc, suffix('.root'), '.mc_cut.root')
def select_mc(infile, outfile):
    select(infile, outfile, plots=None)

@transform(reduce, suffix('.root'), '.cut.root')
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
        'Kplus_ProbNNk > 0.5',
        'piminus_ProbNNk < 0.5',
    ]

    trigger_cut = '(' + ' || '.join(trigger_lines) + ')'

    arr = root2array(infile, selection=prepare_sel(selection) + ' && ' + trigger_cut)

    #if plots:
    #    import matplotlib.pyplot as plt
    #    from matplotlib.colors import LogNorm
    #    import seaborn as sns
    #    sns.set_palette("deep", desat=.6)
    #    sns.set_context('talk')

    #    logging.info('Plotting Kplus_PIDK vs. piminus_PIDK')

    #    mask = (arr['Kplus_PIDK'] > -999) & (arr['piminus_PIDK'] > -999)
    #    plt.hist2d(arr['Kplus_PIDK'][mask], arr['piminus_PIDK'][mask], bins=50, norm=LogNorm())
    #    #jp = sns.jointplot(arr['Kplus_PIDK'], arr['piminus_PIDK'], kind='hex', joint_kws={'norm': LogNorm()})
    #    #jp.set_axis_labels('$K^{+}_\\mathrm{DLLk}}$', '$\\pi^{-}_\\mathrm{DLLk}}$')
    #    plt.xlabel('$K^{+}_{\\mathrm{DLL}K\\pi}}$')
    #    plt.ylabel('$\\pi^{-}_{\\mathrm{DLL}K\\pi}$')
    #    plt.axhline(0, color='r', ls='-')
    #    plt.axvline(0, color='r', ls='-')
    #    # idea: diagonal cut
    #    #plt.plot([-16, 80], [-10, 50], 'r--', alpha=0.6)
    #    plt.xlim(-140, 140)
    #    plt.tight_layout()
    #    plt.text(50, 120, 'True $K$ & False $\\pi$', color='r', fontsize=18)
    #    plt.text(50, -140, 'True $K$ & True $\\pi$', color='r', fontsize=18)
    #    plt.text(-130, 120, 'False $K$ & False $\\pi$', color='r', fontsize=18)
    #    plt.text(-130, -140, 'False $K$ & True $\\pi$', color='r', fontsize=18)
    #    plt.savefig('plots/pid_plot.pdf')
    #    plt.clf()

    arr = root2array(infile, selection=prepare_sel(selection + pid_cuts) + ' && ' + trigger_cut)
    logging.info('Events after selection: ' + str(len(arr)))
    array2root(arr, outfile, mode='recreate')

    if plots:
        plotting.plot(outfile, plots, bins=100)

#@transform(select, suffix('.root'), '.misid.root')
def add_misid(infile, outfile):
    from root_numpy import root2array, array2root
    import numpy as np
    from itertools import product

    arr = root2array(infile)
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
    array2root(arr, outfile, mode='recreate')

@transform(select, suffix('.root'), add_inputs(select_mc), '.classified.root')
def classify(inputs, output):
    from root_numpy import root2array, array2root
    import numpy as np
    fname = inputs[0]
    mcname = inputs[1]

    select_sidebands = [
            'B_M > 5300'
    ]

    step = 100

    mcname_new = mcname.replace('.root', '.classified.root')
    bkg = read_root(fname, columns=bdt_variables, where=prepare_sel(select_sidebands), step=step)
    sig = read_root(mcname, columns=bdt_variables, step=step)

    from sklearn.tree import DecisionTreeClassifier
    from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier

    X = np.vstack([bkg.values, sig.values])
    y = np.append(np.zeros(len(bkg)), np.ones(len(sig)))

    clf = AdaBoostClassifier(
            DecisionTreeClassifier(max_depth=3),
            algorithm='SAMME.R',
            n_estimators=1000,
            learning_rate=0.5
    )

    logging.info('Validating classifier...')
    from analysis.classification import validate_classifier
    validate_classifier(clf, X, y, 'plots')

    logging.info('Fit classifier to data...')
    clf.fit(X, y)

    logging.info('Dump classifier to disk...')
    import pickle
    s = pickle.dumps(clf)
    with open('.classifier.pkl', 'wb') as f:
        f.write(s)

    for inp, outp in zip(inputs, [output, inputs[1].replace('.root', '.classified.root')]):
        df_bdt = read_root(inp, columns=bdt_variables)
        df_data = read_root(inp)

        logging.info('Apply classifier...')
        df_data['classifier'] = clf.decision_function(df_bdt.values)
        df_data.to_root(outp, mode='recreate')

@transform(classify, suffix('.root'), '.bdt_cut.root')
def bdt_cut(infile, outfile):
    read_root(infile, where='classifier > -0.1').to_root(outfile, mode='recreate')

@originate('toy.root')
def generate_toy(output):
    from analysis.fit_model import model, variables
    import ROOT
    ff = ROOT.TFile(output, 'recreate')
    ROOT.RooAbsData.setDefaultStorageType(ROOT.RooAbsData.Tree)
    dset = model.generate(variables, 1e5)
    dset.tree().Write('default')
    ff.Write()

@transform(generate_toy, formatter(), 'plots/fit.pdf')
def fit(infile, outfile):
    import ROOT
    from ROOT import RooDataSet, RooFit, RooArgSet, RooArgList, RooFormulaVar
    from analysis.fit_model import model, variables

    tf = ROOT.TFile(infile)
    tree = tf.Get('default')
    dset = RooDataSet('data', 'data', tree, variables)
    #results = model.fitTo(dset, RooFit.Range('leftband,rightband'), RooFit.Save())
    results = model.fitTo(dset, RooFit.Save())
    results.Print('v')
    model.getParameters(dset).writeToFile(outfile)

    final = results.floatParsFinal()
    lamb = final.iterator().Next().getVal()

    import matplotlib.pyplot as plt
    import numpy as np
    from missing_hep import histpoints
    from scipy.stats import expon
    import seaborn as sns

    sns.set_style('darkgrid')

    df = read_root(infile, columns=['B_M'])
    B_M = df['B_M']
    x, y, norm = histpoints(df.query('B_M > 5200 & B_M < 5450')['B_M'], errorbar={'markersize': 0}, xerr='binwidth')
    xl = np.linspace(5200, 5450, 200)
    #plt.plot(xl, norm * (lamb / (np.exp(lamb * 29) - np.exp(lamb * 0) + np.exp(lamb * 250) - np.exp(lamb * 129))) * np.exp(lamb * (xl - 5200)))
    #plt.plot(xl, norm * (lamb / (np.exp(lamb * 29) - np.exp(lamb * 0) + np.exp(lamb * 250) - np.exp(lamb * 129))) * np.exp(lamb * (xl - 5200)))
    from analysis.plotting import plot_roofit
    xlabel = '$m(B=K\\pi\\mu\\mu)\\ /\\ \\mathrm{MeV}$'
    ax, width = plot_roofit(variables['B_M'], dset, model, components=['sig', 'bkg'], xlabel=xlabel)
    plt.ylabel('Candidates$\\ (1\\ /\\ {:.2f}\\ \\mathrm{{MeV}})$'.format(width[0]), fontsize=14)
    #plt.tight_layout()
    plt.savefig(outfile)

@transform(classify, formatter(), add_inputs(select_mc), 'plots/final.pdf')
def plot_final(infile, outfile):
    cuts = [
        'B_M > 5200 && B_M < 5450',
        'D~0_M > 1800 && D~0_M < 1940',
        'classifier > -0.1',
    ]
    plotting.plot(infile[0], outfile, mcfile=infile[1], cuts=cuts, variables=variables)

@transform(classify, formatter(), 'plots/classifier_mass.pdf')
def plot_bdt_mass(infile, outfile):
    import matplotlib.pyplot as plt
    import seaborn as sns
    from pandas import qcut
    sns.set_style('whitegrid')

    df = read_root(infile, columns=['B_M', 'classifier'])

    N = 10

    cuts = qcut(df.classifier, N, retbins=True)[1][:-1]

    left = min(df.B_M)

    n, bins = np.histogram(df.B_M, bins=50)
    upper = max(n)

    for (c, color) in zip(cuts, sns.color_palette('Blues', N + 2)[2:]):
        data = df[df.classifier > c]['B_M']
        plt.hist(data.ravel(), bins=bins, histtype='stepfilled', color=color, lw=0, label='bdt > {:1.3f}'.format(c))
    plt.xlim(left=left)
    plt.ylim(0, upper)
    plt.ylabel('Candidates')
    plt.xlabel('$m_{B=K\\pi\\mu\\mu}$', fontsize=16)
    plt.legend()
    ax = plt.gca()
    ax.grid(False)
    plt.savefig(outfile)
    plt.clf()

@transform(classify, formatter(), 'plots/correlation.pdf')
def plot_correlation(infile, outfile):
    import matplotlib.pyplot as plt
    import seaborn as sns
    plt.figure(figsize=(10,10))
    df = read_root(infile, columns=bdt_variables + ['B_M'])
    sns.corrplot(df.corr(), diag_names=False)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.clf()

if __name__ == '__main__':
    import sys
    #pipeline_printout_graph("flow.pdf", forcedtorun_tasks = [reduce, reduce_mc], no_key_legend = True)
    pipeline_run(forcedtorun_tasks=sys.argv[1:], logger=logging)

