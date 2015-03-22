# encoding: utf-8

"""
B -> D mu mu main analysis pipeline
"""


from ruffus import *

import os
import sys
import numpy as np
from pandas import DataFrame, Series
from root_pandas import read_root
import matplotlib
matplotlib.use('Agg')
from analysis.log import setup_logging, get_logger, setup_roofit
logger = get_logger()
from analysis import plotting
from analysis.util import *
from analysis.push import notify
from plot_tasks import *


try:
    DATASTORE = os.environ['DATASTORE']
except KeyError:
    DATASTORE='/fhgfs/users/ibabuschkin/DataStore'

BLINDED = True

# Simulated data fetched from the grid (used for training selection and calculating efficiencies)
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
        'B_{OWNPV,ENDVERTEX}_NDOF',
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

        'mu{plus,minus}_PID*',
        'mu{plus,minus}_ProbNN*',
        'mu{plus,minus}_TRACK_GhostProb',
        'mu{plus,minus}_TRACK_CHI2NDOF',
        'mu{plus,minus}_isMuon',

        'nTracks',
]

# Variables used in the multivariate classification
bdt_variables = [
        'B_DiraAngle',
        'B_TAU',
        'B_P',
        'B_PT',
        'B_ISOLATION_BDT_Soft',
        'Kplus_PIDK',
        'piminus_PIDK',
        'muplus_PIDmu',
        'muminus_PIDmu',
        'B_ENDVERTEX_CHI2_NDOF',
]

mc_variables = [
        'B_BKGCAT',
        '*_TRUEID',
]


@transform(input_data, suffix('.root'), '.reduced.root')
def reduce(infile, outfile):

    if BLINDED:
        B_mass = 5279
        width = 50
        blindcut = '(B_M < {}) || (B_M > {})'.format(B_mass - width, B_mass + width)
    else:
        blindcut = None

    df = read_root(infile, 'B2XMuMu_Line_TupleDST/DecayTree', columns=variables, where=blindcut)

    df['B_TAU'] = Series(calc_tau(df), index=df.index)
    logger.info('Initial events: {}'.format(len(df)))
    df['B_DiraAngle'] = np.arccos(df['B_DIRA_OWNPV'])
    df['B_ENDVERTEX_CHI2_NDOF'] = df['B_ENDVERTEX_CHI2'] / df['B_ENDVERTEX_NDOF']
    df.to_root(outfile, mode='recreate')

@transform(input_mc, suffix('.root'), '.reduced_mc.root')
def reduce_mc(infile, outfile):
    df = read_root(infile, 'B2XMuMu_Line_TupleMC/DecayTree', columns=variables + mc_variables, where='B_BKGCAT == 10')

    df['B_TAU'] = Series(calc_tau(df), index=df.index)
    df['B_DiraAngle'] = np.arccos(df['B_DIRA_OWNPV'])
    df['B_ENDVERTEX_CHI2_NDOF'] = df['B_ENDVERTEX_CHI2'] / df['B_ENDVERTEX_NDOF']
    df.to_root(outfile, mode='recreate')

@transform(reduce_mc, suffix('.root'), '.mc_cut.root')
def select_mc(infile, outfile):
    select(infile, outfile, plots=None)

@transform(reduce, suffix('.root'), '.cut.root')
def select(infile, outfile, plots='plots/select.pdf'):
    from root_numpy import root2array, array2root

    selection = [
        # Exclude J/psi
        'Psi_M < 2860 || Psi_M > 3200',
        # Kinematic range ends below this
        'Psi_M < 3500',
    ]

    trigger_selection = [
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
        'Kplus_PIDK > -5',
        'piminus_PIDK < 5',
        'Kplus_isMuonLoose == 0',
        'piminus_isMuonLoose == 0',
    ]

    trigger_cut = '(' + ' || '.join(trigger_selection) + ')'

    arr = read_root(infile, where=prepare_sel(selection) + ' && ' + trigger_cut)

    if plots:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        #import seaborn as sns
        #sns.set_palette("deep", desat=.6)
        #sns.set_context('talk')
        #sns.set_style('white')

        logger.info('Plotting Kplus_PIDK vs. piminus_PIDK')

        plt.figure(figsize=(12,8))
        mask = (arr['Kplus_PIDK'] > -999) & (arr['piminus_PIDK'] > -999)
        plt.hist2d(arr['Kplus_PIDK'][mask], arr['piminus_PIDK'][mask], bins=50, norm=LogNorm(), cmap='afmhot_r')
        plt.colorbar()
        plt.xlim(-140, 140)
        #jp = sns.jointplot(arr['Kplus_PIDK'], arr['piminus_PIDK'], kind='hex', joint_kws={'norm': LogNorm()})
        #jp.set_axis_labels('$K^{+}_\\mathrm{DLLk}}$', '$\\pi^{-}_\\mathrm{DLLk}}$')
        plt.xlabel(u'$\mathrm{Likelihood}(K/\\pi)$ für Kaon', ha='right', x=1)
        plt.ylabel(u'$\mathrm{Likelihood}(K/\\pi)$ für Pion', ha='right', y=1)
        #plt.axhline(0, color='r', ls='-')
        #plt.axvline(0, color='r', ls='-')
        # idea: diagonal cut
        #plt.plot([-16, 80], [-10, 50], 'r--', alpha=0.6)
        plt.tight_layout()
        plt.text(30, 120, 'Echtes $K$ & Falsches $\\pi$', color='r', fontsize=24)
        plt.text(30, -140, 'Echtes $K$ & Echtes $\\pi$', color='r', fontsize=24)
        plt.text(-130, 120, 'Falsches $K$ & Falsches $\\pi$', color='r', fontsize=24)
        plt.text(-130, -140, 'Falsches $K$ & Echtes $\\pi$', color='r', fontsize=24)
        plt.tight_layout()
        plt.savefig('plots/pid_plot.pdf')
        plt.clf()

    before = len(arr)
    after_pid = arr.query(prepare_sel(pid_cuts))
    after = len(after_pid)

    logger.info('PID cut efficiency is {:.2f}%'.format(100.0 * after / before))

    after_pid.to_root(outfile, mode='recreate')

    if plots:
        plotting.plot(outfile, plots, bins=100)

@transform(select, suffix('.root'), add_inputs(select_mc), '.classified.root')
def classify(inputs, output):
    from root_numpy import root2array, array2root
    import numpy as np
    fname = inputs[0]
    mcname = inputs[1]

    select_sidebands = [
            '(B_M > 5300)',
    ]

    step = 1

    mcname_new = mcname.replace('.root', '.classified.root')
    bkg = read_root(fname, columns=bdt_variables, where=prepare_sel(select_sidebands), step=step)
    sig = read_root(mcname, columns=bdt_variables, step=step)

    def set_errors_nan(df):
        for var in df.columns:
            if 'ProbNN' in var:
                df[var][df[var] < 0] = np.nan
            if 'PID' in var:
                df[var][df[var] < -999] = np.nan
            # Surprisingly, the failure case indicates signal with 90% probability
            # This might just be a simulation-specific thing, though
            #if 'ISOLATION' in var:
            #    df[var][df[var] == -2] = np.nan

    set_errors_nan(sig)
    set_errors_nan(bkg)

    from sklearn.tree import DecisionTreeClassifier
    from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import Imputer

    X = np.vstack([bkg.values, sig.values])
    y = np.append(np.zeros(len(bkg)), np.ones(len(sig)))

    imp = Imputer(strategy='mean')
    X = imp.fit_transform(X, y)

    classifiers = {'sklearn_adaboost': AdaBoostClassifier(
                                           DecisionTreeClassifier(max_depth=6),
                                           algorithm='SAMME',
                                           n_estimators=500,
                                           learning_rate=1,
                                       ),
                  }


    for name, clf in classifiers.items():
        logger.info('Validating classifier {}...'.format(name))
        from analysis.classification import validate_classifier
        validate_classifier(clf, X, y, 'plots', name)
        notify('classify', 'x-validation of {} finished.'.format(name))

@originate('toy.root')
def generate_toy(output):
    from analysis.fit_model import model, variables
    import ROOT
    ROOT.RooRandom.randomGenerator().SetSeed(0)
    ff = ROOT.TFile(output, 'recreate')
    ROOT.RooAbsData.setDefaultStorageType(ROOT.RooAbsData.Tree)
    dset = model.generate(variables, 600)
    dset.tree().Write('default')
    ff.Write()

@transform(generate_toy, formatter(), 'plots/fit.pdf')
def fit(infile, outfile):
    import ROOT
    from ROOT import RooDataSet, RooFit, RooArgSet, RooArgList, RooFormulaVar, RooBinning
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
    from scipy.stats import expon
    #import seaborn as sns

    #sns.set_style('darkgrid')

    labels={
        'sig': 'Signal',
        'bkg': 'Background (comb.)',
        'prompt': r'$B^0\to K^+\!\pi^-\!\mu^+\!\mu^-$',
    }

    df = read_root(infile, columns=['B_M'])
    B_M = df['B_M']
    xl = np.linspace(5200, 5450, 200)
    #plt.plot(xl, norm * (lamb / (np.exp(lamb * 29) - np.exp(lamb * 0) + np.exp(lamb * 250) - np.exp(lamb * 129))) * np.exp(lamb * (xl - 5200)))
    #plt.plot(xl, norm * (lamb / (np.exp(lamb * 29) - np.exp(lamb * 0) + np.exp(lamb * 250) - np.exp(lamb * 129))) * np.exp(lamb * (xl - 5200)))
    from analysis.plotting import plot_roofit
    xlabel = '$m(B=K\\pi\\mu\\mu)\\ /\\ \\mathrm{MeV}$'
    binning = RooBinning(25, 5200, 5450)
    ax, width = plot_roofit(variables['B_M'], dset, model, components=['sig', 'bkg', 'prompt'], xlabel=xlabel, labels=labels, binning=binning)
    ax.text(0.75, 0.5, 'LHCb Unofficial', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=26)
    ax.text(0.75, 0.4, 'Toy Simulation', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=26)
    plt.legend()
    plt.ylabel('Candidates$\\ (1\\ /\\ {:.2f}\\ \\mathrm{{MeV}})$'.format(width[0]), ha='right', y=0.9)
    #plt.tight_layout()
    plt.savefig('plots/fit_bmass.pdf')
    plt.clf()
    from analysis.plotting import plot_roofit
    xlabel = '$m(\\overline{D}=K\\pi)\\ /\\ \\mathrm{MeV}$'
    binning = RooBinning(25, 1750, 2000)
    ax, width = plot_roofit(variables['D~0_M'], dset, model, components=['sig', 'bkg', 'prompt'], xlabel=xlabel, labels=labels, binning=binning)
    ax.text(0.75, 0.5, 'LHCb Unofficial', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=26)
    ax.text(0.75, 0.4, 'Toy Simulation', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=26)
    plt.legend()
    plt.ylabel('Candidates$\\ (1\\ /\\ {:.2f}\\ \\mathrm{{MeV}})$'.format(width[0]), ha='right', y=0.9)
    plt.savefig('plots/fit_dmass.pdf')


if __name__ == '__main__':
    setup_logging()
    setup_roofit()
    pipeline_run(forcedtorun_tasks=sys.argv[1:], logger=logger)

