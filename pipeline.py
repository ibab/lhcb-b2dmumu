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

# Simulated data fetched from the grid
input_mc = [
    './storage/MC/2012/Stripping20/AllStreams/DVBd2MuMuD0_MC/combined/Bd2MuMuD0.root'
]

# Real data fetched from the grid
input_data = [
    './storage/Data/AllYears/Stripping20/Dimuon/DVBd2MuMuD0_data/combined/Bd2MuMuD0.root'
]

# Variables that are used in the analysis
variables = [
        'B_M',
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
        'piminus_PX',
        'piminus_PY',
        'piminus_PZ',
        'Psi_PX',
        'Psi_PY',
        'Psi_PZ',
        'Psi_PE',
        'piminus_hasRich',
        'piminus_ProbNNpi',
        'piminus_ProbNNk',
        'piminus_ProbNNmu',
        'piminus_ProbNNp',
]

bdt_variables = [
        'B_OWNPV_CHI2',
        'B_ENDVERTEX_CHI2',
        'B_P',
        'B_PT',
        'B_ISOLATION_BDT_Soft',
        'B_ISOLATION_BDT_Hard',
        'Psi_FD_ORIVX',
        'Psi_FDCHI2_ORIVX',
        'piminus_ProbNNpi',
        'piminus_ProbNNk',
        'piminus_ProbNNmu',
        'Kplus_ProbNNpi',
        'Kplus_ProbNNk',
        'Kplus_ProbNNmu',
        'muminus_ProbNNpi',
        'muminus_ProbNNmu',
        'muplus_ProbNNpi',
        'muplus_ProbNNmu',
]

selection = [
        #'Psi_M < 2850 || Psi_M > 3200',
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
    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

@transform(input_mc, suffix('.root'), '.reduced.root')
def reduce_mc(infile, outfile):
    from root_numpy import root2array, array2root
    arr = root2array(infile, 'B2XMuMu_Line_TupleMC/DecayTree', variables + mc_variables)
    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

@transform(reduce, suffix('.root'), '.blinded.root')
def blind_signalpeak(infile, outfile):
    from root_numpy import root2array, array2root
    B_mass = 5279
    width = 50

    selstr = '(B_M < {}) | (B_M > {})'.format(B_mass - width, B_mass + width)
    arr = root2array(infile, 'B2dD0MuMu', selection=selstr)
    print('Events in blinded dataset: ' + str(len(arr)))
    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

@transform(reduce_mc, suffix('.root'), '.mccut.root')
def select_mc(infile, outfile):
    from root_numpy import root2array, array2root

    selstr = ' & '.join(map(lambda x: ' ( ' + x + ' ) ', mc_selection))

    arr = root2array(infile, 'B2dD0MuMu', selection=selstr)
    print('Events after selection: ' + str(len(arr)))
    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

@transform(blind_signalpeak, suffix('.root'), '.cut.root')
def select(infile, outfile):
    from root_numpy import root2array, array2root

    selstr = ' & '.join(map(lambda x: ' ( ' + x + ' ) ', selection))

    arr = root2array(infile, 'B2dD0MuMu', selection=selstr)
    print('Events after selection: ' + str(len(arr)))
    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

@transform(select, suffix('.root'), '.misid.root')
def add_mis_id(infile, outfile):
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

    for m_K, m_pi in product(hypotheses, hypotheses):
        print 'Hypothesis: {}-{}'.format(m_K, m_pi)
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

        arr = append_field(arr, 'D_M_HYPO_K={},pi={}'.format(m_K, m_pi), D_m)
        arr = append_field(arr, 'B_M_HYPO_K={},pi={}'.format(m_K, m_pi), B_m)

    # Hypothesis: Missing Kaon
    B_px = K_px + Psi_px
    B_py = K_py + Psi_py
    B_pz = K_pz + Psi_pz
    B_pe = K_pe + Psi_pe
    B_m = np.sqrt(B_pe**2 - B_px**2 - B_py**2 - B_pz**2)
    arr = append_field(arr, 'B_M_HYPO_K=Missing,pi=Pion', B_m)

    # Hypothesis: Missing Pion
    B_px = pi_px + Psi_px
    B_py = pi_py + Psi_py
    B_pz = pi_pz + Psi_pz
    B_pe = pi_pe + Psi_pe
    B_m = np.sqrt(B_pe**2 - B_px**2 - B_py**2 - B_pz**2)
    arr = append_field(arr, 'B_M_HYPO_K=Kaon,pi=Missing', B_m)

    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

#@transform(select, suffix('.root'), add_inputs(select_mc), '.bdt.root')
def classify(inputs, output):
    from root_numpy import root2array, array2root
    import numpy as np
    fname = inputs[0]
    mcname = inputs[1]

    select_sidebands = [
            'B_M > 5300'
    ]

    selstr = ' & '.join(map(lambda x: ' ( ' + x + ' ) ', select_sidebands))
    sidebands = root2array(fname, 'B2dD0MuMu', bdt_variables, selection=selstr, step=100)
    mc = root2array(mcname, 'B2dD0MuMu', bdt_variables, step=100)

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
            n_estimators=500,
            learning_rate=0.5
    )

    print('Running x-validation...')
    scores = cross_val_score(clf, X, y, cv=10, n_jobs=4)
    from scipy.stats import sem
    print('Scores: {} +/- {}'.format(np.mean(scores), sem(scores)))
    
    clf.fit(X, y)

    for n, i in sorted(zip(sidebands.dtype.names, clf.feature_importances_), key=lambda x: -x[1]):
        print n, i

    data_vars = root2array(fname, 'B2dD0MuMu', bdt_variables)
    data_vars = np.array(data_vars.tolist())
    data = root2array(fname, 'B2dD0MuMu')

    pred = clf.decision_function(data_vars)

    data = append_field(data, 'classifier', pred)
    array2root(data, output, 'B2dD0MuMu', 'recreate')

def append_field(rec, name, arr, dtype=None):
    import numpy as np
    arr = np.asarray(arr)
    if dtype is None:
        dtype = arr.dtype
    newdtype = np.dtype(rec.dtype.descr + [(name, dtype)])
    newrec = np.empty(rec.shape, dtype=newdtype)
    for field in rec.dtype.fields:
        newrec[field] = rec[field]
    newrec[name] = arr
    return newrec

@transform(add_mis_id, suffix('.root'), add_inputs(select_mc), '.pdf')
def plot_vars(infile, outfile):
    from root_numpy import root2array
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
    sns.set_palette("deep", desat=.6)

    cuts = [
        #'classifier > 0.1',
        #'muplus_ProbNNmu > 0.80',
        #'muminus_ProbNNmu > 0.80',
        #'muplus_ProbNNpi < 0.40',
        #'muminus_ProbNNpi < 0.40',
        #'Psi_M < 3600 || Psi_M > 3800',
        #'Psi_FDCHI2_ORIVX < 1',
        #'Psi_FDCHI2_ORIVX > 0',
        #'Psi_FD_ORIVX / Psi_FDCHI2_ORIVX < 0.1',
        #'Kplus_ProbNNp < 0.10',
        #'D~0_FD_ORIVX > 40',
        #'D~0_FD_ORIVX < 60',
        #'B_FD_OWNPV > 10',
        #'B_FD_OWNPV < 60',
        #'Kplus_ProbNNp > 0.0',
        #'piminus_ProbNNp > 0.0',
        #'D~0_M > 1845 && D~0_M < 1885',
        #'Psi_M < 2950 || Psi_M > 3200',
        #'Psi_M < 3500',
        #'D~0_FD_ORIVX < 120',
        #'B_M > 5150 && B_M < 5400',
        #'Kplus_ProbNNpi < 0.2',
        #'Kplus_ProbNNk > 0.6',
    ]

    data = infile[0]
    mc = infile[1]

    mc_new = mc.replace('.root', '.mc_misid.root')
    add_mis_id(mc, mc_new)

    with PdfPages(outfile) as pdf:
        arr = root2array(data, 'B2dD0MuMu', selection=' && '.join(map(lambda x: '(' + x + ')', cuts)))
        arr_mc = root2array(mc_new, 'B2dD0MuMu', selection=' && '.join(map(lambda x: '(' + x + ')', cuts)))

        for vname in arr.dtype.names:
            print('Plotting ' + vname)
            plt.hist(arr[vname], histtype='stepfilled', bins=100, alpha=0.7, normed=True)
            plt.hist(arr_mc[vname], histtype='stepfilled', bins=100, alpha=0.7, normed=True)
            #plt.yscale('log')
            plt.xlabel(vname)
            pdf.savefig()
            plt.clf()

        jp = sns.jointplot(arr['B_M'], arr['Psi_M'], kind='hex', joint_kws={'norm': LogNorm()})
        jp.set_axis_labels('$m_{B^0}$', '$q^2_{\\mu\\mu}$')
        plt.tight_layout()
        pdf.savefig()
        plt.clf()

        jp = sns.jointplot(arr['D~0_M'], arr['Psi_M'], kind='hex', joint_kws={'norm': LogNorm()})
        jp.set_axis_labels('$m_{\\bar{D}^0}$', '$q^2_{\\mu\\mu}$')
        plt.tight_layout()
        pdf.savefig()
        plt.clf()


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        pipeline_run(forcedtorun_tasks=sys.argv[1:])
    else:
        pipeline_run()

    pipeline_printout_graph("flow.pdf", "pdf",
                            forcedtorun_tasks = [reduce],
                            no_key_legend = True)

