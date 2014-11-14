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

        'Psi_M',

        'D~0_M',

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

        'piminus_hasRich',
        'piminus_ProbNNpi',
        'piminus_ProbNNk',
        'piminus_ProbNNmu',
]

bdt_variables = [
        'D~0_M',
        'piminus_ProbNNpi',
        'piminus_ProbNNk',
        'piminus_ProbNNmu',
        'Kplus_ProbNNpi',
        'Kplus_ProbNNk',
        'Kplus_ProbNNmu',
        #'muminus_ProbNNpi'
        'muminus_ProbNNmu',
        #'muplus_ProbNNpi',
        'muplus_ProbNNmu',
]

selection = [
        'Psi_M < 2850',
]

mc_variables = [
        'B_BKGCAT',
]

mc_selection = [
        'B_BKGCAT == 10',
]

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

#@subdivide([select, select_mc], formatter(), map(lambda x: '{path[0]}/' + x + '.pdf', variables), '{path[0]}/')
def plot_vars(infile, outdir):
    from root_numpy import root2array
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import seaborn as sns
    sns.set_palette("deep", desat=.6)

    arr = root2array(infile, 'B2dD0MuMu')
    for vname in arr.dtype.names:
        print('Plotting ' + vname)
        x = arr[vname]
        plt.hist(x, histtype='stepfilled', bins=100)
        plt.xlabel(vname)
        plt.savefig(outdir + '/' + vname + '.pdf')
        plt.clf()

    jp = sns.jointplot(arr['B_M'], arr['Psi_M'], kind='hex', joint_kws={'norm': LogNorm()})
    jp.set_axis_labels('$m_{B^0}$', '$q^2_{\\mu\\mu}$')
    plt.tight_layout()
    plt.savefig(outdir + '/BdmumuMasses2d.pdf')
    plt.clf()

    jp = sns.jointplot(arr['D~0_M'], arr['Psi_M'], kind='hex', joint_kws={'norm': LogNorm()})
    jp.set_axis_labels('$m_{\\bar{D}^0}$', '$q^2_{\\mu\\mu}$')
    plt.tight_layout()
    plt.savefig(outdir + '/D0mumuMasses2d.pdf')
    plt.clf()

@transform(select, suffix('.root'), add_inputs(select_mc), '.bdt.root')
def classify(inputs, output):
    from root_numpy import root2array, array2root
    import numpy as np
    fname = inputs[0]
    mcname = inputs[1]

    select_sidebands = [
            'B_M > 5300'
    ]

    selstr = ' & '.join(map(lambda x: ' ( ' + x + ' ) ', select_sidebands))
    sidebands = root2array(fname, 'B2dD0MuMu', bdt_variables, selection=selstr)
    mc = root2array(mcname, 'B2dD0MuMu', bdt_variables)

    from sklearn.ensemble import AdaBoostClassifier

    X = np.append(np.array(sidebands), np.array(mc))
    X = np.array(X.tolist())

    y = np.append(np.zeros(len(sidebands)), np.ones(len(mc)))

    clf = AdaBoostClassifier().fit(X, y)
    print(clf.

    data_vars = root2array(fname, 'B2dD0MuMu', bdt_variables)
    data_vars = np.array(data_vars.tolist())
    data = root2array(fname, 'B2dD0MuMu')

    pred = clf.predict_proba(data_vars).T[0]

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

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        pipeline_run(forcedtorun_tasks=sys.argv[1])
    else:
        pipeline_run()

    pipeline_printout_graph("flow.pdf", "pdf",
                            forcedtorun_tasks = [reduce],
                            no_key_legend = True)

