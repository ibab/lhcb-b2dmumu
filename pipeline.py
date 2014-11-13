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
input_mc = []

# Real data fetched from the grid
input_data = [
    './storage/2011/Stripping20r1/Dimuon/DVBd2MuMuD0_data/combined/Bd2MuMuD0.root',
]

# Variables that are used in the analysis
variables = [
        'B_M',
        'B_IPCHI2_OWNPV',
        'B_PT',
        'Psi_M',
        'Psi_IPCHI2_OWNPV',
        'Psi_PT',
        'D~0_M',
        'D~0_IPCHI2_OWNPV',
        'D~0_PT',
        'muplus_isMuon',
]

from ruffus import transform, split, suffix, regex, pipeline_run 

@transform(input_data, suffix('.root'), '.reduced.root')
def reduce(infile, outfile):
    from root_numpy import root2array, array2root

    arr = root2array(infile, 'B2XMuMu_Line_TupleDST/DecayTree', variables)
    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

@transform(reduce, suffix('.root'), '.blinded.root')
def blind_signalpeak(infile, outfile):
    from root_numpy import root2array, array2root
    B_mass = 5279
    width = 50

    selection = '(B_M < {}) | (B_M > {})'.format(B_mass - width, B_mass + width)
    arr = root2array(infile, 'B2dD0MuMu', selection=selection)
    print('Events in blinded dataset: ' + str(len(arr)))
    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

@transform(blind_signalpeak, suffix('.root'), '.cut.root')
def select(infile, outfile):
    from root_numpy import root2array, array2root

    selection = [
            'B_PT < 30000',
            #'Psi_M < 3000 | Psi_M > 3200',
            #'Psi_M > 3000 & Psi_M < 3200',
    ]

    selection = ' & '.join(map(lambda x: ' ( ' + x + ' ) ', selection))

    arr = root2array(infile, 'B2dD0MuMu', selection=selection)
    print('Events after selection: ' + str(len(arr)))
    array2root(arr, outfile, 'B2dD0MuMu', 'recreate')

@split(select, map(lambda x: 'plots/' + x + '.pdf', variables))
def plot_vars(infile, outfiles):
    from root_numpy import root2array
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_palette("deep", desat=.6)

    arr = root2array(infile, 'B2dD0MuMu')
    for vname in arr.dtype.names:
        print('Plotting ' + vname)
        x = arr[vname]
        plt.hist(x, histtype='stepfilled', bins=100)
        plt.xlabel(vname)
        plt.savefig('plots/' + vname + '.pdf')
        plt.clf()

    jp = sns.jointplot(arr['D~0_M'], arr['Psi_M'], kind='hex')
    jp.set_axis_labels('$m_{\\bar{D}^0}$', '$m_{\\mu\\mu}$')
    plt.tight_layout()
    plt.savefig('plots/masses2d.pdf')
    plt.clf()

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        pipeline_run(forcedtorun_tasks=sys.argv[1])
    else:
        pipeline_run()

