
from ruffus import *
from pipeline import reduce, select, reduce_mc, select_mc, classify

@transform(reduce, formatter(), add_inputs(reduce_mc), 'plots/b_mass.pdf')
def plot_entire(infiles, outfile):

    infile = infiles[0]
    mc = infiles[1]

    import matplotlib.pyplot as plt

    arr = read_root(infile, columns=['B_M', 'D~0_M', 'Psi_M'])
    mc = read_root(mc, columns=['B_M', 'D~0_M', 'Psi_M'])
    n, edges, _ = plt.hist(arr['B_M'].ravel(), histtype='stepfilled', bins=200, alpha=0.7, color='blue')
    ax = plt.gca()
    width = np.diff(edges)[0]
    plt.text(0.75, 0.8, 'LHCb Unofficial', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=30)
    plt.xlabel('$m(K\\pi\\mu\\mu)\\,/\\,\\mathrm{MeV}$', ha='right', x=0.9)
    plt.xlim(arr['B_M'].min(), arr['B_M'].max())
    plt.ylabel('Kandidaten $(1 / {:.1f}\\,\\mathrm{{MeV}})$'.format(width), ha='right', y=0.9)
    #plt.tight_layout()
    plt.savefig('plots/b_mass.pdf')
    plt.clf()

    n, edges, _ = plt.hist(arr['Psi_M'].ravel(), histtype='stepfilled', bins=200, color='blue', lw=0, label='Daten')
    plt.hist(mc['Psi_M'].ravel(), histtype='stepfilled', bins=edges, color='green', lw=0, label='Simulation')
    plt.axvspan(2920, 3200, alpha=0.5, color='red')
    plt.axvspan(3500, max(arr['Psi_M']), alpha=0.5, color='red')
    plt.yscale('log')
    ax = plt.gca()
    plt.text(0.25, 0.7, 'LHCb Unofficial', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=30)
    plt.xlim(arr['Psi_M'].min(), arr['Psi_M'].max())
    plt.xlabel('$m(\\mu\\mu)\\,/\\,\\mathrm{MeV}$', ha='right', x=0.9)
    width = np.diff(edges)[0]
    plt.ylabel('Kandidaten $(1 / {:.1f}\\,\\mathrm{{MeV}})$'.format(width), ha='right', y=0.9)
    plt.legend(loc='best')
    #plt.tight_layout()

    plt.savefig('plots/mumu_mass.pdf')
    plt.clf()

    plt.hist(arr['D~0_M'].ravel(), histtype='stepfilled', bins=200, alpha=0.7, color='blue')
    ax = plt.gca()
    plt.text(0.75, 0.8, 'LHCb Unofficial', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=26)
    plt.xlim(arr['D~0_M'].min(), arr['D~0_M'].max())
    plt.xlabel('$m(\\overline{D}^0=K^+\\!\\pi^-\\!)$', ha='right', x=0.9)
    plt.ylabel('Candidates')
    #plt.tight_layout()
    plt.savefig('plots/d_mass.pdf')
    plt.clf()

@transform(select, formatter(), add_inputs(select_mc), 'plots/nice/*.pdf')
def plot_nicevars(infiles, outfiles):
    infile, mc = infiles
    data = read_root(infile)
    mc = read_root(mc)

    import matplotlib.pyplot as plt

    for var in data.columns:
        n, edges, _ = plt.hist(data[var].ravel(), histtype='stepfilled', bins=100, alpha=0.7, color='blue', normed=True, label='Sidebands')
        plt.hist(mc[var].ravel(), histtype='step', bins=edges, color='red', lw=3, label='Simulation', normed=True)
        ax = plt.gca()
        width = np.diff(edges)[0]
        plt.text(0.42, 0.9, 'LHCb Unofficial', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=30)
        plt.xlabel(var, ha='right', x=1)
        plt.xlim(data[var].min(), data[var].max())
        plt.ylabel('Candidates (normalized)', ha='right', y=1)
        plt.legend(loc='upper right')
        plt.tight_layout()
        plt.savefig('plots/nice/{}.pdf'.format(var))
        plt.clf()

@transform(classify, formatter(), 'plots/classifier_mass.pdf')
def plot_classifier_mass(infile, outfile):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import seaborn as sns
    from pandas import qcut
    sns.set_style('whitegrid')
    plt.figure(figsize=(16,10))

    df = read_root(infile, where='B_M < 6000', columns=['B_M', 'classifier'])
    sns.jointplot(df['B_M'], df['classifier'], kind='hex', joint_kws={'norm': LogNorm()}, stat_func=None)
    plt.savefig('plots/partial_mass_hist.pdf')
    plt.clf()

    # Partial dependence
    #df['bin'], bins = qcut(df.B_M, 10, labels=np.arange(10), retbins=True)
    #binned = df.groupby('bin')
    #binned.index = bins[:-1] + np.diff(bins)/2
    #mean = binned.classifier.median()
    #lower = binned.classifier.aggregate(lambda x: np.percentile(x, 5))
    #upper = binned.classifier.aggregate(lambda x: np.percentile(x, 95))
    #plt.errorbar(binned.index, mean, yerr=(lower, upper), fmt='o')
    #plt.savefig('plots/partial_dep_mass.pdf')
    #plt.clf()

    # Response in mass
    N = 10

    cuts = qcut(df.classifier[df.classifier > 0.0], N, retbins=True)[1][:-1]

    left = min(df.B_M)

    n, bins = np.histogram(df.B_M, bins=50)
    upper = 1000

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

