from util import *
import logging
from root_pandas import read_root
from missing_hep import histpoints

def plot(data, plotfile, mcfile=None, cuts=None, variables=None, bins=30):
    import numpy as np
    from root_numpy import root2array
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
    sns.set_palette("deep", desat=.6)
    sns.set_context('talk')

    if cuts is None:
        cuts = []

    if variables is None:
        variables = []

    arr = read_root(data, where=prepare_sel(cuts))

    if mcfile:
        mc = mcfile.replace('.root', '.classified.root')
        arr_mc = read_root(mc, where=prepare_sel(cuts))
        factor = 0.01

    logging.info('Saving plots to {}'.format(plotfile))
    with PdfPages(plotfile) as pdf:
        for col in arr.columns:
            logging.debug('Plotting ' + col)
            x = arr[col]
            n, bine, _ = plt.hist(x.values, histtype='stepfilled', bins=bins, alpha=0.7, color='grey')

            if mcfile:
                x_mc = arr_mc[col]
                if col in arr_mc.columns:
                    n_mc, edges = np.histogram(arr_mc[col], bine)
                    binned_hist(plt.gca(), factor * n_mc, edges, histtype='stepfilled', alpha=0.7)

                    #plt.hist(x_mc, histtype='stepfilled', bins=bins, alpha=0.8, normed=True)
            #plt.yscale('log')
            plt.xlabel(col)
            plt.ylim(0, max(n) * 1.05)
            pdf.savefig()
            plt.clf()

        logging.info('Plotting m_B vs. q^2')
        jp = sns.jointplot(arr['B_M'], arr['Psi_M'], kind='hex', joint_kws={'norm': LogNorm()})
        jp.set_axis_labels('$m(B^0)$', '$m(\\mu^+\\mu^-)$')
        plt.tight_layout()
        pdf.savefig()
        plt.clf()

        logging.info('Plotting m_D vs. q^2')
        jp = sns.jointplot(arr['D~0_M'], arr['Psi_M'], kind='hex', joint_kws={'norm': LogNorm()})
        jp.set_axis_labels('$m_(\\bar{D}^0)$', '$m(\\mu^+\\mu^-)$')
        plt.tight_layout()
        pdf.savefig()
        plt.clf()

        
def plot_roofit(var, data, model, components=None, numcpus=1, xlabel='', extra_params=None, norm=None, log=False, binning=None, labels=None):
    if not extra_params:
        extra_params = []
    if not norm:
        norm = 1 #data.numEntries()

    from numpy import linspace, maximum, array

    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.ticker import MaxNLocator

    cax = plt.gca()
    box = cax.get_position()
    xmin, ymin = box.xmin, box.ymin
    xmax, ymax = box.xmax, box.ymax

    gs = GridSpec(2, 1, height_ratios=[7, 1], left=xmin, right=xmax, bottom=ymin, top=ymax)
    gs.update(hspace=0.12)

    xpos, width, y, yerr = get_binned_data(var, data, extra_params=extra_params, binning=binning)
    x = linspace(var.getMin(), var.getMax(), 200)

    if not components:
        f = get_function(var, model, data, norm=norm, extra_params=extra_params)
    else:
        f, comps = get_function(var, model, data, components=components, norm=norm, extra_params=extra_params)

    ax = plt.subplot(gs[0])

    if components:
        for c, name in zip(comps, components):
            if labels:
                plt.plot(x, c(x), '--', dashes=[8,3], zorder=90, clip_on=False, label=labels[name])
            else:
                plt.plot(x, c(x), '--', dashes=[8,3], zorder=90, clip_on=False)

    plt.plot(x, f(x), 'r-', zorder=95, clip_on=False)
    if log:
        ax.set_yscale('log')
        ax.set_ylim(max(0.1, min(y)), 1.2 * max(y))
        yerr[0] = y - maximum(1e-1, y - yerr[0])
    else:
        pass
        #ax.set_ylim(0, 1.1 * max(y))
        #plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))

    plt.errorbar(xpos, y, yerr=yerr, fmt='o', color='k', zorder=100, markersize=4, linewidth=1.5, clip_on=False, capsize=0)
    plt.setp(ax.get_xticklabels(), visible=False)

    bx = plt.subplot(gs[1], sharex=ax)
    plt.setp(bx, zorder=-10)# Make sure the main histogram can draw over the pull plot
    pull = calc_pull(xpos, f, y, yerr)

    plt.fill_between([var.getMin(), var.getMax()], 2, 1, facecolor='#bbbbbb', linewidth=0, edgecolor='#bbbbbb')
    plt.fill_between([var.getMin(), var.getMax()], -2, -1, facecolor='#bbbbbb', linewidth=0, edgecolor='#bbbbbb')

    color = ['#cc4444' if abs(p) > 3 else '#555555' for p in pull]
    plt.bar(xpos-width/2, pull, width, color=color, linewidth=0.8)

    #plt.axhline(0, color='black')
    #plt.axhline(3, color='black')
    plt.ylabel('Normed\nResiduals')

    if xlabel:
        plt.xlabel(xlabel, fontsize=16)

    plt.xlim(var.getMin(), var.getMax())
    plt.ylim(-3, 3)
    bx.get_yaxis().set_ticks([-3, 0, 3])

    plt.sca(ax)
    return ax, width

def calc_pull(x, f, y, yerr):
    from numpy import zeros, NaN

    try:
        yerr_low, yerr_high = yerr
        yerr = zeros(len(yerr_low))
        lower = y > f(x)
        yerr[lower] = yerr_low[lower]
        yerr[~lower] = yerr_high[~lower]
    except ValueError:
        pass
    yerr[yerr == 0] = NaN
    pull = (y - f(x)) / yerr

    return pull

def get_function(x, model, data, components=None, norm=1, numcpus=1, extra_params=None):
    if not extra_params:
        extra_params = []

    from numpy import vectorize
    from ROOT import RooCurve, Double, RooFit

    frame = x.frame()

    data.plotOn(frame)
    model.plotOn(frame, RooFit.NumCPU(numcpus), *extra_params)

    if components:
        for c in components:
            model.plotOn(frame, RooFit.NumCPU(numcpus), RooFit.Components(c), *extra_params)

    funcs = []
    for idx in range(1, int(frame.numItems())):
        curr = frame.getObject(idx)
        funcs.append(vectorize(lambda x, curr=curr: norm * curr.Eval(x)))
    
    if len(funcs) > 1:
        return funcs[0], funcs[1:]
    else:
        return funcs[0]

def get_binned_data(x, data, extra_params=None, binning=None):
    if not extra_params:
        extra_params = []
    if binning:
        from ROOT import RooFit
        extra_params.append(RooFit.Binning(binning))

    from numpy import array
    from ROOT import RooHist, Double

    frame = x.frame()

    data.plotOn(frame, *extra_params)

    x = []
    y = []
    x_err_low = []
    x_err_high = []
    y_err_low = []
    y_err_high = []

    data = frame.getObject(0)

    d1 = Double()
    d2 = Double()

    ret = 0
    i = 0
    while ret != -1:
        ret = data.GetPoint(i, d1, d2)
        x.append(float(d1))
        y.append(float(d2))
        x_err_low.append(data.GetErrorXlow(i))
        x_err_high.append(data.GetErrorXhigh(i))
        y_err_low.append(data.GetErrorYlow(i))
        y_err_high.append(data.GetErrorYhigh(i))
        i += 1

    width = array(x_err_low) + array(x_err_high)
    return array(x), width, array(y), array([y_err_low, y_err_high])
