from util import *
from analysis.log import get_logger
logger = get_logger()
from root_pandas import read_root
from missing_hep import histpoints

def plot(data, plotfile, mcfile=None, cuts=None, variables=None, bins=30):
    import numpy as np
    from root_numpy import root2array
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.backends.backend_pdf import PdfPages
    #sns.set_palette("deep", desat=.6)
    #sns.set_context('talk')

    if cuts is None:
        cuts = []

    if variables is None:
        variables = []

    arr = read_root(data, where=prepare_sel(cuts))

    if mcfile:
        arr_mc = read_root(mcfile, where=prepare_sel(cuts))

    logger.info('Saving plots to {}'.format(plotfile))
    with PdfPages(plotfile) as pdf:
        for col in arr.columns:
            logger.debug('Plotting ' + col)
            x = arr[col]
            n, bine, _ = plt.hist(x.values, histtype='stepfilled', bins=bins, color='blue')

            if mcfile:
                x_mc = arr_mc[col]
                if col in arr_mc.columns:
                    n_mc, edges = np.histogram(arr_mc[col], bine)
                    binned_hist(plt.gca(), n_mc, edges, histtype='stepfilled', alpha=0.7)

                    #plt.hist(x_mc, histtype='stepfilled', bins=bins, alpha=0.8, normed=True)
            #plt.yscale('log')
            plt.xlabel(col)
            plt.ylim(0, max(n) * 1.05)
            pdf.savefig()
            plt.clf()

        #logger.info('Plotting m_B vs. q^2')
        #jp = sns.jointplot(arr['B_M'], arr['Psi_M'], kind='hex', joint_kws={'norm': LogNorm()}, stat_func=None)
        #jp.set_axis_labels('$m(B^0)$', '$m(\\mu^+\\mu^-)$')
        #plt.tight_layout()
        #pdf.savefig()
        #plt.clf()

        #logger.info('Plotting m_D vs. q^2')
        #jp = sns.jointplot(arr['D~0_M'], arr['Psi_M'], kind='hex', joint_kws={'norm': LogNorm()}, stat_func=None)
        #jp.set_axis_labels('$m_(\\bar{D}^0)$', '$m(\\mu^+\\mu^-)$')
        #plt.tight_layout()
        #pdf.savefig()
        #plt.clf()

        
def plot_roofit(var, data, model, components=None, numcpus=1, xlabel='', extra_params=None, norm=None, log=False, binning=None, labels=None, plot_pull=True):
    if not extra_params:
        extra_params = []
    if not norm:
        norm = 1 #data.numEntries()

    from numpy import linspace, maximum, array

    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.ticker import MaxNLocator
    import matplotlib as mpl

    cax = plt.gca()
    box = cax.get_position()
    xmin, ymin = box.xmin, box.ymin
    xmax, ymax = box.xmax, box.ymax

    if plot_pull:
        gs = GridSpec(2, 1, height_ratios=[7, 1], left=xmin, right=xmax, bottom=ymin, top=ymax, wspace=0., hspace=0.15)
    else:
        gs = GridSpec(1, 1, left=xmin, right=xmax, bottom=ymin, top=ymax)

    frame = var.frame()

    xpos, width, y, yerr = get_binned_data(var, frame, data, extra_params=extra_params, binning=binning)
    x = linspace(var.getMin(), var.getMax(), 200)

    f, comps = get_function(var, frame, model, components=components, norm=norm, extra_params=extra_params)

    ax = plt.subplot(gs[0])

    for c in comps:
        plt.plot(x, c(x), '--', dashes=[8,3], zorder=-10, clip_on=True)

    #if components:
    #    for c, name in zip(comps, components):
    #        if labels:
    #            plt.plot(x, c(x), '--', dashes=[8,3], zorder=90, clip_on=False, label=labels[name])
    #        else:
    #            plt.plot(x, c(x), '--', dashes=[8,3], zorder=90, clip_on=False)

    plt.plot(x, f(x), 'r-', zorder=95, clip_on=False)
    if log:
        ax.set_yscale('log')
        ax.set_ylim(max(0.1, 0.5*min(y)), 1.2 * max(y))
        yerr[0] = y - maximum(1e-1, y - yerr[0])
    else:
        pass
        #ax.set_ylim(0, 1.1 * max(y))
        #plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))

    plt.errorbar(xpos, y, yerr=yerr, fmt='o', color='k', zorder=100, markersize=1, linewidth=0.5, clip_on=False, capsize=0)

    if plot_pull:
        plt.setp(ax.get_xticklabels(), visible=False)
        bx = plt.subplot(gs[1], sharex=ax)
        plt.setp(bx, zorder=-9)# Make sure the main histogram can draw over the pull plot
        pull = calc_pull(xpos, f, y, yerr)

        plt.fill_between([var.getMin(), var.getMax()], 2, 1, facecolor='#bbbbbb', linewidth=0, edgecolor='#bbbbbb')
        plt.fill_between([var.getMin(), var.getMax()], -2, -1, facecolor='#bbbbbb', linewidth=0, edgecolor='#bbbbbb')

        color = ['#cc4444' if abs(p) > 3 else '#555555' for p in pull]
        plt.bar(xpos-width/2, pull, width, color=color, linewidth=0.5)

        #plt.axhline(0, color='black')
        #plt.axhline(3, color='black')
        plt.ylabel('$\\frac{\\hat{n}_i -  n_i}{\\sigma(n_i)}$')

        plt.ylim(-3, 3)
        bx.get_yaxis().set_ticks([-3, 0, 3])

    plt.xlim(var.getMin(), var.getMax())

    if xlabel:
        plt.xlabel(xlabel)

    plt.sca(ax)
    return gs, ax, width

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

def get_function(x, frame, model, components=None, norm=1, numcpus=1, extra_params=None):
    if not extra_params:
        extra_params = []

    from numpy import vectorize, nan
    from ROOT import RooCurve, Double, RooFit

    model.plotOn(frame, RooFit.NumCPU(numcpus), *extra_params)

    if components:
        for c in components:
            model.plotOn(frame, RooFit.NumCPU(numcpus), RooFit.Components(c), *extra_params)

    # Collect all function objects, their names and ranges
    funcs = []
    for idx in range(1, int(frame.numItems())):
        curr = frame.getObject(idx)
        name = curr.GetName()

        logger.debug('Found plot object with name {}'.format(name))
        
        from ROOT import Double
        xmin = Double(0)
        xmax = Double(0)
        ymin = Double(0)
        ymax = Double(0)
        curr.ComputeRange(xmin, ymin, xmax, ymax)

        funcs.append((curr, name, (xmin, xmax)))

    func_results = [[]]

    import re
    # Process the functions, combine them if they belong together
    curr_component = ''
    for fobj, name, (xmin, xmax) in funcs:

        if 'Comp' in name:
            comp_name = re.findall(r'_Comp\[[^\]]*\]', name)[0]
        else:
            comp_name = ''

        if curr_component == comp_name:
            func_results[-1].append((fobj, (xmin, xmax)))
        else:
            logger.debug('Creating new fit object: {} - {}'.format(comp_name, curr_component))
            curr_component = comp_name
            func_results.append([(fobj, (xmin, xmax))])

    logger.debug('Results is: {}'.format(func_results))
    results = []
    for entry in func_results:
        @vectorize
        def myfunc(x, entry=entry):
            for fobj, (xmin, xmax) in entry:
                if xmin < x < xmax:
                    return norm * fobj.Eval(x)
            else:
                return nan
        
        results.append(myfunc)

    return results[0], results[1:]

def get_binned_data(x, frame, data, extra_params=None, binning=None):
    if not extra_params:
        extra_params = []
    if binning:
        from ROOT import RooFit
        extra_params.append(RooFit.Binning(binning))

    from numpy import array
    from ROOT import RooHist, Double

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
    ret = data.GetPoint(i, d1, d2)
    while ret != -1:
        x.append(float(d1))
        y.append(float(d2))
        x_err_low.append(data.GetErrorXlow(i))
        x_err_high.append(data.GetErrorXhigh(i))
        y_err_low.append(data.GetErrorYlow(i))
        y_err_high.append(data.GetErrorYhigh(i))
        i += 1
        ret = data.GetPoint(i, d1, d2)

    width = array(x_err_low) + array(x_err_high)
    return array(x), width, array(y), array([y_err_low, y_err_high])
