from util import *
import logging

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

    arr = root2array(data, 'Bd2D0MuMu', selection=prepare_sel(cuts))

    if mcfile:
        mc = mcfile.replace('.root', '.classified.root')
        arr_mc = root2array(mc, 'Bd2D0MuMu', selection=prepare_sel(cuts))
        total_mc = len(root2array(mc, 'Bd2D0MuMu'))
        print(total_mc)
        #factor = float(50) / total_mc
        #factor = len(arr) / len(arr_mc)
        factor = 0.01

    logging.info('Saving plots to {}'.format(plotfile))
    with PdfPages(plotfile) as pdf:
        for vname in variables:
            logging.debug('Plotting ' + vname)
            #x = arr[vname][arr[vname] > -1000]
            #x_mc = arr_mc[vname][arr_mc[vname] > -1000]
            x = arr[vname]
            n, bine, _ = plt.hist(x, histtype='stepfilled', bins=bins, alpha=0.7)

            if mcfile:
                x_mc = arr_mc[vname]
                if vname in arr_mc.dtype.names:
                    n_mc, edges = np.histogram(arr_mc[vname], bine)
                    binned_hist(plt.gca(), factor * n_mc, edges, histtype='stepfilled', alpha=0.7)

                    #plt.hist(x_mc, histtype='stepfilled', bins=bins, alpha=0.8, normed=True)
            #plt.yscale('log')
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
        jp = sns.jointplot(arr['D~0_M'], arr['Psi_M'], kind='hex', joint_kws={'norm': LogNorm()})
        jp.set_axis_labels('$m_{\\bar{D}^0}$', '$q^2_{\\mu\\mu}$')
        plt.tight_layout()
        pdf.savefig()
        plt.clf()

