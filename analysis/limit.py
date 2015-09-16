from analysis.log import get_logger
import numpy as np
logger = get_logger()

def calc_expected_limit(model_file, data, fix_params=None, set_params=None):
        import ROOT
        from analysis.fit import mle, assemble_model, load_tree

        ROOT.RooAbsData.setDefaultStorageType(ROOT.RooAbsData.Tree)

        ws = assemble_model(model_file)
        model = ws.pdf('model_constrained')

        data = load_tree(ws, data, 'default', '')

        if fix_params:
            for name, val in fix_params.items():
                logger.warn('{}, {}'.format(name, val))
                ws.var(name).setVal(val[0])
                ws.var(name).setConstant(True)

        if set_params:
            for name, val in set_params.items():
                ws.var(name).setVal(val[0])
                ws.var(name).setError(val[1])

        #bkg = ws.var('bkgYield')
        #bkg.setMin(0)
        #bkg.setMax(100)
        #bkg.setVal(100)
        #bkg.setConstant(True)

        # This should be a constraint
        #ws.var('alpha').setConstant()

        ws.Print()

        yld = ws.var('signalBR')

        # These are important
        yld.setMin(0)
        yld.setMax(5e-6)
        yld.setVal(0)

        mc = ROOT.RooStats.ModelConfig(ws)
        mc.SetPdf(model)
        interest = ROOT.RooArgSet(yld)
        mc.SetParametersOfInterest(interest)
        nuisance = ROOT.RooArgSet(ws.var('bkgYield'), ws.var('alpha'))
        mc.SetNuisanceParameters(nuisance)
        variables = ROOT.RooArgSet(yld, ws.var('bkgYield'), ws.var('alpha'))

        ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
        #ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.DEBUG)
        #ROOT.RooMsgService.instance().addStream(ROOT.RooFit.DEBUG,ROOT.RooFit.Topic(ROOT.RooFit.Tracing),ROOT.RooFit.ClassName("RooGaussian")) ;

        datasets = []
        for i in range(10000):
            data_this = data.Clone()
            # This doesn't take into account the fit error
            expected = ws.var('bkgYield').getVal()
            actual = np.random.poisson(expected)
            args = ROOT.RooArgSet(ws.var('B_M'), ws.var('Meson_M'))
            ws.var('B_M').setMin(5279 - 50)
            ws.var('B_M').setMax(5279 + 50)
            data_gen = model.generate(args, actual)
            ws.var('B_M').setMin(4900)
            ws.var('B_M').setMax(7000)
            logger.debug('numEntries: {}'.format(data_gen.numEntries()))
            data_this.append(data_gen)
            if i == 999:
                f = ROOT.TFile('test.root', 'recreate')
                data_this.tree().Write('default')
                f.Write()
            datasets.append(data_this)

        limits = []

        ws.saveSnapshot('snap', variables)

        for i, data in enumerate(datasets):
            logger.info('{}'.format(i))
            ws.loadSnapshot('snap')
            logger.info('{}'.format(ws.var('alpha').getVal()))
            mc.SetName('config_{}'.format(i))
            #data.Print()
            plc = ROOT.RooStats.ProfileLikelihoodCalculator(data, mc)
            plc.SetTestSize(.10)
            plc_interval = plc.GetInterval()

            #fc = ROOT.RooStats.FeldmanCousins(data, mc)
            #fc.UseAdaptiveSampling(True)
            #fc.FluctuateNumDataEntries(False)
            #fc.SetNBins(200)
            #fc.SetTestSize(.05)
            #fc_interval = fc.GetInterval()

            lim = plc_interval.UpperLimit(yld)
            logger.info('UPPER LIMIT LR - {}'.format(lim))
            #lim2 = fc_interval.UpperLimit(yld)
            #logger.info('UPPER LIMIT FC - {}'.format(lim2))
            plc_interval.Delete()

            #fcul = fc_interval.UpperLimit(yld)

            limits.append(lim)

        limits = filter(lambda x: 0 < x <= 5e-7, limits)

        import matplotlib.pyplot as plt
        plt.hist(limits, bins=40, histtype='step', color='k')
        plt.axvline(np.median(limits), c='r')
        plt.axvline(np.percentile(limits, (100-68)/2.), c='r', ls='--')
        plt.axvline(np.percentile(limits, 100 - (100 - 68)/2.), c='r', ls='--')
        plt.xlabel('$90\%\ \mathrm{upper}\ \\mathrm{CL}\ \mathrm{on}\ \mathrm{BR}(B^0\\to \\overline{D}^0\\mu^+\\mu^-)$')
        plt.tight_layout()
        plt.savefig('store/tmp/expected.pdf')
        plt.savefig('store/tmp/expected.pgf')
        return limits
