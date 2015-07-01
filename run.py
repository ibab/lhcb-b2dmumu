
from datapipe import *

import sys
import os
import logging
import numpy as np
import pandas as pd
import joblib
from root_pandas import read_root
from analysis.log import setup_logging
setup_logging()
logger = logging.getLogger('analysis')
from analysis.log import setup_roofit
setup_roofit()

DATASTORE='./store/tmp/'

variables_b2kstmumu = [
        '{B,Kstar,Psi}_M',
        '{B,Kstar,Psi}_P',
        '{B,Kstar,Psi}_PT',
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

        'Kstar_FD_ORIVX',
        'Kstar_CosTheta',
        'Kstar_DIRA_OWNPV',

        '{Kplus,piminus,muplus,muminus}_ProbNN*',
        '{Kplus,piminus,muplus,muminus}_PID*',
        '{Kplus,piminus,muplus,muminus}_hasRich',
        '{Kplus,piminus,muplus,muminus}_TRACK_GhostProb',
        '{Kplus,piminus,muplus,muminus}_TRACK_CHI2NDOF',
        '{Kplus,piminus,muplus,muminus}_isMuonLoose',
        '{Kplus,piminus,muplus,muminus}_isMuon',
        '{Kplus,piminus,muplus,muminus}_CosTheta',
        '{Kplus,piminus,muplus,muminus}_P',
        '{Kplus,piminus,muplus,muminus}_PZ',

        'nTracks',
]

variables_b2dmumu = [
        '{B,D~0,Psi}_M',
        '{B,D~0,Psi}_P',
        '{B,D~0,Psi}_PT',
        '{B,D~0}_TAU',
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

        '{Kplus,piminus,muplus,muminus}_ProbNN*',
        '{Kplus,piminus,muplus,muminus}_PID*',
        '{Kplus,piminus,muplus,muminus}_TRACK_GhostProb',
        '{Kplus,piminus,muplus,muminus}_TRACK_CHI2NDOF',
        '{Kplus,piminus,muplus,muminus}_isMuonLoose',
        '{Kplus,piminus,muplus,muminus}_isMuon',
        '{Kplus,piminus,muplus,muminus}_CosTheta',
        '{Kplus,piminus,muplus,muminus}_P',
        '{Kplus,piminus,muplus,muminus}_PZ',

        'nTracks',
]

mc_variables = [
        'B_BKGCAT',
        '*_TRUEID',
]

class Cut:
    def __init__(self):
        self.cutstring = None

    def add(self, other):
        if not self.cutstring is None:
            self.cutstring = '(' + self.cutstring + ') && (' + other + ')'
        else:
            self.cutstring = other

    def get(self):
        return self.cutstring

class RootAppend(Task):
    infiles = Input()
    outname = Input()

    def outputs(self):
        return LocalFile(DATASTORE + self.outname)

    def run(self):
        from sh import hadd
        out = hadd(['-f'] + [self.outputs().path()] + map(lambda x: x.path(), self.infiles))
        print(out)

class Reduce(Task):
    infile = Input()
    columns = Input()
    treename = Input(default='DecayTree')
    # Remove signal region
    blinded = Input(default=False)

    def outputs(self):
        name = os.path.basename(self.infile.path())
        outfile = LocalFile(DATASTORE + name.replace('.root', '.' + self.__class__.__name__ + '.root'))
        return outfile

    def run(self):
        from analysis.util import calc_tau

        cut = Cut()

        if self.blinded:
            B_mass = 5279
            D_mass = 1864.84
            B_width = 50
            D_width = 50
            #cut.add('((B_M < {}) || (B_M > {})) || ((D~0_M < {}) || (D~0_M > {}))'.format(B_mass - B_width, B_mass + B_width, D_mass - D_width, D_mass + D_width))
            #cut.add('((D~0_M < {}) || (D~0_M > {}))'.format(D_mass - D_width, D_mass + D_width))
            cut.add('((B_M < {}) || (B_M > {}))'.format(B_mass - B_width, B_mass + B_width))

        if 'B_BKGCAT' in self.columns:
            cut.add('B_BKGCAT <= 10')

        df = read_root(self.infile.path(), self.treename, columns=self.columns, where=cut.get())

        df['B_TAU'] = pd.Series(calc_tau(df), index=df.index)
        logger.info('Initial events: {}'.format(len(df)))

        df['B_DiraAngle'] = np.arccos(df['B_DIRA_OWNPV'])
        df['B_ENDVERTEX_CHI2_NDOF'] = df['B_ENDVERTEX_CHI2'] / df['B_ENDVERTEX_NDOF']
        for var in df.columns:
            if 'PZ' in var:
                df[var.replace('PZ', 'ETA')] = np.arctanh(df[var] / df[var.replace('PZ', 'P')])
        df.to_root(self.outputs().path())

class ResamplePID(Task):
    infile = Input()

    def outputs(self):
        name = os.path.basename(self.infile.path())
        outfile = LocalFile(DATASTORE + name.replace('.root', '.' + self.__class__.__name__ + '.root'))
        return outfile

    def run(self):
        import pickle
        from analysis.pid_resample import Resampler
        __import__('__main__').Resampler = Resampler # for pickle

        resamplers = {
            'Kplus': './store/resamplers/Kaon_Stripping20_MagnetUp.pkl',
            'piminus': './store/resamplers/Pi_Stripping20_MagnetUp.pkl',
            'muplus':  './store/resamplers/Mu_Stripping20_MagnetUp.pkl',
            'muminus': './store/resamplers/Mu_Stripping20_MagnetUp.pkl',
        }

        nametrans_pid = {'PIDK': 'CombDLLK',
                         'PIDmu': 'CombDLLmu'}

        nametrans_particle = {'Kplus': 'K',
                              'piminus': 'Pi',
                              'muplus': 'Mu',
                              'muminus': 'Mu',
                              }

        df = read_root(self.infile.path())

        for particle, path in resamplers.items():
            resampler = pickle.load(open(path))
            part = nametrans_particle[particle]
            for pid in ['PIDK', 'PIDmu']:
                key = '{particle}_{pid}'.format(particle=particle, pid=pid)
                df[key + '_OLD'] = df[key]
                res = resampler[part + '_' + nametrans_pid[pid]].sample(df[[particle + '_P', particle + '_ETA', 'nTracks']].values.T)
                df[key] = res
                logger.info('Resampled {} for {}'.format(pid, particle))

        # TODO reuse these dropped samples by resampling them
        df = df.query('Kplus_PIDK > -5')
        df = df.query('muplus_PIDmu > -3')
        df = df.query('muminus_PIDmu > -3')
        df.to_root(self.outputs().path(), 'default')

class ApplyTrigger(Task):
    infile = Input()

    def outputs(self):
        name = os.path.basename(self.infile.path())
        outfile = LocalFile(DATASTORE + name.replace('.root', '.' + self.__class__.__name__ + '.root'))
        return outfile

    def run(self):
        # The signal candidate has to be triggered by one of these strategies
        trigger_selection = [
            'B_L0MuonDecision_TOS              == 1',
            'B_Hlt1TrackAllL0Decision_TOS      == 1',
            'B_Hlt1TrackMuonDecision_TOS       == 1',
            'B_Hlt2Topo2BodyBBDTDecision_TOS   == 1',
            'B_Hlt2Topo3BodyBBDTDecision_TOS   == 1',
            'B_Hlt2Topo4BodyBBDTDecision_TOS   == 1',
            'B_Hlt2TopoMu2BodyBBDTDecision_TOS == 1',
            'B_Hlt2TopoMu3BodyBBDTDecision_TOS == 1',
            'B_Hlt2TopoMu4BodyBBDTDecision_TOS == 1',
            'B_Hlt2SingleMuonDecision_TOS      == 1',
            'B_Hlt2DiMuonDetachedDecision_TOS  == 1',
        ]
        trigger_cut = '(' + ' || '.join(trigger_selection) + ')'
        df = read_root(self.infile.path(), where=trigger_cut)
        df.to_root(self.outputs().path())

class Select(Task):
    infile = Input()
    jpsi_inside = Input(default=False)

    def outputs(self):
        name = os.path.basename(self.infile.path())
        outfile = LocalFile(DATASTORE + name.replace('.root', '.' + self.__class__.__name__ + '.root'))
        efficiency = PyTarget('efficiency')
        return (outfile, efficiency)

    def run(self):
        from analysis.util import prepare_sel

        if self.jpsi_inside:
            selection = [
                'Psi_M > 2850 && Psi_M < 3200',
            ]
        else:
            selection = [
                # Exclude J/psi
                'Psi_M < 2850 | Psi_M > 3200',
                # Kinematic range ends below this
                'Psi_M < 3500',
            ]

        df = read_root(self.infile.path())
        initial = len(df)
        df = df.query(prepare_sel(selection))
        after = len(df)
        df = read_root(self.infile.path(), where=prepare_sel(selection))

        df.to_root(self.outputs()[0].path())

        eff = float(after) / initial

        logger.info('Selection efficiency: {}'.format(eff))

        self.outputs()[1].set(eff)

classifier_variables = [
        'B_DiraAngle',
        'B_TAU',
        'B_ENDVERTEX_CHI2_NDOF',
        'B_P',
        'B_PT',
        'B_ISOLATION_BDT_Soft',
        '{Kplus,piminus}_PIDK',
        '{muplus,muminus}_PIDmu',
        #'{Kplus,piminus,muplus,muminus}_PID{K,mu}',
        '{Kplus,piminus,muplus,muminus}_isMuon',
        #'{Kplus,piminus,muplus,muminus}_TRACK_CHI2NDOF',
        #'D~0_CosTheta',
        # New ideas:
        'B_TAU',
        'D~0_TAU',
]

class ApplyCut(Task):
    infile = Input()
    cuts = Input()
    key = Input(default='')
    insert = Input(default=[])

    def outputs(self):
        if self.key is '':
            keystr = ''
        else:
            keystr = '_{}'.format(self.key)
        return LocalFile(self.infile.path().replace('.root', '.ApplyCut{}.root'.format(keystr)))

    def run(self):
        from analysis.util import prepare_sel

        inserts = []
        for ins in self.insert:
            if isinstance(ins, PyTarget):
                ins = ins.get()
            inserts.insert(ins)

        cuts = self.cuts.format(inserts)

        df = read_root(self.infile.path(), where=prepare_sel(cuts))
        df.to_root(self.outputs().path())

class KFoldTrainAndApply(Task):
    signal = Input()
    background = Input()
    clf = Input()

    def outputs(self):
        return LocalFile(self.signal.path().replace('.root', '.KFoldTrainAndApply.root')), LocalFile(self.signal.path().replace('.root', '.TrainTestSet.root'))

    def run(self):
        clf = self.clf
        step = 1

        select_sidebands = 'B_M > 5800 & B_M < 6300'
        sig = read_root(self.signal.path(), columns=classifier_variables, step=step).dropna()
        bkg = read_root(self.background.path(), columns=classifier_variables, step=step, where=select_sidebands).dropna()
        data = pd.concat([sig, bkg], keys=['sig', 'bkg'])

        logger.info('Using {} events from signal sample'.format(len(sig)))
        logger.info('Using {} events from background sample'.format(len(bkg)))

        X = data.values.astype('float32')
        y = np.append(np.ones(len(sig)), np.zeros(len(bkg)))

        from rep.metaml.folding import FoldingClassifier
        skf = FoldingClassifier(clf, n_folds=5, random_state=0)

        skf.fit(X, y)

        train_data = read_root(self.background.path(), step=step, where=select_sidebands).dropna()
        full_data = read_root(self.background.path(), columns=classifier_variables, where='!(' +  select_sidebands + ')').dropna()
        full_data_allvars = read_root(self.background.path(), where='!(' + select_sidebands + ')').dropna()

        # Get unbiased prediction for train set
        train_probs = skf.predict_proba(X)[:,1]
        logger.debug('{} - {}'.format(len(train_data), len(train_probs[y == 0])))
        train_data['proba'] = train_probs[y == 0]

        # Get max prediction for rest of data
        XX = full_data.values.astype('float32')
        other_probs = skf.predict_proba(full_data.values.astype('float32'), vote_function=lambda xs: np.max(xs[:,:,1], axis=0))
        full_data_allvars['proba'] = other_probs

        # Put them together
        ret = pd.concat([train_data, full_data_allvars], keys=['train', 'other'])

        from scipy.special import logit
        ret['clf'] = logit(ret['proba'])

        ret.to_root(self.outputs()[0].path())

        ret2_vars = dict()
        ret2_vars['y_true'] = y
        ret2_vars['proba'] = skf.predict_proba(X)[:,1]

        ret2 = pd.DataFrame(ret2_vars)

        ret2.to_root(self.outputs()[1].path())


class RooFit(Task):
    infile = Input()
    model = Input()
    model_name = Input(default='model')
    params = Input(default='')
    key = Input(default=0)
    fix_params = Input(default='')
    censor = Input(default='')
    range = Input(default='')

    def outputs(self):
        return [LocalFile(DATASTORE + 'results_{}.params'.format(self.key)),
                PyTarget('workspace_{}'.format(self.key)),
                PyTarget('fitresults_{}'.format(self.key)),
                PyTarget('yield_{}'.format(self.key))]

    def run(self):
        out_params, out_ws, out_results, out_yield = self.outputs()
        out_params = out_params.path()

        import ROOT
        import ROOT.RooFit as RF
        from analysis.fit import mle, assemble_model, load_tree

        ws = assemble_model(self.model.path())
        model = ws.pdf(self.model_name)
        data = load_tree(ws, self.infile.path(), 'default', '')

        if self.fix_params:
            for name, results in self.fix_params.items():
                if isinstance(results, PyTarget):
                    res = results.get().floatParsFinal()
                    var = res.find(name)
                    val = var.getVal()
                else:
                    val = results

                ws.var(name).setVal(val)
                ws.var(name).setConstant(True)

        ROOT.SetOwnership(ws, False)
        if self.params:
            start_params = self.params.path()
        else:
            start_params = None

        # Implement fitting on sub-ranges for censored data
        extra_params = []
        if self.censor:
            ranges = []
            for k, rng in self.censor.items():
                vv = ws.var(k)
                left_name = '{}_leftrange'.format(k)
                right_name = '{}_rightrange'.format(k)
                vv.setRange(left_name, vv.getMin(), rng[0])
                vv.setRange(right_name, rng[1], vv.getMax())
                ranges.append(left_name)
                ranges.append(right_name)

            logger.debug("RANGES: {}".format(ranges))
            rng = ROOT.RooFit.Range(','.join(ranges))
            extra_params.append(rng)

        if self.range:
            ranges = []
            for k, rng in self.range.items():
                vv = ws.var(k)
                thisrange = '{}_thisrange'.format(k)
                vv.setRange(thisrange, rng[0], rng[1])
                ranges.append(thisrange)
            rng = ROOT.RooFit.Range(','.join(ranges))
            extra_params.append(rng)

        results = mle(model, data, out_params=out_params, numcpus=20, extra_params=extra_params)

        ws.var('B_M').setRange('signal', 5279 - 50, 5279 + 50)
        args = ROOT.RooArgSet(ws.var('B_M'), ws.var('D~0_M'))
        integ = data.numEntries() * model.createIntegral(args, ROOT.RooFit.NormSet(args), ROOT.RooFit.Range('signal')).getVal()
        logger.debug('integral: {}'.format(integ))

        #results.Print()
        out_ws.set(ws)
        out_results.set(results)
        out_yield.set(integ)

class CalcExpectedLimit(Task):
    model = Input()
    data = Input()
    fix_params = Input(default=[])
    set_params = Input(default=dict())

    def outputs(self):
        return LocalFile(DATASTORE + 'expected.pdf')

    def run(self):
        from analysis.limit import calc_expected_limit
        import numpy as np

        fix_params = dict()
        set_params = dict()

        for params, args in zip([fix_params, set_params], [self.fix_params, self.set_params]):
            for k, v in args.items():
                if isinstance(v, PyTarget):
                    try:
                        res = v.get().floatParsFinal()
                        var = res.find(k)
                        ret = (var.getVal(), var.getError())
                    except AttributeError:
                        ret = v.get()
                elif isinstance(v, tuple):
                    a, b = v
                    logger.warn('{} - {}'.format(a, b))
                    if isinstance(a, PyTarget):
                        a = a.get()
                    if isinstance(b, PyTarget):
                        b = b.get()
                    ret = (a, b)
                else:
                    ret = v
                params[k] = ret

        limits = calc_expected_limit(self.model.path(), self.data.path(), fix_params, set_params)
        logger.info('{1} |-- {0} --| {2}'.format(np.median(limits), np.percentile(limits, 10), np.percentile(limits, 90)))

class PlotFit(Task):
    infile = Input()
    inws = Input()
    path = Input()
    model_name = Input(default='model')
    plot_var = Input(default='B_M')
    components = Input(default=[])
    binning = Input(default=[])
    range = Input(default=[])
    log = Input(default=False)

    def outputs(self):
        return LocalFile(self.path)

    def run(self):
        import ROOT
        from analysis.plotting import plot_roofit
        from analysis.fit import load_tree
        import matplotlib.pyplot as plt
        ws = self.inws.get()
        model = ws.pdf(self.model_name)
        data = load_tree(ws, self.infile.path(), 'default', '')
        v = ws.var(self.plot_var)
        plt.figure(figsize=(12, 8))

        extra_params = []
        if self.plot_var == 'B_M':
            pass
            #extra_params.append(ROOT.RooFit.Range('B_M_leftrange,B_M_rightrange'))
            #extra_params.append(ROOT.RooFit.NormRange('B_M_leftrange,B_M_rightrange'))
        elif self.plot_var == 'D~0_M':
            pass
            #extra_params.append(ROOT.RooFit.Range('B_M_leftrange,B_M_rightrange'))
            extra_params.append(ROOT.RooFit.NormRange('B_M_leftrange,B_M_rightrange'))

        if self.range:
            v.setMin(self.range[0])
            v.setMax(self.range[1])

        gs, ax, width = plot_roofit(
                            v, data, model,
                            components=self.components,
                            numcpus=20,
                            xlabel='$m(K^+\\!\\pi^-\\!\\mu^+\\!\\mu^-)$',
                            binning=self.binning,
                            log=self.log,
                            #extra_params=extra_params,
                        )

        plt.ylabel('Candidates', ha='right', y=1)
        gs.tight_layout(plt.gcf())
        plt.savefig(self.outputs().path())
        plt.clf()

        c1 = ROOT.TCanvas()
        frame = v.frame()
        data.plotOn(frame)
        model.plotOn(frame)
        frame.Draw()
        c1.SetLogy();
        c1.SaveAs(self.outputs().path().replace('.pdf', '_ROOT.pdf'))

class CalcSWeights(Task):
    infile = Input()
    inws = Input()

    def outputs(self):
        return LocalFile(self.infile.path().replace('.root', '.' + self.__class__.__name__ + '.root'))

    def run(self):
        from analysis.fit import add_weights, load_tree
        from root_numpy import tree2rec
        import ROOT
        ROOT.RooAbsData.setDefaultStorageType(ROOT.RooAbsData.Tree)
        ws = self.inws.get()
        model = ws.pdf('model')
        data = load_tree(ws, self.infile.path(), 'default', '')
        sdata = add_weights(model, data, ['sigYield', 'bkgYield'])
        tf = ROOT.TFile(self.outputs().path(), 'recreate')
        tt = data.tree()
        tt.Write('default')
        tf.Write()
        ROOT.SetOwnership(ws, False)

class RunNotebook(Task):
    notebook = Input()
    dependencies = Input()

    def outputs(self):
        return LocalFile(DATASTORE + os.path.basename(self.notebook.path()))

    def run(self):
        from sh import runipy
        nbpath = self.notebook.path()

        runipy([nbpath, self.outputs().path()], _out='/dev/stdout', _err='/dev/stderr')

class CalculateOptimalMetric(Task):
    signal = Input()
    background = Input()
    traintest = Input()

    def outputs(self):
        return PyTarget('OptimalThreshold')

    def run(self):
        if isinstance(self.signal, PyTarget):
            s = self.signal.get()
        else:
            s = self.signal

        if isinstance(self.background, PyTarget):
            b = self.background.get()
        else:
            b = self.background

        def punzi(s, b, sigma=5):
            return s / (np.sqrt(b) + sigma / 2)

        from rep.report.metrics import OptimalMetric
        metric = OptimalMetric(punzi, s, b)

        from root_pandas import read_root
        df = read_root(self.traintest.path())

        p1 = df.proba.ravel()

        proba = np.zeros((p1.shape[0], 2))
        proba[:,1] = p1

        thresh, m_values = metric.compute(df.y_true, proba)

        from scipy.special import logit

        x = logit(thresh)

        import matplotlib.pyplot as plt
        plt.plot(x, m_values)
        plt.savefig('test.pdf')

        val = x[np.argmax(m_values)]

        logger.info('Optimal FOM threshold: {}'.format(val))

        self.outputs().set(val)

if __name__ == '__main__':

    b2dmumu = {
            'name': 'Bd_D0mumu',
            'contains_jpsi': False,
    }
    b2djpsi = {
            'name': 'Bd_D0Jpsi',
            'contains_jpsi': True,
    }

    # PROCESS SIGNAL
    for decay in [b2dmumu]:
        decay['inputs'] = [
                # Same data files used for mumu and Jpsi
                LocalFile('./store/DATA_Bd_D0mumu_MU11.root'),
                LocalFile('./store/DATA_Bd_D0mumu_MD11.root'),
                LocalFile('./store/DATA_Bd_D0mumu_MU12.root'),
                LocalFile('./store/DATA_Bd_D0mumu_MD12.root'),
        ]

        decay['mc_inputs'] = [
                LocalFile('./store/SIM_{}_MD12.root'.format(decay['name'])),
                LocalFile('./store/SIM_{}_MU12.root'.format(decay['name'])),
        ]

        # Prepare data
        decay['input']     = RootAppend(decay['inputs'], 'DATA_B2D0mumu_ALL.root').outputs()
        decay['reduced']   = Reduce(decay['input'], variables_b2dmumu, treename='B2XMuMu_Line_TupleDST/DecayTree', blinded=True).outputs()
        decay['triggered'] = ApplyTrigger(decay['reduced']).outputs()
        decay['selected'], decay['selected_eff'] = Select(decay['triggered'], jpsi_inside=decay['contains_jpsi']).outputs()

        # Prepare simulation
        decay['mc_input']     = RootAppend(decay['mc_inputs'], 'SIM_Bd_D0mumu_ALL.root').outputs()
        decay['mc_reduced']   = Reduce(decay['mc_input'], variables_b2dmumu + mc_variables, treename='B2XMuMu_Line_TupleMC/DecayTree').outputs()
        decay['mc_resampled'] = ResamplePID(decay['mc_reduced']).outputs()
        decay['mc_triggered'] = ApplyTrigger(decay['mc_resampled']).outputs()
        decay['mc_selected'], decay['mc_selected_eff'] = Select(decay['mc_triggered'], jpsi_inside=decay['contains_jpsi']).outputs()

        # Train and apply classifier
        from rep.estimators.xgboost import XGBoostClassifier
        clf = XGBoostClassifier(n_estimators=150, gamma=12, max_depth=10, verbose=1, nthreads=4)

        #classified_b2dmumu_debug = KFoldCrossValidation(signal=selected_b2dmumu_mc, background=selected_b2dmumu, clf=clf).outputs()
        decay['classified'], decay['traintest'] = KFoldTrainAndApply(signal=decay['mc_selected'], background=decay['selected'], clf=clf).outputs()

        decay['model'] = LocalFile('models/Bd_D0mumu.model')

        bkg_only_fit_precut = RooFit(
                          decay['classified'],
                          decay['model'],
                          model_name='fullBkgMassPdf',
                          key=3,
                       ).outputs()

        bkg_yield_precut = bkg_only_fit_precut[3]
        decay['fom'] = CalculateOptimalMetric(1., bkg_yield_precut, decay['traintest']).outputs()

        decay['classified_cut'] = ApplyCut(decay['classified'], ['clf > {}'], insert=[decay['fom']]).outputs()

        # Perform fits to get parameters for expected limit
        sig_only_fit = RooFit(
                decay['mc_selected'],
                decay['model'],
                model_name='sigMassPdf',
                #range={'B_M': (5210, 5350)},
                key=1,
            ).outputs()
        plot_sig_only_fit = PlotFit(
                decay['mc_selected'],
                sig_only_fit[1],
                model_name='sigMassPdf',
                components=['sigMassPdf1', 'sigMassPdf2'],
                path=DATASTORE + 'b2dmumu_sig_only_fit.pdf',
                range=(5200, 5350)
            ).outputs()
        plot_sig_only_fit_d = PlotFit(
                decay['mc_selected'],
                sig_only_fit[1],
                plot_var='D~0_M',
                model_name='sigMassPdf',
                components=['sigMassPdf1', 'sigMassPdf2'],
                path=DATASTORE + 'b2dmumu_sig_only_fit_d.pdf',
            ).outputs()

        bkg_only_fit = RooFit(
                          decay['classified_cut'],
                          decay['model'],
                          model_name='fullBkgMassPdf',
                          key=2,
                       ).outputs()

        plot_bkg_only_fit = PlotFit(
                                decay['classified_cut'],
                                bkg_only_fit[1],
                                model_name='fullBkgMassPdf',
                                path=DATASTORE + 'b2dmumu_bkg_only_fit.pdf',
                                binning=100,
                                log=False,
                            ).outputs()

        plot_bkg_only_fit_d = PlotFit(
                                decay['classified_cut'],
                                bkg_only_fit[1],
                                plot_var='D~0_M',
                                model_name='fullBkgMassPdf',
                                path=DATASTORE + 'b2dmumu_bkg_only_fit_d.pdf',
                                binning=100,
                                log=False,
                            ).outputs()

        # Calculate the expected limit

        sig_only_fitresults = sig_only_fit[2]
        bkg_only_fitresults = bkg_only_fit[2]
        bkg_only_yield = bkg_only_fit[3]

        decay['expected'] = CalcExpectedLimit(
                decay['model'],
                decay['classified_cut'],
                fix_params={
                    'sigFracB':       sig_only_fitresults,
                    'sigFracD':       sig_only_fitresults,
                    'sigMassMean':    sig_only_fitresults,
                    'sigMassSigma1':  sig_only_fitresults,
                    'sigMassSigma2':  sig_only_fitresults,
                    'sigMassMeanD':   sig_only_fitresults,
                    'sigMassSigmaD1': sig_only_fitresults,
                    'sigMassSigmaD2': sig_only_fitresults,

                    'bkgFrac':        bkg_only_fitresults,
                    'bkgMassSlopeB':  bkg_only_fitresults,
                    'bkgMassSlopeD':  bkg_only_fitresults,
                    'lbgMassSlopeB':  bkg_only_fitresults,
                    'lbgMassMeanD':   bkg_only_fitresults,
                    'lbgMassSigmaD':  bkg_only_fitresults,
                },
                set_params={
                    'bkgYield':      (bkg_only_yield, 1),
                },
                ).outputs()

    """
    # Control channel: B0 -> K* mu mu
    inputs_b2kstjpsi_mc = [
            LocalFile('./store/SIM_Bd_KstJpsi_MD12.root'),
            LocalFile('./store/SIM_Bd_KstJpsi_MU12.root'),
    ]
    inputs_b2kstmumu = [
            LocalFile('./store/DATA_Bd_Kst0mumu_MD11.root'),
            LocalFile('./store/DATA_Bd_Kst0mumu_MU11.root'),
            LocalFile('./store/DATA_Bd_Kst0mumu_MD12.root'),
            LocalFile('./store/DATA_Bd_Kst0mumu_MU12.root'),
    ]
    input_b2kstjpsi_mc    = RootAppend(inputs_b2kstjpsi_mc, 'SIM_Bd_KstJpsi_ALL.root').outputs()
    input_b2kstmumu       = RootAppend(inputs_b2kstmumu, 'DATA_B2Kstmumu_ALL.root').outputs()
    model_b2kstmumu       = LocalFile('models/Bd_KstJpsi_CBall.model')
    init_params_b2kstmumu = LocalFile('models/Bd_KstJpsi_CBall.params')
    control_channel       = LocalFile('control-channel.ipynb')

    reduced_b2kstmumu    = Reduce(input_b2kstmumu, variables_b2kstmumu).outputs()
    triggered_b2kstmumu  = ApplyTrigger(reduced_b2kstmumu).outputs()
    cut_b2kstmumu        = ApplyCut(triggered_b2kstmumu, ['B_M > 5100', 'B_M < 5500', 'Kstar_M > 896 - 150', 'Kstar_M < 896 + 150', 'Psi_M > 3000', 'Psi_M < 3200', 'Kplus_PIDK > -5']).outputs()
    classified_b2kstmumu = ApplyClassifier(cut_b2kstmumu, clf).outputs()
    fit_b2kstmumu        = RooFit(classified_b2kstmumu, model_b2kstmumu, params=init_params_b2kstmumu, key='test').outputs()
    plot_b2kstmumu       = PlotFit(cut_b2kstmumu,
                                   fit_b2kstmumu[1],
                                   path=DATASTORE + 'b2kstmumu_data_fit.pdf',
                                   components=['sigMassPdf1', 'sigMassPdf2', 'bkgMassPdf']).outputs()
    weighted_b2kstmumu   = CalcSWeights(cut_b2kstmumu, fit_b2kstmumu[1]).outputs()
    control_channel      = RunNotebook(control_channel, [weighted_b2kstmumu]).outputs()
    """

    #require([b2dmumu['fom']])
    require([plot_bkg_only_fit, plot_bkg_only_fit_d, plot_sig_only_fit, plot_sig_only_fit_d])
