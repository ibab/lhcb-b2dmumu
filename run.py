
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
logger.setLevel(logging.DEBUG)
from analysis.log import setup_roofit
setup_roofit()

from ROOT import gSystem

gSystem.Load("models/RooExpAndGauss.so")

DATASTORE='./store/tmp/'

def filter_root(infile, outfile, selections, intree='default', outtree='default'):
    import ROOT
    tf = ROOT.TFile(infile)
    intree = tf.Get(intree)

    numbers = []

    numbers.append(intree.GetEntries())
    tmptree = intree
    of = ROOT.TFile(outfile, "RECREATE")
    for sel in selections:
        tmptree = tmptree.CopyTree(sel)
        numbers.append(tmptree.GetEntries())

    tmptree.Write(outtree)
    tf.Close()
    of.Close()

    return numbers

variables_b2kstmumu = [
        '{B,Meson,Jpsi}_M',
        'B_ConstrainJpsi_M[0]',
        'B_TAU',
        'B_P',
        'B_PZ',
        'B_PT',
        'B_DIRA_OWNPV',
        'B_ENDVERTEX_CHI2',
        'B_ENDVERTEX_NDOF',
        'B_FD_OWNPV',
        'B_ISOLATION_BDT_Soft',

        '{Kplus,piminus,muplus,muminus}_ProbNN*',
        '{Kplus,piminus,muplus,muminus}_PID*',
        #'{Kplus,piminus,muplus,muminus}_TRACK_GhostProb',
        '{Kplus,piminus,muplus,muminus}_TRACK_CHI2NDOF',
        '{Kplus,piminus,muplus,muminus}_isMuon',
        #'{Kplus,piminus,muplus,muminus}_CosTheta',
        '{Kplus,piminus,muplus,muminus}_P',
        '{Kplus,piminus,muplus,muminus}_PT',
        '{Kplus,piminus,muplus,muminus}_P{X,Y,Z}',

        'Jpsi_P{,X,Y,Z}',
        'Meson_P{,X,Y,Z}',

        'B_L0MuonDecision_TOS',
        'B_L0DiMuonDecision_TOS',
        'B_Hlt1TrackAllL0Decision_TOS',
        'B_Hlt1TrackMuonDecision_TOS',
        'B_Hlt2Topo{2,3,4}BodyBBDTDecision_TOS',
        'B_Hlt2TopoMu{2,3,4}BodyBBDTDecision_TOS',
        'B_Hlt2SingleMuonDecision_TOS',
        'B_Hlt2DiMuonDetachedDecision_TOS',

        'nTracks',
        'nSPDHits',
]

classifier_variables = [
        'B_DiraAngle',
        'B_TAU',
        'B_ENDVERTEX_CHI2_NDOF',
        'B_P',
        'B_PT',
        'B_ISOLATION_BDT_Soft',
        '{Kplus,piminus,muplus,muminus}_P',
        '{Kplus,piminus,muplus,muminus}_TRACK_CHI2NDOF',
        '{Kplus,piminus}_PID{K,mu}',
        '{muplus,muminus}_PIDmu',
        #'Meson_CosTheta',
        # New ideas:
        #'Meson_TAU',
]

variables_b2dmumu = [
        '{B,Meson,Jpsi}_M',
        '{B,Meson,Jpsi}_P',
        '{B,Meson,Jpsi}_PT',
        '{B,Meson}_TAU',
        'B_ConstrainJpsi_M[0]',
        'B_DIRA_OWNPV',
        'B_FD_OWNPV',
        'B_{OWNPV,ENDVERTEX}_CHI2',
        'B_{OWNPV,ENDVERTEX}_NDOF',
        'B_ISOLATION_BDT_{Hard,Soft}',

        'B_L0MuonDecision_TOS',
        'B_L0DiMuonDecision_TOS',
        'B_Hlt1TrackAllL0Decision_TOS',
        'B_Hlt1TrackMuonDecision_TOS',
        'B_Hlt2Topo{2,3,4}BodyBBDTDecision_TOS',
        'B_Hlt2TopoMu{2,3,4}BodyBBDTDecision_TOS',
        'B_Hlt2SingleMuonDecision_TOS',
        'B_Hlt2DiMuonDetachedDecision_TOS',

        'Jpsi_FD_ORIVX',
        'Jpsi_FDCHI2_ORIVX',
        'Jpsi_P{,X,Y,Z}',

        'Meson_FD_ORIVX',
        'Meson_CosTheta',
        'Meson_DIRA_OWNPV',
        'Meson_P{X,Y,Z}',

        '{Kplus,piminus,muplus,muminus}_ProbNN*',
        '{Kplus,piminus,muplus,muminus}_PID*',
        '{Kplus,piminus,muplus,muminus}_TRACK_GhostProb',
        '{Kplus,piminus,muplus,muminus}_TRACK_CHI2NDOF',
        #'{Kplus,piminus,muplus,muminus}_isMuonLoose',
        '{Kplus,piminus,muplus,muminus}_isMuon',
        '{Kplus,piminus,muplus,muminus}_CosTheta',
        '{Kplus,piminus,muplus,muminus}_P',
        '{Kplus,piminus,muplus,muminus}_PT',
        '{Kplus,piminus,muplus,muminus}_P{X,Y,Z}',

        'nTracks',
        'nSPDHits',
]

mc_variables = [
        'B_BKGCAT',
        '*_TRUEID',
]

class Cuts:
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
    outname = Input(default=False)
    cuts = Input(default=False)

    def outputs(self):
        if self.outname:
            return LocalFile(self.outname)
        else:
            name = os.path.basename(self.infile.path())
            return LocalFile(DATASTORE + name.replace('.root', '.' + self.__class__.__name__ + '.root'))

    def run(self):
        from analysis.util import calc_tau

        cut = Cuts()

        if self.cuts:
            cut.add(self.cuts)

        if self.blinded:
            B_mass = 5279
            D_mass = 1864.84
            B_width = 50
            D_width = 50
            #cut.add('((B_M < {}) || (B_M > {})) || ((Meson_M < {}) || (Meson_M > {}))'.format(B_mass - B_width, B_mass + B_width, D_mass - D_width, D_mass + D_width))
            #cut.add('((Meson_M < {}) || (Meson_M > {}))'.format(D_mass - D_width, D_mass + D_width))
            cut.add('((B_M < {}) || (B_M > {}))'.format(B_mass - B_width, B_mass + B_width))

        if 'B_BKGCAT' in self.columns:
            cut.add('B_BKGCAT <= 10')

        df = read_root(self.infile.path(), self.treename, columns=self.columns, where=cut.get())
        with_nan = len(df)
        df.dropna(inplace=True)
        without_nan = len(df)

        if without_nan < with_nan:
            logger.warn('Dropped {} NaN samples in Reduce'.format(with_nan - without_nan))

        df['B_TAU'] = pd.Series(calc_tau(df), index=df.index)
        logger.info('Initial events: {}'.format(len(df)))

        df['B_DiraAngle'] = np.arccos(df['B_DIRA_OWNPV'])
        df['B_ENDVERTEX_CHI2_NDOF'] = df['B_ENDVERTEX_CHI2'] / df['B_ENDVERTEX_NDOF']
        for var in df.columns:
            if 'PZ' in var:
                df[var.replace('PZ', 'ETA')] = np.arctanh(df[var] / df[var.replace('PZ', 'P')])
        df.to_root(self.outputs().path())

class AddHypos(Task):
    infile = Input()

    def outputs(self):
        name = os.path.basename(self.infile.path())
        outfile = LocalFile(DATASTORE + name.replace('.root', '.' + self.__class__.__name__ + '.root'))
        return outfile

    @staticmethod
    def calc_hypo_mass(a, b, df, a_mass=None, b_mass=None):
        A = a
        B = b

        if a_mass is None:
            A_M = df[a + '_M']
        else:
            A_M = a_mass

        if b_mass is None:
            B_M = df[b + '_M']
        else:
            B_M = b_mass

        X1 = df[A + '_PX']
        Y1 = df[A + '_PY']
        Z1 = df[A + '_PZ']

        X2 = df[B + '_PX']
        Y2 = df[B + '_PY']
        Z2 = df[B + '_PZ']

        E = np.sqrt(A_M**2 + X1**2 + Y1**2 + Z1**2) + np.sqrt(B_M**2 + X2**2 + Y2**2 + Z2**2)
        X = X1 + X2
        Y = Y1 + Y2
        Z = Z1 + Z2

        M = np.sqrt(E**2 - X**2 - Y**2 - Z**2)

        return M

    @staticmethod
    def calc_hypo_mass3(df, a, b, c, a_mass=None, b_mass=None, c_mass=None):
        A = a
        B = b
        C = c

        if a_mass is None:
            A_M = df[a + '_M']
        else:
            A_M = a_mass

        if b_mass is None:
            B_M = df[b + '_M']
        else:
            B_M = b_mass

        if b_mass is None:
            C_M = df[c + '_M']
        else:
            C_M = c_mass

        X1 = df[A + '_PX']
        Y1 = df[A + '_PY']
        Z1 = df[A + '_PZ']

        X2 = df[B + '_PX']
        Y2 = df[B + '_PY']
        Z2 = df[B + '_PZ']

        X3 = df[C + '_PX']
        Y3 = df[C + '_PY']
        Z3 = df[C + '_PZ']

        E = np.sqrt(A_M**2 + X1**2 + Y1**2 + Z1**2) + np.sqrt(B_M**2 + X2**2 + Y2**2 + Z2**2) + np.sqrt(C_M**2 + X3**2 + Y3**2 + Z3**2)
        X = X1 + X2 + X3
        Y = Y1 + Y2 + Y3
        Z = Z1 + Z2 + Z3

        M = np.sqrt(E**2 - X**2 - Y**2 - Z**2)

        return M

    def run(self):

        pi_M = 139.57018
        mu_M = 105.65837
        jpsi_M = 3096.916
        K_M = 493.667
        P_M = 938.272

        try:
            os.remove(self.outputs().path())
        except OSError:
            pass
        for df in read_root(self.infile.path(), chunksize=1000000):
            logger.info('Processing chunk of size {}'.format(len(df)))

            df['Jpsi_SWAPPED_M'] = self.calc_hypo_mass('muplus', 'piminus', df, a_mass=mu_M, b_mass=mu_M)
            df['Dstar_M'] = self.calc_hypo_mass('Meson', 'muminus', df, b_mass=pi_M)
            df['B_FROM_JpsiK_M'] = self.calc_hypo_mass('Jpsi', 'Kplus', df, b_mass=K_M)
            df['Phi_FROM_KK_M'] = self.calc_hypo_mass('Kplus', 'piminus', df, a_mass=K_M, b_mass=K_M)
            df['B_FROM_JpsiPK_M'] = self.calc_hypo_mass3(df, 'Jpsi', 'Kplus', 'piminus', b_mass=K_M, c_mass=P_M)
            df['B_FROM_JpsiPpi_M'] = self.calc_hypo_mass3(df, 'Jpsi', 'Kplus', 'piminus', b_mass=P_M, c_mass=K_M)
            df['Meson_KpiSWAP_M'] = self.calc_hypo_mass('Kplus', 'piminus', df, a_mass=pi_M, b_mass=K_M)

            df.to_root(self.outputs().path(), mode='a')

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
        
        np.random.seed(0)

        df = read_root(self.infile.path())

        for particle, path in resamplers.items():
            resampler = pickle.load(open(path))
            part = nametrans_particle[particle]
            df[particle + '_logP'] = np.log(df[particle + '_P'])
            for pid in ['PIDK', 'PIDmu']:
                key = '{particle}_{pid}'.format(particle=particle, pid=pid)
                df[key + '_OLD'] = df[key]
                res = resampler[part + '_' + nametrans_pid[pid]].sample(df[[particle + '_logP', particle + '_ETA', 'nTracks']].values.T)
                df[key] = res
                logger.info('Resampled {} for {}'.format(pid, particle))
                logger.info('{}'.format(res))

        # TODO reuse these dropped samples by resampling them
        df = df.query('Kplus_PIDK > -5')
        df = df.query('muplus_PIDmu > -3')
        df = df.query('muminus_PIDmu > -3')
        df = df.query('piminus_PIDK < 30')
        df.to_root(self.outputs().path(), 'default')

class Cut(Task):
    infile = Input()
    cuts = Input()
    key = Input(default='')
    insert = Input(default=[])

    def outputs(self):
        if self.key is '':
            keystr = ''
        else:
            keystr = '_{}'.format(self.key)
        return (LocalFile(self.infile.path().replace('.root', '.Cut{}.root'.format(keystr))),
                PyTarget('Cut_eff'),
                PyTarget('all_efficiencies'))

    @staticmethod
    def perform_cut(df, cut):
        initial = len(df)
        df = df.query(cut)
        after = len(df)

        return df, float(after) / initial, initial

    def run(self):
        from analysis.util import prepare_sel

        inserts = []
        for ins in self.insert:
            if isinstance(ins, PyTarget):
                ins = ins.get()
            inserts.append(ins)

        if inserts:
            cuts = map(lambda (cut, ins): cut.format(ins), zip(self.cuts, inserts))
        else:
            cuts = self.cuts

        df = read_root(self.infile.path())
        efficiencies = []
        for cut in cuts:
            logger.info('Performing cut: {}'.format(repr(cut)))
            df, eff, initial = self.perform_cut(df, cut)
            logger.info('Efficiency is: {}'.format(eff))
            efficiencies.append((cut, eff, initial))

        df.to_root(self.outputs()[0].path())

        #numbers = filter_root(self.infile.path(), self.outputs()[0].path(), cuts)

        #efficiencies = []
        #current = numbers[0]
        #for name, n in zip(cuts, numbers[1:]):
        #    efficiencies.append((name, float(n) / current, current))
        #    current = n

        total_eff = np.product(map(lambda x: x[1], efficiencies))
        
        logger.info('Cut total efficiency: {}'.format(total_eff))
        self.outputs()[1].set(total_eff)
        self.outputs()[2].set(efficiencies)

class KFoldTrainAndApply(Task):
    signal = Input()
    background = Input()
    variables = Input()
    clf = Input()
    extra = Input(default=False)

    def outputs(self):
        ret =  [LocalFile(self.signal.path().replace('.root', '.KFoldTrainAndApply.root')),
                LocalFile(self.signal.path().replace('.root', '.TrainTestSet.root'))]
        if self.extra:
            for p in self.extra:
                ret.append(LocalFile(p.path().replace('.root', '.KFoldTrainAndApplyExtra.root')))

        return ret

    def run(self):
        clf = self.clf
        step = 1

        select_sidebands = 'B_M > 5200 & B_M < 5600'
        sig = read_root(self.signal.path(), columns=self.variables, step=step).dropna()
        bkg = read_root(self.background.path(), columns=self.variables, step=step, where=select_sidebands).dropna()
        data = pd.concat([sig, bkg], keys=['sig', 'bkg'])

        weights = read_root(self.signal.path(), columns=['SimpleWeight'])
        #sample_weight = np.append(weights, np.ones(len(bkg)))
        sample_weight = np.append(np.ones(len(sig)), np.ones(len(bkg)))

        logger.info('Using {} events from signal sample'.format(len(sig)))
        logger.info('Using {} events from background sample'.format(len(bkg)))

        X = data.values.astype('float32')
        y = np.append(np.ones(len(sig)), np.zeros(len(bkg)))

        from rep.metaml.folding import FoldingClassifier
        skf = FoldingClassifier(clf, n_folds=5, random_state=0)
        skf.fit(X, y, sample_weight=sample_weight)

        train_data = read_root(self.background.path(), step=step, where=select_sidebands).dropna()
        full_data = read_root(self.background.path(), columns=self.variables, where='!(' +  select_sidebands + ')').dropna()
        full_data_allvars = read_root(self.background.path(), where='!(' + select_sidebands + ')').dropna()

        # Get unbiased prediction for train set
        train_probs = skf.predict_proba(X)[:,1]
        logger.debug('{} - {}'.format(len(train_data), len(train_probs[y == 0])))
        print('DATA:', train_data.shape)
        print('TRAIN:', train_probs[y==0].shape)
        train_data['proba'] = train_probs[y == 0]

        # Get mean prediction for rest of data
        XX = full_data.values.astype('float32')
        other_probs = skf.predict_proba(full_data.values.astype('float32'), vote_function=lambda xs: np.mean(xs[:,:,1], axis=0))
        full_data_allvars['proba'] = other_probs

        # Put them together
        ret = pd.concat([train_data, full_data_allvars], keys=['train', 'other'])

        from scipy.special import logit
        ret['clf'] = logit(ret['proba'])

        ret.to_root(self.outputs()[0].path())

        ret2_vars = dict()
        ret2_vars['y_true'] = y
        ret2_vars['proba'] = skf.predict_proba(X)[:,1]
        ret2_vars['clf'] = logit(ret2_vars['proba'])
        ret2 = pd.DataFrame(ret2_vars)

        if self.extra:
            for p in self.extra:
                data_full = read_root(p.path()).dropna()
                data_clf = read_root(p.path(), columns=self.variables).dropna()
                probs = skf.predict_proba(data_clf.values.astype('float32'), vote_function=lambda xs: np.mean(xs[:,:,1], axis=0))
                data_full['clf'] = logit(probs)
                data_full.to_root(p.path().replace('.root', '.KFoldTrainAndApplyExtra.root'))

        import matplotlib.pyplot as plt
        #import seaborn as sns
        ret2.query('(y_true == 1)')['clf'].hist(bins=30, histtype='stepfilled', lw=2, alpha=0.5, normed=True, grid=False, color='blue', label='signal')
        ret2.query('(y_true == 0)')['clf'].hist(bins=30, histtype='step', normed=True, lw=2, grid=False, color='red', label='background')
        plt.xlabel('classifier response', fontsize=18)
        plt.legend(loc='best')
        plt.savefig(DATASTORE + 'classifier.pdf')
        plt.clf()

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
    rangeX = Input(default=False)
    integral = Input(default=False)

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
            for i, (k, rng) in enumerate(self.range.items()):
                vv = ws.var(k)
                thisrange = '{}_thisrange_{}'.format(k, i)
                vv.setRange(thisrange, rng[0], rng[1])
                ranges.append(thisrange)
            rng = ROOT.RooFit.Range(','.join(ranges))
            extra_params.append(rng)

        if self.rangeX:
            for i, (k, rng) in enumerate(self.rangeX.items()):
                vv = ws.var(k)
                vv.setMin(rng[0])
                vv.setMax(rng[1])

        data = load_tree(ws, self.infile.path(), 'default', '')
        results = mle(model, data, out_params=out_params, numcpus=10, extra_params=extra_params)

        #try:
        #    integ = results.floatParsFinal().find('sigYield').getVal()
        #except:
        #    integ = 0

        if self.integral:
            try:
                ws.var('B_M').setRange('signal', self.integral[0], self.integral[1])
                args = ROOT.RooArgSet(ws.var('B_M'), ws.var('Meson_M'))
                args = ROOT.RooArgSet(ws.var('B_M'))
                integ = data.numEntries() * model.createIntegral(args, ROOT.RooFit.NormSet(args), ROOT.RooFit.Range('signal')).getVal()
                logger.debug('integral: {}'.format(integ))
            except:
                integ = -1
        else:
            integ = -1

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
        import ROOT

        fix_params = dict()
        set_params = dict()

        for params, args in zip([fix_params, set_params], [self.fix_params, self.set_params]):
            for k, v in args.items():
                if isinstance(v, PyTarget):
                    if isinstance(v.get(), ROOT.RooFitResult):
                        res = v.get().floatParsFinal()
                        var = res.find(k)
                        logger.debug('{}'.format(var))
                        logger.debug('{}'.format(v.get()))
                        ret = (var.getVal(), var.getError())
                    else:
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
        expected = np.median(limits)
        logger.info('{0} + {1} - {2}'.format(expected, np.percentile(limits, (100-68)/2.) - expected, expected - np.percentile(limits, 100-(100-68)/2.)))

class PlotFit(Task):
    infile = Input()
    inws = Input()
    path = Input()
    model_name = Input(default='model')
    plot_var = Input(default='B_M')
    components = Input(default=[])
    binning = Input(default=[])
    range = Input(default=[])
    named_range = Input(default=False)
    log = Input(default=False)
    fix_norm = Input(default=[])
    xlabel = Input(default='')
    plot_pull = Input(default=True)
    scale = Input(default=False)

    def outputs(self):
        return LocalFile(self.path)

    def run(self):
        import ROOT
        from analysis.plotting import plot_roofit
        from analysis.fit import load_tree
        import matplotlib.pyplot as plt
        from scale import mksize
        ws = self.inws.get()
        model = ws.pdf(self.model_name)
        data = load_tree(ws, self.infile.path(), 'default', '')
        v = ws.var(self.plot_var)
        plt.figure(figsize=(12, 8))

        extra_params = []

        #if self.plot_var == 'B_M':
        #    extra_params.append(ROOT.RooFit.Range('B_M_leftrange,B_M_rightrange'))
        #    extra_params.append(ROOT.RooFit.NormRange('B_M_leftrange,B_M_rightrange'))
        #elif self.plot_var == 'Meson_M':
        #    #extra_params.append(ROOT.RooFit.Range('B_M_leftrange,B_M_rightrange'))
        #    extra_params.append(ROOT.RooFit.NormRange('B_M_leftrange,B_M_rightrange'))

        if self.range:
            v.setMin(self.range[0])
            v.setMax(self.range[1])

        if self.fix_norm:
            extra_params.append(ROOT.RooFit.Normalization(self.fix_norm))

        if self.scale:
            scale = self.scale
        else:
            scale = 0.8

        plt.figure(figsize=mksize(scale, 0.7))
        gs, ax, width = plot_roofit(
                            v, data, model,
                            components=self.components,
                            numcpus=20,
                            xlabel=self.xlabel,
                            binning=self.binning,
                            log=self.log,
                            extra_params=extra_params,
                            plot_pull=self.plot_pull,
                        )

        plt.ylabel('Candidates')
        #gs.tight_layout(plt.gcf())

        #plt.text(0.75, 0.8,'LHCb',
        #         horizontalalignment='center',
        #         verticalalignment='center',
        #         transform = ax.transAxes,
        #         fontsize=30)

        #import Image
        #im = Image.open('/home/igor/b2dmumu/lhcb.png')
        #height = im.size[1]
        #im = np.array(im).astype(np.float) / 255
        #fig = plt.gcf()
        #fig.figimage(im, 0.7 * fig.bbox.xmax, 0.8 * fig.bbox.ymax)

        plt.savefig(self.outputs().path())
        plt.savefig(self.outputs().path().replace('.pdf', '.png'))
        plt.savefig(self.outputs().path().replace('.pdf', '.pgf'))
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
        tmpname = self.infile.path().replace('.root', '.OnlySWeights.root')
        tf = ROOT.TFile(tmpname, 'recreate')
        tt = data.tree()
        tt.Write('default')
        tf.Write()
        ROOT.SetOwnership(ws, False)

        all_data = read_root(self.infile.path())
        weighted = read_root(tmpname)
        for col in weighted.columns:
            if col.endswith('_sw'):
                all_data[col] = weighted[col].ravel()

        all_data.to_root(self.outputs().path())

class RunNotebook(Task):
    notebook = Input()
    dependencies = Input()

    def outputs(self):
        return LocalFile(DATASTORE + os.path.basename(self.notebook.path()))

    def run(self):
        from sh import runipy
        nbpath = self.notebook.path()

        runipy([nbpath, self.outputs().path()], _out='/dev/stdout', _err='/dev/stderr')

class CalcSimpleWeights(Task):
    mc = Input()
    data = Input()
    data_weights_column = Input(default='sigYield_sw')
    apply_to = Input(default=False)
    variables = Input(default=False)

    def outputs(self):
        if self.apply_to:
            rets = []
            rets.append(LocalFile(self.mc.path().replace('.root', '.' + self.__class__.__name__ + '.root')))
            for p in self.apply_to:
                rets.append(LocalFile(p.path().replace('.root', '.ApplySimpleWeights.root')))
            return rets
        else:
            return LocalFile(self.mc.path().replace('.root', '.' + self.__class__.__name__ + '.root'))

    def run(self):
        import numpy as np
        from root_pandas import read_root
        from hep_ml.reweight import BinsReweighter

        mc = read_root(self.mc.path(), columns=self.variables).values
        data = read_root(self.data.path(), columns=self.variables).values

        logger.info('{} entries in MC dataset'.format(len(mc)))
        logger.info('{} entries in DATA dataset'.format(len(data)))

        weights = read_root(self.data.path(), columns=[self.data_weights_column])[self.data_weights_column].ravel()
        rw = BinsReweighter()
        rw.fit(mc, data, target_weight=weights)
        mc_weights = rw.predict_weights(mc)
        mc_weights *= len(mc) / np.sum(mc_weights)

        all_mc = read_root(self.mc.path())

        #mask = mc_weights < 5
        #logger.info('Lost {} events due to overly large weights'.format(len(mask) - np.sum(mask)))
        #mc_weights *= len(all_mc[mask]) / np.sum(mc_weights[mask])
        # Delete events with overly large weights
        all_mc['SimpleWeight'] = mc_weights
        all_mc.to_root(self.outputs()[0].path())

        if self.apply_to:
            for i, o in zip(self.apply_to, self.outputs()[1:]):
                logger.info('Applying to {}'.format(i.path()))
                app = read_root(i.path(), columns=self.variables).values
                app_all = read_root(i.path())
                weights = rw.predict_weights(app)
                weights *= len(app) / np.sum(weights)
                app_all['SimpleWeight'] = weights
                app_all.to_root(o.path())

class CalcSuperWeights(Task):
    mc = Input()
    data = Input()
    data_weights_column = Input(default='sigYield_sw')
    apply_to = Input(default=False)

    def outputs(self):
        if self.apply_to:
            rets = []
            rets.append(LocalFile(self.mc.path().replace('.root', '.' + self.__class__.__name__ + '.root')))
            for p in self.apply_to:
                rets.append(LocalFile(p.path().replace('.root', '.ApplySuperWeights.root')))
            return rets
        else:
            return LocalFile(self.mc.path().replace('.root', '.' + self.__class__.__name__ + '.root'))

    def run(self):
        import numpy as np
        from root_pandas import read_root

        from rep.estimators.xgboost import XGBoostClassifier
        clf = XGBoostClassifier(n_estimators=300, eta=0.3, gamma=12, max_depth=10, verbose=1, nthreads=20)
        from rep.metaml.folding import FoldingClassifier
        skf = FoldingClassifier(clf, n_folds=5, random_state=0)

        mc = read_root(self.mc.path(), columns=classifier_variables)
        data = read_root(self.data.path(), columns=classifier_variables)
        X = pd.concat([mc, data], keys=['mc', 'data']).values
        y = np.append(np.zeros(len(mc)), np.ones(len(data)))

        logger.info('{} entries in MC dataset'.format(len(mc)))
        logger.info('{} entries in DATA dataset'.format(len(data)))

        weights = read_root(self.data.path(), columns=[self.data_weights_column])[self.data_weights_column].ravel()
        weights = np.append(np.ones(len(mc)), weights)

        logger.info('Training first fold')
        skf.fit(X, y, sample_weight=weights)
        p = skf.predict_proba(X)[:,1]
        calc_weights = p / (1 - p)
        mc_weights = calc_weights[y == 0]
        print(mc_weights)

        all_mc = read_root(self.mc.path())
        mask = mc_weights < 5
        logger.info('Lost {} events due to overly large weights'.format(len(mask) - np.sum(mask)))
        mc_weights *= len(all_mc[mask]) / np.sum(mc_weights[mask])
        # Delete events with overly large weights
        all_mc.ix[mask,'SuperWeight'] = mc_weights[mask]
        all_mc.iloc[mask].to_root(self.outputs()[0].path())

        if self.apply_to:
            for i, o in zip(self.apply_to, self.outputs()[1:]):
                logger.info('Applying to {}'.format(i.path()))
                app = read_root(i.path(), columns=classifier_variables).values
                app_all = read_root(i.path())
                def vote(xs):
                    return np.mean(xs[:,:,1], axis=0)
                app_probs = skf.predict_proba(app, vote_function=vote)
                weights = app_probs / (1 - app_probs)
                mask = weights < 5
                logger.info('Lost {} events due to overly large weights'.format(len(mask) - np.sum(mask)))
                weights *= len(app_all[mask]) / np.sum(weights[mask])
                app_all.ix[mask, 'SuperWeight'] = weights[mask]
                app_all.iloc[mask].to_root(o.path())

class CalculateOptimalMetric(Task):
    signal = Input()
    background = Input()
    traintest = Input()
    weight_column = Input(default=False)

    def outputs(self):
        return PyTarget('OptimalThreshold'), PyTarget('OptimalThreshold_eff')

    def run(self):
        if isinstance(self.signal, PyTarget):
            s = self.signal.get()
        else:
            s = self.signal

        if isinstance(self.background, PyTarget):
            b = self.background.get()
        else:
            b = self.background

        def punzi(s, b, sigma=3):
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
        from scale import mksize
        plt.figure(figsize=mksize(0.8))
        plt.plot(x, m_values, color='cornflowerblue')

        val = x[np.argmax(m_values)]
        plt.xlim(0, 6)
        plt.axvline(val, color='k')
        plt.ylabel('Punzi figure of merit (a=3)')
        plt.xlabel('classifier threshold')
        plt.savefig(DATASTORE + 'fom.pdf')
        plt.savefig(DATASTORE + 'fom.pgf')
        plt.clf()

        logger.info('Optimal FOM threshold: {}'.format(val))

        self.outputs()[0].set(val)

        sig_x = df.clf[df.y_true == 1]
        initial = len(sig_x)
        after = sum(sig_x > val)
        eff = float(after) / initial
        logger.info('CalculateOptimalMetric efficiency: {}'.format(eff))

        self.outputs()[1].set(eff)

class CalculateNormalization(Task):
    signal_effs = Input()
    signal_base = Input()
    norm_effs = Input()
    norm_base = Input()
    norm_BR = Input()
    norm_yield = Input()

    def outputs(self):
        return [PyTarget('alpha'), PyTarget('alpha_err')]

    def run(self):
        from uncertainties import ufloat

        def unpack(x):
            if isinstance(x, PyTarget):
                return x.get()
            else:
                return x

        signal_effs = map(unpack, self.signal_effs)
        signal_base = unpack(self.signal_base)
        norm_effs = map(unpack, self.norm_effs)
        norm_base = unpack(self.norm_base)
        norm_BR = ufloat(self.norm_BR)
        logger.info('Signal efficiencies: {}'.format(signal_effs))
        logger.info('Norm efficiencies: {}'.format(norm_effs))

        res = self.norm_yield.get().floatParsFinal()
        var = res.find('sigYield')
        norm_yield = ufloat(var.getVal(), var.getError())

        logger.info('Normalization yield: {}'.format(norm_yield))

        norm_total_eff = np.product(norm_effs)
        sig_total_eff = np.product(signal_effs) 

        norm_total_eff = ufloat(norm_total_eff, np.sqrt(norm_total_eff * (1 - norm_total_eff) / norm_base))
        sig_total_eff = ufloat(sig_total_eff, np.sqrt(sig_total_eff * (1 - sig_total_eff) / signal_base))

        logger.info('Norm total eff: {}'.format(norm_total_eff))
        logger.info('Signal total eff: {}'.format(sig_total_eff))

        alpha = norm_BR * norm_total_eff / (sig_total_eff * norm_yield)

        logger.info('ALPHA is {}'.format(alpha))
        self.outputs()[0].set(alpha.n)
        self.outputs()[1].set(alpha.s)

if __name__ == '__main__':

    b2dmumu = {
            'name': 'B_Dmumu',
            'contains_jpsi': False,
    }
    b2djpsi = {
            'name': 'Bd_D0Jpsi',
            'contains_jpsi': True,
    }

    inputs = [
            LocalFile('./store/DATA_B_Dmumu_MD11.root'),
            LocalFile('./store/DATA_B_Dmumu_MU11.root'),
            LocalFile('./store/DATA_B_Dmumu_MD12.root'),
            LocalFile('./store/DATA_B_Dmumu_MU12.root'),
    ]
    input = RootAppend(inputs, 'DATA_B_Dmumu_ALL.root').outputs()

    # Control channel: B0 -> Jpsi Kst
    inputs_b2kstjpsi_mc = [
            LocalFile('./store/SIM_B_JpsiKst_MD12.root'),
            LocalFile('./store/SIM_B_JpsiKst_MU12.root'),
    ]
    input_b2kstjpsi_mc    = RootAppend(inputs_b2kstjpsi_mc, 'SIM_B_JpsiKst_ALL.root').outputs()

    b2jpsikst_selection = [
                          '(B_M > 5100 & B_M < 5700)',
                          'Kplus_PIDK > -5',
                          '(Meson_M > 800 & Meson_M < 995.9)',
                          'Jpsi_M > (3070 - 100) & Jpsi_M < (3070 + 100)',
                          # Remove B+->JpsiK+
                          '(B_FROM_JpsiK_M < 5220 | B_FROM_JpsiK_M > 5340)',
                          # Remove B0->JpsiPhi
                          '(Phi_FROM_KK_M < 1020-15 | Phi_FROM_KK_M > 1020+15)',
                          # Remove B0->JpsiK* with K<->pi swap
                          '(Meson_KpiSWAP_M < 895 - 100 | Meson_KpiSWAP_M > 895 + 100 )',
                          #'piminus_ProbNNpi > 0.4',
                          #'Kplus_ProbNNk > 0.4',
                          ]
    
    trigger_selection = [
        '(B_L0MuonDecision_TOS == 1 | B_L0DiMuonDecision_TOS == 1)',
        '(B_Hlt1TrackAllL0Decision_TOS == 1 | B_Hlt1TrackMuonDecision_TOS == 1)',
        '(B_Hlt2Topo2BodyBBDTDecision_TOS == 1 | B_Hlt2Topo3BodyBBDTDecision_TOS == 1 | B_Hlt2Topo4BodyBBDTDecision_TOS == 1 | B_Hlt2TopoMu2BodyBBDTDecision_TOS == 1 | B_Hlt2TopoMu3BodyBBDTDecision_TOS == 1 | B_Hlt2TopoMu4BodyBBDTDecision_TOS == 1 | B_Hlt2SingleMuonDecision_TOS == 1 | B_Hlt2DiMuonDetachedDecision_TOS == 1)',
    ]

    model_b2kstmumu       = LocalFile('models/Bd_JpsiKst_Constr.model')
    model_b2kstmumu_mc    = LocalFile('models/Bd_JpsiKst_Constr_MC.model')
    init_params_b2kstmumu = LocalFile('models/Bd_KstJpsi_CBall.params')
    control_channel       = LocalFile('control-channel.ipynb')

    reduced_b2kstmumu_mc  = Reduce(input_b2kstjpsi_mc, variables_b2kstmumu + mc_variables, treename='MC_B_JpsiKst/default').outputs()
    triggered_b2jpsikst_mc, b2jpsikst_trigger_eff, b2jpsikst_trigger_all_eff = Cut(reduced_b2kstmumu_mc, trigger_selection).outputs()
    hypos_b2jpsikst_mc = AddHypos(triggered_b2jpsikst_mc).outputs()
    cut_b2jpsikst_mc, b2jpsikst_cut_eff, b2jpsikst_all_eff = Cut(hypos_b2jpsikst_mc, b2jpsikst_selection).outputs()
    #input_b2kstmumu       = RootAppend(inputs_b2kstmumu, 'DATA_B2Kstmumu_ALL.root').outputs()
    fit_b2jpsikst_mc      = RooFit(cut_b2jpsikst_mc,
                                   model_b2kstmumu_mc,
                                   #params=init_params_b2kstmumu,
                                   model_name='sigMassPdf',
                                   key='B2JpsiKst_MC_fit').outputs()
    b2jpsikst_mc_fitresults = fit_b2jpsikst_mc[2]

    plot_b2jpsikst_mc     = PlotFit(cut_b2jpsikst_mc,
                                   fit_b2jpsikst_mc[1],
                                   plot_var='B_ConstrainJpsi_M_0',
                                   path=DATASTORE + 'b2kstmumu_mc_fit.pdf',
                                   model_name='sigMassPdf',
                                   #log=True,
                                   components=['sigMassPdf1', 'sigMassPdf2'],
                                   log=True,
                                   plot_pull=False,
                                   xlabel='$m(K^+\\pi^-\\mu^+\\mu^-)$',
                                   #range=(5220, 5340)
                                   ).outputs()

    resampled_b2jpsikst_mc = ResamplePID(cut_b2jpsikst_mc).outputs()


    reduced_b2kstmumu    = Reduce(input, variables_b2kstmumu, treename='DATA_B_JpsiKst/default', outname=DATASTORE+'DATA_B_JpsiKst.Reduce.root', cuts='B_M > 5000 & B_M < 5700').outputs()
    triggered_b2kstmumu, _, _  = Cut(reduced_b2kstmumu, trigger_selection).outputs()
    hypos_b2kstmumu = AddHypos(triggered_b2kstmumu).outputs()
    cut_b2kstmumu, _, _ = Cut(hypos_b2kstmumu, b2jpsikst_selection + ['B_ConstrainJpsi_M_0 > 5220 & B_ConstrainJpsi_M_0 < 5340']).outputs()
    fit_b2kstmumu        = RooFit(cut_b2kstmumu,
                                  model_b2kstmumu,
                                  #model_name='bkgMassPdf',
                                  #params=init_params_b2kstmumu,
                                  fix_params={
                                    #'bkgMassSlope':  0.0044466,
                                    #'bkgYield':      115000,
                                    #'sigMassMean':   b2jpsikst_mc_fitresults,
                                    #'sigMassSigma1': b2jpsikst_mc_fitresults,
                                    #'sigMassSigma2': b2jpsikst_mc_fitresults,
                                    'sigFrac':       b2jpsikst_mc_fitresults,
                                    #'bsYield':       0,
                                    #'deltaM':        0,
                                    'alpha1':        b2jpsikst_mc_fitresults,
                                    'alpha2':        b2jpsikst_mc_fitresults,
                                    'n1':            b2jpsikst_mc_fitresults,
                                    'n2':            b2jpsikst_mc_fitresults,
                                  },
                                  #range={'B_ConstrainJpsi_M_0': (5450, 5600)},
                                  range={'B_ConstrainJpsi_M_0': (5220, 5340)},
                                  key='B2JpsiKst_DATA_fit').outputs()
    plot_b2kstmumu       = PlotFit(cut_b2kstmumu,
                                   fit_b2kstmumu[1],
                                   plot_var='B_ConstrainJpsi_M_0',
                                   path=DATASTORE + 'b2kstmumu_data_fit.pdf',
                                   components=['sigMassPdf', 'bkgMassPdf'],
                                   log=True,
                                   #range=(5220, 5340)
                                   xlabel='$m(K^+\\pi^-\\mu^+\\mu^-)$',
                                   plot_pull=False,
                                   ).outputs()
    weighted_b2kstmumu   = CalcSWeights(cut_b2kstmumu, fit_b2kstmumu[1]).outputs()

    # PROCESS SIGNAL
    for decay in [b2dmumu]:
        decay['mc_inputs'] = [
                LocalFile('./store/SIM_{}_MD12.root'.format(decay['name'])),
                LocalFile('./store/SIM_{}_MU12.root'.format(decay['name'])),
        ]


        b2dmumu_selection = [
          # Exclude J/psi
          'Jpsi_M < 2900 | Jpsi_M > 3200',
          # Kinematic range ends below this
          'Jpsi_M < 3500',
          # Kill mis-id mu- -> pi-
          'Jpsi_SWAPPED_M < 2900 | Jpsi_SWAPPED_M > 3200',
          'Jpsi_SWAPPED_M < 3500 | Jpsi_SWAPPED_M > 3800',
          # Exclude B->D*munu
          'Dstar_M < 1990 | Dstar_M > 2030',
        ]

        # Prepare data
        decay['input'] = input 
        decay['reduced']  = Reduce(decay['input'], variables_b2dmumu, treename='DATA_{}/default'.format(decay['name']), blinded=True).outputs()
        decay['triggered'], decay['trigger_eff'], decay['trigger_all_eff'] = Cut(decay['reduced'], trigger_selection).outputs()
        decay['withhypos'] = AddHypos(decay['triggered']).outputs()
        decay['selected'], decay['selected_eff'], decay['selected_all_eff'] = Cut(decay['withhypos'], cuts=b2dmumu_selection).outputs()

        # Prepare simulation
        decay['mc_input']     = RootAppend(decay['mc_inputs'], 'SIM_B_Dmumu_ALL.root').outputs()
        decay['mc_reduced']   = Reduce(decay['mc_input'], variables_b2dmumu + mc_variables, treename='MC_B_Dmumu/default').outputs()
        decay['mc_triggered'], decay['mc_trigger_eff'], decay['mc_trigger_all_eff'] = Cut(decay['mc_reduced'], trigger_selection).outputs()
        decay['mc_withhypos'] = AddHypos(decay['mc_triggered']).outputs()
        decay['mc_selected'], decay['mc_selected_eff'], decay['mc_selected_all_eff'] = Cut(decay['mc_withhypos'], b2dmumu_selection).outputs()
        decay['mc_resampled'] = ResamplePID(b2dmumu['mc_selected']).outputs()

        weighted_b2jpsikst_mc, weighted_b2dmumu_mc = CalcSuperWeights(cut_b2jpsikst_mc, weighted_b2kstmumu, apply_to=[decay['mc_selected']]).outputs()

        simpleweighted_b2jpsikst_mc, simpleweighted_b2dmumu_mc = CalcSimpleWeights(resampled_b2jpsikst_mc, weighted_b2kstmumu, apply_to=[decay['mc_resampled']], variables=['nSPDHits', 'nTracks']).outputs()

        # Train and apply classifier
        from rep.estimators.xgboost import XGBoostClassifier
        clf = XGBoostClassifier(n_estimators=150, eta=0.3, gamma=12, max_depth=10, verbose=1, nthreads=20)

        #classified_b2dmumu_debug = KFoldCrossValidation(signal=selected_b2dmumu_mc, background=selected_b2dmumu, clf=clf).outputs()
        decay['classified'], decay['traintest'] = KFoldTrainAndApply(signal=simpleweighted_b2dmumu_mc,
                                                                     background=decay['selected'],
                                                                     extra=[weighted_b2kstmumu,
                                                                            cut_b2jpsikst_mc,
                                                                            simpleweighted_b2dmumu_mc,
                                                                            decay['selected'],
                                                                            simpleweighted_b2jpsikst_mc,
                                                                            weighted_b2jpsikst_mc,
                                                                            weighted_b2dmumu_mc,
                                                                           ],
                                                                     variables=classifier_variables,
                                                                     clf=clf).outputs()[:2]

        decay['model'] = LocalFile('models/Bd_D0mumu.model')

        bkg_only_fit_precut = RooFit(
                          decay['classified'],
                          decay['model'],
                          model_name='fullBkgMassPdf',
                          key=3,
                          censor={'B_M':(5229, 5329)},
                          integral=(5279-20, 5279+20),
                       ).outputs()

        plot_bkg_only_fit_precut = PlotFit(
                            decay['classified'],
                            bkg_only_fit_precut[1],
                            model_name='fullBkgMassPdf',
                            fix_norm=0.5,
                            xlabel='$m(K^+\\!\\pi^-\\!\\mu^+\\!\\mu^-)$',
                            path=DATASTORE + 'b2dmumu_bkg_only_fit_precut.pdf').outputs()

        bkg_yield_precut = bkg_only_fit_precut[3]
        decay['fom'], decay['mc_bdt_eff'] = CalculateOptimalMetric(1., bkg_yield_precut, decay['traintest']).outputs()

        decay['classified_cut'], decay['classified_cut_eff'], decay['classified_cut_eff_all'] = Cut(decay['classified'], ['clf > {}'], insert=[decay['fom']]).outputs()
        #decay['classified_cut'], decay['classified_cut_eff'] = Cut(decay['classified'], 'clf > 4.5').outputs()

        # Perform fits to get parameters for expected limit
        sig_only_fit = RooFit(
                decay['mc_selected'],
                decay['model'],
                model_name='sigMassPdf',
                rangeX={'B_M': (5225, 5335), 'Meson_M': (1840, 1890)},
                key=1,
            ).outputs()
        plot_sig_only_fit = PlotFit(
                decay['mc_selected'],
                sig_only_fit[1],
                model_name='sigMassPdf',
                components=['sigMassPdf1', 'sigMassPdf2'],
                path=DATASTORE + 'b2dmumu_sig_only_fit.pdf',
                range=(5225, 5335),
                xlabel='$m(K^+\\!\\pi^-\\!\\mu^+\\!\\mu^-)$',
                scale=0.48
            ).outputs()
        plot_sig_only_fit_d = PlotFit(
                decay['mc_selected'],
                sig_only_fit[1],
                plot_var='Meson_M',
                model_name='sigMassPdf',
                components=['sigMassPdfD1', 'sigMassPdfD2'],
                path=DATASTORE + 'b2dmumu_sig_only_fit_d.pdf',
                range=(1840, 1890),
                xlabel='$m(K^+\\!\\pi^-)$',
                scale=0.48
            ).outputs()

        bkg_only_fit = RooFit(
                          decay['classified_cut'],
                          decay['model'],
                          model_name='fullBkgMassPdf',
                          censor={'B_M': (5229, 5329)},
                          integral=(5279-50, 5279+50),
                          key=2,
                       ).outputs()

        plot_bkg_only_fit = PlotFit(
                                decay['classified_cut'],
                                bkg_only_fit[1],
                                model_name='fullBkgMassPdf',
                                path=DATASTORE + 'b2dmumu_bkg_only_fit.pdf',
                                binning=80,
                                components=['bkgMassPdf', 'lbgMassPdf'],
                                fix_norm=0.5,
                                log=False,
                                xlabel='$m(K^+\\!\\pi^-\\!\\mu^+\\!\\mu^-)$',
                            ).outputs()

        plot_bkg_only_fit_d = PlotFit(
                                decay['classified_cut'],
                                bkg_only_fit[1],
                                plot_var='Meson_M',
                                model_name='fullBkgMassPdf',
                                path=DATASTORE + 'b2dmumu_bkg_only_fit_d.pdf',
                                components=['bkgMassPdf', 'lbgMassPdf'],
                                binning=80,
                                fix_norm=0.25,
                                log=False,
                                xlabel='$m(K^+\\!\\pi^-)$',
                            ).outputs()

        decay['alpha'], decay['alpha_err'] = CalculateNormalization([
            # Values for signal channel
                # Generator level efficiency
                0.15791,
                # Stripping efficiency
                0.109938,
                # Trigger efficiency
                decay['mc_trigger_eff'],
                # Selection efficiency preselection
                decay['mc_selected_eff'],
                # Selection efficiency classifier
                decay['mc_bdt_eff'],
            ],
            # Total number of generated events
            507614 + 510169,
            # Values for control channel
            [
                # Generator level efficiency
                0.160500,
                # Stripping efficiency
                0.096866,
                # Trigger efficiency
                b2jpsikst_trigger_eff,
                # Selection efficiency
                b2jpsikst_cut_eff,
            ],
            # Total number of generated events
            4435958 + 4425822,
            # Branching fraction from PDG
            (1.32e-3, 0.06e-3),
            fit_b2kstmumu[2],
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

                    'mu_alpha':       (decay['alpha'], 0),
                    'sigma_alpha':    (decay['alpha_err'], 0),
                },
                set_params={
                    'bkgYield':      (bkg_only_yield, 1),
                    'alpha':          (decay['alpha'], 1)
                },
                ).outputs()

    def show(val):
        require(val)
        def rec(val):
            if isinstance(val, PyTarget):
                rec(val.get())
            if isinstance(val, list):
                for x in val:
                    rec(x)
            else:
                print(val)
        rec(val)

    #require([b2dmumu['alpha']])
    #require([weighted_b2kstmumu,plot_bkg_only_fit, plot_bkg_only_fit_d, plot_sig_only_fit, plot_sig_only_fit_d])
    #require([b2dmumu['reduced']])
    #require([plot_b2kstmumu, plot_b2jpsikst_mc])
    #require([plot_bkg_only_fit, plot_bkg_only_fit_d])
    #require(b2dmumu['classified'])
    #require([plot_sig_only_fit, plot_sig_only_fit_d])
    #require([plot_bkg_only_fit_precut])
    show(b2dmumu['mc_bdt_eff'])
    #show(b2dmumu['mc_selected_all_eff'])
    #show([b2dmumu['mc_selected_all_eff'], b2dmumu['mc_trigger_all_eff']])
    #require(decay['classified'])
    #require(decay['expected'])
    #require(simpleweighted_b2dmumu_mc)
