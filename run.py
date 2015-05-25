
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

    def outputs(self):
        name = os.path.basename(self.infile.path())
        outfile = LocalFile(DATASTORE + name.replace('.root', '.' + self.__class__.__name__ + '.root'))
        return outfile

    def run(self):
        from analysis.util import prepare_sel

        selection = [
            # Exclude J/psi
            'Psi_M < 2860 || Psi_M > 3200',
            # Kinematic range ends below this
            'Psi_M < 3500',
        ]

        # Filter particle ID variables for true Kaon and true Pion
        pid_cuts = [
            'Kplus_PIDK > 0',
            'piminus_PIDK < 0',
            'Kplus_isMuonLoose == 0',
            'piminus_isMuonLoose == 0',
        ]

        df = read_root(self.infile.path(), where=prepare_sel(selection))
        df.to_root(self.outputs().path())

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
]

class ApplyCut(Task):
    infile = Input()
    cuts = Input()
    key = Input(default='')

    def outputs(self):
        if self.key is '':
            keystr = ''
        else:
            keystr = '_{}'.format(self.key)
        return LocalFile(self.infile.path().replace('.root', '.ApplyCut{}.root'.format(keystr)))

    def run(self):
        from analysis.util import prepare_sel
        df = read_root(self.infile.path(), where=prepare_sel(self.cuts))
        df.to_root(self.outputs().path())

class TrainClassifier(Task):
    signal = Input()
    background = Input()
    name = Input(default='xgboost')
    folds = Input(default=5)

    def outputs(self):
        return PyTarget('clf')

    def run(self):

        from analysis.util import prepare_sel
        from sklearn.tree import DecisionTreeClassifier
        from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
        from rep.classifiers.xgboost import XGBoostClassifier
        from rep.classifiers.tmva import TMVAClassifier

        classifiers = {
          'xgboost':  XGBoostClassifier(n_estimators=300, gamma=2, max_depth=7, verbose=1, nthreads=4),
        }

        clf = classifiers[self.name]

        select_sidebands = [
                '(B_M > 5300)',
        ]

        step = 1
        bkg = read_root(self.background.path(), columns=classifier_variables, where=prepare_sel(select_sidebands), step=4*step).dropna()
        sig = read_root(self.signal.path(), columns=classifier_variables, step=step).dropna()

        logger.info('Signal events: {}'.format(len(sig)))
        logger.info('Backg. events: {}'.format(len(bkg)))

        X = np.vstack([bkg.values, sig.values])
        y = np.append(np.zeros(len(bkg)), np.ones(len(sig)))

        clf.fit(X, y)

        self.outputs().set(clf)

class ApplyClassifier(Task):
    infile = Input()
    clf = Input()
    name = Input(default='classifier')

    def outputs(self):
        return LocalFile(self.infile.path().replace('.root', '.ApplyClassifier_{}.root'.format(self.name)))
    
    def run(self):
        clf = self.clf.get()
        data = read_root(self.infile.path(), columns=classifier_variables)
        def logit(x):
            return np.log(x / (1 - x))

        pred = logit(clf.predict_proba(data.values)[:,1])
        data = read_root(self.infile.path())
        data[self.name] = pred
        data.to_root(self.outputs().path())

class KFoldCrossValidation(Task):
    signal = Input()
    background = Input()
    clf = Input()

    def outputs(self):
        return LocalFile(self.signal.path().replace('.root', '.KFoldCrossValidation.root'))
    
    def run(self):
        from analysis.classification import evaluate_classifier

        clf = self.clf
        step = 1

        sig = read_root(self.signal.path(), columns=classifier_variables, step=step).dropna()
        bkg = read_root(self.background.path(), columns=classifier_variables, step=step * 5).dropna()
        data = pd.concat([sig, bkg], keys=['sig', 'bkg'])

        X = data.values.astype('float32')
        y = np.append(np.ones(len(sig)), np.zeros(len(bkg)))

        from sklearn.cross_validation import StratifiedKFold
        from sklearn.metrics import roc_auc_score
        from sklearn.base import clone
        from scipy.special import logit

        classifiers = []
        data['signal'] = y
        data['fold'] = 0
        data['decision'] = 0.
        data['proba'] = 0.

        skf = StratifiedKFold(y, n_folds=5, shuffle=True, random_state=0)
        for i, (train, test) in enumerate(skf):
            classifier = clone(clf)
            logger.info('Training fold {}'.format(i))
            classifier.fit(X[train], y[train])
            logger.info('Predicting fold {}'.format(i))
            proba = classifier.predict_proba(X[test])[:,1]
            data['fold'].values[test] = i + 1
            data['proba'].values[test] = proba

            data['decision'].values[test] = logit(proba)
            classifiers.append(classifier)

        score = roc_auc_score(data['signal'], data['decision'])
        logger.info('ROC AUC: {}'.format(score))
        data.to_root(self.outputs().path())

class RooFit(Task):
    infile = Input()
    model = Input()
    model_name = Input(default='model')
    params = Input(default='')
    key = Input(default=0)
    fix_params = Input(default='')
    censor = Input(default='')

    def outputs(self):
        return [LocalFile(DATASTORE + 'results_{}.params'.format(self.key)),
                PyTarget('workspace_{}'.format(self.key)),
                PyTarget('fitresults_{}'.format(self.key))]

    def run(self):
        out_params, out_ws, out_results = self.outputs()
        out_params = out_params.path()

        import ROOT
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

            rng = ROOT.RooFit.Range(','.join(ranges))
            extra_params.append(rng)

        results = mle(model, data, out_params=out_params, numcpus=20, extra_params=extra_params)
        results.Print()
        out_ws.set(ws)
        out_results.set(results)

class GenToy(Task):
    model = Input()
    variables = Input()
    model_name = Input(default='model')
    outpath = Input(default='toy.root')
    fix_params = Input(default='')

    def outputs(self):
        return LocalFile(DATASTORE + self.outpath)
    
    def run(self):
        import ROOT
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

        vrs = [ws.var(x) for x in self.variables]
        print(vrs)
        args = ROOT.RooArgSet(*vrs)
        ff = ROOT.TFile(self.outputs().path(), 'recreate')
        ROOT.RooAbsData.setDefaultStorageType(ROOT.RooAbsData.Tree)
        data = model.generate(args)
        data.tree().Write('tree')
        ff.Close()

class PlotFit(Task):
    infile = Input()
    inws = Input()
    path = Input()
    model_name = Input(default='model')
    plot_var = Input(default='B_M')
    components = Input(default=[])

    def outputs(self):
        return LocalFile(self.path)

    def run(self):
        from analysis.plotting import plot_roofit
        from analysis.fit import load_tree
        import matplotlib.pyplot as plt
        ws = self.inws.get()
        model = ws.pdf(self.model_name)
        data = load_tree(ws, self.infile.path(), 'default', '')
        v = ws.var(self.plot_var)
        plt.figure(figsize=(12, 8))
        gs, ax, width = plot_roofit(v, data, model, components=self.components, numcpus=20, xlabel='$m(K^+\\!\\pi^-\\!\\mu^+\\!\\mu^-)$')

        plt.ylabel('Candidates', ha='right', y=1)
        gs.tight_layout(plt.gcf())
        plt.savefig(self.outputs().path())
        plt.clf()

        import ROOT
        c1 = ROOT.TCanvas()
        frame = v.frame()
        data.plotOn(frame)
        model.plotOn(frame)
        frame.Draw()
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

if __name__ == '__main__':

    inputs_b2dmumu = [
            LocalFile('./store/DATA_Bd_D0mumu_MU11.root'),
            LocalFile('./store/DATA_Bd_D0mumu_MD11.root'),
            LocalFile('./store/DATA_Bd_D0mumu_MU12.root'),
            LocalFile('./store/DATA_Bd_D0mumu_MD12.root'),
    ]

    inputs_b2dmumu_mc = [
            LocalFile('./store/SIM_Bd_D0mumu_MD12.root'),
            LocalFile('./store/SIM_Bd_D0mumu_MU12.root'),
    ]

    # B0 -> D~0 mu mu
    input_b2dmumu        = RootAppend(inputs_b2dmumu, 'DATA_B2D0mumu_ALL.root').outputs()
    reduced_b2dmumu      = Reduce(input_b2dmumu, variables_b2dmumu, treename='B2XMuMu_Line_TupleDST/DecayTree', blinded=True).outputs()
    triggered_b2dmumu    = ApplyTrigger(reduced_b2dmumu).outputs()
    selected_b2dmumu     = Select(triggered_b2dmumu).outputs()

    input_b2dmumu_mc     = RootAppend(inputs_b2dmumu_mc, 'SIM_Bd_D0mumu_ALL.root').outputs()
    reduced_b2dmumu_mc   = Reduce(input_b2dmumu_mc, variables_b2dmumu + mc_variables, treename='B2XMuMu_Line_TupleMC/DecayTree').outputs()
    resampled_b2dmumu_mc = ResamplePID(reduced_b2dmumu_mc).outputs()
    triggered_b2dmumu_mc = ApplyTrigger(resampled_b2dmumu_mc).outputs()
    selected_b2dmumu_mc  = Select(triggered_b2dmumu_mc).outputs()

    TrainClassifier(signal=selected_b2dmumu_mc, background=selected_b2dmumu).outputs()

    from rep.classifiers.xgboost import XGBoostClassifier
    clf = XGBoostClassifier(n_estimators=150, gamma=12, max_depth=10, verbose=1, nthreads=4)
    classified_b2dmumu = KFoldCrossValidation(signal=selected_b2dmumu_mc, background=selected_b2dmumu, clf=clf).outputs()

    model_b2dmumu = LocalFile('models/Bd_D0mumu.model')

    sig_only_fit = RooFit(
            selected_b2dmumu_mc,
            model_b2dmumu,
            model_name='sigMassPdf',
        ).outputs()
    plot_sig_only_fit = PlotFit(
            selected_b2dmumu_mc,
            sig_only_fit[1],
            model_name='sigMassPdf',
            path=DATASTORE + 'b2dmumu_sig_only_fit.pdf',
        ).outputs()

    bkg_only_fit = RooFit(
                      selected_b2dmumu,
                      model_b2dmumu,
                      model_name='bkgMassPdf',
                      censor={'B_M': (5279 - 50, 5279 + 50)},
                   ).outputs()
    plot_bkg_only_fit = PlotFit(
                            selected_b2dmumu,
                            bkg_only_fit[1],
                            model_name='bkgMassPdf',
                            path=DATASTORE + 'b2dmumu_bkg_only_fit.pdf'
                        ).outputs()
    
    sig_only_fitresults = sig_only_fit[2]
    bkg_only_fitresults = bkg_only_fit[2]

    #toys = []
    #for sig_yield in np.linspace(0, 100, 10):
    #    toy = GenToy(
    #            model_b2dmumu,
    #            variables=['B_M'],
    #            fix_params={
    #                'frac':          sig_only_fitresults,
    #                'sigMassMean':   sig_only_fitresults,
    #                'sigMassSigma1': sig_only_fitresults,
    #                'sigMassSigma2': sig_only_fitresults,
    #                'bkgMassSlope':  bkg_only_fitresults,
    #                'sigYield':      sig_yield,
    #                'bkgYield':      100,
    #            }).outputs()

    #    toys.append(toy)

    # B0 -> K* mu mu
    inputs_b2kstmumu = [
            LocalFile('./store/DATA_Bd_Kst0mumu_MD11.root'),
            LocalFile('./store/DATA_Bd_Kst0mumu_MU11.root'),
            LocalFile('./store/DATA_Bd_Kst0mumu_MD12.root'),
            LocalFile('./store/DATA_Bd_Kst0mumu_MU12.root'),
    ]
    input_b2kstmumu       = RootAppend(inputs_b2kstmumu, 'DATA_B2Kstmumu_ALL.root').outputs()
    model_b2kstmumu       = LocalFile('models/Bd_KstJpsi_CBall.model')
    init_params_b2kstmumu = LocalFile('models/Bd_KstJpsi_CBall.params')
    control_channel       = LocalFile('control-channel.ipynb')

    reduced_b2kstmumu    = Reduce(input_b2kstmumu, variables_b2kstmumu).outputs()
    triggered_b2kstmumu  = ApplyTrigger(reduced_b2kstmumu).outputs()
    cut_b2kstmumu        = ApplyCut(triggered_b2kstmumu, ['B_M > 5100', 'B_M < 5500', 'Kstar_M > 896 - 150', 'Kstar_M < 896 + 150', 'Psi_M > 3000', 'Psi_M < 3200', 'Kplus_PIDK > -5']).outputs()
    classified_b2kstmumu = ApplyClassifier(cut_b2kstmumu, clf).outputs()
    fit_b2kstmumu        = RooFit(classified_b2kstmumu, model_b2kstmumu, params=init_params_b2kstmumu).outputs()
    plot_b2kstmumu       = PlotFit(cut_b2kstmumu,
                                   fit_b2kstmumu[1],
                                   path=DATASTORE + 'b2kstmumu_data_fit.pdf',
                                   components=['sigMassPdf1', 'sigMassPdf2', 'bkgMassPdf']).outputs()
    weighted_b2kstmumu   = CalcSWeights(cut_b2kstmumu, fit_b2kstmumu[1]).outputs()
    control_channel      = RunNotebook(control_channel, [weighted_b2kstmumu]).outputs()

    require([plot_bkg_only_fit])

