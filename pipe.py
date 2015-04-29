
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

DATASTORE='/fhgfs/users/ibabuschkin/DataStore/tmp/'

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
            width = 50
            cut.add('(B_M < {}) || (B_M > {})'.format(B_mass - width, B_mass + width))

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
            'Kplus': '/fhgfs/groups/e5/lhcb/PIDCalib/Kaon_Stripping20_MagnetUp.pkl',
            'piminus': '/fhgfs/groups/e5/lhcb/PIDCalib/Pi_Stripping20_MagnetUp.pkl',
            'muplus': '/fhgfs/groups/e5/lhcb/PIDCalib/Mu_Stripping20_MagnetUp.pkl',
            'muminus': '/fhgfs/groups/e5/lhcb/PIDCalib/Mu_Stripping20_MagnetUp.pkl',
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
                df[key] = resampler[part + '_' + nametrans_pid[pid]].sample(df[[particle + '_P', particle + '_ETA', 'nTracks']].values.T)
                logger.info('Resampled {} for {}'.format(pid, particle))
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
        'Kplus_PIDK',
        'piminus_PIDK',
        'muplus_PIDmu',
        'muminus_PIDmu',
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
          'xgboost':  XGBoostClassifier(n_estimators=200, max_depth=12, verbose=1),
        }

        clf = classifiers[self.name]

        select_sidebands = [
                '(B_M > 5300)',
        ]

        step = 1

        bkg = read_root(self.background.path(), columns=classifier_variables, where=prepare_sel(select_sidebands), step=step).dropna()
        sig = read_root(self.signal.path(), columns=classifier_variables, step=step).dropna()

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

class RooFit(Task):
    infile = Input()
    model = Input()
    plot_var = Input(default='B_M')
    key = Input(default=0)

    def outputs(self):
        return [LocalFile(DATASTORE + 'results_{}.params'.format(self.key)), PyTarget('workspace_{}'.format(self.key))]

    def run(self):
        out_params, out_ws = self.outputs()
        out_params = out_params.path()

        import ROOT
        from analysis.fit import mle, assemble_model, load_tree

        ws = assemble_model(self.model.path())
        model = ws.pdf('model')
        data = load_tree(ws, self.infile.path(), 'default', '')

        ROOT.SetOwnership(ws, False)
        mle(model, data, start_params="models/Bd_KstJpsi.params", out_params=out_params, numcpus=20)
        out_ws.set(ws)

class PlotFit(Task):
    infile = Input()
    inws = Input()
    plot_var = Input(default='B_M')

    def outputs(self):
        return LocalFile(DATASTORE + 'plot.pdf')

    def run(self):

        from analysis.plotting import plot_roofit
        from analysis.fit import load_tree
        import matplotlib.pyplot as plt
        ws = self.inws.get()
        model = ws.pdf('model')
        data = load_tree(ws, self.infile.path(), 'default', '')
        v = ws.var(self.plot_var)
        ax, width = plot_roofit(v, data, model, components=['sigGaussian1', 'sigGaussian2', 'sigGaussian3', 'bkgMassPdf'], numcpus=20)
        plt.savefig(self.outputs().path())
        plt.clf()

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

class CalcROCFromWeights(Task):

    infile = Input()
    inweights = Input()

    def outputs(self):
        return LocalFile(DATASTORE + 'weighted_roc.pdf'), LocalFile(DATASTORE + 'sig_eff.pdf'), LocalFile(DATASTORE + 'bkg_eff.pdf')
    
    def run(self):
        import pandas as pd
        df = read_root(self.infile.path())
        weights = read_root(self.inweights.path())

        assert len(df) == len(weights), 'length of data and weights must match'

        logger.warn('df: {}, weights: {}'.format(len(df), len(weights)))

        df['sigweights'] = pd.Series(weights['sigYield_sw'].ravel(), index=df.index)
        df['bkgweights'] = pd.Series(weights['bkgYield_sw'].ravel(), index=df.index)

        total_sig = df.sigweights.sum()
        total_bkg = df.bkgweights.sum()

        logger.warn('total_sig: {}'.format(total_sig))
        logger.warn('total_bkg: {}'.format(total_bkg))

        x = []
        y = []
        cuts = np.linspace(-15, 10, 200)

        for c in cuts:
            x.append(df['bkgweights'][df.classifier > c].sum())
            y.append(df['sigweights'][df.classifier > c].sum())

        import matplotlib.pyplot as plt
        plt.plot(1 - x / total_bkg, y / total_sig, 'b-')
        plt.savefig(self.outputs()[0].path())
        plt.clf()
        plt.plot(cuts, y, 'b-')
        plt.savefig(self.outputs()[1].path())
        plt.clf()
        plt.plot(cuts, x, 'b-')
        plt.savefig(self.outputs()[2].path())
        plt.clf()

class MergeROOT(Task):
    inputfiles = Input()
    outputname = Input()

    def outputs(self):
        return LocalFile(DATASTORE + self.outputname)

    def run(self):
        from sh import hadd
        print(hadd('-f', self.outputs().path(), *list(map(lambda x: x.path(), self.inputfiles)), _out=sys.stdout, _err=sys.stderr))

if __name__ == '__main__':
    # B0->D~0mumu
    input_b2dmumu        = LocalFile('/fhgfs/users/ibabuschkin/DataStore/Data/AllYears/Stripping20/Dimuon/DVBd2MuMuD0_data/combined/DATA_Bd2D0mumu.root')
    reduced_b2dmumu      = Reduce(input_b2dmumu, variables_b2dmumu, treename='B2XMuMu_Line_TupleDST/DecayTree', blinded=True).outputs()
    triggered_b2dmumu    = ApplyTrigger(reduced_b2dmumu).outputs()
    selected_b2dmumu     = Select(triggered_b2dmumu).outputs()

    input_b2dmumu_mc     = LocalFile('/fhgfs/users/ibabuschkin/DataStore/MC/2012/Stripping20/AllStreams/DVBd2MuMuD0_MC/combined/SIM_Bd2D0mumu.root')
    reduced_b2dmumu_mc   = Reduce(input_b2dmumu_mc, variables_b2dmumu + mc_variables, treename='B2XMuMu_Line_TupleMC/DecayTree').outputs()
    resampled_b2dmumu_mc = ResamplePID(reduced_b2dmumu_mc).outputs()
    triggered_b2dmumu_mc = ApplyTrigger(resampled_b2dmumu_mc).outputs()
    selected_b2dmumu_mc  = Select(triggered_b2dmumu_mc).outputs()

    clf = TrainClassifier(signal=selected_b2dmumu_mc, background=selected_b2dmumu).outputs()

    # B->K*mumu
    input_b2kstmumu = LocalFile('/fhgfs/users/ibabuschkin/DataStore/tmp/DATA_Kstmumu.root')

    reduced_b2kstmumu = Reduce(input_b2kstmumu, variables_b2kstmumu).outputs()
    triggered_b2kstmumu = ApplyTrigger(reduced_b2kstmumu).outputs()
    cut_b2kstmumu = ApplyCut(triggered_b2kstmumu, ['B_M > 5180', 'B_M < 5380', 'Kstar_M > 896 - 150', 'Kstar_M < 896 + 150', 'Psi_M > 3000', 'Psi_M < 3200']).outputs()
    classified_b2kstmumu = ApplyClassifier(cut_b2kstmumu, clf).outputs()
    model = LocalFile('models/Bd_KstJpsi.model')
    fit_b2kstmumu = RooFit(classified_b2kstmumu, model).outputs()
    plot_b2kstmumu = PlotFit(cut_b2kstmumu, fit_b2kstmumu[1]).outputs()
    weighted_b2kstmumu = CalcSWeights(cut_b2kstmumu, fit_b2kstmumu[1]).outputs()
    roc_plot  = CalcROCFromWeights(classified_b2kstmumu, weighted_b2kstmumu).outputs()

    require(roc_plot)
    #require(roc_plot)
    #ret = []
    #for i, c in enumerate(np.linspace(0, 5, 10)):
    #    cut = ApplyCut(classified_b2kstmumu, cuts=['classifier > {}'.format(c)], key=i).outputs()
    #    fit = RooFit(cut, model, key=i).outputs()
    #    ret.append(fit)
    #require(ret)

