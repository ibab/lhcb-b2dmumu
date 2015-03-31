
from analysis.log import get_logger
logger = get_logger()
import luigi
from luigi import Task

# Variables that are used in the analysis
variables = [
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

# These are MC-only variables
mc_variables = [
        'B_BKGCAT',
        '*_TRUEID',
]

class Reduce(Task):
    """
    Read in dataset and apply preparation:
        - Select only certain variables
        - Apply blinding by removing signal region from data
        - Calculate additional variables
    """
    in_task = luigi.Parameter()
    extra_vars = luigi.Parameter(default=[])
    blinded = luigi.BoolParameter(default=False)

    def requires(self):
        return self.in_task

    def run(self):
        logger.info('Running {}'.format(self))

        import numpy as np
        from pandas import Series
        from root_pandas import read_root

        if self.blinded:
            B_mass = 5279
            width = 50
            blindcut = '(B_M < {}) || (B_M > {})'.format(B_mass - width, B_mass + width)
        else:
            blindcut = None

        df = read_root(self.input().path, self.input().tree, columns=variables + list(self.extra_vars), where=blindcut)
        logger.info('Initial events: {}'.format(len(df)))

        # Calculate some extra variables
        from scipy.constants import c
        df['B_TAU'] = Series(df['B_FD_OWNPV'] * df['B_M'] / (df['B_P'] * c * 10**3) * 10**12, index=df.index)
        df['B_DiraAngle'] = Series(np.arccos(df['B_DIRA_OWNPV']), index=df.index)
        df['B_ENDVERTEX_CHI2_NDOF'] = Series(df['B_ENDVERTEX_CHI2'] / df['B_ENDVERTEX_NDOF'], index=df.index)
        for var in df.columns:
            if 'PZ' in var and var.replace('PZ', 'P') in df.columns:
                df[var.replace('PZ', 'ETA')] = np.arctanh(df[var] / df[var.replace('PZ', 'P')])

        df.to_root(self.output().path)

    def output(self):
        return self.input().add_step(self.__class__.__name__)

