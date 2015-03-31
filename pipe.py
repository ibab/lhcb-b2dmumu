
import luigi
from luigi import Task, LocalTarget, run

import os

from analysis.log import setup_logging
logger = setup_logging()
from analysis.util import *
from analysis.reduce import *

try:
    DATASTORE = os.environ['DATASTORE']
except KeyError:
    DATASTORE='/fhgfs/users/ibabuschkin/DataStore'

BLINDED = True

bdt_variables = [
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
        'D~0_CosTheta',
]

class DataFiles(Task):
    def output(self):
        path = os.path.join(DATASTORE, 'Data/AllYears/Stripping20/Dimuon/DVBd2MuMuD0_data/combined/Bd2MuMuD0.root')
        return RootTarget(path, 'B2XMuMu_Line_TupleDST/DecayTree')

class MCFiles(Task):
    def output(self):
        path = os.path.join(DATASTORE, 'MC/2012/Stripping20/AllStreams/DVBd2MuMuD0_MC/combined/Bd2MuMuD0.root')
        return RootTarget(path, 'B2XMuMu_Line_TupleMC/DecayTree')

class Analysis(Task):
    def requires(self):
        reduce_data = Reduce(DataFiles(), blinded=BLINDED)
        reduce_mc = Reduce(MCFiles(), extra_vars=mc_variables)
        return [reduce_data, reduce_mc]

if __name__ == '__main__':
    luigi.run(main_task_cls=Analysis, local_scheduler=True)
