
from Gaudi.Configuration import *

from Configurables import PrintHeader, PrintDecayTree, DaVinci, StoreExplorerAlg, FilterDesktop, DecodeVeloRawBuffer, DeterministicPrescaler, CheckSelResult, CombineParticles

DaVinci().appendToMainSequence([DecodeVeloRawBuffer("DecodeVeloClusters")])

from DiLeptonTuple.DiLeptonTuple import addTuple, Bd2MuMuD0

DaVinci().appendToMainSequence(addTuple(
    name="B2XMuMu_Line",
    head="/Event/AllStreams/Phys",
    decay=Bd2MuMuD0,
    dtf=False,
    resonant=False,
    addendum="MC",
    verbose=["Pid", "iso"]
))

DaVinci().EvtMax = -1
DaVinci().TupleFile = "Bd2MuMuD0.root"

