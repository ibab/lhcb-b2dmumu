#
#   Stripping selections test job
#
#   @date 2009-09-29
#

from Gaudi.Configuration import *

MessageSvc().Format = "% F%80W%S%7W%R%T %0W%M"
from Configurables import PrintHeader, PrintDecayTree, DaVinci, StoreExplorerAlg, FilterDesktop, DecodeVeloRawBuffer, DeterministicPrescaler, CheckSelResult, CombineParticles

"""
# mumuK*
seq1 = GaudiSequencer("Bd2MuMuKstarSeqMain")
VX =  "(((VFASPF ( VCHI2/VDOF) ) < 9) & (BPVVDZ > 0))"
JVX = "(INTREE( ('J/psi(1S)'==ABSID) & "+VX+"))"
KVX = "(INTREE( ('K*(892)0'==ABSID) & "+VX+"))"
IPS = "(NINTREE( (ISBASIC) & (MIPCHI2DV(PRIMARY)>9))==4)"
TRC = "(NINTREE( (ISBASIC) & (TRCHI2DOF < 5))==4)"
PT  = "(INTREE( ('K+'==ABSID) & (PT>300))) & (INTREE( ('pi+'==ABSID) & (PT>300)))"
Code = VX+" & "+JVX+" & "+KVX+" & "+IPS+" & "+TRC+" & "+PT
seq1.Members = [
     DeterministicPrescaler("OfficialBd2MuMuKstarPrescaler", AcceptFraction = 0.5 ) # NEVER CHANGE ANYTHING IN THIS LINE
   , FilterDesktop("RefineSignalMain", Inputs = [ '/Event/Dimuon/Phys/Bd2KstarMuMuLT' ], Code = Code )
     ]

seq2 = GaudiSequencer("Bd2MuMuKstarSeq2")
seq2.Members = [
    CheckSelResult(Notmode=True, Algorithms=["OfficialBd2MuMuKstarPrescaler"])
   , FilterDesktop("RefineSignalHidden", Inputs = [ '/Event/Dimuon/Phys/Bd2KstarMuMuLT' ], Code = Code )
]
"""

# velo clusters
DaVinci().appendToMainSequence( [ #PrintHeader(),
                                  # PrintDecayTree(Inputs = ["/Event/Dimuon/Phys/Bd2KstarMuMu_Signal"]),
                                  DecodeVeloRawBuffer("DecodeVeloClusters") ] ) # , StoreExplorerAlg()


# ADDTUPLE
from DiLeptonTuple.DiLeptonTuple import addTuple, Bd2MuMuD0
DaVinci().appendToMainSequence( addTuple(name="Bd2D0MuMu_Signal", head="/Event/Dimuon/Phys/", decay=Bd2MuMuD0, dtf=True, resonant=False) )
DaVinci().appendToMainSequence( addTuple(name="Bd2D0MuMuLT", head="/Event/Dimuon/Phys/", decay=Bd2MuMuD0, dtf=True, resonant=False) )
DaVinci().appendToMainSequence( addTuple(name="Bd2D0MuMu_BdToKstarMuMuLine", head="/Event/Dimuon/Phys/", decay=Bd2MuMuD0, dtf=True, resonant=False) )
DaVinci().appendToMainSequence( addTuple(name="StrippingB2D0MuMu_B2MuMuXLine", head="/Event/Dimuon/Phys/", decay=Bd2MuMuD0, dtf=True, resonant=False) )

DaVinci().EvtMax = 10000
DaVinci().TupleFile = "Bd2MuMuD0.root"               
DaVinci().DataType = "2012"
DaVinci().Simulation = True
    
DaVinci().Lumi = True
# DaVinci().CondDBtag = "head-20111111"
# DaVinci().DDDBtag = "head-20111102"

"""
TTree* T = _file0->Get("Bd2KstarMuMu_BdToKstarMuMuLine_Tuple/DecayTree")
T->Draw("muminus_CosTheta")
T->Draw("muplus_CosTheta")
T->Draw("Kst_892_0_CosTheta")
T->Draw("Kplus_CosTheta")
T->Draw("piminus_CosTheta")
T->Draw("B_ThetaL")
T->Draw("B_ThetaK") 
T->Draw("B_Phi") 
T->Draw("B_ThetaTr")
T->Draw("B_PhiTr") 
T->Draw("B_ThetaVtr")

T->Draw("acos(muminus_CosTheta):B_ThetaL")
T->Draw("acos(muminus_CosTheta):B_ThetaTr")
T->Draw("acos(muminus_CosTheta):B_PhiTr")
T->Draw("acos(muminus_CosTheta):B_ThetaTr")
T->Draw("acos(Kst_892_0_CosTheta):B_ThetaL")
T->Draw("acos(Kst_892_0_CosTheta):B_ThetaTr")
T->Draw("acos(Kst_892_0_CosTheta):B_PhiTr")
T->Draw("acos(Kst_892_0_CosTheta):B_ThetaTr")
T->Draw("acos(Kst_892_0_CosTheta):B_Phi")
T->Draw("acos(Kst_892_0_CosTheta):B_ThetaK")
T->Draw("acos(Kplus_CosTheta):B_ThetaL")
T->Draw("acos(Kplus_CosTheta):B_ThetaTr")
T->Draw("acos(Kplus_CosTheta):B_PhiTr")
T->Draw("acos(Kplus_CosTheta):B_ThetaTr")
T->Draw("acos(Kplus_CosTheta):B_Phi")
T->Draw("acos(Kplus_CosTheta):B_ThetaK")
"""
