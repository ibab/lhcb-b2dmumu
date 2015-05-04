

from Gaudi.Configuration import *
from Configurables import  DaVinci
from DiLeptonTuple.DiLeptonTuple import addTuple

def create_options(decay, outname, year, name="B2XMuMu_Line", mc=False):

    if mc:
        head = "/Event/AllStreams/Phys"
        addendum = "MC"
    else:
        head = "/Event/Dimuon/Phys"
        addendum = "DST"

    DaVinci().appendToMainSequence(addTuple(
        name=name,
        head=head,
        decay=decay,
        dtf=True,
        resonant=False,
        addendum=addendum,
        verbose=["Pid", "iso", "Kinematic"]
    ))

    DaVinci().EvtMax = -1
    DaVinci().TupleFile = outname
    DaVinci().Lumi = True
    DaVinci().DataType = str(year)

    if year == 2012:
        DaVinci().CondDBtag = "cond-20120831"
        DaVinci().DDDBtag = "dddb-20120831"
    elif year == 2011:
        DaVinci().CondDBtag = "cond-20130114"
        DaVinci().DDDBtag = "dddb-20130111"
    else:
        raise ValueError("unknown year: {}".format(year))

