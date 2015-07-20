

from Gaudi.Configuration import *
from Configurables import  DaVinci, DecayTreeTuple, TupleToolTISTOS
from DiLeptonTuple.AllTriggerLines import muonLines

def create_tuples(mc):

    if mc:
        kind='MC'
        path = "/Event/AllStreams/Phys/B2XMuMu_Line/Particles"
    else:
        kind='DATA'
        path = "/Event/Dimuon/Phys/B2XMuMu_Line/Particles"

    mesons = [
        {
            'particle': 'D~0',
            'decay': 'B_Dmumu',
        },
        {
            'particle': 'K*(892)',
            'decay': 'B_JpsiKst',
        }]

    sequence = GaudiSequencer('DataSequence')

    for meson in mesons:

        data = DecayTreeTuple('{}_{}'.format(kind, meson['decay']))
        data.TupleName = 'default'
        
        geometry = data.addTupleTool('TupleToolGeometry')
        geometry.FillMultiPV = True

        triggerLines = muonLines()
        tistos = TupleToolTISTOS(TriggerList=triggerLines,
                                 VerboseHlt1=True,
                                 VerboseHlt2=True,
                                 VerboseL0=True)
        
        data.Branches = {
                'B' :       '[B0 ->  (J/psi(1S) ->  mu+  mu-)  ({} ->  K+  pi-) ]CC'.format(meson['particle']),
                'mumu' :    '[B0 -> ^(J/psi(1S) ->  mu+  mu-)  ({} ->  K+  pi-) ]CC'.format(meson['particle']),
                'AntiD0':   '[B0 ->  (J/psi(1S) ->  mu+  mu-) ^({} ->  K+  pi-) ]CC'.format(meson['particle']),
                'Kplus':    '[B0 ->  (J/psi(1S) ->  mu+  mu-)  ({} -> ^K+  pi-) ]CC'.format(meson['particle']),
                'piminus':  '[B0 ->  (J/psi(1S) ->  mu+  mu-)  ({} ->  K+ ^pi-) ]CC'.format(meson['particle']),
                'muplus':   '[B0 ->  (J/psi(1S) -> ^mu+  mu-)  ({} ->  K+  pi-) ]CC'.format(meson['particle']),
                'muminus':  '[B0 ->  (J/psi(1S) ->  mu+ ^mu-)  ({} ->  K+  pi-) ]CC'.format(meson['particle']),
        }

        mumu = data.addTupleTool('TupleToolDecay/mumu')
        mumu.ToolList += ['TupleToolTISTOS']

        B = data.addTupleTool('TupleToolDecay/B')

        p2vv = B.addTupleTool('TupleToolP2VV')
        p2vv.Calculator = 'Bd2KstarMuMuAngleCalculator'

        angles = data.addTupleTool('TupleToolAngles')
        angles.WRTMother = False

        constrainJpsi = B.addTupleTool('TupleToolDecayTreeFitter/ConstrainJpsi')
        constrainJpsi.Verbose = True
        constrainJpsi.constrainToOriginVertex = True
        constrainJpsi.daughtersToConstrain + ['J/psi(1S)']

        constrainMeson = B.addTupleTool('TupleToolDecayTreeFitter/ConstrainMeson')
        constrainMeson.Verbose = True
        constrainJpsi.constrainToOriginVertex = True
        constrainJpsi.daughtersToConstrain + ['D~0']

        recoStats = data.addTupleTool('TupleToolRecoStats')
        eventInfo = data.addTupleTool('TupleToolEventInfo')
        kinematic = data.addTupleTool('TupleToolKinematic')
        primares = data.addTupleTool('TupleToolPrimaries')
        trackInfo = data.addTupleTool('TupleToolTrackInfo')
        pid = data.addTupleTool('TupleToolPid')
        pid.Verbose = True
        data.addTupleTool('TupleToolL0Calo')

        from configIso import configIso
        configIso()

        from Configurables import TupleToolApplyIsolation
        data.B.addTool(TupleToolApplyIsolation, name="TupleToolApplyIsolationHard")
        data.B.TupleToolApplyIsolationHard.OutputSuffix="_Hard"
        data.B.TupleToolApplyIsolationHard.WeightsFile="weightsHard.xml"
        data.B.ToolList+=["TupleToolApplyIsolation/TupleToolApplyIsolationHard"]
        data.B.addTool(TupleToolApplyIsolation, name="TupleToolApplyIsolationSoft")
        data.B.TupleToolApplyIsolationSoft.OutputSuffix="_Soft"
        data.B.TupleToolApplyIsolationSoft.WeightsFile="weightsSoft.xml"
        data.B.ToolList+=["TupleToolApplyIsolation/TupleToolApplyIsolationSoft"]  
        vtxiso = data.B.addTupleTool("TupleToolVtxIsoln")
        data.B.TupleToolApplyIsolationHard.OutputLevel = 3 
        data.B.TupleToolApplyIsolationSoft.OutputLevel = 3 

        data.addTupleTool('TupleToolPropertime')
        tts = data.addTupleTool('TupleToolStripping')
        data.addTupleTool('TupleToolPropertime')

        if mc:
            mcTruth = data.addTupleTool('TupleToolMCTruth')
            mcTruth.addTupleTool('MCTupleToolHierarchy')
            data.addTupleTool('TupleToolMCBackgroundInfo')

        sequence.Members.append(data)

    return [sequence]

def create_options(year, magnet, mc):
    DaVinci().appendToMainSequence(create_tuples(mc))

    DaVinci().EvtMax = -1
    DaVinci().TupleFile = 'out.root'

    DaVinci().DataType = str(year)

    if mc and magnet == 'up':
        DaVinci.Simulation() = True
        DaVinci().CondDBtag = "sim-20130522-1-vc-mu100"
        DaVinci().DDDBtag = "dddb-20130929-1"
    elif mc and magnet == 'down':
        DaVinci().Simulation = True
        DaVinci().CondDBtag = "sim-20130522-1-vc-md100"
        DaVinci().DDDBtag = "dddb-20130929-1"
    if year == 2012:
        DaVinci().Lumi = True
        DaVinci().CondDBtag = "cond-20120831"
        DaVinci().DDDBtag = "dddb-20120831"
    elif year == 2011:
        DaVinci().Lumi = True
        DaVinci().CondDBtag = "cond-20130114"
        DaVinci().DDDBtag = "dddb-20130111"
    else:
        raise ValueError("unknown year: {}".format(year))


