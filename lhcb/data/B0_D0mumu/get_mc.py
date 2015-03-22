
extraopts = {
'up': '''
DaVinci().DataType = "2012"
DaVinci().Simulation = True
DaVinci().CondDBtag = "sim-20130522-1-vc-mu100"
DaVinci().DDDBtag = "dddb-20130929-1"''',
'down': '''
DaVinci().DataType = "2012"
DaVinci().Simulation = True
DaVinci().CondDBtag = "sim-20130522-1-vc-md100"
DaVinci().DDDBtag = "dddb-20130929-1"'''
}

datasets = { 'up': './bookkeeping/BK_MC_2012_Beam4000GeV-2012-MagDown-Nu2.5-Pythia8_Sim08f_Digi13_Trig0x409f0045_Reco14a_Stripping20NoPrescalingFlagged_11174011_ALLSTREAMS.DST.py',
             'down': './bookkeeping/BK_MC_2012_Beam4000GeV-2012-MagDown-Nu2.5-Pythia8_Sim08f_Digi13_Trig0x409f0045_Reco14a_Stripping20NoPrescalingFlagged_11174011_ALLSTREAMS.DST.py'
           }

mc_jobs = {}

for kind in ['up', 'down']:

    j = Job()
    j.application = DaVinci()

    data = j.application.readInputData(datasets[kind])
    j.inputdata = data

    j.application.make()
    j.application.optsfile = './DVBd2MuMuD0_MC.py'
    j.application.extraopts = extraopts[kind]

    j.backend = Dirac()
    j.splitter = SplitByFiles(filesPerJob=1)
    j.merger = RootMerger(files=['Bd2MuMuD0.root'])

    # Reuse data vertex isolation classifier
    j.inputsandbox = ['./weightsHard.xml', './weightsSoft.xml']

    mc_jobs[kind] = j

