
extraopts = {
'up2012':
'''
DaVinci().DataType = "2012"
DaVinci().Lumi = True
DaVinci().CondDBtag = "cond-20120831"
DaVinci().DDDBtag = "dddb-20120831"
''',
'down2012':
'''
DaVinci().DataType = "2012"
DaVinci().Lumi = True
DaVinci().CondDBtag = "cond-20120831"
DaVinci().DDDBtag = "dddb-20120831"
''',
'up2011':
'''
DaVinci().DataType = "2011"
DaVinci().Lumi = True
DaVinci().CondDBtag = "cond-20130114"
DaVinci().DDDBtag = "dddb-20130111"
''',
'down2011':
'''
DaVinci().DataType = "2011"
DaVinci().Lumi = True
DaVinci().CondDBtag = "cond-20130114"
DaVinci().DDDBtag = "dddb-20130111"
'''
}

datasets = { 'up2011':   'bookkeeping/BK_LHCb_Collision11_Beam3500GeV-VeloClosed-MagUp_RealData_Reco14_Stripping20r1_FullStream.py',
             'down2011': 'bookkeeping/BK_LHCb_Collision11_Beam3500GeV-VeloClosed-MagDown_RealData_Reco14_Stripping20r1_FullStream.py',
             'up2012':   'bookkeeping/BK_LHCb_Collision12_Beam4000GeV-VeloClosed-MagUp_RealData_Reco14_Stripping20_FullStream.py',
             'down2012': 'bookkeeping/BK_LHCb_Collision12_Beam4000GeV-VeloClosed-MagDown_RealData_Reco14_Stripping20_FullStream.py',
           }

data_jobs = {}

for kind in ['up2012', 'down2012', 'up2011', 'down2011']:

    j = Job(name='Bd2D0MuMu_' + kind)
    j.application = DaVinci()
    data = j.application.readInputData(datasets[kind])
    j.inputdata = data

    j.application.make()
    j.application.optsfile = './DVBd2MuMuD0_data.py'
    j.application.extraopts = extraopts[kind]

    j.backend = Dirac()
    j.splitter = SplitByFiles(filesPerJob=20)
    j.merger = RootMerger(files=['Bd2MuMuD0.root'])

    j.inputsandbox = ['./weightsHard.xml', './weightsSoft.xml']

    data_jobs[kind] = j

