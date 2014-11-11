
extraopts = { 'up2012':
                 '''
                 DaVinci().DataType = "2012"
                 DaVinci().Lumi = True
                 DaVinci().CondDBtag = ?
                 DaVinci().DDDBtag = ?
                 ''',
              'down2012':
                 '''
                 DaVinci().DataType = "2012"
                 DaVinci().Lumi = True
                 DaVinci().CondDBtag = ?
                 DaVinci().DDDBtag = ?
                 ''',
              'up2011':
                 '''
                 DaVinci().DataType = "2011"
                 DaVinci().Lumi = True
                 DaVinci().CondDBtag = ?
                 DaVinci().DDDBtag = ?
                 ''',
              'down2011':
                 '''
                 DaVinci().DataType = "2011"
                 DaVinci().Lumi = True
                 DaVinci().CondDBtag = ?
                 DaVinci().DDDBtag = ?
                 '''
            }

datasets = { 'up2012': ?,
             'down2012': ?,
             'up2011': ?,
             'down2011': ?,
           }

data_jobs = {}

for kind in ['up2012', 'down2012', 'up2011', 'down2011']:

    j = Job()
    j.application = DaVinci()
    data = j.application.readInputData(datasets[kind])
    j.inputdata = data

    j.application.make()
    j.application.prepare()
    j.application.optsfile = './DVBd2MuMuD0_data.py'
    j.application.extraopts = extraopts[kind]

    j.backend = Dirac()
    j.splitter = SplitByFiles(filesPerJob=1)
    j.merger = RootMerger(files=['Bd2MuMuD0.root'])

    j.inputsandbox = ['./weightsHard.xml', './weightsSoft.xml']

    data_jobs[kind] = j

