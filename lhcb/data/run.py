
decays = {
    'Bd_D0mumu':   '[B0   -> ^(D~0 -> ^K+ ^pi-) ^(J/psi(1S) -> ^mu+ ^mu-) ]CC',
}

stream = "/LHCb/Collision{year}/Beam{energy}GeV-VeloClosed-{polarity}/Real Data/Reco14/Stripping{stripping}/90000000/DIMUON.DST"

for decay in decays.keys():
    for year in [2011, 2012]:
        for polarity in ['MagUp', 'MagDown']:
            if year == 2011:
                energy = '3500'
                stripping = '20r1'
            if year == 2012:
                energy = '4000'
                stripping = '20'

            filepath = stream.format(year=str(year)[-2:], energy=energy, polarity=polarity, stripping=stripping)
            data = BKQuery(filepath, dqflag=['OK']).getDataset()

            name = 'DATA_B_Dmumu_' + str(year) + '_' + polarity
            j = Job(name=name)
            print(j.name)
            j.application = DaVinci()
            j.application.version = 'v36r5'
            j.inputdata = data
            j.application.make()
            j.application.optsfile = './options.py'
            j.application.extraopts = 'create_options(year={year}, )'.format(year=repr(year))
            j.inputsandbox = ['./weightsHard.xml', './weightsSoft.xml']
            j.backend = Dirac()
            j.splitter = SplitByFiles(filesPerJob=5)
            j.merger = RootMerger(files=['out.root'])

