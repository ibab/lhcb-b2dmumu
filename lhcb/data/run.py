
decays = {
    #'Bd_D0mumu':   '[B0   -> ^(D~0 -> ^K+ ^pi-) ^(J/psi(1S) -> ^mu+ ^mu-) ]CC',
    #'Bc_Ds+mumu':  '[B_c+ -> ^(D_s+ -> ^K+ ^K- ^pi+) ^(J/psi(1S) -> ^mu+ ^mu-)]CC',
    #'Bu_Ds+mumu':  '[B+   -> ^(D_s+ -> ^K+ ^K- ^pi+) ^(J/psi(1S) -> ^mu+ ^mu-]CC',
    #'Bu_Dst+mumu': '[B+   -> ^(D*(2010)+ -> ^(D0 -> ^K- ^pi+) ^pi+) ^(J/psi(1S) -> ^mu+ ^mu-)]CC',
    'Bd_KstJpsi':   '[B0   -> ^(K*(892)0 -> ^K+ ^pi-) ^(J/psi(1S) -> ^mu+ ^mu-) ]CC',
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

            name = 'DATA_' + decay + "_" + str(year) + "_" + polarity
            outname = decay + '.root'
            j = Job(name=name)
            print(j.name)
            j.application = DaVinci()
            j.application.version = 'v36r5'
            j.inputdata = data
            j.application.make()
            j.application.optsfile = './options.py'
            j.application.extraopts = 'create_options(decay="{decay}", year={year}, outname="{outname}")'.format(decay=decays[decay], year=year, outname=outname)
            j.inputsandbox = ['./weightsHard.xml', './weightsSoft.xml']
            j.backend = Dirac()
            j.splitter = SplitByFiles(filesPerJob=5)
            j.merger = RootMerger(files=[outname])

