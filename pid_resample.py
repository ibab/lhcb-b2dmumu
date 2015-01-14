# Adapted from Konstantin Schubert's c253dde3b31fa27e8081514cd37e10076a291f7b

def create_resampler(use_probnn=True):
    import sys
    sys.path.append("/home/igor/Mount/eve/fhgfs/users/ashires/")
    import os
    import time
    from math import isnan, sqrt, log
    from ROOT import TFile, TRandom, TRandom3, TMath, gROOT

    particles = {
        211  : "Pion",
        321  : "Kaon",
        13   : "Muon",
        2212 : "Proton",
    }

    Particle = {
            'Pion': {},
            'Kaon': {},
            'Muon': {},
            'Proton': {},
    }

    def perform_resample(which_pid, true_pid, p, pt, numTracks, muons=False):
        pl = sqrt(p**2 - pt**2)
        eta = 0.5 * log(p + pl) - log(p - pl)

        try:
            name = particles[abs(true_pid)]
        except KeyError:
            return -1

        if not name in Particle:
            # FIXME 
            return -1

        if name == 'Muon':
            which_pid = which_pid.replace('Down', '').replace('Up', '')
        res = Particle[name][which_pid].GetDLLhisto(p, eta, numTracks).GetRandom()
        if isnan(res):
            return -1
        else:
            return res

    basedir = "/home/igor/Mount/eve/fhgfs/users/ashires/"
    if use_probnn:
        pid_table_path = os.path.join(basedir,"PIDEffTables_v20r1/forProbNN/" )
    else :
        pid_table_path = os.path.join(basedir,"PIDEffTables_v20r1/forDLL/" )

    from reader.memtest import stacksize, status, change #####* #stacksize
    from reader.ReadEffTables import DLLclass

    for part in ['Pion', 'Kaon', 'Muon']:
        for mag in ['Up', 'Down']:
            for target in ['DLLpi', 'DLLK', 'DLLmu']:
                if part == 'Muon':
                    # the muon dataset does not have polarity for some reason
                    dllClass = DLLclass(Particle=part, whichDLL=target, Mag=None, options='normal, both', tablepath=pid_table_path)
                    dllClass.Initialize()
                    Particle[part][target] = dllClass
                else:
                    dllClass = DLLclass(Particle=part, whichDLL=target, Mag=mag, options='normal, both', tablepath=pid_table_path)
                    dllClass.Initialize()
                    Particle[part][target+mag] = dllClass

    return perform_resample
