# lines from Greg Ciezarek -
# this is work in progress 
#will be comitted to SVN properly at a future date  as of 3/4/14
# this needs to go before the tuple tool is run
def configIso():
    from Configurables import ChargedProtoParticleMaker, DaVinci
    veloprotos = ChargedProtoParticleMaker("ProtoPMaker")
    veloprotos.Inputs = ["Rec/Track/Best"]
    veloprotos.Output = "Rec/ProtoP/myProtoPMaker/ProtoParticles"
    DaVinci().appendToMainSequence( [ veloprotos ])

    from Configurables       import ProtoParticleCALOFilter, CombinedParticleMaker,NoPIDsParticleMaker
    from CommonParticles.Utils import trackSelector, updateDoD
    algorithm = NoPIDsParticleMaker('StdNoPIDsVeloPions',  Particle =
            'pion',  )
    algorithm.Input = "Rec/ProtoP/myProtoPMaker/ProtoParticles"
    selector = trackSelector ( algorithm , trackTypes = ['Velo'] )
    locations = updateDoD ( algorithm )
    DaVinci().appendToMainSequence( [ algorithm ]) 
