from Gaudi.Configuration import *
from DecayTreeTuple.Configuration import *
from AllTriggerLines import *

# Decays
Bd2JpsiKS0 = "B0 -> ^(J/psi(1S) -> ^mu+ ^mu-) ^(KS0 -> ^pi+ ^pi-)"
Jpsi2MuMu = "J/psi(1S) -> ^mu+ ^mu-"
Bu2MuMuK = "[B+ -> ^(J/psi(1S) -> ^mu- ^mu+) ^K+ ]CC"
Bu2eeK = "[B+ -> ^(J/psi(1S) -> ^e- ^e+) ^K+ ]CC"
Bu2MueK = "[B+ -> ^(J/psi(1S) -> ^mu- ^e+) ^K+ ]CC"
Bu2eMuK = "[B+ -> ^(J/psi(1S) -> ^e- ^mu+) ^K+ ]CC"
Bu2JpsiPi = "[B+ -> ^(J/psi(1S) -> ^mu+ ^mu-) ^pi+]CC"
Bu2JpsiK = "[B+ -> ^(J/psi(1S) -> ^mu+ ^mu-) ^K+]CC"
Bs2JpsiPhi = "B_s0 -> ^(J/psi(1S) -> ^mu+ ^mu-) ^(phi(1020) -> ^K+ ^K-)"
Bd2JpsiPi0 = "B0 -> ^(J/psi(1S) -> ^mu+ ^mu-) ^pi0"
Bd2JpsiKstar = "[B0 -> ^(J/psi(1S) -> ^mu+ ^mu-) ^(K*(892)0 -> ^K+ ^pi-)]CC"

Bd2eeKst = "[B0 -> ^(J/psi(1S) -> ^e+ ^e-) ^(K*(892)0 -> ^K+ ^pi-) ]CC"
Bd2MuMuKst = "[B0 -> ^(J/psi(1S) -> ^mu+ ^mu-) ^(K*(892)0 -> ^K+ ^pi-) ]CC"
Bd2MuMuKstSS = "[B0 -> ^(J/psi(1S) -> ^mu+ ^mu+) ^(K*(892)0 -> ^K+ ^pi-) ]CC"
Bd2MueKst = "[B0 -> ^(J/psi(1S) -> ^mu+ ^e-) ^(K*(892)0 -> ^K+ ^pi-) ]CC"
Bd2eMuKst = "[B0 -> ^(J/psi(1S) -> ^e+ ^mu-) ^(K*(892)0 -> ^K+ ^pi-) ]CC"
Lambdab2MuMupK = "[Lambda_b0 -> ^(J/psi(1S) -> ^mu+ ^mu-) ^(Lambda(1520)0 -> ^p+ ^K-) ]CC"

Bd2MuMuD0 = "[B0 -> ^(J/psi(1S) -> ^mu+ ^mu-) ^(D~0 -> ^K+ ^pi-) ]CC"

Lambdab2Jpsippi = "[Lambda_b0 -> ^(J/psi(1S) -> ^mu+ ^mu-) ^p+ ^pi-]CC"
Lambdab2JpsiLambda = "[Lambda_b0 -> ^(J/psi(1S) -> ^mu+ ^mu-) ^(Lambda0 -> ^p+ ^pi-)]CC"
Xib2Lambdabpi = "[Xi_b- -> ^(Lambda_b0 -> ( ^J/psi(1S) -> ^mu+ ^mu-) ^(Lambda0 -> ^p+ ^pi-)) ^pi-  ]CC"
Xib2LambdabpiWS = "[Xi_b~+ -> ^(Lambda_b0 -> ( ^J/psi(1S) -> ^mu+ ^mu-) ^(Lambda0 -> ^p+ ^pi-)) ^pi+ ]CC"
Xib2PsiLpi = "[Xi_b- -> ( ^J/psi(1S) -> ^mu+ ^mu-) ^(Lambda0 -> ^p+ ^pi-) ^pi-  ]CC"
Xib2PsiLpiWS = "[Xi_b~+ -> ( ^J/psi(1S) -> ^mu+ ^mu-) ^(Lambda0 -> ^p+ ^pi-) ^pi+ ]CC"

##############################################################################
#
# Tuple maker
#
def addTuple(name="", decay="", addendum="", head="/Event/Phys/", dtf=True, resonant=True, shortname="", verbose=[] ):
    from Configurables import DecayTreeTuple, PrintDecayTree, FilterDesktop, GaudiSequencer, PrintHeader, TESCheck
    
    if shortname == "": shortname = name
    shortname = shortname+"_Tuple"+addendum
    shortname = shortname.replace("MCMC","MC")
    seq = GaudiSequencer("Seq"+shortname)
    if ( not "/"==head[-1] ): head = head+'/'

    location = head+name+"/Particles" 
    from Configurables import LoKi__HDRFilter
    if ( "/Event/Phys/" == head ):
        filter = TESCheck("Check"+shortname,Inputs = [ location ], Stop = False)
    else :              # I am not running the selection, hence the stripping decision must be here     
        filter = LoKi__HDRFilter( "Check"+shortname,
                                  Code = "HLT_PASS('Stripping"+name+"Decision')",
                                  Location="/Event/Strip/Phys/DecReports" )

    #filter.OutputLevel = 1
    seq.Members += [  filter ] # PrintHeader(),
    tuple = DecayTreeTuple(shortname)
    isMDST = (addendum.upper()=="MDST")        
    if (isMDST):
        RIT = head.replace("/Phys","")
        print "RootInTES set to", RIT
        tuple.RootInTES = RIT
        tuple.Inputs = [ "Phys/"+name+"/Particles" ]
    else :
        tuple.Inputs = [ location ]

#    tuple.OutputLevel = 1
    tuple.ToolList =  []
        
    tuple.Decay = decay

    tg = tuple.addTupleTool("TupleToolGeometry")
    if not isMDST: tg.FillMultiPV = True


    tlist = []
    if ("e+" in decay):
        tlist = electronLines()            
    elif ("mu+" in decay):
        tlist = muonLines()

    if ( False ):
        tlist = allLines()
    print tlist
    
    if ( Jpsi2MuMu != decay ): bpsi = (decay.replace("^","")).replace("(J/psi(1S)","^(J/psi(1S)")
    else : bpsi = "^("+decay.replace("^","")+")"

    print "J/psi branch is `` ", bpsi, "''" 
    tuple.Branches["Psi"] = bpsi

    # sort out kstars
    if "892" in decay:
        bkstar = (decay.replace("^","")).replace("(K*(892)","^(K*(892)")
        tuple.Branches["Kstar"] = bkstar
        Kstar = tuple.addTupleTool("TupleToolDecay/Kstar")

    from Configurables import TupleToolTISTOS
    tistos = TupleToolTISTOS(TriggerList = tlist
            , VerboseHlt1 = True, VerboseHlt2 = True, VerboseL0 = True)
    
    Psi = tuple.addTupleTool("TupleToolDecay/Psi")
    Psi.addTool(tistos)
    Psi.ToolList += [ "TupleToolTISTOS" ]
#    if (not isMDST):
#        vi = tuple.Psi.addTupleTool("TupleToolVtxIsoln")
#        vi.InputParticles = [ "/Event/Phys/MyGoodPions" ]

    if ( Jpsi2MuMu == decay ):
        if (dtf):
            pvfit = tuple.Psi.addTupleTool("TupleToolDecayTreeFitter/PVFit")        # fit with all constraints I can think of
            pvfit.Verbose = True
            pvfit.constrainToOriginVertex = True
        
    else:       
        B = tuple.addTupleTool("TupleToolDecay/B")
        if ( Bs2JpsiPhi==decay ):
            p2vv = B.addTupleTool("TupleToolP2VV") 
            p2vv.Calculator = "Bs2JpsiPhiAngleCalculator"
        elif ( "K*(892)0" in decay and not Bd2MuMuKstSS==decay ):
            p2vv = B.addTupleTool("TupleToolP2VV") 
            p2vv.Calculator = "Bd2KstarMuMuAngleCalculator"
        if (Lambdab2Jpsippi==decay ): B.addTupleTool("TupleToolDalitz") 
 
        if ('Xi_b-' in decay ): bh = 'Xi_b-'
        elif ('Xi_b~+' in decay ): bh = 'Xi_b~+'
        elif ('Lambda_b0' in decay ): bh = 'Lambda_b0'
        elif ('B0' in decay): bh = 'B0'
        elif ('B+' in decay ): bh = 'B+'
        elif ('B_s0' in decay ): bh = 'B_s0'
        if ('CC' in decay): bh = '['+bh+']cc'
        print "Branch will be ``", bh+" : "+decay.replace("^",""), "''"

        tuple.Branches["B"] = "^("+decay.replace("^","")+")"
        
        # This is needed for ConstB
        bhh = bh.replace('[','').replace(']cc','')

        if (('pi0' not in decay) and ('KS0' not in decay) and ('Lambda0' not in decay)):
#            B.ToolList +=  [ "TupleToolAngles"] 
            tta = tuple.addTupleTool("TupleToolAngles") 
            if (not resonant): tta.WRTMother = False 
#            if (( 'K*(892)0' in decay)  or ('phi' in decay)): tuple.B.ToolList +=  [ "TupleToolP2VV" ]  # really wants a VV
        if (dtf):
            FullFit = B.addTupleTool("TupleToolDecayTreeFitter/FullFit")
            FullFit.Verbose = True                          # fills daughters
            FullFit.constrainToOriginVertex = True
            FullFit.daughtersToConstrain += [ "J/psi(1S)" ] # even non resonant
            if ("ConstBFit" in verbose):
                ConstBFit = B.addTupleTool("TupleToolDecayTreeFitter/ConstBFit")
                ConstBFit.Verbose = True
                ConstBFit.constrainToOriginVertex = True
                ConstBFit.daughtersToConstrain += [ "J/psi(1S)", bhh ] # constrain B as well
                ConstBFit.UpdateDaughters = True                       # Maurizio's daughters

            for d in [ 'KS0', 'pi0', 'Lambda0' ]:
                if d in decay:
                    FullFit.daughtersToConstrain += [ d ]
                    if ( "ConstBFit" in verbose ): ConstBFit.daughtersToConstrain += [ d ]
            if (not resonant):
                NoPsiFit = B.addTupleTool("TupleToolDecayTreeFitter/NoPsiFit")
                NoPsiFit.Verbose = True
                NoPsiFit.constrainToOriginVertex = True
            if (not resonant): #b constraint for calculating qsq
                ConstBNoPsiFit = B.addTupleTool("TupleToolDecayTreeFitter/ConstBNoPsiFit")
                ConstBNoPsiFit.Verbose = True
                ConstBNoPsiFit.constrainToOriginVertex = True
                ConstBNoPsiFit.daughtersToConstrain += [ bhh ] # constrain B as well
                ConstBNoPsiFit.UpdateDaughters = True                       # Maurizio's daughters

            #if (resonant and ('KS0' in decay or "pi0" in decay or 'Lambda0' in decay)):
            #    OnlyPsiFit = B.addTupleTool("TupleToolDecayTreeFitter/OnlyPsiFit")
            #    OnlyPsiFit.Verbose = True
            #    OnlyPsiFit.constrainToOriginVertex = True
            #    OnlyPsiFit.daughtersToConstrain += [ "J/psi(1S)" ]
            if ('phi' in decay):
                substitute1 = { 'Beauty -> Meson (Strange -> ^K+ K-)': 'pi+'  }
                from Configurables import TupleToolDecayTreeFitter
                subDTF = TupleToolDecayTreeFitter("SubPipKm", Verbose=False,
                                                  daughtersToConstrain = [ "J/psi(1S)" ],
                                                  constrainToOriginVertex=True,
                                                  Substitutions=substitute1)
                B.addTool(subDTF)
                B.ToolList += [ "TupleToolDecayTreeFitter/SubPipKm" ]
                
                B.addTool(subDTF)
                substitute2 = { 'Beauty -> Meson (Strange -> K+ ^K-)': 'pi-' }
                B.addTool(subDTF.clone("SubKpPim",Substitutions=substitute2))
                B.ToolList += [ "TupleToolDecayTreeFitter/SubKpPim" ]

            if ( Lambdab2Jpsippi==decay ):
                substitute1 = { 'Beauty -> Meson p+ ^pi-': 'K-', 'Beauty -> Meson p~- ^pi+': 'K+'  }
                from Configurables import TupleToolDecayTreeFitter
                B.addTupleTool(TupleToolDecayTreeFitter("SubpK", Verbose=False,
                                                        daughtersToConstrain = [ "J/psi(1S)" ],
                                                        constrainToOriginVertex=True,
                                                        Substitutions=substitute1))
                        
                substitute2 = { 'Beauty -> Meson ^p+ pi-': 'K+', 'Beauty -> Meson ^p~- pi+': 'K-'  }
                B.addTupleTool(TupleToolDecayTreeFitter("SubKpi", Verbose=False,
                                                        daughtersToConstrain = [ "J/psi(1S)" ],
                                                        constrainToOriginVertex=True,
                                                        Substitutions=substitute2))
                if ("ConstBFit" in verbose):
                    ConstBSubFit = B.addTupleTool("TupleToolDecayTreeFitter/ConstBSubFit")
                    ConstBSubFit.Substitutions=substitute1   # substitute and then constrain
                    ConstBSubFit.Verbose = True
                    ConstBSubFit.constrainToOriginVertex = True
                    ConstBSubFit.daughtersToConstrain += [ "J/psi(1S)", bhh ] # constrain B as well
                    ConstBSubFit.UpdateDaughters = True                       # Maurizio's daughters

            if ( Bd2JpsiKS0==decay ):
                substitute1 = { 'Beauty -> Meson (Strange -> ^pi+ pi-)' : 'K+' }
                substitute2 = { 'Beauty -> Meson (Strange -> pi+ ^pi-)' : 'K-' }
                from Configurables import TupleToolDecayTreeFitter
                B.addTupleTool(TupleToolDecayTreeFitter("PsiKppim", Verbose=False,
                                                        daughtersToConstrain = [ "J/psi(1S)" ],
                                                        constrainToOriginVertex=True,
                                                        Substitutions=substitute1))
                        
                B.addTupleTool(TupleToolDecayTreeFitter("PsiKmpip", Verbose=False,
                                                        daughtersToConstrain = [ "J/psi(1S)" ],
                                                        constrainToOriginVertex=True,
                                                        Substitutions=substitute2))

            if ( Bu2JpsiK == decay):
            
                substitute2 = { 'Beauty -> Meson ^K+': 'pi+', 'Beauty -> Meson ^K-': 'pi-' }
                subDTF2 = B.addTupleTool(TupleToolDecayTreeFitter("Subpi", Verbose=False,
                                                                  daughtersToConstrain = [ "J/psi(1S)" ],
                                                                  constrainToOriginVertex=True,
                                                                  Substitutions=substitute2))
            
            if ( Bu2eeK == decay or Bu2MuMuK == decay or Bu2MueK == decay or Bu2eMuK == decay):
                if Bu2eeK == decay : substitute1 = {  "Beauty -> ( J/psi(1S) -> ^e+ X- ) Strange " : "pi+" ,
                                                      "Beauty -> ( J/psi(1S) -> ^e- X+ ) Strange " : "pi-" }
                elif Bu2MuMuK == decay : substitute1 = {  "Beauty -> ( J/psi(1S) -> ^mu+ X- ) Strange " : "pi+" ,
                                                          "Beauty -> ( J/psi(1S) -> ^mu- X+ ) Strange " : "pi-" }
                elif Bu2MueK == decay : substitute1 = {  "Beauty -> ( J/psi(1S) -> ^e+ X- ) Strange " : "pi+" ,
                                                         "Beauty -> ( J/psi(1S) -> ^mu- X+ ) Strange " : "pi-" }
                elif Bu2eMuK == decay : substitute1 = {  "Beauty -> ( J/psi(1S) -> ^mu+ X- ) Strange " : "pi+" ,
                                                         "Beauty -> ( J/psi(1S) -> ^e- X+ ) Strange " : "pi-" }

                from Configurables import TupleToolDecayTreeFitter
                subDTF = B.addTupleTool(TupleToolDecayTreeFitter("Subpipi", Verbose=False,
                                                                 constrainToOriginVertex=True,
                                                                 Substitutions=substitute1,
                                                                 UpdateDaughters = True)) # Maurizio's daughters
                    
        """
        if (('KS0' in decay) or ("K*(892)" in decay) or ("Lambda0" in decay)):
            mh = tuple.addTupleTool("TupleToolMassHypo")
            if ("KS0" in decay):
                mh.PIDReplacements = { "pi+" : "p+" }    
            #if ("K*(892)" in decay):
            #    mh.PIDReplacements = { "pi-" : "K+" }     
            if ("Lambda0" in decay):
                mh.PIDReplacements = { "p+" : "pi+" }
        """
                
        B.addTool(tistos)
        B.ToolList += [ "TupleToolTISTOS"  ]

    rs = tuple.addTupleTool("TupleToolRecoStats")
    rs.Verbose = ("RecoStats" in verbose)
    ei = tuple.addTupleTool("TupleToolEventInfo")
    ei.Verbose = ("EventInfo" in verbose)          # gives GpsYear, month, day, etc
    ki = tuple.addTupleTool("TupleToolKinematic")  # gives _AtVtx_P
    ki.Verbose=  ("Kinematic" in verbose)
    tuple.addTupleTool("TupleToolPrimaries")
    ti = tuple.addTupleTool("TupleToolTrackInfo")
    ti.Verbose = ("TrackInfo" in verbose)          # gives many more chi2
    if ( "TrackPosition" in  verbose):             # extrapolates to TT, IT
        tp1 = tuple.addTupleTool("TupleToolTrackPosition/TupleToolTrackPositionTT")
        tp2 = tuple.addTupleTool("TupleToolTrackPosition/TupleToolTrackPositionIT")
        tp1.ExtraName = "TT"
        tp2.ExtraName = "IT"
        tp2.Z = 7800    
    pid = tuple.addTupleTool("TupleToolPid")
    pid.Verbose = ("Pid" in verbose)               # many more vars in RICH
    if ("Pid" in verbose) :
        tuple.addTupleTool("TupleToolANNPID")
    tuple.addTupleTool("TupleToolL0Calo")
    if ('e+' in decay):
        brem = tuple.addTupleTool("TupleToolBremInfo")
    #isolation
    if "iso" in verbose :
        # lines from Greg Ciezarek
        from configIso import configIso
        configIso()
        # this is work in progress 
        #will be comitted to SVN properly at a future date as of 3/4/14
        from Configurables import TupleToolApplyIsolation
        tuple.B.addTool(TupleToolApplyIsolation,    name="TupleToolApplyIsolationHard")
        tuple.B.TupleToolApplyIsolationHard.OutputSuffix="_Hard"
        tuple.B.TupleToolApplyIsolationHard.WeightsFile="weightsHard.xml"
        tuple.B.ToolList+=["TupleToolApplyIsolation/TupleToolApplyIsolationHard"]
        tuple.B.addTool(TupleToolApplyIsolation,    name="TupleToolApplyIsolationSoft")
        tuple.B.TupleToolApplyIsolationSoft.OutputSuffix="_Soft"
        tuple.B.TupleToolApplyIsolationSoft.WeightsFile="weightsSoft.xml"
        tuple.B.ToolList+=["TupleToolApplyIsolation/TupleToolApplyIsolationSoft"]  
        vtxiso = tuple.B.addTupleTool("TupleToolVtxIsoln")
        tuple.B.TupleToolApplyIsolationHard.OutputLevel = 3 
        tuple.B.TupleToolApplyIsolationSoft.OutputLevel = 3 
#    tuple.addTupleTool("TupleToolPropertime")
    if ( "/Event/Phys/" == head): # not reading stripping output
        tts = tuple.addTupleTool("TupleToolStripping")
        tts.StrippingList = [ 
            "StrippingBetaSBd2JpsiKsDetachedLineDecision", 
            "StrippingBetaSBd2JpsiKsPrescaledLineDecision", 
            "StrippingBetaSBd2JpsiKstarDetachedLineDecision", 
            "StrippingBetaSBd2JpsiKstarPrescaledLineDecision", 
            "StrippingBetaSBs2ChicPhi_Chic2KKPiPiNominalLineDecision", 
            "StrippingBetaSBs2ChicPhi_Chic2PiPiPiPiNominalLineDecision", 
            "StrippingBetaSBs2EtacPhiNominalLineDecision", 
            "StrippingBetaSBs2JpsiEtaDetachedLineDecision", 
            "StrippingBetaSBs2JpsiEtaPrescaledLineDecision", 
            "StrippingBetaSBs2JpsiKstarLineDecision", 
            "StrippingBetaSBs2JpsiPhiDetachedLineDecision", 
            "StrippingBetaSBs2JpsiPhiPrescaledLineDecision", 
            "StrippingBetaSBs2JpsieePhiDetachedLineDecision", 
            "StrippingBetaSBs2JpsieePhiLineDecision", 
            "StrippingBetaSBs2Jpsif0LineDecision", 
            "StrippingBetaSBs2K0stK0stNominalLineDecision", 
            "StrippingBetaSBs2KstKstNominalLineDecision", 
            "StrippingBetaSBs2KstKstSameChargeLineDecision", 
            "StrippingBetaSBs2PhiKstNominalLineDecision", 
            "StrippingBetaSBs2PhiPhiLineDecision", 
            "StrippingBetaSBs2PhiPhiWideLineDecision", 
            "StrippingBetaSBs2Q2Body4piLineDecision", 
            "StrippingBetaSBu2JpsiKDetachedLineDecision", 
            "StrippingBetaSBu2JpsiKNoPIDDetachedLineDecision", 
            "StrippingBetaSBu2JpsiKPrescaledLineDecision", 
            "StrippingBetaSJpsi2MuMuDetachedLineDecision", 
            "StrippingBetaSJpsi2MuMuLineDecision", 
            "StrippingBetaSLambdab2JpsiLambdaUnbiasedLineDecision", 
            "StrippingBu2LLK_eeLineDecision", 
            "StrippingBu2LLK_mmLineDecision", 
            "StrippingBd2KstarMuMu_BdToKstarMuMuLineDecision", 
            "StrippingBd2KstarMuMu_BdToKstarMuMuLowPLineDecision", 
            "StrippingBd2KstarMuMu_BdToKstarMuMuSSLineDecision", 
            "StrippingBd2KstarMuMu_BdToKstarMuMuSSLowPLineDecision", 
            "StrippingBd2KstarMuMu_BuToKMuMuLineDecision", 
            "StrippingBd2KstarMuMu_BuToKMuMuSSLineDecision" ]

    tuple.addTupleTool("TupleToolPhotonInfo")
    tuple.addTupleTool("TupleToolPi0Info")
#    tuple.OutputLevel = 1
    if ( addendum == "MC" ):   # it's MC!
        tmc = tuple.addTupleTool("TupleToolMCTruth")
        #tmc.addTupleTool("MCTupleToolKinematic")
        tmc.addTupleTool("MCTupleToolHierarchy")
        tuple.addTupleTool("TupleToolMCBackgroundInfo")

    if (isMDST):
        seq.Members += [ #PrintDecayTree( Inputs = [ location]),
            EventNodeKiller(Nodes=["DAQ"] ) ]

    if ( addendum == "MC" ):   # it's not MC!
        from Configurables import TrackSmearState as SMEAR
        smear = SMEAR()
#        seq.Members += [ smear ]
    elif (not isMDST):   # it's not mdst
        from Configurables import TrackScaleState as SCALER
        scaler = SCALER('StateScale')   
#        seq.Members += [ scaler ]

    seq.Members += [ # PrintDecayTree( Inputs = [ location]),
                     tuple ]
    return [ seq ]
