

from ROOT import * # WHY!?
import dircache
import os
import sys
from array import array

#PIDTABLEPATH = "/Users/nserra/LHCb/Kstmumu/PIDEffTables_v17/"
DEFAULTPIDTABLEPATH = "/vols/lhcbdisk03/ashires/EvtEffData/PIDEffTables_v17/"

class DLLclass:

    def __init__(self, Mag="Up", Particle="Pion", whichDLL="DLLK", 
                 options = None, tablepath=DEFAULTPIDTABLEPATH):

        ### Setting if I want PiDLL or MuonDLL 
        self.whichDLL = whichDLL
        self.Mag = Mag
        self.hFluctuation = None
        self.Particle = Particle
        self.ErrorType = "Poisson"
        self.options = None
        if options is not None:
            self.options = options.lower()

        if self.Mag!=None:
            self.fileName1D = tablepath+self.Particle+self.whichDLL+"Mag"+self.Mag+"TH1s.root"
            self.fileName3D = tablepath+self.Particle+self.whichDLL+"Mag"+self.Mag+".root"

        else:
            self.fileName1D = tablepath+self.Particle+self.whichDLL+"TH1s.root"
            self.fileName3D = tablepath+self.Particle+self.whichDLL+".root"


        self.File1D = None
        self.File3D = None
        self.TH3histos = {}
        self.TH1histos = {}

        self.sys_n = False
        self.sys_p = False
        self.sys_e = False
        self.sys_level = 0.1
        from ROOT import TRandom3
        self.r = TRandom3(0) 



    def setSystematic( self, _sys_ntracks, _sys_p, _sys_e, _syslevel=0  ) :
        self.sys_n = _sys_ntracks
        self.sys_p = _sys_p
        self.sys_e = _sys_e
        self.sys_level = _syslevel


    def FillHistos1D(self):
        """Gets the relevant 1D histograms from the file"""
        for ik in self.File1D.GetListOfKeys():

            obj = self.File1D.Get(ik.GetName())            
            if type(obj) == TH1F :
                num = int(obj.GetName().replace("h", ""))
                self.TH1histos[num] = obj

            if  type(obj) == TH3F or type(obj) == TH3D:
                self.TH1histos['aux'] = obj
        
    def FillHistos3D(self):
        """Gets the relevant 3D histograms from the file"""
        for ik in self.File3D.GetListOfKeys():

            obj = self.File3D.Get(ik.GetName())            

            if type(obj) == TH3F or type(obj) == TH3D:
                num = obj.GetName().replace("HistoCut", "")
                if num.find(">")!=-1:
                    num = num.replace(">", "")
                if num.find("<")!=-1:
                    print "The histograms in input are suppose to imply always > sign"
                    assert(False)
                    
                num = float(num)

                
                self.TH3histos[num] = obj

    
    def Initialize(self):
        """Opens the files and gets the histograms for this instance"""

        if self.options.find('efficiency') or self.options.find('both'):
            if self.options.find('verbose')!=-1:
                print "Initializing 3D hstos to compute Efficiency"

            # open and check 3D histogram file    
            self.File3D = TFile(self.fileName3D)
            if (not self.File3D) or (not self.File3D.IsOpen()):
                print "DLLclass::Initialize ERROR! problem with the file:", self.fileName3D
                sys.exit(-1)

            self.FillHistos3D()

        if self.options.find('pid') or self.options.find('both'):
            if self.options.find('verbose')!=-1:
                print "Initializing 1D hstos to compute PIDs"
                
            # open and check 1D histogram file
            self.File1D = TFile(self.fileName1D)
            if (not self.File1D) or (not self.File1D.IsOpen()):
                print "DLLclass::Initialize ERROR! problem with the file:", self.fileName1D
                sys.exit(-1)

            self.FillHistos1D()            


    def GetOrderXYZ(self, h, p, eta, ntrack):
        """Gets (p, eta, ntrack) variables in the correct (x, y, z) axis order
        for sampling the histogram"""
        xyz = {'p':p, 'eta':eta, 'ntrack':ntrack}

        xvar = h.GetXaxis().GetTitle().lower()
        x = xyz[xvar]

        yvar = h.GetYaxis().GetTitle().lower()
        y = xyz[yvar]

        zvar = h.GetZaxis().GetTitle().lower()
        z = xyz[zvar]

        return {'x':x, 'y':y, 'z':z}


    def FindElement(self, DLL, listCut):
        
        for il, i in zip(listCut, xrange(0, len(listCut))):
            if il == DLL:
                return i-1
        
        if DLL > listCut[len(listCut)-1]:
            return len(listCut)-1
        
        return -1
            
                

    def FindBin3D(self, h, x, y, z):

        ix = h.GetXaxis().FindBin(x)
        iy = h.GetYaxis().FindBin(y)
        iz = h.GetZaxis().FindBin(z)

        ### If I am outside the range I set the range extreme
        if ix ==0 :
            ix = 1
        if ix>h.GetXaxis().GetNbins():
            ix = h.GetXaxis().GetNbins()

        if iy ==0 :
            iy = 1
        if iy>h.GetYaxis().GetNbins():
            iy = h.GetYaxis().GetNbins()

        if iz ==0 :
            iz = 1
        if iz>h.GetZaxis().GetNbins():
            iz = h.GetZaxis().GetNbins()

        return h.GetBin(ix, iy, iz)


    def AdjustBin3D(self, h, ibin , xyz, adjxyz ) :
        import ROOT
        newbin = ibin
        #print ivals
        if adjxyz['x']  or adjxyz['y'] or adjxyz['z'] :
            ivals = { 'x': ROOT.Long(0), 'y' : ROOT.Long(0) , 'z': ROOT.Long(0) }
            h.GetBinXYZ( ibin, ivals['x'], ivals['y'], ivals['z'] )
            ###adjust for systematic studies
            for d in adjxyz.iterkeys() :
                if adjxyz[d] :
                    #gettop and low edges
                    highedge = 0.0
                    lowedge = 0.0
                    if d == 'x' :
                        lowedge = h.GetXaxis().GetBinLowEdge( ivals[d] )  
                        highedge = h.GetXaxis().GetBinLowEdge( ivals[d] + 1 )  
                    elif d == 'y' :
                        lowedge = h.GetYaxis().GetBinLowEdge( ivals[d] )  
                        highedge = h.GetYaxis().GetBinLowEdge( ivals[d] + 1 )  
                    elif d == 'z' :
                        lowedge = h.GetZaxis().GetBinLowEdge( ivals[d] )  
                        highedge = h.GetZaxis().GetBinLowEdge( ivals[d] + 1 )  
                    else : 
                            print "unkown dimension"
                            return ibin
                    # calc min / max 10 per cent hard coded 
                    lval = lowedge + ( 1.0 * ( highedge - lowedge ) ) 
                    hval = highedge - ( 1.0 * ( highedge - lowedge ) ) 
                    
                    #val = self.r.Uniform(0,1) 
                    if self.sys_level < 0  and xyz[d] < lval :
                        ivals[d] -= 1
                    elif self.sys_level > 0 and xyz[d] > hval : 
                        ivals[d] +=1
                    ###set minumum
                    if ivals[d] < 1 :
                        ivals[d] = 1
            ### If I am above the range I set the range extreme
            if ivals['x']>h.GetXaxis().GetNbins():
                ivals['x'] = h.GetXaxis().GetNbins()
            if ivals['y']>h.GetYaxis().GetNbins():
                ivals['y'] = h.GetYaxis().GetNbins()
            if ivals['z']>h.GetZaxis().GetNbins():
                ivals['z'] = h.GetZaxis().GetNbins()
            newbin = h.GetBin(ivals['x'], ivals['y'], ivals['z'])                       
            #print adjxyz, ivals, ibin, newbin
        return newbin

    def GetDLLhisto(self, P, ETA, nTrack):

        if self.TH1histos == None :
            print 'Please initialize efficiency histograms!'
            self.PrintInstructions()
            return
        
        ######## HERE #########  
        # print self.TH1histos
        xyz = self.GetOrderXYZ(self.TH1histos['aux'], P, ETA, nTrack)
        ibin = self.FindBin3D(self.TH1histos['aux'], xyz['x'], xyz['y'], xyz['z'])

        adjxyz = self.GetOrderXYZ( self.TH1histos['aux'] , self.sys_p, self.sys_e, self.sys_n )
        ibin = self.AdjustBin3D( self.TH1histos['aux'] ,ibin,xyz,  adjxyz ) 
        #print xyz, ibin, adjxyz
        return self.TH1histos[ibin]

                
    def ComputeEff(self, DLLK, P, ETA, nTrack, CUT="Efficiency"):

        if self.options.find('verbose')!=-1:
            print "The default cut is DLL>value"
            print "To have the cut DLL<value, specify the option CUT='Invert' "

        if self.TH3histos == None :
            print 'Please initialize efficiency histograms!'
            self.PrintInstructions()
            return

        listCut = self.TH3histos.keys()
        listCut.sort()
        found = 0

        if self.FindElement(DLLK, listCut) !=-1:
            found =1
            
        if found ==0:
            listCut.append(DLLK)
            
        listCut.sort()
        if self.options.find('verbose')!=-1:
            print listCut
        
        
        i = self.FindElement(DLLK, listCut)

        if self.options.find('verbose')!=-1:
            print "The cut is ",DLLK
            print "The next closer cut is ",listCut[i]
           

        xyz = self.GetOrderXYZ(self.TH3histos[listCut[i]], P, ETA, nTrack)

        if self.options.find('verbose')!=-1:
            print "This is the ordered coordinate"
            print xyz

        ibin = self.FindBin3D(self.TH3histos[listCut[i]], xyz['x'], xyz['y'], xyz['z'])

        if self.options.find('verbose')!=-1:
            print "Bin Number ",ibin


        val1 = self.TH3histos[listCut[i]].GetBinContent(ibin)

        err1 = self.TH3histos[listCut[i]].GetBinError(ibin)

        if found == 1:
            if CUT == 'efficiency':
                return {'eff':val1, 'err':err1}
            if CUT.lower() == 'invert':
                return {'eff':1-val1, 'err':err1}

        if self.options.find('verbose')!=-1:
            print "The efficiency value is ",val1

        ibin = self.FindBin3D(self.TH3histos[listCut[i+2]], xyz['x'], xyz['y'], xyz['z'])
        val2 = self.TH3histos[listCut[i+2]].GetBinContent(ibin)
        err2 = self.TH3histos[listCut[i+2]].GetBinError(ibin)

        Dval = val2-val1
        Derr = err2-err1

        val = val1+Dval*(float(DLLK-listCut[i])/float(listCut[i+2]-listCut[i]))
        err = err1+Derr*(float(DLLK-listCut[i])/float(listCut[i+2]-listCut[i]))


        if CUT.lower() == 'efficiency':
            return {'eff':val, 'err':err}
        if CUT.lower() == 'invert':
            return {'eff':1-val, 'err':err}



        
    def PrintInstructions(self):
        print "************************************************"
        print "Possible options are:"
        print "PID: It simulate the PID distributions for each phase space bin"
        print "Efficiency: It computes the efficiency for a particular PID cut"
        print "Both: You can do both things but it uses more RAM"
        print "The default is None, to initialize it chooce the option and the do initialize()"
        print "************************************************"

        return 

    def CalcN(self, sigma, epsilon):
        if sigma == 0:
            return 0
        
        N = epsilon*(1-epsilon)/(sigma*sigma)
        return N


    def ProduceTH1histos(self, error="Piosson"):

        f = TFile(self.fileName3D)

        
        keys = f.GetListOfKeys()

        histos3D = {}
        binning = {}
        cutList = []

        for ik in keys:
            obj = f.Get(ik.GetName())


            if type(obj) == TH3F or type(obj) == TH3D:
                cut = ik.GetName().replace("HistoCut", "")
                if cut.find(">")!=-1:
                    cut = cut.replace(">", "")
                if cut.find("<")!=-1:
                    print "The histograms in input are suppose to imply always > sign"
                    assert(False)
                    
                cut = float(cut)
                histos3D[cut]=obj
                cutList.append(cut)

            if type(obj) == RooBinning:
                binning[ik]=obj

        cutList.sort()
        if self.options.find('verbose')!=-1:
            print cutList
        cutList.append(cutList[-1]+1.0)

        binning = array('f', cutList)
        histos1D = []

        print histos3D
        print cutList
        #assert(False)

        aux = cutList[0]
        nameX = histos3D[aux].GetXaxis().GetTitle()
        nameY = histos3D[aux].GetYaxis().GetTitle()
        nameZ = histos3D[aux].GetZaxis().GetTitle()


        ###########################################################
        ##### I save in TH1 Format!!! 
        ###########################################################
        nameF =  self.fileName3D.replace(".root", "")
        fOutput = TFile(nameF+"TH1s.root", "RECREATE")

        histos3D[aux].SetName('h3D_example')
        histos3D[aux].Write()

        for ix in xrange(1, histos3D[aux].GetXaxis().GetNbins()+1):
            for iy in xrange(1, histos3D[aux].GetYaxis().GetNbins()+1):
                for iz in xrange(1, histos3D[aux].GetZaxis().GetNbins()+1):

                    lowX =  histos3D[aux].GetXaxis().GetBinLowEdge(ix)
                    lowY = histos3D[aux].GetYaxis().GetBinLowEdge(iy)
                    lowZ = histos3D[aux].GetZaxis().GetBinLowEdge(iz)
                    highX = histos3D[aux].GetXaxis().GetBinUpEdge(ix)
                    highY = histos3D[aux].GetYaxis().GetBinUpEdge(iy)
                    highZ = histos3D[aux].GetZaxis().GetBinUpEdge(iz)

                    titleHist = '%.1f< %s <%.1f,  %.1f< %s < %.1f, %.1f< %s < %.1f'%(lowX, nameX, highX, lowY, nameY, highY, lowZ, nameZ, highZ)

                    binTot = histos3D[aux].GetBin(ix, iy, iz)

                    h1 = TH1F('h%s'%(binTot), titleHist, len(cutList)-1, binning)

                    for ibin in xrange(0, (len(cutList)-1)):

                        if(ibin<(len(cutList)-2)):

                            ikey1 = cutList[ibin]
                            ikey2 = cutList[ibin+1]

                            val1 = float(histos3D[ikey1].GetBinContent(ix, iy, iz))
                            val2 = float(histos3D[ikey2].GetBinContent(ix, iy, iz))
                            err1 = float(histos3D[ikey1].GetBinError(ix, iy, iz))
                            err2 = float(histos3D[ikey2].GetBinError(ix, iy, iz))
                            prob = TMath.fabs(val2-val1)

                            if self.ErrorType.lower() == 'poisson':
                                N = self.CalcN(err1, val1)
                                if N == 0:
                                    err = 0.0
                                else:
                                    err = TMath.sqrt(prob*(1-prob))/TMath.sqrt(N)
                            else:
                                err = TMath.sqrt(err1*err1+err2*err2)


                        else:
                            ikey1 = cutList[ibin]
                            val1 = float(histos3D[ikey1].GetBinContent(ix, iy, iz))
                            val2 = 0.0
                            prob = TMath.fabs(val2-val1)
                            err = float(histos3D[ikey1].GetBinError(ix, iy, iz))
                            if error.lower() == 'poisson':
                                N = self.CalcN(err, val1)
                                if N == 0:
                                    err = 0.0
                                else:
                                    err = TMath.sqrt(prob*(1-prob))/TMath.sqrt(N)
                            

                        h1.SetBinContent(ibin+1, prob)
                        h1.SetBinError(ibin+1, err)
                    h1.Write()


        return 

