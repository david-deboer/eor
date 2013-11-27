###Imports
import matplotlib.pyplot as plt
import math
import sys
import string

import eorview
evu = eorview.view()
import eorio
sys.path.append('/Users/daviddeboer1/Documents/Code/cosmo')
import nedclass as ned
sys.path.append('/Users/daviddeboer1/Documents/Code/lib')
import astrolib
al = astrolib.astrolib()

###Various constants
edgeLeft=0
ctr=1
edgeRight=2
plotSens = 1
plotBoost = 2
plotOther = 3
Units = {'Hz':1.0, 'kHz':1.0E3, 'MHz':1.0E6, 'GHz':1.0E9, 'sec':1.0, 'min':60.0, 'hr':3600.0, 'm':1.0, 'cm':0.01, 'km':1000.0}

###Everything else is in the class bleor
class eor:
    def __init__(self, config='13ele.dat', rotate=True, Nplat=1, Tsys=50.0, Npol=1.0, dlogk=0.2, tobs=1000.0, tobsUnit='hr', z=6.0, restFreq=230.0, freqUnit='GHz',\
                 BW=1.0, bwUnit='GHz', Nch=128, HAmin=-5.0, HAmax=5.0, HAunit='hr', tint=0.5, tintUnit='hr', dec=25.0, lat=34.0, boxcar=10.,\
                 modelType='auto:file', H0=72.0):
        """eorpy runs the baseline stuff for EoR power spectrum measurement.  Needs nedclass, eorview, pwrspec, astrolib.  Removes PGAFF f-o-m stuff
        December 7, 2012 version to do more clean-up and include modeling in addition to the data files.
        March 3, 2012 version makes it more general and does some clean-up.
        March 6, 2012 version combines sens.py (except the 'special' functions that will be in sensUtil.py)
        September 5, 2012 version 'getting back into it' and making tweaks.  There is no sens.py nor sensUtil.py.
        config='Npp,Nplat,Dant[,bmin]' or 'filename' (usual)
        For 'array' consider that earth is the platform.
        This version handles a full track, as opposed to individual integrations within a track."""

        ###-----This is all superceded by the config files and the full calculation - but remains to verify the calculations and used in the "MMQboost"-----###
        #########First set arrays for the "classic" DACOTA configuration
        self.B2coeff = [1.0, 3.0, 4.0, 7.0, 9.0, 12.0, 13.0, 16.0, 19.0, 21.0, 25.0, 27.0, 28.0, 31.0, 36.0, 37.0, 39.0, 43.0, 48.0, 49.0, 49.0, 52.0, 57.0, 61.0, 63.0, 64.0, 67.0, 73.0, 75.0, 76.0, 79.0, 81.0, 84.0, 91.0, 100.0]
        self.nocoeff = [3.0, 3.0, 3.0, 6.0, 3.0,  3.0,  6.0,  3.0,  6.0,  6.0,  3.0,  3.0,  6.0,  6.0,  3.0,  6.0,  6.0,  6.0,  3.0,  6.0,  3.0,  6.0,  6.0,  6.0,  6.0,  3.0,  6.0,  6.0,  3.0,  6.0,  6.0,  3.0,  6.0,  6.0,   3.0]
        # The "arr" prefix lists the array config parameter for 'shortcut' dacota configurations
        self.arrNpp =   [ 7, 13, 19, 31, 37, 55, 61, 85, 91]       #Number of elements per platform    
        self.arrNplat = [37, 20, 14,  8,  7,  4,  4,  3,  3]       #Number of platforms in "full" system
        self.arrNuv = [ [4,2,1],
                        [8,6,5,2,2,1],
                        [14,10,9,6,4,3,2,1],
                        [24,20,17,14,12,9,8,7,4,4,2,2,1],
                        [30,24,23,18,16,13,12,9,8,6,4,4,3,2,1],
                        [46,40,37,32,28,25,24,21,18,16,14,12,11,10, 7, 6, 6, 4, 3, 2, 2, 1],
                        [52,44,43,36,34,29,28,25,22,20,16,16,15,12, 9,10, 8, 6, 5, 4, 4, 3, 2, 0, 0, 1],
                        [74,66,63,56,52,47,46,41,38,36,32,30,29,26,23,22,20,18,15,14,14,13,10, 8, 8,7,6,4,4,3,2,2,1],
                        [80,70,69,60,58,51,50,47,42,40,36,34,33,30,25,26,24,20,19,18,16,15,12,12,10,9,8,6,6,5,4,4,3,2,1] ] #Number of baselines per scale
        self.arrNk = []   # Nk refers to the last Bcoeff used to compute baseline spacings (so it uses that index number)
        for i in range(len(self.arrNuv)):
            self.arrNk.append(len(self.arrNuv[i]))
        ###-----End of MMQboost

        ##########initialize any variable as needed (but be careful vis a vis defaults etc)
        self.kbin=[ [], [], [] ]
        self.kmin=-1
        self.kmax=-1
        self.DantUnit = 'm'
        self.bUnit = 'm'
        self.arrayType = 'NA'
        self.elLimit = 10.0
        self.dkcoher = None
        self.dkcoherScale = 1.0
        self.cohMode = 'normal'  # vs full or un
        self.shadow = False
        self.dacotaShortcut = False
        ##########Now set system defaults (do this order to not screw up updateWeb)
        self.setModel(modelType=modelType,H0=H0)
        self.setData(boxcar=boxcar,coherentTracks=True,outputType='sq')
        self.setTrack(HAmin=HAmin,HAmax=HAmax,HAunit=HAunit,tobs=tobs,tobsUnit=tobsUnit,tint=tint,tintUnit=tintUnit,Npt=1,dec=dec,lat=lat)
        self.setPlot(plotModel=False,plotSteps=False,plotBins=True,plotVersion='errbar',plotMMQ=False,plotLines=True,showLegend=True,ltype='-')
        self.setSys(config=config,rotate=rotate,Nplat=Nplat,Npol=Npol,dlogk=dlogk,z=z,Tsys=Tsys,restFreq=restFreq,freqUnit=freqUnit,BW=BW,bwUnit=bwUnit,Nch=Nch)
        evu.updateWeb(self)
    

    def arraySensitivity(self,color=None,label=None,linewidth=2,plotExtraStuff=False,kukvdot='k.'):
        """Does the sensitivity analysis - used to be self.array, and also subsumed self.dacota
             - arrayFile - the config file, None is self.arrayFile
             - plotMMQ - plot Matt's expression along with calculated value
             - dacotaShortcut - this sidesteps the iterations by assuming isoplanactic coherent integration
             it has a whole bunch of extraneous arguments to make the figures per the proposal..."""
        # setup and initialize
        hu = Units[self.HAunit]/Units['hr']    # makes sure HA in hours
        tu = Units[self.tintUnit]/Units['hr']  # makes sure tint in hours
        HA = self.HAmin*hu
        tint = self.tint*tu
        HAmax = self.HAmax*hu
        oldValues = []
        dacotaShortcut = self.dacotaShortcut
        if dacotaShortcut and self.arrayType == 'platform':
            oldValues = [self.HAmin,self.HAmax,self.tint,self.tintUnit]
            HA=0.0
            HAmax = 0.01*hu
            tint = self.tobs*Units[self.tobsUnit]/Units['hr']
            self.setTrack(HAmin=0.0, HAmax=HAmax, tint=tint, tintUnit=self.tobsUnit)
        if self.plotSteps or self.plotBins:
            oldBoxcar = self.boxcar
            self.setData(boxcar=1)
        self.dkcoher=None
        self.setkbins()
        self.setCoherenceGrid()
        self.HA = []
        # loop over HA to set the coherence grid and count points
        i = 0
        while HA <= 1.01*HAmax:
            self.HA.append(HA)
            self.computeK(HA)
            if plotExtraStuff:
                plt.figure(100)
                plt.plot(self.ku,self.kv,kukvdot)
                plt.figure(101)
                for kperp in self.kperp:
                    kperpvalarray = []
                    for i in range(len(self.kpar)):
                        kperpvalarray.append(kperp)
                    plt.plot(kperpvalarray,self.kpar,kukvdot)
            self.coherenceGrid()
            HA+=tint
            i+=1
            if (i % 5) == 0:
                print ' '
        print ' '
        # take the populated coherence grid and calculate the boost per k bin for a track
        self.calcBoost()
        print ' '
        # compute the sensitivity for a track
        self.d = self.computeDValues(self.kbin[ctr],self.boost)
        # compute for total number of tracks
        self.d = self.computeTracks(self.d)
        evu.updateWeb(self)

        # plot the data
        color = 'b'
        eorio.plotBoost(self,plotBoost)
        eorio.plotSensitivity(self,plotSens,label=label,color=color)
                
        if dacotaShortcut and self.arrayType == 'platform':  #Go ahead and reset the values back
            self.setTrack(HAmin=oldValues[0], HAmax=oldValues[1], tint=oldValues[2], tintUnit=oldValues[3])
        if self.plotSteps or self.plotBins:
            self.setData(boxcar=oldBoxcar)
        self.HA = []

        return i

    def set2paper(self):
        self.setSys(config='paper9-12.dat',BW=0.01,Tsys=800.0,restFreq=1.42)
        self.setTrack(lat=-30.,dec=-25,elLimit=40.0)
        self.setData(dkcoherScale=100.0)
    def set2sza(self,attempt='now'):
        """now vs best vs karto"""
        print 'Setting parameters for sza for '+attempt

        if self.h<0.9:
            kh = 1.0
        else:
            kh = 0.7

        if attempt == 'now':
            self.setSys(config='sza_L.ant', Tsys=66., BW=1.0, Nch=30, Npol=2., restFreq=230., z=6.)
            self.setTrack(tobs=1000.)
            self.setData(coherentTracks=False)
        elif attempt == 'best':
            self.setSys(config='sza_L.ant', Tsys=66., BW=1.2, Nch=33, Npol=2., restFreq=230., z=6.)
            self.setTrack(tobs=2000.)
            self.setData(coherentTracks=True)
        elif attempt == 'karto':
            if self.z>5.0:
                fp = open('infiles/szaObs6.dat','r')
            else:
                fp = open('infiles/szaObs3.dat','r')
            x = []
            y = []
            for line in fp:
                dat = line.split()
                k = pow(10.0,float(dat[0]))/kh
                d2= 1.0E12*pow(10.0,float(dat[1]))*math.sqrt(11.0)
                x.append(k)
                y.append(d2)
            kstart = pow(10.0,-0.1)*x[0]
            kstop = pow(10.0,-0.1)*x[-1]
            self.setData(boxcar=1.0)
            self.setkbins(kstart=kstart,kstop=kstop)
            #errx,erry=self.errorBars()
            #plt.errorbar(x,y,yerr=None,xerr=errx,fmt='_',color='b')
            label = 'SZA z='+str(self.z)+' (4 fields)'
            plt.loglog(x,y,color='b',label=label)
            for i in range(len(y)):
                y[i] = y[i]/math.sqrt(11.0)
            #plt.errorbar(x,y,yerr=None,xerr=errx,fmt='_',color='c')
            label = 'SZA z='+str(self.z)+' (44 fields)'
            plt.loglog(x,y,color='c',label=label)            
            self.setData(boxcar=10.0)
        else:
            print 'Not a valid attempt...'
        return attempt
    def set2amiba(self,z=6):
        """6 vs 3"""
        if z<4.:
            self.setSys(config='amiba.dat', Tsys=50., BW=2.0, Nch=256, Npol=2., restFreq=115., z=z)
            self.setTrack(tobs=1500.0, Npt=1)
        else:
            self.setSys(config='amiba.dat', Tsys=50., BW=4.0, Nch=256, Npol=2., restFreq=230., z=z)
            self.setTrack(tobs=1500.0, Npt=1.)
        self.setData(coherentTracks=True)

    def computeTracks(self,d2):
        dd = list(d2)
        if self.outputType!='sq':
            for i in range(self.nbin):
                dd[i]=dd[i]*dd[i]
                
        if self.coherentTracks==True:
            for i in range(self.nbin):
                dd[i]/=self.nTrack
        else:
            for i in range(self.nbin):
                dd[i]/=math.sqrt(self.nTrack)
    
        if self.outputType!='sq':
            for i in range(nbin):
                dd[i] = math.sqrt(dd[i])
        return dd


    def setTrack(self,HAmin=None,HAmax=None,HAunit=None,tobs=None,tobsUnit=None,tint=None,tintUnit=None,Npt=None,\
                 dec=None,lat=None,elLimit=None,shadow=None):
        """Sets track parameters:  (HAmin, HAmax, HAstep) [hrs], dec [deg]"""

        if HAunit==None:
            HAunit=self.HAunit
        else:
            self.HAunit = HAunit
        if tobsUnit==None:
            tobsUnit=self.tobsUnit
        else:
            self.tobsUnit=tobsUnit
        if tintUnit==None:
            tintUnit=self.tintUnit
        else:
            self.tintUnit=tintUnit

        if HAmin!=None:
            print 'TRACK: Setting HAmin to '+str(HAmin)+' '+self.HAunit
            self.HAmin=float(HAmin)
        if HAmax!=None:
            print 'TRACK: Setting HAmax to '+str(HAmax)+' '+self.HAunit
            self.HAmax=float(HAmax)
        if tobs!=None:
            print 'TRACK: Setting tobs to '+str(tobs)+' '+self.tobsUnit
            self.tobs=float(tobs)
        if tint!=None:
            if tint=='auto':
                tint = (self.Dant*Units[self.DantUnit])/(self.bmax*Units[self.bUnit])/2.0 # this is 1/4 the dwell time
                self.tintUnit = 'hr'
            print 'TRACK: Setting tint to '+str(tint)+' '+self.tintUnit
            self.tint=float(tint)
        if Npt!=None:
            print 'TRACK: Setting Npt to '+str(Npt)
            self.Npt=Npt
        if dec!=None:
            print 'TRACK: Setting dec to '+str(dec)+' degrees'
            self.dec = dec
        if lat!=None:
            print 'TRACK: Setting lat to '+str(lat)+' degrees'
            self.lat = lat

        self.trackTime = (self.HAmax-self.HAmin)*Units[self.HAunit]
        self.nint = int(math.ceil( self.trackTime / (self.tint*Units[self.tintUnit])))
        if self.tint==self.tobs:
            self.nTrack=1
        else:
            self.nTrack = int(math.ceil( (self.tobs*Units[self.tobsUnit])/self.trackTime))

        if elLimit!=None:
            print 'TRACK: Setting elLimit to '+str(lat)+' degrees'
            self.elLimit = elLimit

        if shadow!=None:
            print 'TRACK:  Setting shadow to '+str(shadow)
            print '        (Note that shadow==True ignores shadowed antennas'
            self.shadow = shadow

        Az, El = al.eq2hor(15.0*self.HAmin,self.dec,self.lat)
        if El < self.elLimit:
            while El < self.elLimit:
                self.HAmin+=0.001
                if self.HAmin > 12.0:
                    print "Error resetting HAmin.  Using %.2f" % self.HAmin
                    return
                Az, El = al.eq2hor(15.0*self.HAmin,self.dec,self.lat)
            print "Reseting HAmin to "+str(self.HAmin)
        Az, El = al.eq2hor(15.0*self.HAmax,self.dec,self.lat)
        if El < self.elLimit:
            while El < self.elLimit:
                self.HAmax-=0.001
                if self.HAmax < -12.0:
                    print "Error resetting HAmax.  Using %.2f" % self.HAmax
                    return
                Az, El = al.eq2hor(15.0*self.HAmax,self.dec,self.lat)
            print "Reseting HAmax to "+str(self.HAmax)
        evu.updateWeb(self)

    def setData(self,boxcar=None,dkcoher=None,dkcoherScale=None,cohMode=None,\
                coherentTracks=None,outputType=None,useCoherentGridding=None):
        """Sets parameters related to data generation/display"""
        if boxcar!=None:
            print 'DATA:  Setting boxcar to ',boxcar
            self.boxcar = float(boxcar)
        if dkcoher!=None:
            print 'DATA:  Setting dkcoher to ',dkcoher
            self.dkcoher = float(dkcoher)
        if dkcoherScale!=None:
            print 'DATA:  Setting dkcoherScale to ',dkcoherScale
            self.dkcoherScale = float(dkcoherScale)
            if isinstance(self.dkcoher,float):
                self.dkcoher*=self.dkcoherScale
                print '\treset dkcoher by new scale:  ',self.dkcoher
        if coherentTracks!=None:
            print 'DATA:  Setting coherentTracks to ',coherentTracks
            self.coherentTracks=coherentTracks
        if outputType!=None:
            print 'DATA: Setting outputType to ',outputType
            self.outputType=outputType
        if cohMode!=None:
            print 'DATA: Setting cohMode to ',cohMode
            self.cohMode = cohMode
        evu.updateWeb(self)

    def setSys(self,config=None,rotate=None,Npp=None,Nplat=None,Npol=None,dlogk=None,dkpk=None,Dant=None,bmin=None,z=None,Tsys=None,restFreq=None,obsFreq=None,freqUnit=None,\
               BW=None,bwUnit=None,Nch=None,Npt=None,lat=None,dacotaShortcut=None):
        """sets config parameters.  List in dispSys()"""
        if config!=None:
            print 'SYS: Setting config to ',config
            self.config = config
            self.dkcoher = None
            v = config.split(',')
            if len(v) == 1:
                self.arrayFile = config
                self.readArrayFile(arrayFile=config)  #gets Npp, Nplat, Dant, bmin
            else:
                self.arrayFile = 'generate'
                self.Npp = int(v[0])
                self.Nplat = int(v[1])
                self.Dant = float(v[2])
                if len(v)==4:
                    self.bmin = float(v[3])
                else:
                    self.bmin = self.Dant
        if rotate!=None:
            print 'SYS: Setting rotate to '+str(rotate)
            if self.arrayType!='platform':
                print 'SYS: note rotate only for platforms, not '+self.arrayType
            self.rotate = rotate
        if dacotaShortcut!=None:
            print 'SYS: Setting dacotaShortcut to ',str(dacotaShortcut)
            if self.arrayType!='platform':
                print 'SYS: note dacotaShortcut only for platforms, not '+self.arrayType
            self.rotate = rotate            
        if Npp!=None:
            print 'SYS: Setting Npp to ',Npp
            self.Npp=float(Npp)
        if Nplat!=None:
            print 'SYS: Setting Nplat to ',Nplat
            self.Nplat=float(Nplat)
        if Npol!=None:
            print 'Sys: Setting Npol to ',Npol
            self.Npol = float(Npol)
        try:  # if these are all defined compute Nant and nBaseline
            self.Nant = self.Npp*self.Nplat
            self.nBaseline = int((self.Npp)*(self.Npp-1.0)/2.0)
        except:
            print 'Nant and nBaseline still undefined (check Npp, Nplat)'
        if Dant!=None:
            print 'SYS: Setting Dant to ',Dant
            self.Dant = float(Dant)
        if bmin!=None:
            print 'SYS: Setting bmin to ',bmin
            self.bmin = float(bmin)
        if Tsys!=None:
            print 'SYS: Setting Tsys to '+str(Tsys)
            self.Tsys=float(Tsys)
        if dkpk!=None:
            print 'SYS: Setting dkpk to '+str(dkpk),
            self.dkpk = float(dkpk)
            self.dlogk = math.log10(dkpk+1.0)
            print '(dlogk = '+str(self.dlogk)+')'
        if dlogk!=None:
            print 'SYS: Setting dlogk to '+str(dlogk),
            self.dlogk = float(dlogk)
            self.dkpk = math.pow(10.0,dlogk) - 1.0
            print '(dkpk = '+str(self.dkpk)+')'
            if dkpk!=None:
                print "BTW, you should just set dlogk or dkpk, not both"
        if z!=None:
            print 'SYS: Setting z to ',z
            self.z=float(z)
            if restFreq==None and obsFreq==None:
                self.obsFreq = self.restFreq/(1.+self.z)
                print '\tSetting obsFreq to '+str(self.obsFreq)+' '+self.freqUnit
        if freqUnit!=None:
            print 'SYS: Setting freqUnit to '+freqUnit
            self.freqUnit=freqUnit
        if restFreq!=None:
            print 'SYS: Setting restFreq to '+str(restFreq)+' '+self.freqUnit,
            self.restFreq=float(restFreq)
            self.obsFreq = self.restFreq/(1.+self.z)
            print '(obsFreq = '+str(self.obsFreq)+' '+self.freqUnit+')'
        if obsFreq!=None:
            print 'SYS: Setting obsFreq to '+str(obsFreq)+' '+self.freqUnit,
            self.obsFreq=float(obsFreq)
            self.restFreq = self.obsFreq*(1.+self.z)
            print '(restFreq = '+str(self.restFreq)+' '+self.freqUnit+')'
            if restFreq!=None:
                print "BTW, you should just set restFreq or obsFreq, not both"
        self.obsWavelength = 3.0E8/(self.obsFreq*Units[self.freqUnit])  # in meters
        if bwUnit!=None:
            print 'SYS: Setting bwUnit to '+bwUnit
            self.bwUnit = bwUnit
        if BW!=None:
            print 'SYS: Setting BW to '+str(BW)+' '+self.bwUnit
            self.BW=float(BW)
        if Nch!=None:
            print 'SYS: Setting Nch to ',Nch
            self.Nch=int(Nch)
        self.nKline = int(self.nBaseline*self.Nch)
        evu.updateWeb(self)

    def setPlot(self,plotModel=None, plotSteps=None, plotBins=None, plotVersion=None, plotMMQ=None, plotLines=None, showLegend=None, ltype=None):
        if plotModel!=None:
            print 'PLOT:  Setting plotModel to ',plotModel
            self.plotModel = plotModel
        if plotSteps!=None:
            print 'PLOT:  Setting plotSteps to ',plotSteps
            self.plotSteps = plotSteps
        if plotBins!=None:
            print 'PLOT:  Setting plotBins to ',plotBins
            self.plotBins = plotBins
        if plotMMQ!=None:
            print 'PLOT:  Setting plotMMQ to ',plotMMQ
            self.plotMMQ = plotMMQ
        if plotLines!=None:
            print 'PLOT:  Setting plotLines to ',plotLines
            self.plotLines = plotLines
        if ltype!=None:
            print 'PLOT:  Setting linetype to '+ltype
            self.ltype = ltype
        if plotVersion!=None:
            print 'PLOT:  Setting plotVersion to '+plotVersion
            print '       options are errbar or senscurve'
            self.plotVersion = plotVersion
        if showLegend!=None:
            print 'PLOT:  Setting showLegend to '+str(showLegend)
            self.showLegend = showLegend
        evu.updateWeb(self)

    def setModel(self,modelType=None,H0=None):
        if modelType!=None:
            print "MODEL:  Setting modelType to ",modelType
            self.modelType = modelType
        if H0!=None:
            print "MODEL:  Setting H0 to ",H0
            self.H0 = H0
            self.h = H0/100.0
            self.ned = ned.ned(H0=H0)
        evu.updateWeb(self)

    def scaledPerformancePlot(self,scaleFactor=1.5,dsub=0.1,taper=11.0):
        # initialize metric parameters
        Tsys_start = self.Tsys
        E = pow(10.0,taper/10.0)
        Cb = -1.0*math.log(E)/(1.0-math.sqrt(E))
        print 'Cb = ',Cb
        self.rescaleArray(0.5/scaleFactor)
        self.Dant_scale = []
        self.kmin_scale = []
        self.dmin_scale = []
        self.eff_scale = []
        scaleFactor = 1.5
        while self.Dant < 6.0:
            self.rescaleArray(scaleFactor)
            Dratio = dsub/self.Dant
            eff = (1.0 - Cb*(1.0 + 4.0*math.sqrt(1.0-Dratio))*Dratio**2)**2
            self.setSys(Tsys=Tsys_start/eff)
            s = 'D=%.1f' % (self.Dant)
            self.arraySensitivity(label=s)
            self.Dant_scale.append(self.Dant)
            self.kmin_scale.append(self.kmin)
            self.dmin_scale.append(self.d[0])
            self.eff_scale.append(eff)
        print 'end'
    def rescaleArray(self,scale=1.0):
        self.scale = scale
        Dnew = self.Dant*scale
        self.setSys(Dant=Dnew)
        for i in range(self.Npp):
            self.N[i]*=scale
            self.E[i]*=scale
            self.U[i]*=scale
        self.calcStats()

    def readArrayFile(self,arrayFile=None, scale=None,  nameIncluded=1, path='/Users/daviddeboer1/Documents/Code/antutil/telescopeData/configs/',verbose=False):
        """Reads a configuration file of form:  [name] N E [U]
           For type 'platform'
               Npp=Npp, Nplat=Nplat
           For type 'array', the earth is the (single) platform
               Npp=number of antennas, Nplat=1
           Nant = Npp*Nplat
           nBaseline = Npp*(Npp-1)/2"""
        if arrayFile=='generate':
            print "Can't read file in 'generate' mode."
            return 0
        try:
            fp = open(path+arrayFile,"r")
        except IOError:
            print path+arrayFile+' not found'
            return 0
        print '    Reading '+arrayFile
        self.elementName = []
        self.N = []
        self.E = []
        self.U = []
        Npp = 0
        foundExcl=False
        for line in fp:
            if line[0] == '!':
                foundExcl=True
                line = line.strip('!')
                v=line.split()
                if v[0]=='name':
                    nameIncluded=1
                else:
                    nameIncluded=0
                self.arrayType=v[1]
                self.Dant = float(v[2])
                self.DantUnit = v[3]
                print '    arrayType = '+self.arrayType
                print '    Dant = '+str(self.Dant)+' '+self.DantUnit
                if len(v)==5:
                    self.scale=float(v[4])
                    print '    Scale by '+str(self.scale)
                else:
                    self.scale=1.0
                if scale != None:
                    self.scale = scale # override if fed this in argument list
                continue
            if line[0] == '#':
                continue
            if foundExcl!=True:   #past the header
                print arrayFile+"didn't have an antenna size or scale factor as first line"
                self.Dant=1.2
                self.DantUnit='m'
                self.scale=1.0
                self.arrayType = 'array'
                print "\tassuming Dant = "+str(self.Dant)+" "+self.DantUnit
                print "\tassuming scale = "+str(self.scale)
                foundExcl=True   #so we don't do this everytime
            v = line.split()
            if len(v) == 0:
                continue
            if nameIncluded:
                self.elementName.append(v[0])
            else:
                self.elementName.append('ant'+str(Npp))
            self.N.append(float(v[0+nameIncluded])*self.scale)
            self.E.append(float(v[1+nameIncluded])*self.scale)
            if len(v) > 2+nameIncluded:
                self.U.append(float(v[2+nameIncluded]))
            else:
                self.U.append(0.0)
            Npp+=1
        self.bUnit = self.DantUnit
        if self.arrayType=='platform':
            print "    Platform"
            self.Npp = Npp
        else:
            print "    Array"
            self.Npp = Npp
            self.Nplat = 1
        fp.close()
        if verbose:
            for i in range(Nant):
                print self.elementName[i],self.N[i],self.E[i],self.U[i].self.bUnit
        self.calcStats()
        return Npp

    def computeUV(self,coord='uv'):
        """Computes the uv distribution.  Can select uv coordinates or kperp. (uv or k)
           Note that this does the symmetric version...have j in range(i) for unique.  This is not used in the standard
           progression (computeK() is).  Append 'd' if want differential.  Append 'w' if want all 'Hermitian' points.
           Normals options are then:  k, kd, kw, kdw, uv, uvd, uvw, uvdw"""
        X = self.calcX(self.z)
        scale = 1.0
        if 'k' in coord:
            scale = 2.0*math.pi/X
        else:
            scale = 1.0

        if 'd' in coord:
            kind = 'differential'
        else:
            kind = 'absolute'

        if 'w' in coord:
            allPts = True
            mall = range(self.Npp)
        else:
            allPts = False
        print kind+'  '+str(scale)+' use allPts = '+str(allPts)
        hu = Units[self.HAunit]/Units['hr']    # make sure HA in hours
        tu = Units[self.tintUnit]/Units['hr']  # make sure tint in hours
        self.u = []
        self.v = []
        self.w = []
        self.uvperp = []
        self.HA = []
        HA = self.HAmin*hu
        nBaseline = 0
        while HA <= self.HAmax*hu:
            print "HA = %.4f" % HA
            nBaseline = 0
            self.HA.append(HA)
            for i in range(self.Npp):
                if allPts==False:
                    mall = range(i)
                for j in mall:
                    if kind == 'absolute':
                        u, v, w = self.calcuv(i,j,HA)
                    else:
                        u, v, w = self.calcduv(i,j,HA)
                    # kperp, kpar = ditherData(kperp,kpar)
                    self.u.append(u*scale)
                    self.v.append(v*scale)
                    self.w.append(w*scale)
                    self.uvperp.append( math.sqrt(u*u + v*v)*scale)
                    nBaseline += 1
            HA+= self.tint*tu
        print "Compare "+str(nBaseline)+" with "+str(self.nBaseline)
        self.nBaseline = nBaseline
	return nBaseline

    def computeDishUVsize(self,coord='uv'):
        """Draws a circle representing the dish diameter in wavelengths"""
        X = self.calcX(self.z)
        scale = 1.0
        if 'k' in coord:
            scale = 2.0*math.pi/X
        self.dishSizeX = []
        self.dishSizeY = []
        for th in range(361):
            dsw = scale*self.Dant/self.obsWavelength
            thr = th*math.pi/180.0
            self.dishSizeX.append(dsw*math.cos(thr))
            self.dishSizeY.append(dsw*math.sin(thr))
        return th

    def computeK(self,HA=0.0):
        """Computes the k distribution (in u, v) at an HA."""
        print "HA = %.4f\t" % HA,
        X = self.calcX(self.z)
        scale = 2.0*math.pi/X
        f = self.restFreq*Units[self.freqUnit]
        Y = self.calcY(self.z,f)
        B = self.BW*Units[self.bwUnit]
        Nch = self.Nch
        self.ku    = []  # only has u,v plane ones (nBaseline)
        self.kv    = []  # only has u,v plane ones (nBaseline)
        self.kperp = []  # keep this just for plotting purposes
        nKline = 0
        for i in range(self.Npp):
            for j in range(i):
                u, v, w = self.calcuv(i,j,HA)
                if self.shadow==True and self.obsWavelength*math.sqrt(u*u + v*v) < self.Dant:
                    print 'Shadowed baseline not used'
                    continue
                nKline+=1
                ku = u*scale
                kv = v*scale
                self.ku.append(ku)
                self.kv.append(kv)
                a = 1.0
                if ku < 0.0:
                    a = -1.0
                self.kperp.append(a*math.sqrt(ku**2 + kv**2))
        nKline*=self.Nch
	return nKline

    def calcKextremes(self):
        """Computes the k extremes over a full track: kmin, kmax, kperpmax."""

        X = self.calcX(self.z)
        scale = 2.0*math.pi/X
        f = self.restFreq*Units[self.freqUnit]
        Y = self.calcY(self.z,f)
        B = self.BW*Units[self.bwUnit]
        Nch = self.Nch
        self.kmin  = 1.0E9
        self.kmax  = 0.0
        self.kperpmax = 0.0
        hu = Units[self.HAunit]/Units['hr']    # make sure HA in hours
        tu = Units[self.tintUnit]/Units['hr']  # make sure tint in hours
        if self.arrayType == 'platform':
            HA = 0.0
            HAmax = 0.01*hu
            tint = 1.0*tu
        else:
            HA = self.HAmin*hu
            HAmax = self.HAmax*hu
            tint = self.tint*tu
        nKline = 0
        nBaseline = 0
        while HA <= HAmax:
            nKline = 0
            nBaseline = 0
            for i in range(self.Npp):
                for j in range(i):
                    u, v, w = self.calcuv(i,j,HA)
                    if u is False:
                        continue
                    ku = u*scale
                    kv = v*scale
                    kperp = math.sqrt(ku**2. + kv**2.)
                    if kperp > self.kperpmax:
                        self.kperpmax = kperp
                    nBaseline += 1
                    for n in range(int(Nch)):
                        kpar  = n*2.0*math.pi/(Y*B)
                        k = math.sqrt(kperp**2. + kpar**2.)
                        if k > self.kmax:
                            self.kmax = k
                        if k < self.kmin:
                            self.kmin = k
                        nKline += 1
            HA+=tint
        print "nKline = "+str(nKline)+" ("+str(self.nKline)+")"
        self.nKline = nKline
        print "nBaseline = "+str(nBaseline)+" ("+str(self.nBaseline)+")"
        self.nBaseline = nBaseline
        print "kmin = ",self.kmin
        print "kmax = ",self.kmax
        print "kperpmax = ",self.kperpmax
        #print "---from sens"
        #self.calckmin()  This assume kmin set by Dant
        #print "---"
	return nKline

    def setkbins(self,kstart=None,kstop=None):
        """Sets up the kbins to use.  Calculates the natural ranges for non-specified kstart/kstop.
           self.kbin[edgeLeft] self.kbin[edgeRight] are the edges and self.kbin[ctr] are the centers.
           So, plot against kbin[ctr].
           boxcar does a boxcar version.  boxcar=1 is 1 per interval, etc"""
        boxcar = self.boxcar
        if kstart==None or kstop==None:
            self.calcKextremes()
        if kstart==None:
            kstart = 0.98*self.kmin
        if kstop==None:
            kstop = self.kmax
        print 'kbin from '+str(kstart)+' to '+str(kstop)+' with dlogk = '+str(self.dlogk)
        self.kbin[edgeLeft] = []
        self.kbin[ctr] = []
        self.kbin[edgeRight] = []
        kLeft = kstart
        i = 0
        while kLeft<kstop:
            self.kbin[edgeLeft].append(kLeft)
            #k = kstart*math.pow(10.0,(i+1.0)*self.dlogk)
            widthk = kLeft*(math.pow(10.0,self.dlogk)-1.0)
            kRight=kLeft+widthk
            self.kbin[edgeRight].append(kRight)
            deltak = kLeft*(math.pow(10.0,self.dlogk/boxcar) - 1.0)
            kLeft+=deltak
            i+=1
        self.nbin = i
        for i in range(self.nbin):
            k = (self.kbin[edgeLeft][i] + self.kbin[edgeRight][i])/2.0
            self.kbin[ctr].append(k)
            #print "Bin "+str(i)+" :  "+str(self.kbin[edgeLeft][i])+" < "+str(k)+" < "+str(self.kbin[edgeRight][i])
        return i

    def setCoherenceGrid(self):
        """initialize the coherent sky grid"""        
        self.kCohBins = []      # 3-D matrix of mode counts
        self.kCohPerpVal = []   # Linear array of kperp-values of bins
        self.kCohParVal = []    # Linear array of kpar-values of bins
        self.boost = []         # Reset the boost array
        # set dkcoher, assuming Gaussian beam
        if self.dkcoher == 0.0 or not isinstance(self.dkcoher,float):
            print 'dkcoher:  was',self.dkcoher,
            X = self.calcX(self.z)
            coeff = self.dkcoherScale * 4.0 * math.log(1.2)   # this is 1/1.2=80% point for Gaussian beam
            Dant = self.Dant*Units[self.DantUnit]/Units['m']  # convert to meters since obsWavelength always in meters
            self.dkcoher = coeff*( Dant/self.obsWavelength ) / X
            print '  set to %.4f (dkcoherScale=%.4f)' % (self.dkcoher, self.dkcoherScale)
        dkcoher=self.dkcoher
        f = self.restFreq*Units[self.freqUnit]
        Y = self.calcY(self.z,f)
        B = self.BW*Units[self.bwUnit]
        self.nCoherBin = int(2.01*math.ceil(self.kperpmax/dkcoher))
        nCoherBin = self.nCoherBin
        Nch = self.Nch
        for i in range(nCoherBin):
            iv = i - nCoherBin/2
            kcohbin = (iv+0.5)*dkcoher
            self.kCohPerpVal.append(kcohbin)
            self.kCohBins.append([])
            for j in range(nCoherBin):
                self.kCohBins[i].append([])
                for n in range(int(Nch)):
                    self.kCohBins[i][j].append([])
        for n in range(int(Nch)):
            kpar  = n*2.0*math.pi/(Y*B)
            self.kCohParVal.append(kpar)
        evu.updateWeb(self)
        self.kpar = self.kCohParVal  # rename and save just for plotting purposes
        return i

    def calcBoost(self):
        """This counts the number of measurements and puts them in k bins"""
        self.sp=[]
        dkcoher = self.dkcoher
        nCoherBin = self.nCoherBin
        Nch = self.Nch
        if self.cohMode == 'un':
            cohExp = 1.0
        else:
            cohExp = 2.0
        print '\nCalculating boost factor in cohMode= '+self.cohMode+' cohExp=',cohExp
        for ibin in range(self.nbin):
            print str(ibin)+', ',
            if (ibin+1)%15 == 0:
                print ' '
            klo = self.kbin[edgeLeft][ibin]
            khi = self.kbin[edgeRight][ibin]
            ss = 0.0
            for i in range(nCoherBin):
                for j in range(nCoherBin):
                    if len(self.kCohBins[i][j][0]) == 0 :  # don't bother - nothing here...
                        continue
                    nPerpOcc = len(self.kCohBins[i][j][0])
                    #print '%d/%d (%d,%d): %d' % (ibin,self.nbin,i,j,nPerpOcc)
                    # pre-compute decorrelation factor (so you don't do it Nch times)
                    #    this computes the decorrelation within a cell if 'normal' (cohExp=2)
                    #   otherwise spread = 0 and cohExp = 1
                    spread = 0.0
                    if self.cohMode == 'normal' and nPerpOcc > 1:
                        nctr = 0
                        for io in range(nPerpOcc):
                            inuv = self.kCohBins[i][j][0][io]
                            for jo in range(io):
                                jnuv = self.kCohBins[i][j][0][jo]
                                dist = math.sqrt((inuv[0] - jnuv[0])**2. + (inuv[1] - jnuv[1])**2.)
                                spread+=dist
                                nctr+=1
                        spread = (spread/nctr)/(dkcoher)
                        self.sp.append(math.exp(-spread**2.))
                    #check if in big bin
                    for n in range(int(Nch)):
                        kcoh = math.sqrt(self.kCohPerpVal[i]**2. + self.kCohPerpVal[j]**2. + self.kCohParVal[n]**2.)
                        if kcoh > klo and kcoh <= khi:
                            if nPerpOcc - len(self.kCohBins[i][j][n]) != 0:
                                print 'Incorrect nPerpOcc value:  ',nPerpOcc
                            ss += math.pow(nPerpOcc*math.exp(-(spread**2.)),cohExp)
            self.boost.append(ss)
        return i
        
    def coherenceGrid(self):
        """Goes through all ku, kv to find occupancy.  Saves all points into coherence grid"""
        dkcoher = self.dkcoher
        Nch = self.Nch
        nCoherBin = self.nCoherBin
        nOcc = 0
        # Check all of the perp bins for occupancy...
        for i in range(nCoherBin):
            kubin = self.kCohPerpVal[i]
            for j in range(nCoherBin):
                kvbin = self.kCohPerpVal[j]
                kubinmin = kubin - dkcoher/2.0
                kubinmax = kubin + dkcoher/2.0
                kvbinmin = kvbin - dkcoher/2.0
                kvbinmax = kvbin + dkcoher/2.0
                kperpcoh = math.sqrt( (kubinmax)**2 + (kvbinmax)**2)
                if kperpcoh < self.kmin:
                    continue
                Occupied = False
                perpOcc = []
                for nuv in range(self.nKline/self.Nch):
                    if (self.ku[nuv] >= kubinmin and self.ku[nuv] < kubinmax) and (self.kv[nuv] >= kvbinmin and self.kv[nuv] < kvbinmax):
                        Occupied = True
                        perpOcc.append(nuv)
                if Occupied:   # then append the ku,kv values
                    for n in range(int(Nch)):
                        nOcc+=1
                        for inum in perpOcc:
                            self.kCohBins[i][j][n].append([self.ku[inum],self.kv[inum]])
        return nOcc

    def calcCentroid(self):
        """compute centroid, etc"""
        self.centroidE = 0.0
        self.centroidN = 0.0
        for i in range(self.Nant):
            self.centroidE+=s[i].E;
	    self.centroidN+=s[i].N;
	self.centroidE/=(self.Nant);
	self.centroidN/=(self.Nant);
        return i

    def calcuv(self,i,j,HAhr,decdeg=None,latdeg=None):
        """calculates a uvw point.  Assumes HA in hours, lat in deg and dec in deg"""

        bu = Units[self.bUnit]/Units['m']   # make sure baseline in meters to match obsWavelength, which always is
	dE = bu*(self.E[i] - self.E[j])
	dN = bu*(self.N[i] - self.N[j])
	dU = bu*(self.U[i] - self.U[j])
        if self.arrayType == 'platform':
            u = dE/self.obsWavelength
            v = dN/self.obsWavelength
            w = dU/self.obsWavelength
            if self.rotate == False:  # we don't rotate the platform, so rotate uv
                p = al.parallacticAngle(HAhr*15.0, self.dec, self.lat)*(math.pi/180.0)
                up = u*math.cos(p) - v*math.sin(p)
                v  = u*math.sin(p) + v*math.cos(p)
                u = up
        else:
            HA  = HAhr*(15.0*math.pi/180.0)
            if decdeg==None:
                dec = self.dec*math.pi/180.0
            else:
                dec = decdeg*math.pi/180.0
            if latdeg==None:
                lat = self.lat*math.pi/180.0
            else:
                lat = latdeg*math.pi/180.0

            X = (-math.sin(lat)*dN + math.cos(lat)*dU)/self.obsWavelength
            Y = dE/self.obsWavelength
            Z = (math.cos(lat)*dN + math.sin(lat)*dU)/self.obsWavelength
            u =  X*math.sin(HA)               + Y*math.cos(HA)
            v = -X*math.sin(dec)*math.cos(HA) + Y*math.sin(dec)*math.sin(HA) + Z*math.cos(dec)
            w =  X*math.cos(dec)*math.cos(HA) - Y*math.cos(dec)*math.sin(HA) + Z*math.sin(dec)

	return u,v,w

    def calcduv(self,i,j,HAhr,decdeg=None,latdeg=None):
        """calculates a uvw velocity.  Assumes HA in hours, lat in deg and dec in deg"""
	HA  = HAhr*(15.0*math.pi/180.0)
	if decdeg==None:
            dec = self.dec*math.pi/180.0
        else:
            dec = decdeg*math.pi/180.0
        if latdeg==None:
            lat = self.lat*math.pi/180.0
        else:
            lat = latdeg*math.pi/180.0
        dHdt = 2.0*math.pi/(23.0*3600.0 + 56.0*60.0)
	HA  = HAhr*(15.0*math.pi/180.0)
	lat = self.lat*math.pi/180.0
	dec = self.dec*math.pi/180.0
        bu = Units[self.bUnit]/Units['m']   # make sure baseline in meters to match obsWavelength, which always is
	dE = bu*(self.E[i] - self.E[j])
	dN = bu*(self.N[i] - self.N[j])
	dU = bu*(self.U[i] - self.U[j])
	X = (-math.sin(lat)*dN + math.cos(lat)*dU)/self.obsWavelength
	Y = dE/self.obsWavelength
	Z = (math.cos(lat)*dN + math.sin(lat)*dU)/self.obsWavelength
	u = dHdt*( X*math.cos(HA)               - Y*math.sin(HA))
	v = dHdt*( X*math.sin(dec)*math.sin(HA) + Y*math.sin(dec)*math.cos(HA) + Z*math.cos(dec))
	w = dHdt*(-X*math.cos(dec)*math.sin(HA) - Y*math.cos(dec)*math.cos(HA) + Z*math.sin(dec))
	return u,v,w

    def calcStats(self):
        """Calculates min max etc"""
	self.bmax = -1E9;
	self.bmin = 1E9;
	self.Emax = -1E9;
	self.Emin = 1E9;
	self.Nmax = -1E9;
	self.Nmin = 1E9;
	self.Umax = -1E9;
	self.Umin = 1E9;
	for i in range(len(self.E)):
	    self.bmax = 0.0
	    self.bimax = 0
	    self.bmin = 1E6
	    self.bimin = 0
	    if self.E[i] > self.Emax:
                self.Emax = self.E[i]
	    if self.E[i] < self.Emin:
                self.Emin = self.E[i]
	    if self.N[i] > self.Nmax:
                self.Nmax = self.N[i]
	    if self.N[i] < self.Nmin:
                self.Nmin = self.N[i]
	    if self.U[i] > self.Umax:
                self.Umax = self.U[i]
	    if self.U[i] < self.Umin:
                self.Umin = self.U[i]
	    for j in range(i):
		dtmp = math.sqrt((self.E[j] - self.E[i])**2.0 + (self.N[j] - self.N[i])**2)
		if dtmp > self.bmax:
                    self.bmax = dtmp
		if dtmp < self.bmin:
                    self.bmin = dtmp
        print "    E-W limits: "+str(self.Emin)+" - "+str(self.Emax)+' '+self.bUnit
        print "    N-S limits: "+str(self.Nmin)+" - "+str(self.Nmax)+' '+self.bUnit
        print "    baseline-min - baseline-max: "+str(self.bmin)+" - "+str(self.bmax)+' '+self.bUnit
        return i

    def calckmin(self,z=None,bmin=None,restFreq=None):
        """Computes kmin for bmin in meters, freq in GHz and z.  Returns kmin - does NOT set self.kmin"""
        if z==None:
            z = self.z
        if bmin==None:
            bmin = self.bmin*Units[self.bUnit]/Units['m']
        if restFreq==None:
            restFreq = self.restFreq
            obsFreq = self.obsFreq
            wavelength = self.obsWavelength  # always in meters
        else:
            obsFreq = restFreq/(1.+z)
            wavelength = 0.3/obsFreq  # wavelength in m            
        twopi = 2.0*math.pi
        if bmin == 0.0:
            print "Must have non-zero minimum baseline - resetting to 1.4"
            bmin = 1.4
        u = bmin/wavelength
        kmin = twopi*u/self.calcX(z)
        print 'kmin='+str(self.kmin)
        print 'z=',z,'bmin=',bmin,'restFreq=',restFreq,'obsFreq=',obsFreq
        return kmin

    def calcX(self,z,units="MpcRad"):
        """Computes the spatial conversion factor X for co-moving distances.  Units:  MpcRad (default) and MpcArcmin"""
        DA_Mpc = self.ned.calcUniverse(z)
        X = (1.+z)*DA_Mpc
        if units=='MpcArcmin':
            X = X*math.pi/10800.0
        return X

    def calcY(self,z,restFreq):
        """Computes the frequency conversion factor Y for co-moving distances.  Units Mpc/[nu] where nu is the rest frequency"""
        DA_Mpc = self.ned.calcUniverse(z)
        Y = ( self.ned.Ynu )/restFreq
        return Y

    def writeData(self,x,y,fileName='data.out'):
        """writes all x,y data to file"""

        fp = open(fileName,"w")
        for i,v in enumerate(x):
            s = str(v)+'\t'+str(y[i])+'\n'
            fp.write(s)

        fp.close()
        return i


    def computeDValues(self,kbin=None,boost=None):
        """Computes array of D values.  If no arguments are given uses default with MMQboost
           If kbin is given it calculates on those bin centers using MMQboost
           if kbin and boost are given, it uses those values."""

        Delta = []

        if self.bmin < self.Dant:
            print "Your minimum baseline ("+str(self.bmin)+") is less than the antenna diameter ("+str(self.Dant)+")!!!"
        if kbin == None:
            kbin = self.k
        if boost == None:
            boost = self.R
            boostType = 'MMQ'
        else:
            boostType = 'count'
        if len(kbin) != len(boost):
            print "The bins and boosts have different numbers of elements (k: "+str(len(kbin))+" vs R: "+str(len(boost))+")"
            return Delta

        #print "Using "+boostType
        for i,k in enumerate(kbin):
            R = boost[i]
            f = self.restFreq*Units[self.freqUnit]/Units['GHz']    # make sure in GHz
            bw= self.BW*Units[self.bwUnit]/Units['GHz']            # make sure in GHz
            tint = self.tint*Units[self.tintUnit]/Units['hr']      # make sure in hours
            Dant = self.Dant*Units[self.DantUnit]/Units['m']            # make sure in meters
            Delta.append( self.calcD(k,self.Tsys,self.z,R,self.Npt,self.Npol,f,Dant,bw,tint,self.dkpk,boostType) )
            i+=1
        return Delta

    def calcD(self,k,Tsys,z,R,Npt,Npol,restFreq,Dant,B,tobs,dkpk,boostType):
        """Calculate a single D value, given all parameters no default values.  Returns microKelvin or uK^2"""
        """\tk in Mpc-1"""
        """\tTsys/efficiency in K"""
        """\tz"""
        """\tR value/boost"""
        """\tNpt number of independent pointings"""
        """\tNpol number of independent polarizations"""
        """\tnu rest frequency in GHz"""
        """\tDant in meters"""
        """\tB bandwidth in GHz"""
        """\ttobs in hours"""
        """\tdkpk delta k_parallel over k"""
        """\tboostType is either the MMQ approximation or the counted ones"""

        # First check for all denominator terms
        if restFreq==0.0:
            print "Frequency = 0!!! resetting to 230 GHz"
            restFreq = 230.
        if Dant==0.0:
            print "Dant = 0!!! Not valid - reset to 1.2 m"
            Dant = 1.2
        if tobs==0.0:
            print "tobs = 0!!! Not valid - reset to 1000 hours"
            tobs = 1000.0
        if R == 0.0:
            R = 1E-20
        if R < 0.0:
            print 'R<0, reset to 0.001 to just get on with it'
            R = 0.001
        if Npt == 0:
            Npt = 1
        if boostType == 'MMQ':
            if B==0.0:
                print "Bandwidth = 0!!! resetting to 1 GHz"
                B = 1.0
            if dkpk==0.0:
                print "dkpk = 0!!! resetting to 1"
                dkpk=1.0

        #Do a single baseline:
        tobs = tobs*3600.0            #Convert to seconds
        restFreq = restFreq*1E9       #Convert to Hertz
        obsFreq = restFreq/(1.+z)
        obsLamb = 3E8/(obsFreq)       #in meters
        FOV = 1.18*(obsLamb/Dant)**2.0
        kv = (k**3.0)/(2.0*(math.pi**2.0))
        X = self.calcX(z)
        Y = self.calcY(z,restFreq)
        singleBaselineSensSq = (X**2.0)*Y*kv*FOV*(Tsys**2.0)/(2.0*tobs)
        
        #Calculate "boost"
        if boostType == 'MMQ':
            B = 1E9*B #Convert to Hz
            kpar = 2.0*math.pi/(Y*B)
            Nkbin = k*dkpk/kpar
        else:
            Nkbin = 1.0
        boost = math.sqrt(Npt*Npol*R*Nkbin)*self.Nplat
        if boost==0.0:
            print 'boost = 0, reset to 0.001 to just go ahead with things'
            boost=0.001

        #Calculate sensitivity
        D2 = singleBaselineSensSq/boost
        D = math.sqrt(D2)*1.0E6    # convert to microKelvin
        if self.outputType=='sq':
            D*=D

        return D


    def computeMMQboost(self):
        """Generates the equivalent boost to compare"""
        nbin = self.nbin
        self.computeRMMQValues(klo=self.kbin[ctr][0],khi=self.kbin[ctr][nbin-1],Npp=self.Npp,kstep=nbin)
        B = self.BW*Units[self.bwUnit]
        f = self.restFreq*Units[self.freqUnit]
        Y = self.calcY(self.z, f)
        kpar = 2.0*math.pi/(Y*B)
        self.MMQboost=[]
        for i,k in enumerate(self.kbin[ctr]):
            Nkbin = k*self.dkpk/kpar
            R = self.R[i]
            boost = math.sqrt(self.Npol)*R*Nkbin*(self.nint**2.0)
            self.MMQboost.append(boost)
        return self.MMQboost

    def calcRMMQ(self,k,kmin,Npp):
        """Calculate a single R value (MMQ) given k, kmin, array"""
        ind = self.arrNpp.index(Npp)

        R = 0.0
        for i in range(self.arrNk[ind]):
            if k >= math.sqrt(self.B2coeff[i])*kmin:
                R+=self.nocoeff[i]*(self.arrNuv[ind][i])**2.0
        return R
    
    def computeRMMQValues(self,klo=None,khi=None,kstep=None,Npp=None):
        """computes a range of R values - uses selected config
           if kstep>1 then number of steps in log (e.g. 100 vs 0.01)"""
        kmin = self.kmin
        kmax = self.kmax
        if klo == None:
            klo = 0.95*kmin
        if khi == None:
            khi = 1.05*kmax
        if kstep == None:
            kstep = self.nbin
        if Npp == None:
            Npp = self.Npp
        # Assume using non-log steps
        N = int(math.ceil((khi-klo)/kstep))
        useLogStep = False
        if kstep>1.0:
            N = int(kstep)
            useLogStep=True
            kstep = (math.log10(khi)-math.log10(klo))/N #compute just in case...
        self.k = []
        self.R = []
        for i in range(N):
            if useLogStep:
                self.k.append(klo*math.pow(10.0,i*kstep))
            else:
                self.k.append(klo + i*kstep)
            self.R.append(self.calcRMMQ(self.k[i],kmin,Npp))
        return self.R

    
