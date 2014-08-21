###Imports
import matplotlib.pyplot as plt
import math
import sys
import string
import eorpy
import eormodel
import numpy as np
em = eormodel.modelData()

print 'loading eorio'

def plotSensitivity(eor,fignum=1,color=None,label='auto',linewidth=2):
    """plot sensitivity"""

    if label == 'auto':
        label = '%s z=%.0f' % (eor.config,eor.z)

    print 'EORIO.PY:  Note that to get both bins and lines, you need to run with plotBins first'

    little_h = eor.h
    edgeLeft = eorpy.edgeLeft
    ctr = eorpy.ctr
    edgeRight = eorpy.edgeRight

    plt.figure(fignum)
    panels=False
    if panels:
        if eor.z < 4.5:
            plt.subplot(121)
        else:
            plt.subplot(122)
    
    xdat=np.array(eor.kbin[ctr])*little_h
    ydat=np.array(eor.d)
    if eor.plotVersion == 'errbar':
        em.plotModelData(eor,[xdat,ydat],fignum,color,linewidth)
    else:
        if eor.plotModel:
            em.plotModelData(eor,None,fignum,color,linewidth)
        if eor.plotSteps:
            xdat,ydat=stepWise(eor,little_h)
        if eor.plotBins:
            errx,erry=errorBars(eor,little_h)
##            if eor.z < 4.5:
##                plt.subplot(121,xscale='log',yscale='log')
##            else:
##                plt.subplot(122,xscale='log',yscale='log')
            if color==None:
                plt.errorbar(xdat,ydat,yerr=None,xerr=errx,fmt='_')
            else:
                color = 'b'
                plt.errorbar(xdat,ydat,yerr=None,xerr=errx,fmt='_')#,ecolor=color)
        if eor.plotLines:
            if color==None:
                plt.loglog(xdat,ydat,ls=eor.ltype,label=label,linewidth=linewidth)
            else:
                plt.loglog(xdat,ydat,color=color,ls=eor.ltype,label=label,linewidth=linewidth)
        else:
            plt.loglog(xdat,ydat,'+')
        #if eor.plotModel:
        #    em.plotModelData(eor,None,fignum,color,linewidth)
    plt.xlabel(r'$k$ [h$\cdot$Mpc$^{-1}$]')
    #plt.ylabel(r'$\Delta^2_{CO}$ or $\Delta^2_N$ [$\mu$K$^2$]')
    plt.ylabel(r'$\Delta^2_{CO}$ [$\mu$K$^2$]')
    if eor.showLegend:
        plt.legend(loc='upper left') #lower right')
        
    if eor.plotMMQ:     # compute approx
        bmmq = eor.computeMMQboost()
        dmmq = eor.computeDValues(eor.kbin[ctr],bmmq)
        dmmq = eor.computeTracks(dmmq)
        plt.figure(fignum)
        if eor.plotLines:
            plt.loglog(eor.kbin[ctr],dmmq)
        else:
            plt.loglog(eor.kbin[ctr],dmmq,'x')

def plotBoost(eor,fignum=2):
    plt.figure(fignum)
    plt.title('Boost [k vs boost]')
    plt.loglog(eor.kbin[eorpy.ctr],eor.boost)
    if eor.plotMMQ:     # compute approx
        bmmq = eor.computeMMQboost()
        plt.figure(fignum)
        if eor.plotLines:
            plt.loglog(eor.kbin[eorpy.ctr],bmmq)
        else:
            plt.loglog(eor.kbin[eorpy.ctr],bmmq,'x')    

def stepWise(eor,h):
    stepK=[]
    stepD=[]
    for i in range(eor.nbin-1):
        stepK.append(h*eor.kbin[eorpy.edgeLeft][i])
        stepK.append(h*eor.kbin[eorpy.edgeRight][i])
        stepK.append(h*eor.kbin[eorpy.edgeLeft][i+1])
        stepD.append(eor.d[i])
        stepD.append(eor.d[i])
        stepD.append(eor.d[i+1])
    return stepK,stepD
def errorBars(eor,h):
    errx = []
    erry = []
    for i in range(eor.nbin):
        errx.append(h*(eor.kbin[eorpy.edgeRight][i] - eor.kbin[eorpy.edgeLeft][i])/2.0)
        erry.append(0.0)
    return errx, erry

def writeData(eor,x,y,fileName='data.out'):
    """writes all x,y data to file"""

    fp = open(fileName,"w")
    for i,v in enumerate(x):
        s = str(v)+'\t'+str(y[i])+'\n'
        fp.write(s)

    fp.close()
    return i


def drawBins(eor,thetaStep=1,fignum=100,plotBins=True,writeBins=False,span=360.0,color='k'):
    """Generates the circles corresponding to the kedge bins: bins[ [x,y] ]
       thetaStep in degrees"""

    bins = []
    ntheta = int(math.ceil(span/thetaStep)) + 1
    
    for k in eor.kbin[eorpy.edgeLeft]:
        x = []
        y = []
        for i in range(ntheta):
            x.append(k*math.cos(i*thetaStep*math.pi/180.0))
            y.append(k*math.sin(i*thetaStep*math.pi/180.0))
        d = [x,y]
        bins.append(d)

    print "plot(eor.bins[nedge][0],eor.bins[nedge][1])"
    if writeBins:
        for i in range(len(eor.kbin[eorpy.edgeLeft])):
            filename = 'circ'+str(i)+'.dat'
            print "Writing "+filename
            eor.writeData(bins[i][0],bins[i][1],filename)
    if plotBins:
        plt.figure(fignum)
        for i in range(len(eor.kbin[eorpy.edgeLeft])):
            plt.plot(bins[i][0],bins[i][1],color)
        
    return ntheta


