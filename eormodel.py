###Imports
import math
import sys
import string
import os.path
import eorpy
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

print 'loading eormodel'

class modelData:
    def __init__(self, fileName=None, path='infiles'):
        if fileName!=None:
            self.fileName = os.path.join(path,fileName)
        
    def readModelData(self,fileName='cosmoz6.txt',path='infiles',z=None):
        """Reads in file so I can plot with other data"""
        fileName = os.path.join(path,fileName)
        fp = open(fileName,"r")
        h = self.h
        if z == None:
            z = eor.z
        if z < 4.5:
            zref = 2.0
        else:
            zref = 7.0
            
        print "\nReading %s and scaling from %.1f to %.1f" % (fileName,zref,z)
        self.kModel = []
        self.co = []
        self.kKISS = []
        self.dKISS = []
        i=1
        for line in fp:
            line=line.strip()
            t = line.split()
            if fileName == 'fpaSensRead.out':
                k = float(t[0])
                self.kKISS.append(k)
                d = float(t[1])
                self.dKISS.append(d)
            else:
                k = float(t[0])
                co = float(t[1])*(1.0+zref)/(1.0+z)
                if self.outputType!='sq':
                    print 'changing co#########'
                    co = math.sqrt(co)
                self.kModel.append(k)
                self.co.append(co)
            i+=1
        return i


    def plotModelData(self,eor,yerrbar=None,fignum=1,color=None,linewidth=2,splitOut=1.0):
        """Plots Model data (including reading in the file)"""

        self.h = eor.H0/100.0
        h = self.h
        self.outputType = eor.outputType
        plt.figure(fignum)

        mType = eor.modelType.split(':')[0]
        try:
            mPar  = eor.modelType.split(':')[1]
        except:
            mPar = ''

        if mType=='auto' and mPar=='file':
            if eor.z < 4.0:
                fileName1 = 'cosmoz2_13Sep_modelB.txt'
                fileName2 = 'cosmoz2_13Sep_modelA.txt'
                cl1 = 'c'
                cl2 = 'm'
            elif eor.z >= 4.0:
                fileName1 = 'cosmoz7_13Sep_10.txt'
                fileName2 = 'cosmoz7_13Sep_8.txt'
                cl1 = 'k'
                cl2 = 'r'

        self.readModelData(fileName=fileName1,z=eor.z)
        k1 = np.copy(self.kModel)
        co1 = np.copy(self.co)
        self.readModelData(fileName=fileName2,z=eor.z)
        k2 = np.copy(self.kModel)
        co2 = np.copy(self.co)
            
        if eor.z < 4.0:
            plt.subplot(121)
            plt.loglog(k1,co1,cl1, label='z=3, Model B') #r'z=3, SFR$_{min}$=0.01 M$_{\odot}$/yr')
            plt.loglog(k2,co2,cl2, label='z=3, Model A') #r'z=3, $10^{9}$ M$_{\odot}$')            
        else:
            plt.subplot(122)
            plt.loglog(k2,co2,cl2, label=r'z=6, $10^{8}$ M$_{\odot}$')
            plt.loglog(k1,co1,cl1, label=r'z=6, $10^{10}$ M$_{\odot}$')

        if eor.z > splitOut:
            a = co1[-1]/( (k1[-1])**3.0)
            co_p = a*(np.array(k1))**3.0
            co_c = co1 - co_p
            sm = interpolate.UnivariateSpline(k1,co_c,s=40)
            sm_co_c = sm(k1)
            st = cl1+'--'
            plt.loglog(k1,co_p,st)
            if eor.z > 4.5:
                st = cl1+':'
                plt.loglog(k1,sm_co_c,st)
            a = co2[-1]/( (k2[-1])**3.0)
            co_p = a*(np.array(k2))**3.0
            co_c = co2 - co_p
            sm = interpolate.UnivariateSpline(k2,co_c,s=40)
            sm_co_c = sm(k2)
            st = cl2+'--'
            plt.loglog(k2,co_p,st)
            if eor.z > 4.5:
                st = cl2+':'
                plt.loglog(k2,sm_co_c,st)

        if yerrbar:
            #plt.figure(1)
            #ax = plt.subplot(111)
            #ax.set_xscale("log",nonposx='clip')
            #ax.set_yscale("log",nonpoxy='clip')
            dither = 1.0
            if eor.z < 1.0:  # these are the dacota points/error bars
                cld = 'g'
            else:
                cld = 'b'
            scale = 1.0
            if cl1=='k':
                dither = 0.98
                scale = 2.5
                print 'SCALING ERRBAR BY '+cl1+'  ',scale
            x = np.copy(yerrbar[0])*dither
            f = interpolate.interp1d(k1,co1)
            d = f(x)
            yerr = fix_yerr(d,yerrbar[1],1e-6,scale)
            plt.errorbar(x,d,yerr=yerr,xerr=None,fmt='_',color=cld)
            plt.loglog(x,d,'o',color=cld)

            if cl2=='r':
                dither = 1.02
                scale = 2.5
                print 'SCALING ERRBAR BY '+cl2+'  ',scale
            x = np.array(yerrbar[0],copy=True)*dither
            f2 = interpolate.interp1d(k2,co2)
            d2 = f2(x)
            yerr = fix_yerr(d2,yerrbar[1],1e-6,scale)
            plt.errorbar(x,d2,yerr=yerr,xerr=None,fmt='_',color=cld)
            if cl2=='m':
                plt.loglog(x,d2,'o',color=cld)
                #plt.loglog(x,d2,'o',label='DACOTA',color=cld)
            else:
                plt.loglog(x,d2,'o',color=cld)

            # calc snr
            x = np.copy(yerrbar[0])
            f = interpolate.interp1d(k1,co1)
            d = f(x)
            outfile = 'snr_'+fileName1
            fp = open(outfile,'w')
            for i in range(len(x)):
                snr = d[i]/yerrbar[1][i]
                sout = '%.4f\t%.4f\n' % (x[i],snr)
                fp.write(sout)
            fp.close()
            x = np.copy(yerrbar[0])
            f = interpolate.interp1d(k2,co2)
            d = f(x)
            outfile = 'snr_'+fileName2
            fp = open(outfile,'w')
            for i in range(len(x)):
                snr = d[i]/yerrbar[1][i]
                sout = '%.4f\t%.4f\n' % (x[i],snr)
                fp.write(sout)
            fp.close()

        
        return 1


def fix_yerr(a,b,c,scale):
    if len(a) != len(b):
        print 'lengths not the same'
        return b
    yer = [[],[]]
    for i, v in enumerate(b):
        if a[i]-v*scale<0.0:
            yer[1].append(v*scale)
            yer[0].append(a[i]-c)
        else:
            yer[1].append(v*scale)
            yer[0].append(v*scale)
    return np.array(yer)
