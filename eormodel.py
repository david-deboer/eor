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
        print "\nReading "+fileName+" and scaling for z="+str(z)
        self.kModel = []
        self.kModel2= []
        self.kKISS = []
        self.dKISS = []
        self.co1 = []
        self.co2 = []
        self.co3 = []
        dfn = ' '
        i=1
        for line in fp:
            line=line.strip()
            t = line.split()
            if z<4.:
                k = float(t[0])*h
                self.kModel.append(k)
                co3 = float(t[1])*3.0/(1.0+z)
                if self.outputType!='sq':
                    co3=math.sqrt(co3)
                self.co3.append(co3)
                dfn = 'co3'
            elif fileName == 'fpaSensRead.out':
                k = float(t[0])
                self.kKISS.append(k)
                d = float(t[1])
                self.dKISS.append(d)
            else:
                k = h*float(t[0])
                self.kModel.append(k)
                k = h*float(t[2])
                self.kModel2.append(k)
                co1 = math.sqrt(7.0/(1.0+z))*float(t[1])
                co2 = math.sqrt(7.0/(1.0+z))*float(t[3])
                if self.outputType=='sq':
                    co1*=co1
                    co2*=co2
                self.co1.append(co1)
                self.co2.append(co2)
                dfn = 'co1, co2'
            i+=1
        print "Note arrays are:  kModel, "+dfn
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
                fileName = 'cosmoz3_12Oct_modelB.txt'
            elif eor.z >= 4.0:
                fileName = 'cosmoz6.txt'
            else:
                fileName = 'cosmoz6.txt'
        elif mType=='file':
            fileName = mPar
        else:
            print mType+' '+mPar+' not found'
            return 0

        self.readModelData(fileName=fileName,z=eor.z)
            
        if eor.z < 4.0:
            cl1 = 'c'
            plt.loglog(self.kModel,self.co3,cl1, label='z=3, Model B') #r'z=3, SFR$_{min}$=0.01 M$_{\odot}$/yr')
            k1 = np.copy(self.kModel)
            c1 = np.copy(self.co3)
        else:
            cl1 = 'k'
            cl2 = 'r'
            plt.loglog(self.kModel,self.co1,cl1, label=r'z=6, $10^{10}$ M$_{\odot}$')
            plt.loglog(self.kModel2,self.co2,cl2, label=r'z=6, $10^{8}$ M$_{\odot}$')
            k1 = self.kModel
            c1 = self.co1
            k2 = self.kModel2
            c2 = self.co2
        if (eor.outputType=='sq'):
            plt.axis([0.1,5.0,0.01,40000.0])
        else:
            plt.axis([0.06,5.0,.1,40])

        if mType=='auto' and mPar == 'file' and eor.z < 4.0:
            fileName = 'cosmoz3_12Oct_modelA.txt'
            self.readModelData(fileName=fileName,z=eor.z)
            k2 = np.copy(self.kModel)
            c2 = np.copy(self.co3)
            cl2 = 'm'
            plt.loglog(k2,c2,cl2, label='z=3, Model A') #r'z=3, $10^{9}$ M$_{\odot}$')

        if eor.z > splitOut:
            if eor.z > 4.5:
                a = self.co1[-1]/( (self.kModel[-1]/h)**3.0)
                co1_p = a*(np.array(self.kModel)/h)**3.0
                co1_c = np.array(self.co1) - co1_p
                sm = interpolate.UnivariateSpline(k1,co1_c,s=40)
                sm_co1_c = sm(k1)
                plt.loglog(self.kModel,co1_p,'k:')
                plt.loglog(self.kModel,sm_co1_c,'k--')
                
                b = self.co2[-1]/((self.kModel2[-1]/h)**3.0)
                co2_p = b*(np.array(self.kModel2)/h)**3.0
                co2_c = np.array(self.co2) - co2_p
                sm = interpolate.UnivariateSpline(k2,co2_c,s=40)
                sm_co2_c = sm(k2)        
                plt.loglog(self.kModel2,co2_p,'r:')
                plt.loglog(self.kModel2,sm_co2_c,'r--')
            else:
                a = c1[-1]/( (k1[-1]/h)**3.0)
                co1_p = a*(k1/h)**3.0
                co1_c = c1 - co1_p
                sm = interpolate.UnivariateSpline(k1,co1_c)
                sm_co1_c = sm(k1)
                hfil = np.where(k1>1.0)
                sm_co1_c[hfil] = 1e-6
                plt.loglog(k1,co1_p,'c:')
                #plt.loglog(k1,sm_co1_c,'c--')

                b = c2[-1]/( (k2[-1]/h)**3.0)
                co2_p = b*(k2/h)**3.0
                co2_c = c2 - co2_p
                sm = interpolate.UnivariateSpline(k2,co2_c)
                sm_co2_c = sm(k2)
                hfil = np.where(k2 > 0.7)
                sm_co2_c[hfil] = 1e-6
                plt.loglog(k2,co2_p,'m:')
                #plt.loglog(k2,sm_co2_c,'m--')

        if yerrbar:
            #plt.figure(1)
            #ax = plt.subplot(111)
            #ax.set_xscale("log",nonposx='clip')
            #ax.set_yscale("log",nonpoxy='clip')
            dither = 1.0
            scale = 1.0
            if cl1=='k':
                dither = 0.98
                scale = 2.5
                print 'SCALING ERRBAR BY '+cl1+'  ',scale
            x = np.copy(yerrbar[0])*dither
            f = interpolate.interp1d(k1,c1)
            d = f(x)
            yerr = fix_yerr(d,yerrbar[1],1e-6,scale)
            plt.errorbar(x,d,yerr=yerr,xerr=None,fmt='_',color=cl1)
            plt.loglog(x,d,'o',color=cl1)

            if cl2=='r':
                dither = 1.02
                scale = 2.5
                print 'SCALING ERRBAR BY '+cl2+'  ',scale
            x = np.array(yerrbar[0],copy=True)*dither
            f2 = interpolate.interp1d(k2,c2)
            d2 = f2(x)
            yerr = fix_yerr(d2,yerrbar[1],1e-6,scale)
            plt.errorbar(x,d2,yerr=yerr,xerr=None,fmt='_',color=cl2)
            if cl2=='m':
                plt.loglog(x,d2,'o',label='DACOTA',color=cl2)
            else:
                plt.loglog(x,d2,'o',color=cl2)
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
