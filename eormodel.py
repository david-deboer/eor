###Imports
import math
import sys
import string
import os.path
import eorpy
import matplotlib.pyplot as plt

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
                co3 = float(t[1]) #*3.0/(1.0+z)
                #if eor.outputType=='sq':
                #    co3=co3
                self.co3.append(co3)
                dfn = 'co3'
            elif fileName == 'fpaSensRead.out':
                k = float(t[0])
                self.kKISS.append(k)
                d = float(t[1])
                self.dKISS.append(d)
            else:
                k = h*(float(t[0]) + float(t[2]))/2.0
                self.kModel.append(k)
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


    def plotModelData(self,eor,fignum=1,color=None,linewidth=2):
        """Plots Model data (including reading in the file)"""

        self.h = eor.H0/100.0
        self.outputType = eor.outputType

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
            plt.loglog(self.kModel,self.co3,'r--', label='z=3, Model B') #r'z=3, SFR$_{min}$=0.01 M$_{\odot}$/yr')
        else:
            plt.loglog(self.kModel,self.co1,'g', label=r'z=6, $10^{10}$ M$_{\odot}$')
            plt.loglog(self.kModel,self.co2,'r', label=r'z=6, $10^{8}$ M$_{\odot}$')
        if (eor.outputType=='sq'):
            plt.axis([0.1,5.0,0.01,40000.0])
        else:
            plt.axis([0.06,5.0,.1,40])

        if mType=='auto' and mPar == 'file' and eor.z < 4.0:
            fileName = 'cosmoz3_12Oct_modelA.txt'
            self.readModelData(fileName=fileName,z=eor.z)
            plt.loglog(self.kModel,self.co3,'g--', label='z=3, Model A') #r'z=3, $10^{9}$ M$_{\odot}$')

            
        return 1


