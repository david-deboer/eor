import datetime

Units = {'Hz':1.0, 'kHz':1.0E3, 'MHz':1.0E6, 'GHz':1.0E9, 'sec':1.0, 'min':60.0, 'hr':3600.0, 'm':1.0, 'cm':0.01, 'km':1000.0}

print 'loading eorview'

class view:
    def __init__(self):
        """This just prints out the parameters to an html page
        as well as any other outputs"""
        

    def updateWeb(self,bl):
        """Write web-page with system, track and data parameters"""

        if bl.arrayType == 'NA':
            fp = open("index.html","w")
            s = "<html>\n<body>\n<table border cellspacing=0 cellpadding=5>\n"
            fp.write(s)
            s = "<tr>\n\t<th colspan=2 align=center>arrayType</th><td colspan=6 align=center>"+bl.arrayType+"</td></tr>"
            fp.write(s)
            s = "</table>\n</body>\n</html>\n"
            fp.write(s)
            fp.close()
            return 0
        elif bl.arrayType=='platform':
            rotopt = bl.rotate
            dacsh = bl.dacotaShortcut
        else:
            rotopt = 'N/A'
            dacsh = 'N/A'
            
        fp = open("index.html","w")
        dt = str(datetime.datetime.today())
        s = "<html>\n<body>\n<p>"+dt+"\n<p>\n"
        fp.write(s)
        s = "<table border cellspacing=0 cellpadding=5>\n"
        fp.write(s)
        s = "<tr>\n\t<th colspan=2 align=center>arrayType</th><td colspan=6 align=center>"+bl.arrayType+"</td></tr>"
        fp.write(s)
        s = "<tr>\n\t<th colspan=8 bgColor='black'> </th>\n</tr>"
        fp.write(s)
        s = "<tr>\n\t<th colspan=2 align=center>System</th><td> </td><th colspan=2 align=center>Track</th><td> </td><th colspan=2 align=center>Data</th>\n</tr>"
        fp.write(s)
        #####Row 1
        s = "<tr>\n\t<th align=left>Config</th><td>"+bl.config+\
            "</td><td> </td><th align=left>HAmin</th><td>"+str(bl.HAmin)+" "+bl.HAunit+\
            "</td><td> </td><th align=left>boxcar</th><td>"+str(bl.boxcar)+"</td>\n</tr>\n"
        fp.write(s)
        #####Row 2
        s = "<tr>\n\t<th align=left>Npp</th><td>"+str(bl.Npp)+\
            "</td><td> </td><th align=left>HAmax</th><td>"+str(bl.HAmax)+" "+bl.HAunit+\
            "</td><td> </td><th align=left>outputType</th><td>"+bl.outputType+"</td>\n</tr>\n"
        fp.write(s)
        #####Row 3
        a = '%.4f %s' % (bl.tint, bl.tintUnit)
        s = "<tr>\n\t<th align=left>Nplat</th><td>"+str(bl.Nplat)+\
            "</td><td> </td><th align=left>tint</th><td>"+a+\
            "</td><td> </td><th align=left>coherentTracks</th><td>"+str(bl.coherentTracks)+"</td>\n</tr>\n"
        fp.write(s)
        #####Row 4
        if isinstance(bl.dkcoher,float):
            a = '%.4f' % bl.dkcoher
        else:
            a = str(bl.dkcoher)        
        s = "<tr>\n\t<th align=left>Dant</th><td>"+str(bl.Dant)+" "+bl.DantUnit+\
            "</td><td> </td><th align=left>tobs</th><td>"+str(bl.tobs)+" "+bl.tobsUnit+\
            "</td><td> </td><th align=left>dkcoher</th><td>"+a+"</td>\n</tr>\n"
        fp.write(s)
        #####Row 5
        a = '%.4f / %.4f ' % (bl.bmin,bl.bmax)
        b = '%.4f' % bl.dkcoherScale
        s = "<tr>\n\t<th align=left>bmin / bmax</th><td> "+a+bl.bUnit+\
            "</td><td> </td><th align=left>[nint]</th><td>"+str(bl.nint)+\
            "</td><td> </td><th align=left>dkcoherScale</th><td>"+b+"</td></tr>\n"
        fp.write(s)
        #####Row 6
       
        s = "<tr>\n\t<th align=left>Tsys/eta</th><td>"+str(bl.Tsys)+" K"+\
            "</td><td> </td><th align=left>[nTrack]</th><td>"+str(bl.nTrack)+\
            "</td><td> </td><th align=left>cohMode</th><td>"+bl.cohMode+" </td>\n</tr>\n"
        fp.write(s)
        #####Row 7
        a = '%.4f / %.4f' % (bl.dlogk,bl.dkpk)
        s = "<tr>\n\t<th align=left>dlogk/[dkpk]</th><td>"+a+\
            "</td><td> </td><th align=left>Npt</th><td>"+str(bl.Npt)+\
            "</td><td> </td><td align=left> </td><td> </td>\n</tr>\n"
        fp.write(s)
        #####Row 8
        s = "<tr>\n\t<th align=left>z</th><td>"+str(bl.z)+\
            "</td><td> </td><th align=left>dec</th><td>"+str(bl.dec)+\
            "</td><td> </td><th align=left> </th><td> </td>\n</tr>\n"
        fp.write(s)
        #####Row 9
        a = '%.4f / %.4f' % (bl.restFreq, bl.obsFreq) + ' '+bl.freqUnit
        s = "<tr>\n\t<th align=left>restFreq/[obsFreq]</th><td>"+a+\
            "</td><td> </td><th align=left>lat</th><td>"+str(bl.lat)+\
            "</td><td> </td><th align=center colspan=2>Plot</th>\n</tr>\n"
        fp.write(s)
        #####Row 10
        a = '%.2f' % bl.elLimit
        s = "<tr>\n\t<th align=left>[Nant]</th><td>"+str(bl.Nant)+\
            "</td><td> </td><th align=left>elLimit</th><td> "+a+\
            "</td><td> </td><th align=left>plotModel </th><td>"+repr(bl.plotModel)+"</td>\n</tr>\n"
        fp.write(s)
        #####Row 11
        s = "<tr>\n\t<th align=left>[nBaseline]</th><td>"+str(bl.nBaseline)+\
            "</td><td> </td><th align=left>shadow </th><td> "+str(bl.shadow)+\
            "</td><td> </td><th align=left>plotSteps </th><td>"+repr(bl.plotSteps)+"</td>\n</tr>\n"
        fp.write(s)
        #####Row 12
        s = "<tr>\n\t<th align=left>BW</th><td>"+str(bl.BW)+" "+bl.bwUnit+\
            "</td><td> </td><th align=left> </th><td> "+\
            "<td> </td><th align=left>plotBins</th><td>"+repr(bl.plotBins)+"</td>\n</tr>\n"
        fp.write(s)
        #####Row 13
        s = "<tr>\n\t<th align=left>Nch</th><td>"+str(bl.Nch)+\
            "</td><td> </td><th align=center colspan=2>Model</th>"+\
            "<td> </td><th align=left>plotMMQ</th><td>"+repr(bl.plotMMQ)+"</td>\n</tr>\n"
        fp.write(s)
        #####Row 14
        s = "<tr>\n\t<th align=left>Npol</th><td>"+str(bl.Npol)+\
            "</td><td> </td><th align=left>modelType</th><td>"+bl.modelType+" </td>"+\
            "</td><td> </td><th align=left>plotLines</th><td>"+repr(bl.plotLines)+" </td>\n</tr>\n"
        fp.write(s)
        #####Row 15
        s = "<tr>\n\t<th align=left>rotate</th><td>"+str(rotopt)+\
            "</td><td> </td><th align=left> H0 </th><td> "+str(bl.H0)+\
            "</td><td> </td><th align=left>linetype</th><td>"+bl.ltype+" </td>\n</tr>\n"
        fp.write(s)
        #####Row 16
        s = "<tr>\n\t<th align=left>dacotaShortcut </th><td>"+str(dacsh)+\
            "</td><td> </td><th align=left> </th><td> "+\
            "</td><td> </td><th align=left>showLegend</th><td>"+str(bl.showLegend)+" </td>\n</tr>\n"
        fp.write(s)
        #####Row 17
        s = "<tr>\n\t<th align=left></th><td>"+" "+\
            "</td><td> </td><th align=left> </th><td> "+\
            "</td><td> </td><th align=left>plotVersion</th><td>"+str(bl.plotVersion)+" </td>\n</tr>\n"
        fp.write(s)
        
        s = "</table>\n"
        fp.write(s)
        s = "<p><p>\n"
        fp.write(s)
        if bl.kmin != -1. or bl.kmax != -1:
            s = "<table border cellspacing=0 cellpadding=5>\n"
            fp.write(s)
            s = "<tr>\n\t<th>kmin</th><td>"+str(bl.kmin)+"</td><th>kmax</th><td>"+str(bl.kmax)+"</td>\n</tr>\n"
            fp.write(s)
            s = "<tr>\n\t<th>nbin</th><td>"+str(bl.nbin)+"</td><td> </td><td> </td>\n</tr>\n"
            fp.write(s)
            s = "<tr>\n\t<th>#</th><th>left</th><th>center</th><th>right</th>\n</tr>\n"
            fp.write(s)
            for i,kctr in enumerate(bl.kbin[1]):
                s = "<tr>\n\t<th>bin "+str(i)+"</th>"+"<td>"+str(bl.kbin[0][i])+"</td><td>"+str(kctr)+"</td><td>"+str(bl.kbin[2][i])+"</td>\n</tr>\n"
                fp.write(s)
            s = "</table>"
            fp.write(s)
        
        s = "</body>\n</html>"
        fp.write(s)
        fp.close()

    def dispTrack(self,bl):
        """Displays the track parameters"""
        print str(bl.nTrack)+' tracks from '+str(bl.HAmin)+' to '+str(bl.HAmax)+' '+bl.HAunit
        print '\tIntegration '+str(bl.tint)+' '+bl.tintUnit+'  ('+str(bl.nint)+' steps)'
        print '\tDeclination = '+str(bl.dec)+' degrees and Latitute = '+str(bl.lat)
        
    def dispData(self,bl):
        """Displays parameters set by setData"""

        print 'DATA Parameters:'
        print '\tDATA: boxcar = ',bl.boxcar
        print '\tDATA: plotLines = ',bl.plotLines
        print '\tDATA: coherentTracks = ',bl.coherentTracks
        print '\tDATA: outputType:  '+bl.outputType
        print '\tDATA: dkcoher:  ',bl.dkcoher
        print '\tDATA: dkcoherScale:  ',bl.dkcoherScale

    def dispSys(self,bl):
        """Prints sys parameters.  Set in sysSys"""
        print 'System parameters:'
        print '\tconfig:  ',bl.config
        print '\tNplat:  ',bl.Nplat,
        print '=> Nant:  ',bl.Nant,
        print ' nBaseline: ',bl.nBaseline
        print '\tdlogk:  ',bl.dlogk
        print '\tz:  ',bl.z
        print '\trestFreq:  '+str(bl.restFreq)+' '+bl.freqUnit,
        print '=>obsFreq:  '+str(bl.obsFreq)+' '+bl.freqUnit,
        print ' obsWavelength:  '+str(bl.obsWavelength)+' m'
        print '\tBW:  '+str(bl.BW)+' '+bl.bwUnit
        print '\tNch:  ',bl.Nch
        print '\tNpt:  ',bl.Npt
        print '\t\t=>nKline ',bl.nKline

    def dispArray(self,bl):
        """Types out the antenna locations"""
        for i,name in enumerate(bl.elementName):
            print name+'\t('+str(bl.N[i])+','+str(bl.E[i])+') '+bl.bUnit

        return i

    def dispAll(self,bl):
        """Prints all variables"""
        g = bl.__dict__
        keys = g.keys()

        print '---------------\n'
        for k in sorted(g.iterkeys()):
            s = str(g[k])
            if len(s) > 100:
                s = '---many---'
            print k+' = '+s
        print '--------------\n\n'

        print '---float---'
        for k in sorted(g.iterkeys()):
            if isinstance(g[k],float):
                s = str(g[k])
                if len(s) > 100:
                    s = '---many---'
                print '['+k+': ('+s+')]',
        print '\n\n---int---'
        for k in sorted(g.iterkeys()):
            if isinstance(g[k],int):
                s = str(g[k])
                if len(s) > 100:
                    s = '---many---'
                print '['+k+': ('+s+')]',
        print '\n\n---list---'
        for k in sorted(g.iterkeys()):
            if isinstance(g[k],list):
                s = str(g[k])
                if len(s) > 100:
                    s = '---many---'
                print '['+k+': ('+s+')]',
        print '\n'
