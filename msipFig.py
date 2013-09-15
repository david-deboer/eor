import eorpy
import matplotlib.pyplot as plt
pt = 'senscurve'
#pt = 'errb'
plt.rcParams.update({'font.size':20})
#e.setPlot(showLegend=False)
e = eorpy.eor()
e.set2amiba()
e.setTrack(tobs=1500.,Npt=4.0)
dlogk_z6 = 0.18
dlogk_z3 = 0.18

#plot something to get DACOTA as first entry in legend
#plt.loglog([1e-6],[1e-6],'bo',label='DACOTA')

if pt == 'senscurve':
        e.setPlot(plotVersion = 'sens')
else:
        e.setPlot(plotVersion = 'errbar')

###z=3
e.setSys(Tsys=50.0,z=3,restFreq=115.0,BW=4.0,dlogk=dlogk_z3)
e.setPlot(plotModel=True,plotBins=True,plotLines=False)
e.arraySensitivity(color='b')
if pt == 'senscurve':
	e.setPlot(plotBins=False,plotModel=False,plotLines=True)
	e.arraySensitivity(color='b',label=r'DACOTA z=3')

###z=6, trusting all defaults...
e.setSys(Tsys=50.0,z=6,restFreq=230.0,BW=8.0,dlogk=dlogk_z6)
e.setPlot(plotModel=True,plotBins=True,plotLines=False)
e.arraySensitivity(color='k')
if pt == 'senscurve':
	e.setPlot(plotBins=False,plotModel=False,plotLines=True)
	e.arraySensitivity(color='k',label=r'DACOTA z=6')


plt.subplot(121)
plt.axis([0.1,3.0,0.001,20000])
plt.subplot(122)
plt.axis([0.1,3.0,0.001,20000])
