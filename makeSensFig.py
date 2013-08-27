from bleor import *
b = bleor()

b.set2amiba()
b.arraySensitivity(color='k')
b.arraySensitivity(color='k',plotLidz=False,plotBins=False)

b.set2sza('karto')
##b.arraySensitivity(color='b',plotLidz=False)
##b.arraySensitivity(color='b',plotLidz=False,plotBins=False)

#b.set2amiba(z=3)
#b.arraySensitivity(color='0.75')
#b.arraySensitivity(color='k',plotLidz=False,plotBins=False)
b.z = 3.
b.plotLidzData(z=3.)

plt.rc('legend',**{'fontsize':12})
plt.legend(loc='lower right',labelspacing=0)
