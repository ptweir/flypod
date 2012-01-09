import numpy as np
import pylab
pylab.ion()

COLORS = dict(N=[.9,.9,.9],E=[.8,.8,.8],S=[.7,.7,.7],W=[1,1,1])
ROTATIONS = dict(N=0,E=-np.pi/2,S=np.pi,W=np.pi/2) #CHECK

def plotIt(ax,flyToPlot,skyToPlot,changebuffer):
    startTime = skyToPlot['changeTimes'][0]
    
    posToPlot = np.ma.masked_invalid(np.mod(-flyToPlot['orientations']*np.pi/180.0 + np.pi/2 + np.pi,2*np.pi),copy=True)
    time = np.ma.masked_invalid(flyToPlot['times'],copy=True)
    timesToPlot = np.ma.masked_invalid((flyToPlot['times']-startTime)/60.0,copy=True)
    rotations = np.ma.masked_all(posToPlot.shape)
    
    for changeTimeInd, changeTime in enumerate(skyToPlot['changeTimes'][:-1]):
        patchHandle = ax.axvspan((changeTime-changebuffer-startTime)/60.0, (skyToPlot['changeTimes'][changeTimeInd+1]-changebuffer-startTime)/60.0)
        patchHandle.set_facecolor(COLORS[skyToPlot['directions'][changeTimeInd]])
        patchHandle.set_edgecolor('none')
        
        changeInds = (time > changeTime-changebuffer) & (time <= skyToPlot['changeTimes'][changeTimeInd+1]-changebuffer)
        rotations[changeInds] = ROTATIONS[skyToPlot['directions'][changeTimeInd]]
    
    worldOrientations = np.mod(posToPlot+rotations,2*np.pi)
    worldOrientations[abs(np.diff(worldOrientations))>np.pi] = np.nan

    posToPlot[abs(np.diff(posToPlot))>np.pi] = np.nan
    ax.plot(timesToPlot,worldOrientations,'k',linewidth=.5)
    #ax.plot(timesToPlot,posToPlot,'k',linewidth=.5)
    ax.set_xlim((0, (skyToPlot['changeTimes'][8]-changebuffer-startTime)/60.0))
    ax.set_ylim((0, 2*np.pi))
    ax.spines['left'].set_position(('outward',8))
    ax.spines['bottom'].set_position(('outward',8))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_yticks(np.linspace(0,2*np.pi,5))
    #ax.set_yticklabels(np.array(ax.get_yticks()*180/np.pi,dtype=int))
    ax.set_yticklabels(['E','N','W','S','E'])
    ax.set_xticks(np.arange(0,25,3))
    for tickLine in ax.get_xticklines() + ax.get_yticklines():
        tickLine.set_markeredgewidth(1)
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Compass \nHeading')
    
