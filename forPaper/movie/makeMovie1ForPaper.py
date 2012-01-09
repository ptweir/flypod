#make_enhanced_movie.py
#PTW 10/12/2011 based on make_movie.py
# to make movie, run at command prompt:
# ffmpeg -b 8000000 -r 20 -i frame%05d.png supplementalMovie3Weir.mpeg
# this will result in movie 10x actual speed

#import flypod, sky_times
import motmot.FlyMovieFormat.FlyMovieFormat as FMF
import numpy as np
import os
import pylab, time, imp
from PIL import Image
import matplotlib as mpl
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['svg.embed_char_paths'] = False
pylab.ion()

flypod = imp.load_source('flypod','/home/peter/src/flypod/flypod.py')
sky_times = imp.load_source('sky_times','/home/peter/src/flypod/sky_times.py')

inDirName = '/home/peter/data/grayFilter/fly03/'
outDirName = './frames'

def circle_fit(dataX,dataY):
    """fit a circle to data, returns x,y,radius
    
    arguments:      
    dataX       numpy array containing x data
    dataY       numpy array containing y data (must be same size as dataX)
    
    example:
    cx, cy, r = circle_fit(x,y)
    """ 
    n = sum(~np.isnan(dataX))
    a = np.ones((n,3))
    a[:,0] = dataX[~np.isnan(dataX)]
    a[:,1] = dataY[~np.isnan(dataY)]
    b = -dataX[~np.isnan(dataX)]**2 - dataY[~np.isnan(dataY)]**2
    ai = np.linalg.pinv(a)
    out = np.dot(ai,b)
    circCenterX = -.5*out[0]
    circCenterY = -.5*out[1]
    circR  =  ((out[0]**2+out[1]**2)/4-out[2])**.5;
    return circCenterX, circCenterY, circR
def rose(ax,data,line_color,fill_color):
    """makes polar histogram plot, returns wrapped data, n, bins, binCenters, axes
    
    arguments:
    ax              axes to plot in 
    data            numpy array containing data, must be in radians
    plotArgs        dict of plot arguments (color, etc.)
    
    example:
    orw,n,b,bc,ax = rose(orientations)
    """ 
    NUMBINS = 36

    n, bins, patches = pylab.hist(data,NUMBINS,range=(0,2*np.pi),normed=True,visible=False)
    #n = n*2*np.pi/sum(n)
    binCenters = bins[:-1] + (bins[1:] - bins[:-1])/2
    n = np.append(n,n[0])
    binCenters = np.append(binCenters,binCenters[0])
    lineHandle = ax.plot(binCenters,n,color=line_color)
    fillHandle = ax.fill(binCenters,n,color=fill_color)
    ax.set_rmax(.25*NUMBINS/(2*np.pi))
    ax.set_rgrids([2],'')
    ax.set_thetagrids([0,90,180,270],['','','',''])
    return lineHandle, fillHandle

def make_frames(inDirName,outDirName):
    """saves frames in movie with orientation superimposed
    
    arguments:      
    inDirName       directory with fly, sky, and .fmf files
    outDirName      directory to save frames to
    
    example:
    make_frames('/home/peter/data/grayFilter/fly03/','./frames/')
    """ 
    CHANGEBUFFER = 9.5
    
    fly = flypod.analyze_directory(inDirName)
    sky = sky_times.analyze_directory(inDirName)
    sky['changeTimes'] = sky['changeTimes'] - CHANGEBUFFER
    
    startTime = sky['changeTimes'][0]
    experimentEndTime = sky['changeTimes'][8] - startTime
    
    OFFSET = np.pi/2
    orientationInFrame = fly['orientations'].copy() + 180

    posToPlot = np.mod(-fly['orientations'].copy()*np.pi/180.0 + np.pi/2 + np.pi,2*np.pi)
    times = fly['times'].copy() - startTime
    experimentEndInd = np.argmin(abs(times-experimentEndTime))
    
    arenaHeadingEachFrame = -np.ones(times.shape) # we will store frame-based rotator state in here
    
    COLORS = dict(N=[.9,.9,.9],E=[.8,.8,.8],S=[.7,.7,.7],W=[1,1,1])
    ROTATIONS = dict(N=0,E=-np.pi/2,S=np.pi,W=np.pi/2)
    ROTATIONSDEG = dict(N=0,E=-90,S=180,W=90)
    FRAMESTEP = 145
    
    circleCenterX, circleCenterY, circleRadius = circle_fit(fly['x'],fly['y']) #NOTE:should have been saving this info from start.
    
    fmf = FMF.FlyMovie(os.path.join(inDirName,str(fly['fileName'])))

    nFrames = fmf.get_n_frames()
    timestamps = np.ma.masked_all(len(range(0,int(nFrames),FRAMESTEP)))
    
    fig = pylab.figure(figsize=(5.5,5.5))
    fig.set_facecolor('w')
    fig.text(.05,.96,'10x normal speed, natural sky polarization')
    fig.text(.6,.92,'West up')
    fig.text(.8,.92,'North up')
    fig.text(.6,.62,'East up')
    fig.text(.8,.62,'South up')
    fig.text(.55,.375,'West',color='r')
    frameAxes = fig.add_axes([0.05, 0.4, 0.5, 0.5])
    timePlotAxes = fig.add_axes([0.17, 0.1, 0.8, 0.25])
    polarWAxes = fig.add_axes([0.57, 0.7, 0.2, 0.2],polar=True)
    polarNAxes = fig.add_axes([0.78, 0.7, 0.2, 0.2],polar=True)
    polarEAxes = fig.add_axes([0.57, 0.4, 0.2, 0.2],polar=True)
    polarSAxes = fig.add_axes([0.78, 0.4, 0.2, 0.2],polar=True)
    #polarWAxes.set_rmax(.25*36/(2*np.pi))
    #polarEAxes.set_rmax(.25*36/(2*np.pi))
    
    posToPlotWithoutJumps = posToPlot.copy()
    posToPlotWithoutJumps[abs(np.diff(posToPlotWithoutJumps))>np.pi] = np.nan
    timePlotAxes.plot(times,posToPlotWithoutJumps,'k',linewidth=.5,zorder=10)
    for i, cT in enumerate(np.array(sky['changeTimes'][:-1])-startTime):
        bgpatches = timePlotAxes.axvspan(cT, sky['changeTimes'][i+1]-startTime, facecolor=COLORS[sky['directions'][i]],edgecolor='none')
        timePlotAxes.plot([cT, sky['changeTimes'][i+1]-startTime],np.pi-ROTATIONS[sky['directions'][i]]*np.ones(2),'r',linewidth=3)
        startInd = np.argmin(abs(times-cT))
        endInd = np.argmin(abs(times-(sky['changeTimes'][i+1]-startTime)))
        arenaHeadingEachFrame[startInd:endInd] = ROTATIONSDEG[sky['directions'][i]]
    NPosToPlot = np.ma.masked_where(arenaHeadingEachFrame != 0,posToPlot,copy=True)
    EPosToPlot = np.ma.masked_where(arenaHeadingEachFrame != -90,posToPlot,copy=True)
    WPosToPlot = np.ma.masked_where(arenaHeadingEachFrame != 90,posToPlot,copy=True)
    SPosToPlot = np.ma.masked_where(arenaHeadingEachFrame != 180,posToPlot,copy=True)
    if sum(np.isnan(orientationInFrame)) > 0:
        for trackingErrorTime in times[np.isnan(orientationInFrame)]:
            timePlotAxes.axvline(trackingErrorTime,linewidth=3,color='g')
    if fly.has_key('stopTimes'):
        for i, sT in enumerate(np.array(fly['stopTimes'])-startTime):
            timePlotAxes.axvspan(sT, fly['startTimes'][i]-startTime, facecolor='g')

    timePlotAxes.set_ylim((0,2*np.pi))
    timePlotAxes.set_xlim((0,experimentEndTime))
    timePlotAxes.set_xlabel('time')
    timePlotAxes.spines['left'].set_position(('outward',8))
    timePlotAxes.spines['bottom'].set_position(('outward',8))
    timePlotAxes.spines['top'].set_visible(False)
    timePlotAxes.spines['right'].set_visible(False)
    timePlotAxes.spines['bottom'].set_bounds(0,24*60)
    timePlotAxes.yaxis.set_ticks_position('left')
    timePlotAxes.xaxis.set_ticks_position('bottom')
    timePlotAxes.set_yticks(np.linspace(0,2*np.pi,5))
    timePlotAxes.set_yticklabels([str(int(round(yt*180/np.pi)))+u'\u00b0' for yt in timePlotAxes.get_yticks()])
    timePlotAxes.set_xticks(np.linspace(0,24*60,9))
    timePlotAxes.set_xticklabels([str(int(round(xt/60))) for xt in timePlotAxes.get_xticks()])
    timePlotAxes.set_ylabel('Heading in arena')
    for tickLine in timePlotAxes.get_xticklines() + timePlotAxes.get_yticklines():
        tickLine.set_markeredgewidth(1)
    timePlotAxes.set_xlabel('Time (min)')
    
    for fn,frameNumber in enumerate(range(0,int(experimentEndInd),FRAMESTEP)):
        frame,timestamps[fn] = fmf.get_frame(frameNumber)
        frameAxes.imshow(frame,cmap='gray',origin='lower')

        lineX = [circleCenterX+fly['ROI'][2], 2*(circleCenterX+fly['ROI'][2])-(fly['x'][frameNumber]+fly['ROI'][2])]
        lineY = [circleCenterY+fly['ROI'][1], 2*(circleCenterY+fly['ROI'][1])-(fly['y'][frameNumber]+fly['ROI'][1])]

        frameAxes.plot(lineX,lineY,'b', linewidth=2)
        frameAxes.set_ylim((0,64))
        frameAxes.set_xlim((0,64))
        frameAxes.set_axis_off()
        
        movingDotList = timePlotAxes.plot(times[frameNumber],posToPlot[frameNumber],'ob',markeredgecolor='w',zorder=20)
        movingDot = movingDotList[0]
        
        lineHandle, fillHandle = rose(polarWAxes,WPosToPlot[:frameNumber].compressed(),'k',COLORS['W'])
        lineHandle[0].set_zorder(20)
        fillHandle[0].set_zorder(20)
        lh = polarWAxes.plot([0, np.pi/2],[0, .25*36/(2*np.pi)],'r')
        lh[0].set_linewidth(3)
        polarWAxes.set_rmax(.25*36/(2*np.pi))
        
        lineHandle, fillHandle = rose(polarNAxes,NPosToPlot[:frameNumber].compressed(),'k',COLORS['N'])
        lineHandle[0].set_zorder(20)
        fillHandle[0].set_zorder(20)
        lh = polarNAxes.plot([0,np.pi],[0,.25*36/(2*np.pi)],'r')
        lh[0].set_linewidth(1)
        polarNAxes.set_rmax(.25*36/(2*np.pi))
        
        lineHandle, fillHandle = rose(polarEAxes,EPosToPlot[:frameNumber].compressed(),'k',COLORS['E'])
        lineHandle[0].set_zorder(20)
        fillHandle[0].set_zorder(20)
        lh = polarEAxes.plot([0, -np.pi/2],[0, .25*36/(2*np.pi)],'r',zorder=5)
        lh[0].set_linewidth(3)
        polarEAxes.set_rmax(.25*36/(2*np.pi))
        
        lineHandle, fillHandle = rose(polarSAxes,SPosToPlot[:frameNumber].compressed(),'k',COLORS['S'])
        lineHandle[0].set_zorder(20)
        fillHandle[0].set_zorder(20)
        lh = polarSAxes.plot([0, 0],[0, .25*36/(2*np.pi)],'r',zorder=5)
        lh[0].set_linewidth(1)
        polarSAxes.set_rmax(.25*36/(2*np.pi))
        
        if arenaHeadingEachFrame[frameNumber] == 0: # N
            arrowHandle = frameAxes.plot([55,59],[5,5],'r',linewidth=3)
            frameAxes.plot([55],[5],'+r')
            #frameAxes.arrow(56,4,4,0,color='r')
        elif arenaHeadingEachFrame[frameNumber] == -90: # E
            arrowHandle = frameAxes.plot([59,59],[1,5],'r',linewidth=3)
            frameAxes.plot([59],[1],'+r')
        elif arenaHeadingEachFrame[frameNumber] == 180: # S
            arrowHandle = frameAxes.plot([59,63],[5,5],'r',linewidth=3)
            frameAxes.plot([63],[5],'+r')
        elif arenaHeadingEachFrame[frameNumber] == 90: # W
            arrowHandle = frameAxes.plot([59,59],[5,9],'r',linewidth=3)
            frameAxes.plot([59],[9],'+r')
        #pylab.draw()
        fig.savefig(os.path.join(outDirName,'frame')+("%05d" %(fn+1))+'.png')
        polarWAxes.cla()
        polarNAxes.cla()
        polarEAxes.cla()
        polarSAxes.cla()
        frameAxes.cla()
        movingDot.remove()
    print 'frame rate = ' + str(1/np.mean(np.diff(timestamps)))
    
make_frames(inDirName,outDirName)

