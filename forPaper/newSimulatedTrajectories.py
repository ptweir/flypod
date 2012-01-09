import numpy as np
import pylab
pylab.ion()

def get_time_stopped(fly, beginTime, endTime):
    """calculate total time stopped between beginTime and endTime
    """
    if fly.has_key('stopTimes'):
        a = np.ma.masked_all((2,np.size(fly['startTimes'])))
        a[0,:] = np.copy(fly['startTimes'])
        a[1,:] = beginTime*np.ones(np.shape(fly['startTimes']))
        startTimes = np.max(a,0)
        a[0,:] = np.copy(startTimes)
        a[1,:] = endTime*np.ones(np.shape(fly['startTimes']))
        startTimes = np.min(a,0)
        
        a[0,:] = np.copy(fly['stopTimes'])
        a[1,:] = beginTime*np.ones(np.shape(fly['stopTimes']))
        stopTimes = np.max(a,0)
        a[0,:] = np.copy(stopTimes)
        a[1,:] = endTime*np.ones(np.shape(fly['stopTimes']))
        stopTimes = np.min(a,0)
        time_stopped = np.sum(startTimes - stopTimes)
    else:
        time_stopped = 0
        
    return time_stopped
    
def get_num_stops(fly, beginTime, endTime):
    """calculate total number of stops between beginTime and endTime
    """
    #num_stops = sum([(sT <= endTime) & (sT >= beginTime) for sT in stopTimes])
    num_stops = 0
    if fly.has_key('stopTimes'):
        stopTimes = np.copy(fly['stopTimes'])
        startTimes = np.copy(fly['startTimes'])
        for i, stopTime in enumerate(stopTimes):
            if stopTime >= beginTime and stopTime <= endTime:
                num_stops += 1
            elif stopTime < beginTime and startTimes[i] >= beginTime:
                num_stops += 1 # started stopped
  
    return num_stops
def cumprobdist(ax,data,xmax=None,plotArgs={}):
    if xmax is None:
        xmax = np.max(data)
    elif xmax < np.max(data):
        warnings.warn('value of xmax lower than maximum of data')
        xmax = np.max(data)
    num_points = len(data)
    X = np.concatenate(([0.0],data,data,[xmax]))
    X.sort()
    X = X[-1::-1]
    Y = np.concatenate(([0.0],np.arange(num_points),np.arange(num_points)+1,[num_points]))/num_points
    Y.sort()
    line = axIndv.plot(X,Y,**plotArgs)
    return line[0]

def plotIt(flies,skies,axIndv,axAll,maxNumStops,maxTimeStopped,changebuffer):

    #baseDirs = ['indoor/halogen','circularPolarizer','noFilter','blueFilter','grayFilter']
    baseDirs = ['indoor/dark','circularPolarizer','noFilter','blueFilter','grayFilter']
    PLOTCOLORS = [[1,.5,0],[1,0,0],[0,.9,0],[0,0,1],[.2,.2,.2]]
    MIN_EXPERIMENT_TIME = 11
    #maxNumStops = 2
    #maxTimeStopped = 60
    #ROTATIONS = dict(N=0,E=-np.pi/2,S=np.pi,W=np.pi/2) #CHECK
    ROTATIONS = dict(N=0,E=-np.pi/2,S=-np.pi,W=-3*np.pi/2) #CHECK
    SPEED = 0.5
    DMAX = 460
    THETA_FOR_CIRCLE = np.arange(0,2*np.pi,.01)
    RADIUS_FOR_CIRCLE = 100 #400
    CIRCLE_X = RADIUS_FOR_CIRCLE*np.cos(THETA_FOR_CIRCLE)
    CIRCLE_Y = RADIUS_FOR_CIRCLE*np.sin(THETA_FOR_CIRCLE)
    
    allDistances = list()
    #fig = pylab.figure()
    #fig.set_facecolor('w')

    for baseDirInd, baseDir in enumerate(baseDirs):
        #ax = fig.add_subplot(3,6,12+baseDirInd+1)
        #axIndv.set_position([.02+baseDirInd*.05,.1,.3,.3])
        axIndv.fill(CIRCLE_X+DMAX*baseDirInd-200,CIRCLE_Y,color=[.9,.9,.9])
        
        numFlies = len(flies[baseDir])
        distances = np.ma.masked_all((numFlies))
        for flyInd in range(numFlies):
            fly = flies[baseDir][flyInd].copy()
            sky = skies[baseDir][flyInd].copy()
            
            if fly['times'][-1] - fly['times'][0] > MIN_EXPERIMENT_TIME*60:
                alpha = np.ma.masked_invalid(np.mod(-fly['orientations']*np.pi/180.0 + np.pi/2 + np.pi,2*np.pi),copy=True)
                time = np.ma.masked_invalid(fly['times'],copy=True)
                experimentStartTime = sky['changeTimes'][0]
                experimentEndTime = sky['changeTimes'][4] - changebuffer
                timeStopped = get_time_stopped(fly,experimentStartTime,experimentEndTime)
                numStops = get_num_stops(fly,experimentStartTime,experimentEndTime)
                if timeStopped <= maxTimeStopped and numStops <= maxNumStops:
                    rotations = np.ma.masked_all(time.shape)
                    for changeTimeInd, changeTime in enumerate(sky['changeTimes'][:-1]):
                        changeInds = (time > changeTime-changebuffer) & (time <= sky['changeTimes'][changeTimeInd+1]-changebuffer)
                        rotations[changeInds] = ROTATIONS[sky['directions'][changeTimeInd]]
                    
                    worldOrientations = alpha+rotations
                    worldOrientations[time>experimentEndTime] = np.ma.masked
                    time[time>experimentEndTime] = np.ma.masked
                    if fly.has_key('stopTimes'):
                        for stopTimeInd, stopTime in enumerate(fly['stopTimes']):
                            stoppedInds = (time > stopTime) & (time < fly['startTimes'][stopTimeInd])
                            worldOrientations[stoppedInds] = np.ma.masked
                            time[stoppedInds] = np.ma.masked
                    timeSteps=np.diff(time)
                    simVelX = (np.cos(worldOrientations[1:])*SPEED*timeSteps).compressed()
                    simVelY = (np.sin(worldOrientations[1:])*SPEED*timeSteps).compressed()
                    simTrajX = simVelX.cumsum()
                    simTrajY = simVelY.cumsum()
                    distances[flyInd]=np.sqrt(simTrajX[-1]**2 + simTrajY[-1]**2)
                    lines = axIndv.plot(simTrajX+DMAX*baseDirInd-200,simTrajY,color=PLOTCOLORS[baseDirInd])
                    lines[0].set_linewidth(.5)
                    #changeTimeInds = [np.argmin(abs(time.compressed() - (sky['changeTimes'][i]-changebuffer))) for i in range(len(sky['changeTimes']))]
                    #axIndv.plot(simTrajX[changeTimeInds[1:4]],simTrajY[changeTimeInds[1:4]],'.',color=l[0].get_color())
                    #axIndv.scatter(simTrajX[-1],simTrajY[-1],edgecolors=l[0].get_color(),facecolors='k')
                    #dots = axIndv.plot(simTrajX[-1]+DMAX*baseDirInd-200,simTrajY[-1],'o',markeredgecolor=PLOTCOLORS[baseDirInd],markerfacecolor='k')
                    dots = axIndv.plot(simTrajX[-1]+DMAX*baseDirInd-200,simTrajY[-1],'o',markeredgecolor=lines[0].get_color(),markerfacecolor='k')
                    dots[0].set_markersize(3)
        #axIndv.set_xlim((-DMAX,DMAX))
        #axIndv.set_ylim((-DMAX,DMAX))
        axIndv.set_aspect('equal')
        #axIndv.set_title(baseDir)
        axIndv.set_axis_off()
        #lineHandle = cumprobdist(axAll,distances.compressed(),xmax=DMAX)
        #lineHandle.set_color(PLOTCOLORS[baseDirInd])
        boxplotHandle = axAll.boxplot(distances.compressed(),positions=[baseDirInd])
        boxHandle = boxplotHandle['boxes'][0]
        boxHandle.set_color(PLOTCOLORS[baseDirInd])
        for whiskerHandle in boxplotHandle['whiskers']:
            whiskerHandle.set_linestyle('-')
            whiskerHandle.set_color('k')
        for capHandle in boxplotHandle['caps']:
            capHandle.set_visible(False)
        for flierHandle in boxplotHandle['fliers']:
            flierHandle.set_marker('+')
            flierHandle.set_markerfacecolor(PLOTCOLORS[baseDirInd])
            flierHandle.set_markeredgecolor(PLOTCOLORS[baseDirInd])              
        allDistances.append(distances.compressed())
    axIndv.annotate('N', xy=(30,250), xytext=(30,150),arrowprops=dict(width=.5,headwidth=2,facecolor='black',shrink=0.05))
    axIndv.annotate(str(RADIUS_FOR_CIRCLE)+' m', xy=(-200, -RADIUS_FOR_CIRCLE), xytext=(-100,-RADIUS_FOR_CIRCLE - 100),arrowprops=dict(width=.5,headwidth=2,facecolor='black',shrink=0.05))
    axAll.set_xlim((-.5,len(baseDirs))) # axAll.set_xlim((-1,len(baseDirs)-.5))
    axAll.spines['top'].set_visible(False)
    axAll.spines['bottom'].set_visible(False)
    axAll.spines['left'].set_visible(False)
    try:
        #axAll.spines['left'].set_bounds(0,300)
        axAll.spines['right'].set_bounds(0,300)
    except AttributeError:
        print "running old matplotlib version"
    axAll.yaxis.set_ticks_position('right')
    axAll.xaxis.set_ticks_position('bottom')
    XLABELS = ['D','C','N','B','G']
    axAll.set_xticks(range(5)) #axAll.set_xticks([])
    axAll.set_xticklabels(XLABELS) #axAll.set_xticklabels([],visible=False)
    axAll.set_yticks(np.linspace(0,300,4))
    axAll.text(1.65+len(baseDirs),300,'Fictive distance (m)',rotation=-90) # axAll.set_ylabel('Fictive distance (m)')
    for tickLine in axAll.get_xticklines() + axAll.get_yticklines():
        tickLine.set_markeredgewidth(1)
    return allDistances
    
