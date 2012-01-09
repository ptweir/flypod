import numpy as np
import pylab
import imp
pylab.ion()
fisher = imp.load_source('fisher','/home/peter/src/fisher/fisher.py')
circ = imp.load_source('circ','/home/peter/src/circ/circ.py')
    
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
    
def plotIt(flies,skies,fig,maxNumStops,maxTimeStopped,changebuffer):
    #baseDirs = ['indoor/halogen','circularPolarizer','noFilter','blueFilter','grayFilter']
    baseDirs = ['indoor/dark','circularPolarizer','noFilter','blueFilter','grayFilter']
    TITLES = ['Dark','Circular polarizer','No filter','Blue filter','Gray filter']
    FONTDICT = {'fontsize':8}
    PLOTCOLORS = [[1,.5,0],[1,0,0],[0,.9,0],[0,0,1],[.2,.2,.2]]
    SECONDS_BEFORE = 10
    SECONDS_AFTER = 30
    NUM_SAMPLES = 2000
    NUM_ROTATIONS = 3
    MIN_EXPERIMENT_TIME = 11
    NUM_BOOTSTRAPS = 1000
    #maxNumStops = 8
    #maxTimeStopped = 60
    SUBPLOT_HEIGHT = .16
    allAvgOrientationChange = list()
    #fig = pylab.figure()
    ax1 = fig.add_subplot(4,6,18)
    ax1.set_position([.81,.36,.12,SUBPLOT_HEIGHT]) # ax1.set_position([.86,.4,.13,SUBPLOT_HEIGHT])
    for baseDirInd, baseDir in enumerate(baseDirs):
        numFlies = len(flies[baseDir])
        orientations = np.ma.masked_all((numFlies,NUM_SAMPLES,NUM_ROTATIONS))
        times = np.ma.masked_all((numFlies,NUM_SAMPLES,NUM_ROTATIONS))
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
                    for rotationNumber in range(NUM_ROTATIONS):
                        changeTime = sky['changeTimes'][rotationNumber+1] - changebuffer
                        changeInd = np.argmin(abs(time - changeTime))
                        startTime = changeTime - SECONDS_BEFORE
                        startInd = np.argmin(abs(time - startTime))
                        endTime = changeTime + SECONDS_AFTER
                        endInd = np.argmin(abs(time - endTime))
                        indicesToSample = np.array(np.round(np.linspace(startInd,endInd,NUM_SAMPLES)),dtype=int)
                        #timesToSample = np.linspace(changeTime - SECONDS_BEFORE,changeTime + SECONDS_AFTER, NUM_SAMPLES)
                        #numStops = get_num_stops(fly,timesToSample[0],timesToSample[-1])
                        #indicesToSample = np.zeros(NUM_SAMPLES,dtype=int)
                        #for timeToSampleInd, timeToSample in enumerate(timesToSample):
                        #    indicesToSample[timeToSampleInd] = np.argmin(abs(time - timeToSample))
                        numStopsDuring = get_num_stops(fly,startTime,endTime)
                        if numStopsDuring == 0:
                            da = alpha[indicesToSample] - alpha[changeInd]
                            orientations[flyInd,:,rotationNumber] = np.arctan2(np.sin(da),np.cos(da))
                            #orientations[flyInd,:,rotationNumber] = alpha[indicesToSample] - alpha[changeInd]
                            times[flyInd,:,rotationNumber] = time[indicesToSample] - time[changeInd]
                            #orientations[flyInd*NUM_ROTATIONS+rotationNumber,:] = alpha[changeInd - INDS_BEFORE:changeInd + INDS_AFTER] - alpha[changeInd]
                            #times[flyInd*NUM_ROTATIONS+rotationNumber,:] = time[changeInd - INDS_BEFORE:changeInd + INDS_AFTER] - time[changeInd]
                    
        flyOrientations = circ.circmean(orientations,2)
        flyTimes = np.mean(times,2)
        ax = fig.add_subplot(4,6,12+baseDirInd+1)
        ax.set_position([.1+baseDirInd*.145,.36,.095,SUBPLOT_HEIGHT]) # ax.set_position([.1+baseDirInd*.155,.4,.105,SUBPLOT_HEIGHT])
        patchHandle = ax.fill_between(np.mean(flyTimes,0),circ.circmean(flyOrientations,0)-circ.circvar(flyOrientations,0),circ.circmean(flyOrientations,0)+circ.circvar(flyOrientations,0))
        patchHandle.set_facecolor([.8, .8, .8])
        patchHandle.set_edgecolor('none')
        lineHandle = ax.plot(np.mean(flyTimes,0),circ.circmean(flyOrientations,0),color=PLOTCOLORS[baseDirInd])
        #ax.plot(np.mean(flyTimes,0),circ.circmean(flyOrientations,0)+circ.circvar(flyOrientations,0))
        #ax.plot(np.mean(flyTimes,0),circ.circmean(flyOrientations,0)-circ.circvar(flyOrientations,0))
        #ax.plot(flyTimes.transpose(),flyOrientations.transpose())
        ax.set_title(TITLES[baseDirInd],FONTDICT)
        ax.set_xlim((-SECONDS_BEFORE,SECONDS_AFTER))
        ax.spines['left'].set_position(('outward',8))
        try:
            if baseDir == 'indoor/dark':
                ax.set_ylim((-np.pi,np.pi))
                ax.spines['left'].set_bounds(-np.pi/2,np.pi/2)
                ax.set_yticks(np.linspace(-np.pi/2,np.pi/2,3))
                ax.set_ylabel('Heading change')
            else:
                ax.set_ylim((-.7,2.1))
                ax.spines['left'].set_bounds(0,np.pi/2)
                ax.set_yticks(np.linspace(0,np.pi/2,2))
        except AttributeError:
            print "running old matplotlib version"
            ax.set_ylim((-np.pi/2,np.pi))
            ax.set_yticks(np.linspace(-np.pi/2,np.pi,4))
        #ax.set_yticklabels([int(yt) for yt in np.round(ax.get_yticks()*180/np.pi)])
        #ax.set_yticklabels(np.array(ax.get_yticks()*180/np.pi,dtype=int))
        ax.set_yticklabels([str(int(round(yt*180/np.pi)))+u'\u00b0' for yt in ax.get_yticks()])
        ax.spines['bottom'].set_position(('outward',8))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xticks(np.array(np.linspace(-10,30,5),dtype=int))
        for tickLine in ax.get_xticklines() + ax.get_yticklines():
            tickLine.set_markeredgewidth(1)
        ax.set_xlabel('Time (s)')
        flyOrientations[flyTimes<10] = np.ma.masked
        flyOrientations[flyTimes>30] = np.ma.masked
        avgOrientationChange = circ.circmean(flyOrientations,1)

        #boxplotHandle = ax1.boxplot(avgOrientationChange.compressed(),positions=[baseDirInd])
        #boxHandle = boxplotHandle['boxes'][0]
        #boxHandle.set_color(PLOTCOLORS[baseDirInd])
        
        #lower, upper = fisher.confidence_interval_for_mean_direction(avgOrientationChange.compressed())
        lowers, uppers = fisher.confidence_interval_for_mean_direction(avgOrientationChange.compressed(),B=NUM_BOOTSTRAPS,alpha=[.001,.01,.05])
        #P, Rbar = fisher.rayleigh_test(avgOrientationChange.compressed())
        #print baseDir, P
        lineHandle = ax1.plot(baseDirInd*np.ones(2),[lowers[2], uppers[2]])
        meanDotHandle = ax1.plot(baseDirInd*np.ones(1),[circ.circmean(avgOrientationChange.compressed())],'o')
        #dotsHandle = ax1.plot(.2+baseDirInd*np.ones(avgOrientationChange.compressed().shape),avgOrientationChange.compressed(),'.')
        meanDotHandle[0].set_color(PLOTCOLORS[baseDirInd])
        #dotsHandle[0].set_color(PLOTCOLORS[baseDirInd])
        lineHandle[0].set_color(PLOTCOLORS[baseDirInd])
        #lineHandle[0].set_color('k')
        allAvgOrientationChange.append(avgOrientationChange.compressed())
        print baseDir,circ.circmean(avgOrientationChange.compressed()),lowers,uppers
        if lowers[0] > 0 or uppers[0] < 0:
            ax1.text(baseDirInd-.25,2.1,'***')
        elif lowers[1] > 0 or uppers[1] < 0:
            ax1.text(baseDirInd-.25,2.1,'**')        
        elif lowers[2] > 0 or uppers[2] < 0:
            ax1.text(baseDirInd-.25,2.1,'*')
        else:
            ax1.text(baseDirInd-.25,2.1,'NS')
    ax1.set_xlim((-.5,len(baseDirs))) # ax1.set_xlim((-1,len(baseDirs)-.5)) 
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_position(('outward',8))
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.set_ylim(ax.get_ylim())
    try:
        ax1.spines['left'].set_bounds(ax.spines['left'].get_bounds()[0],ax.spines['left'].get_bounds()[1])
    except AttributeError:
        print "running old matplotlib version"
    ax1.set_yticks(ax.get_yticks())
    #ax1.set_yticklabels([int(yt) for yt in np.round(ax1.get_yticks()*180/np.pi)])
    #ax1.set_yticklabels(np.array(ax1.get_yticks()*180/np.pi,dtype=int))
    ax1.set_yticklabels([str(int(round(yt*180/np.pi)))+u'\u00b0' for yt in ax1.get_yticks()])
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_xticks([])
    ax1.set_xticklabels([],visible=False)
    for tickLine in ax1.get_xticklines() + ax1.get_yticklines():
        tickLine.set_markeredgewidth(1)
        
    return allAvgOrientationChange
