import numpy as np
import imp
import pylab
pylab.ion()
circ = imp.load_source('circ','/home/peter/src/circ/circ.py')
fisher = imp.load_source('fisher','/home/peter/src/fisher/fisher.py')
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
    line = ax.plot(X,Y,**plotArgs)
    return line[0]

def plotIt(ax,flies,skies,baseDirs,maxNumStops,maxTimeStopped,changebuffer,numChunks,maxTimeStoppedChunk):

    PLOTCOLORS = [[1,.5,0],[1,0,0],[0,.9,0],[0,0,1],[.2,.2,.2]]
    MIN_EXPERIMENT_TIME = 11
    ROTATIONS = dict(N=0,E=-np.pi/2,S=-np.pi,W=-3*np.pi/2) #CHECK
    
    allDirsTestStatistics = list()
    allDirsPValues = list()

    for baseDirInd, baseDir in enumerate(baseDirs):
        fig = pylab.figure()
        axInd = 0
        
        numFlies = len(flies[baseDir])
        meanAngleChunks = np.ma.masked_all((numFlies,numChunks))
        pValues = np.ma.masked_all(numFlies)
        testStatistics = np.ma.masked_all(numFlies)
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
                    axInd += 1
                    ax0 = fig.add_subplot(5,5,axInd,polar=True)
                    rotations = np.ma.masked_all(time.shape)
                    if fly.has_key('stopTimes'):
                        for stopTimeInd, stopTime in enumerate(fly['stopTimes']):
                            stoppedInds = (time > stopTime) & (time < fly['startTimes'][stopTimeInd])
                            alpha[stoppedInds] = np.ma.masked
                    
                    for changeTimeInd, changeTime in enumerate(sky['changeTimes'][:4]):
                        meanAngleChunks
                        trialStartTime = changeTime-changebuffer
                        trialEndTime = sky['changeTimes'][changeTimeInd+1]-changebuffer
                        trialInds = (time > trialStartTime) & (time <= trialEndTime)
                        rotations[trialInds] = ROTATIONS[sky['directions'][changeTimeInd]]

                    worldOrientations = alpha+rotations
                    chunkDuration = (experimentEndTime - experimentStartTime)/numChunks
                    for chunkInd in range(numChunks):
                        chunkStartTime = experimentStartTime + chunkInd*chunkDuration
                        chunkStartInd = np.argmin(abs(time-chunkStartTime))
                        chunkEndTime = chunkStartTime + chunkDuration
                        chunkEndInd = np.argmin(abs(time-chunkEndTime))
                        
                        timeStoppedChunk = get_time_stopped(fly, chunkStartTime, chunkEndTime)
                        if timeStoppedChunk <= maxTimeStoppedChunk:
                            meanAngleChunks[flyInd,chunkInd] = circ.circmean(worldOrientations[chunkStartInd:chunkEndInd])
                                
                    testStat, p_val = fisher.rayleigh_test(meanAngleChunks[flyInd,:].compressed(),None)
                    testStatistics[flyInd] = testStat
                    pValues[flyInd] = p_val
                    ax0.scatter(meanAngleChunks[flyInd,:].compressed(),1+0*meanAngleChunks[flyInd,:].compressed())
                    ax0.set_title([str(flyInd), ' p=', str(round(p_val,4))])
                    
        allDirsTestStatistics.append(testStatistics.compressed())
        allDirsPValues.append(pValues.compressed())

        plotArgs = dict(color=PLOTCOLORS[baseDirInd])
        cpLine = cumprobdist(ax,1-pValues.compressed(),xmax=1.1,plotArgs=plotArgs)
        ax.axvline(.95,color='k')
    
    pylab.draw()
    return allDirsTestStatistics, allDirsPValues
                                
                        

                                            




