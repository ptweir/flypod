import time, numpy, pylab, os
import flypod, sky_times
from scipy.stats.morestats import circmean, circvar	
pylab.ion()

def compare_dir_names(dn1, dn2):
    n1, n2 = int(dn1[3:]),int(dn2[3:])
    return cmp(n1,n2)

COLORS = dict(N='b',E='y',S='r',W='g')
ROTATIONS = dict(N=0,E=90,S=180,W=270)

MAX_TIME = 12*3*60+20
MIN_TIME = 12*1*60

FPS = 291
SPEED = 0.5

baseDirs = ['/home/cardini/data/grayFilter','/home/cardini/data/circularPolarizer']
circVar, distances, totalTimes = {}, {}, {}
for b, baseDir in enumerate(baseDirs):
    dNames = os.listdir(baseDir)
    flyDirNames = [D for D in dNames if D[:3] == 'fly']
    flyDirNames.sort(compare_dir_names)
    M = numpy.empty(len(flyDirNames))
    M.fill(nan)
    V = numpy.empty(len(flyDirNames))
    V.fill(nan)
    D = numpy.empty(len(flyDirNames))
    D.fill(nan)
    T = numpy.empty(len(flyDirNames))
    T.fill(nan)
    totalTimeRecording = numpy.zeros(len(flyDirNames))
    timeStopped = numpy.zeros(len(flyDirNames))
    #meanAngSp = numpy.zeros(len(flyDirNames))

    pylab.figure()
    for dNum, dName in enumerate(flyDirNames):
        dirName=os.path.join(baseDir,dName)
        fly = flypod.analyze_directory(dirName)
        sky = sky_times.analyze_directory(dirName)

        worldTotals=numpy.array([])
        totals = dict(N=numpy.array([]),E=numpy.array([]),S=numpy.array([]),W=numpy.array([]))

        orientations = numpy.copy(fly['orientations'])
        orientations = orientations + 180
        #orientations = numpy.unwrap(orientations,180)
        
        times = numpy.copy(fly['times'])
        if fly.has_key('stopTimes'):
            for i, sT in enumerate(fly['stopTimes']):
                inds = (times > sT) & (times < fly['startTimes'][i])
                orientations[inds] = numpy.nan
                timeStopped[dNum] = timeStopped[dNum] + numpy.sum(numpy.diff(times[inds]))
                
        totalTimeRecording[dNum] = (times[-1] - times[0])
        inds = (times - times[0]) > MAX_TIME
        orientations[inds] = numpy.nan
        
        for i, cT in enumerate(sky['changeTimes'][:-1]):
            inds = (times > cT) & (times < sky['changeTimes'][i+1])
            ors = orientations[inds]
            ors = ors[~numpy.isnan(ors)]
            if len(ors)>0:
                totals[sky['directions'][i]] = numpy.concatenate((totals[sky['directions'][i]],ors))
        
        for i, d in enumerate(COLORS):
            worldTotals = numpy.concatenate((worldTotals,totals[d]+ROTATIONS[d]))
        
        if totalTimeRecording[dNum] > MIN_TIME and (totalTimeRecording[dNum] - timeStopped[dNum]) > MIN_TIME*.8:
            M[dNum]=circmean(worldTotals,high=180,low=-180)
            V[dNum]=circvar(worldTotals*numpy.pi/180,high=numpy.pi,low=-numpy.pi)
            #meanAngSp[dNum] = numpy.mean(abs(numpy.diff(numpy.unwrap(orientations[~numpy.isnan(orientations)],180))))
            #V[dNum]=numpy.mean(abs(numpy.diff(worldTotals)))
            #fig = pylab.figure()
            pylab.subplot(4,5,dNum+1,polar=True)
            plotArgs = dict(color='k',linewidth=2)
            orw,n,b,bc,ax = flypod.rose(worldTotals,plotArgs=plotArgs)
            pylab.hold('on')
            
            #n={}
            for i, d in enumerate(COLORS):
                plotArgs = dict(color=COLORS[d],linewidth=.5)
                flypod.rose(totals[d]+ROTATIONS[d],plotArgs=plotArgs)
                #NUMBINS = 8
                #n[d], bins, patches = pylab.hist(numpy.mod(totals[d]+ROTATIONS[d],360),bins=numpy.arange(NUMBINS+1)*360/NUMBINS,range=(0,360),normed=True,visible=False)
                #bins[:-1]+numpy.diff(bins)/2.0
            pylab.polar([0,numpy.pi/2-M[dNum]*numpy.pi/180],[0,ax.get_rmax()*(1-V[dNum])],color='k')
            ax.set_rmax(.15)
            ax.set_rgrids([1],'')
            ax.set_thetagrids([0,90,180,270],['E','N','W','S'])
            pylab.title(dName)
            pylab.draw()
            
            #sim traj stuff:
            rotations = numpy.empty(fly['times'].shape)
            rotations.fill(numpy.nan)
            for i, cT in enumerate(sky['changeTimes'][:-1]):
                inds = (times > cT) & (times <= sky['changeTimes'][i+1])
                rotations[inds] = ROTATIONS[sky['directions'][i]]
                
            worldOrientations = orientations+rotations

            simVelX = numpy.sin(worldOrientations*numpy.pi/180)*SPEED/FPS
            simVelY = numpy.cos(worldOrientations*numpy.pi/180)*SPEED/FPS
            simTrajX = simVelX[~numpy.isnan(simVelX)].cumsum()
            simTrajY = simVelY[~numpy.isnan(simVelY)].cumsum()
            D[dNum]=numpy.sqrt(simTrajX[-1]**2 + simTrajY[-1]**2)
            times[isnan(orientations)] = numpy.nan
            T[dNum] = numpy.nansum(numpy.diff(times))

    pylab.suptitle(baseDir)
    circVar[baseDir] = numpy.copy(V)
    distances[baseDir] = numpy.copy(D)
    totalTimes[baseDir] = numpy.copy(T)

pylab.figure()
pylab.hold('on')
ax = pylab.gca()
for b, baseDir in enumerate(baseDirs):
    pylab.scatter(circVar[baseDir]*0+b,circVar[baseDir])
    
ax.set_xticks(range(len(baseDirs)+1))
ax.set_xticklabels(baseDirs)
