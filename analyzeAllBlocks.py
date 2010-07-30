import time, numpy, pylab, os
import flypod, sky_times
pylab.ion()

def compare_dir_names(dn1, dn2):
    n1, n2 = int(dn1[3:]),int(dn2[3:])
    return cmp(n1,n2)
    
def circmean(alpha):
    mean_angle = numpy.arctan2(numpy.mean(numpy.sin(alpha)),numpy.mean(numpy.cos(alpha)))
    return mean_angle

def circvar(alpha):
    R = numpy.sqrt(numpy.sum(numpy.sin(alpha))**2 + numpy.sum(numpy.cos(alpha))**2)/len(alpha)
    V = 1-R
    return V

rootDir = '/home/cardini/data/'
#baseDirs = ['grayFilter']
#baseDirs = ['grayFilter','circularPolarizer']
baseDirs = ['grayFilter','circularPolarizer','noFilter' ]
try:
    flies
except NameError:
    print "loading flies"
    skies, flies = {}, {}
    for b, baseDir in enumerate(baseDirs):
        sks, fls = {}, {}
        dNames = os.listdir(os.path.join(rootDir,baseDir))
        flyDirNames = [D for D in dNames if D[:3] == 'fly']
        flyDirNames.sort(compare_dir_names)

        for dNum, dName in enumerate(flyDirNames):
            dirName=os.path.join(rootDir,baseDir,dName)
            print dirName
            fly = flypod.analyze_directory(dirName)
            fls[dNum] = fly
            sky = sky_times.analyze_directory(dirName)
            sks[dNum] = sky
            
        flies[baseDir] = fls
        skies[baseDir] = sks
    
COLORS = dict(N='b',E='y',S='r',W='g')
ROTATIONS = dict(N=0,E=90,S=180,W=270)

PRECHANGE_BUFFER = 10
POSTCHANGE_BUFFER = 10

circVar = {}

for b, baseDir in enumerate(flies.keys()):
    numFlies = len(flies[baseDir])
    M = numpy.empty(numFlies)
    M.fill(nan)
    V = numpy.empty(numFlies)
    V.fill(nan)
    totalTimeRecording = numpy.zeros(numFlies)
    timeStopped = numpy.zeros(numFlies)
    #meanAngSp = numpy.zeros(numFlies)

    pylab.figure()
    #for fNum in range(3):
    for fNum in range(numFlies):
        print fNum
        dirName=os.path.join(baseDir,dName)
        fly = flies[baseDir][fNum]
        sky = skies[baseDir][fNum]

        worldTotals=numpy.array([])
        totals = dict(N=numpy.array([]),E=numpy.array([]),S=numpy.array([]),W=numpy.array([]))

        orientations = numpy.copy(fly['orientations'])
        orientations = orientations + 180
        #orientations = numpy.unwrap(orientations,180)
        
        times = numpy.copy(fly['times'])
        beginTime = times[0]
        if fly.has_key('stopTimes'):
            for i, sT in enumerate(fly['stopTimes']):
                inds = (times > sT) & (times < fly['startTimes'][i])
                orientations[inds] = numpy.nan
                timeStopped[fNum] = timeStopped[fNum] + numpy.sum(numpy.diff(times[inds]))
                times[inds] = numpy.nan
                
        totalTimeRecording[fNum] = (times[-1] - beginTime)
        #num_blocks = floor((len(sky['changeTimes'])-1)/4)
        num_blocks = 1
        min_time = .9*num_blocks*12*60
        max_time = num_blocks*12*60
        inds = (times - beginTime) > max_time
        orientations[inds] = numpy.nan
        times[inds] = numpy.nan
        while numpy.nansum(numpy.diff(times)) < min_time:
            print fNum, num_blocks
            num_blocks = num_blocks - 1
            min_time = .9*num_blocks*12*60
            max_time = num_blocks*12*60
            inds = (times - beginTime) > max_time
            orientations[inds] = numpy.nan
            times[inds] = numpy.nan
        if num_blocks > 0:
            for i, cT in enumerate(sky['changeTimes'][:-1]):
                inds = (times > cT+POSTCHANGE_BUFFER) & (times < sky['changeTimes'][i+1]-PRECHANGE_BUFFER)
                ors = orientations[inds]
                ors = ors[~numpy.isnan(ors)]
                if len(ors)>0:
                    totals[sky['directions'][i]] = numpy.concatenate((totals[sky['directions'][i]],ors))
            
            for i, d in enumerate(COLORS):
                worldTotals = numpy.concatenate((worldTotals,totals[d]+ROTATIONS[d]))
            
            M[fNum]=circmean(worldTotals*numpy.pi/180)
            V[fNum]=circvar(worldTotals*numpy.pi/180)
            #pylab.figure()
            pylab.subplot(4,5,fNum+1,polar=True)
            plotArgs = dict(color='k',linewidth=2)
            orw,n,b,bc,ax = flypod.rose(worldTotals,plotArgs=plotArgs)
            pylab.hold('on')
            
            for i, d in enumerate(COLORS):
                plotArgs = dict(color=COLORS[d],linewidth=.5)
                flypod.rose(totals[d]+ROTATIONS[d],plotArgs=plotArgs)

            pylab.polar([0,numpy.pi/2-M[fNum]],[0,ax.get_rmax()*(1-V[fNum])],color='k')
            ax.set_rmax(.15)
            ax.set_rgrids([1],'')
            ax.set_thetagrids([0,90,180,270],['E','N','W','S'])
            pylab.title(fly['dirName'])
            

    pylab.suptitle(baseDir)
    circVar[baseDir] = numpy.copy(V)
    pylab.draw()

pylab.figure()
pylab.hold('on')
ax = pylab.gca()
for b, baseDir in enumerate(baseDirs):
    pylab.scatter((numpy.random.rand(len(circVar[baseDir]))-.5)/3+b,circVar[baseDir],marker='+')
    
ax.set_xticks(range(len(baseDirs)+1))
ax.set_xticklabels(baseDirs)

pylab.figure()
pylab.hold('on')
ax = pylab.gca()
for b, baseDir in enumerate(baseDirs):
    n,bins,p=pylab.hist(circVar[baseDir][~numpy.isnan(circVar[baseDir])],arange(11)/10.0,visible=False)
    binCenters = bins[:-1] + (bins[1:] - bins[:-1])/2
    plot(binCenters,n)
    
    
