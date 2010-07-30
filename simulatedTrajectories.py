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
    
def cumprobdist(data,xmax=None):
    if xmax is None:
        xmax = numpy.max(data)
    elif xmax < numpy.max(data):
        warnings.warn('value of xmax lower than maximum of data')
        xmax = numpy.max(data)
    num_points = len(data)
    X = numpy.concatenate(([0.0],data,data,[xmax]))
    X.sort()
    Y = numpy.concatenate(([0.0],arange(num_points),arange(num_points)+1,[num_points]))/num_points
    Y.sort()
    line = plot(X,Y)
    return line[0]

rootDir = '/home/cardini/data/'
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
    
CHANGEBUFFER = 10
COLORS = dict(N='b',E='y',S='r',W='g')
ROTATIONS = dict(N=0,E=90,S=180,W=270)
SPEED = 0.5
DMAX = 330
Dist = {}
angSpeed = {}
pylab.figure()
for b, baseDir in enumerate(flies.keys()):
    numFlies = len(flies[baseDir])
    D = numpy.empty(numFlies)
    D.fill(nan)
    S = numpy.empty(numFlies)
    S.fill(nan)
    pylab.subplot(1,len(flies),b+1)
    pylab.hold('on')
    for fNum in range(numFlies):
        print fNum
        dirName=os.path.join(baseDir,dName)
        fly = flies[baseDir][fNum].copy()
        sky = skies[baseDir][fNum].copy()

        orientations = fly['orientations'].copy()
        orientations = orientations + 180
        
        times = fly['times'].copy()
        beginTime = times[0]
        if fly.has_key('stopTimes'):
            for i, sT in enumerate(fly['stopTimes']):
                inds = (times > sT) & (times < fly['startTimes'][i])
                orientations[inds] = numpy.nan
                #timeStopped[fNum] = timeStopped[fNum] + numpy.sum(numpy.diff(times[inds]))
                times[inds] = numpy.nan
                
        num_blocks = 1
        min_time = .9*num_blocks*12*60
        max_time = num_blocks*12*60
        inds = (times - beginTime) > max_time
        orientations[inds] = numpy.nan
        times[inds] = numpy.nan

        if numpy.nansum(numpy.diff(times))>min_time:
            rotations = numpy.empty(times.shape)
            rotations.fill(numpy.nan)
            for i, cT in enumerate(sky['changeTimes'][:-1]):
                inds = (times > cT-CHANGEBUFFER) & (times <= sky['changeTimes'][i+1]-CHANGEBUFFER)
                rotations[inds] = ROTATIONS[sky['directions'][i]]
                
            worldOrientations = orientations+rotations
            
            t=numpy.diff(times)
            aS = abs(numpy.diff(numpy.unwrap(orientations,180)))
            S[fNum] = numpy.median(aS[~numpy.isnan(aS)])
            simVelX = numpy.sin(worldOrientations[1:]*numpy.pi/180)*SPEED*t
            simVelY = numpy.cos(worldOrientations[1:]*numpy.pi/180)*SPEED*t
            simTrajX = simVelX[~numpy.isnan(simVelX)].cumsum()
            simTrajY = simVelY[~numpy.isnan(simVelY)].cumsum()
            D[fNum]=numpy.sqrt(simTrajX[-1]**2 + simTrajY[-1]**2)
            print numpy.nansum(numpy.diff(times))
            #T[fNum] = numpy.nansum(numpy.diff(times))
            l = pylab.plot(simTrajX,simTrajY)
            pylab.scatter(simTrajX[-1],simTrajY[-1],color=l[0].get_color())
            
            changeInds = [abs(times[~numpy.isnan(times)] - (sky['changeTimes'][i]-CHANGEBUFFER)).argmin() for i in range(len(sky['changeTimes']))]
            pylab.plot(simTrajX[changeInds[1:4]],simTrajY[changeInds[1:4]],'.',color=l[0].get_color())

    ax=gca()
    th = numpy.arange(0,2*numpy.pi,.01)
    R1 = 115
    xcirc = R1*cos(th)
    ycirc = R1*sin(th)
    plot(xcirc,ycirc,color='k')
    ax.set_xlim((-DMAX,DMAX))
    ax.set_ylim((-DMAX,DMAX))
    ax.set_aspect('equal')
    ax.set_title(baseDir)
    pylab.draw()
    Dist[baseDir] = D.copy()
    angSpeed[baseDir] = S.copy()

pylab.figure()
pylab.hold('on')
ax = pylab.gca()
for b, baseDir in enumerate(baseDirs):
    line = cumprobdist(Dist[baseDir][~numpy.isnan(Dist[baseDir])],DMAX)
    #line = cumprobdist(angSpeed[baseDir][~numpy.isnan(angSpeed[baseDir])])

    line.set_label(baseDir)
    
ax.axvline(R1,color='k')
ax.legend(loc='lower right')
ax.set_ylim((-.1,1.1))
ax.set_xlim((-25,350))
