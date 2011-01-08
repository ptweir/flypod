import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['svg.embed_char_paths'] = False

import time, numpy, pylab, os
import flypod, sky_times, tools
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
    
def cumprobdist(ax,data,xmax=None,plotArgs={}):
    if xmax is None:
        xmax = numpy.max(data)
    elif xmax < numpy.max(data):
        warnings.warn('value of xmax lower than maximum of data')
        xmax = numpy.max(data)
    num_points = len(data)
    X = numpy.concatenate(([0.0],data,data,[xmax]))
    X.sort()
    X = X[-1::-1]
    Y = numpy.concatenate(([0.0],arange(num_points),arange(num_points)+1,[num_points]))/num_points
    Y.sort()
    line = ax.plot(X,Y,**plotArgs)
    return line[0]
    
def detect_saccades0(orientations, timestep):
    orientations = orientations[~numpy.isnan(orientations)]

    unwrappedOrientations = numpy.unwrap(orientations,180)
    SMOOTH_TIME = .5
    smoothWindow = int(round(SMOOTH_TIME/timestep))
    smoothedOrientations = tools.smooth(unwrappedOrientations,smoothWindow)
    angSpds = abs(numpy.diff(smoothedOrientations))/timestep
    sacInds = tools.local_minima(-angSpds) & (angSpds > 10)
    saccades = angSpds[sacInds]
    #sacs = tools.zigzag(smoothedOrientations)
    #sacInds = abs(sacs) > 10
    #saccades = sacs[sacInds]
    return saccades

def detect_saccades(orientations, timestep):
    orientations = orientations[~numpy.isnan(orientations)]

    unwrappedOrientations = numpy.unwrap(orientations,180)
    SMOOTH_TIME = 1
    smoothWindow = int(round(SMOOTH_TIME/timestep))
    smoothedOrientations = tools.smooth(unwrappedOrientations,smoothWindow)
    angSpds = abs(numpy.diff(smoothedOrientations))/timestep
    MIN_PEAK_SPD=30
    sacInds = tools.local_minima(-angSpds) & (angSpds > MIN_PEAK_SPD)
    sacSpds = angSpds[sacInds]

    MIN_SAC_SPD=5
    inSac = angSpds > MIN_SAC_SPD
    startInds = pylab.find(numpy.diff(inSac.view(dtype=int8)) == 1)
    stopInds = pylab.find(numpy.diff(inSac.view(dtype=int8)) == -1)
    if stopInds[-1] < startInds[-1]:
        l = stopInds.tolist()
        l.append(len(inSac)-1)
        stopInds = numpy.array(l)
    if startInds[0] > stopInds[0]:
        l = startInds.tolist()
        l.insert(0,0)
        startInds = numpy.array(l)

    angDisp = numpy.zeros(len(angSpds))
    for i, startInd in enumerate(startInds):
        stopInd = stopInds[i]
        angDisp[startInd:stopInd] = smoothedOrientations[stopInd] - smoothedOrientations[startInd]

    sacAmps = angDisp[sacInds]
    return sacAmps

#rootDir = '/home/peter/data/'
#baseDirs = ['noFilter', 'circularPolarizer', 'grayFilter', 'blueFilter' ]
#baseDirs = ['noFilter','circularPolarizer','indoor/halogen','grayFilter','blueFilter','noFilter/cloudy','circularPolarizer/cloudy','grayFilter/cloudy']
#PLOTCOLORS = ['g','r','k','b','c','m','y','k']
#PLOTCOLORS = [[0,.9,0],[1,0,0],[1,.5,0],[.2,.2,.2],[0,0,1],[0,.5,0],[.6,0,0],[.5,.5,.5]]

rootDir = '/'
baseDirs = ['/media/weir05/data/indoor/dark','/home/peter/data/circularPolarizer','/home/peter/data/noFilter','/home/peter/data/blueFilter']
PLOTCOLORS = ['m','r','g','k','b','c','k','y']
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
    
MAX_NUM_STOPS = 6
SAVEFIGS = False
CHANGEBUFFER = 10
STARTBUFFER = 0
STOPBUFFER = 0
COLORS = dict(N='b',E='y',S='r',W='g')
ROTATIONS = dict(N=0,E=90,S=180,W=270)
SPEED = 0.5
DMAX = 330 #800
Dist = {}
sacAmp = {}
sacNum = {}
fig1 = pylab.figure()
#for b, baseDir in enumerate(flies.keys()):
for b, baseDir in enumerate(baseDirs):
    numFlies = len(flies[baseDir])
    D = numpy.empty(numFlies)
    D.fill(nan)
    S = numpy.empty(numFlies)
    S.fill(nan)
    Sn = numpy.empty(numFlies)
    Sn.fill(nan)
    ax = fig1.add_subplot(2,4,b+1)
    th = numpy.arange(0,2*numpy.pi,.01)
    R1 = 100 #400
    xcirc = R1*cos(th)
    ycirc = R1*sin(th)
    #plot(xcirc,ycirc,color='k')
    ax.fill(xcirc,ycirc,color=[.9,.9,.9])
    pylab.hold('on')
    for fNum in range(numFlies):
        dirName=os.path.join(baseDir,dName)
        fly = flies[baseDir][fNum].copy()
        sky = skies[baseDir][fNum].copy()
        
        print fly['dirName']

        orientations = fly['orientations'].copy()
        orientations = orientations + 180
        
        times = fly['times'].copy()
        beginTime = times[0]
        
        num_blocks = 1
        min_time = .9*num_blocks*12*60
        #max_time = num_blocks*12*60 + beginTime

        if len(sky['changeTimes']) >=5:
            max_time = sky['changeTimes'][4] - CHANGEBUFFER
        else:
            max_time = 0

        inds = times > max_time
        orientations[inds] = numpy.nan
        times[inds] = numpy.nan
        
        if fly.has_key('stopTimes'):
            #num_stops = sum(fly['stopTimes'] < max_time)
            num_stops = sum([sT < max_time for sT in fly['stopTimes']])
            for i, sT in enumerate(fly['stopTimes']):
                inds = (times > sT - STOPBUFFER) & (times < fly['startTimes'][i]+STOPBUFFER)
                orientations[inds] = numpy.nan
                #timeStopped[fNum] = timeStopped[fNum] + numpy.sum(numpy.diff(times[inds]))
                times[inds] = numpy.nan
        else:
            num_stops = 0
                
        print 'num stops = ', num_stops
        if numpy.nansum(numpy.diff(times))>min_time and num_stops <= MAX_NUM_STOPS:
            t=numpy.diff(times)
            timestep = numpy.median(t[~numpy.isnan(t)])
            #saccades = detect_saccades(orientations,timestep)
            #S[fNum] = beginTime - 60*60
            #S[fNum], Sn[fNum] = numpy.mean(abs(saccades)), len(saccades)
            S, Sn = numpy.array([]), numpy.array([])
            
            rotations = numpy.empty(times.shape)
            rotations.fill(numpy.nan)
            for i, cT in enumerate(sky['changeTimes'][:-1]):
                inds = (times > cT-CHANGEBUFFER+STARTBUFFER) & (times <= sky['changeTimes'][i+1]-CHANGEBUFFER)
                rotations[inds] = ROTATIONS[sky['directions'][i]]
                
            worldOrientations = orientations+rotations
            
            #aS = abs(numpy.diff(numpy.unwrap(orientations,180)))/t
            #S[fNum] = numpy.median(aS[~numpy.isnan(aS)])
            #aS = numpy.unwrap(orientations,180)
            #aS = aS[~numpy.isnan(aS)]
            #if len(aS) > 0:
            #    S[fNum] = aS[-1]
            #    print aS[-1]
            simVelX = numpy.sin(worldOrientations[1:]*numpy.pi/180)*SPEED*t
            simVelY = numpy.cos(worldOrientations[1:]*numpy.pi/180)*SPEED*t
            simTrajX = simVelX[~numpy.isnan(simVelX)].cumsum()
            simTrajY = simVelY[~numpy.isnan(simVelY)].cumsum()
            D[fNum]=numpy.sqrt(simTrajX[-1]**2 + simTrajY[-1]**2)
            print D[fNum]
            #T[fNum] = numpy.nansum(numpy.diff(times))
            l = ax.plot(simTrajX,simTrajY,color=PLOTCOLORS[b])
            
            changeInds = [abs(times[~numpy.isnan(times)] - (sky['changeTimes'][i]-CHANGEBUFFER)).argmin() for i in range(len(sky['changeTimes']))]
            #ax.plot(simTrajX[changeInds[1:4]],simTrajY[changeInds[1:4]],'.',color=l[0].get_color())
            #ax.scatter(simTrajX[-1],simTrajY[-1],edgecolors=l[0].get_color(),facecolors='k')
            ax.plot(simTrajX[-1],simTrajY[-1],'o',markeredgecolor=l[0].get_color(),markerfacecolor='k')

    ax.set_xlim((-DMAX,DMAX))
    ax.set_ylim((-DMAX,DMAX))
    ax.set_aspect('equal')
    ax.set_title(baseDir)
    ax.set_axis_off()
    fig1.set_facecolor('w')
    Dist[baseDir] = D.copy()
    sacAmp[baseDir] = S.copy()
    sacNum[baseDir] = Sn.copy()

fig2 = pylab.figure()
pylab.hold('on')
ax = fig2.add_subplot(111)
for b, baseDir in enumerate(baseDirs):
    plotArgs = dict(color=PLOTCOLORS[b])
    line = cumprobdist(ax,Dist[baseDir][~numpy.isnan(Dist[baseDir])],DMAX,plotArgs=plotArgs)
    #line = cumprobdist(angSpeed[baseDir][~numpy.isnan(angSpeed[baseDir])])
    ax.axvline(mean(Dist[baseDir][~numpy.isnan(Dist[baseDir])]),color=PLOTCOLORS[b])
    line.set_label(baseDir)
    
ax.axvline(R1,color='k')
ax.legend(loc='upper right')
ax.set_ylim((-.1,1.1))
ax.set_xlim((-25,DMAX))


fig3 = pylab.figure()
ax = fig3.add_subplot(111)
#pylab.boxplot([Dist[k][~numpy.isnan(Dist[k])] for k in Dist.keys()])
ax.boxplot([Dist[k][~numpy.isnan(Dist[k])] for k in baseDirs])
xtickText = [bd[-24:] + ', N=' + str(sum(~numpy.isnan(Dist[bd]))) for bd in baseDirs]
ax.set_xticklabels(xtickText)
pylab.draw()

fig4 = pylab.figure()
ax = fig4.add_subplot(111)
ax.hold('on')
for b, baseDir in enumerate(baseDirs):
    plotArgs = dict(color=PLOTCOLORS[b])
    N = sum(~numpy.isnan(Dist[baseDir]))
    ax.scatter(b*numpy.ones(N)+(pylab.rand(N)-.5)/8.0,Dist[baseDir][~numpy.isnan(Dist[baseDir])],color=PLOTCOLORS[b])
    ax.scatter(b,mean(Dist[baseDir][~numpy.isnan(Dist[baseDir])]),marker='+')
    

#meanDists = [mean(Dist[bd][~numpy.isnan(Dist[bd])]) for bd in baseDirs]
#ax.scatter(range(b+1),meanDists,marker='s')

ax.set_xlim((-.5,b+.5))
ax.set_xticks(range(b+1))
#ax.set_xticklabels(baseDirs)
ax.set_xticklabels(xtickText)

if SAVEFIGS == True:
    fig1.savefig('simulatedTrajectories.svg', dpi=600)
    fig2.savefig('cumulativeSum.svg', dpi=300)
    fig3.savefig('boxPlot.svg', dpi=300)
    fig4.savefig('scatter.svg', dpi=300)


