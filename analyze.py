import time, numpy
import flypod, sky_times
from scipy.stats.morestats import circmean, circvar	
pylab.ion()

#dirName = '/home/cardini/data/indoor/fly16/'
#dirName = '/home/cardini/data/circularPolarizer/fly02'
#dirName = '/home/cardini/data/grayFilter/fly20'
#dirName = '/home/cardini/data/noFilter/fly38'
dirName = '/home/cardini/data/uvFilter/fly03'
fly = flypod.analyze_directory(dirName)
CHANGEBUFFER = 10
if fly is not None:
    sky = sky_times.analyze_directory(dirName)

    COLORS = dict(N='b',E='y',S='r',W='g')
    ROTATIONS = dict(N=0,E=90,S=180,W=270)

    TOTAL_TIME = 12*3*60

    worldTotals=numpy.array([])
    totals = dict(N=numpy.array([]),E=numpy.array([]),S=numpy.array([]),W=numpy.array([]))


    orientations = numpy.copy(fly['orientations']) + 180
    times = numpy.copy(fly['times'])
    #orientations = numpy.unwrap(orientations,180)
    pylab.figure()
    pylab.plot(times,orientations,'k')
    ax = gca()
    
    ax.set_xticklabels([time.ctime(float(ti))[11:19] for ti in ax.get_xticks()])
    ax.set_yticks((0,90,180,270,360))
    for i, cT in enumerate(sky['changeTimes'][:-1]):
        pylab.text(cT,380,sky['directions'][i])
        pylab.axvspan(cT-CHANGEBUFFER, sky['changeTimes'][i+1]-CHANGEBUFFER, facecolor=COLORS[sky['directions'][i]], alpha=0.2)
        #pylab.draw()
    if sum(numpy.isnan(orientations)) > 0:
        for trackingErrorTime in times[numpy.isnan(orientations)]:
            pylab.axvline(trackingErrorTime,linewidth=2,color='k')
    if fly.has_key('stopTimes'):
        for i, sT in enumerate(fly['stopTimes']):
            pylab.axvspan(sT, fly['startTimes'][i], facecolor='r', alpha=0.5)
            #pylab.draw()
            inds = (times > sT) & (times < fly['startTimes'][i])
            pylab.plot(times[inds],orientations[inds],color='0.5')
            orientations[inds] = numpy.nan
    inds = (times > TOTAL_TIME + times[0])
    orientations[inds] = numpy.nan

    pylab.draw()
    fig = pylab.figure()
    fig.set_facecolor('w')
    pylab.suptitle(fly['fileName'])
    for i, cT in enumerate(sky['changeTimes'][:-1]):
        pylab.subplot(5,4,1+i,polar=True)
        ax=pylab.gca()
        inds = (times > cT) & (times < sky['changeTimes'][i+1])
        ors = orientations[inds]
        ors = ors[~numpy.isnan(ors)]
        if len(ors)>0:
            orw,n,b,bc,ax = flypod.rose(ors,360)
            m=circmean(ors,high=180,low=-180)
            v=circvar(ors*numpy.pi/180,high=numpy.pi,low=-numpy.pi)
            hold('on')
            polar([0,numpy.pi/2-m*numpy.pi/180],[0,ax.get_rmax()*(1-v)])
            ax.set_rmax(.4)
            ax.set_rgrids([1],'')
            ax.set_thetagrids([0,90,180,270],['','','',''])
            title(sky['directions'][i])
            #ax.set_axis_bgcolor(COLORS[sky['directions'][i]])
            ax.axesPatch.set_facecolor(COLORS[sky['directions'][i]])
            ax.axesPatch.set_alpha(0.4)
            
            totals[sky['directions'][i]] = numpy.concatenate((totals[sky['directions'][i]],ors))
        
    pylab.draw()
        
    fig = pylab.figure()
    fig.set_facecolor('w')
    pylab.suptitle(fly['dirName'])
    for i, d in enumerate(COLORS):
        #m=circmean(totals[d],high=180,low=-180)
        m=circmean(totals[d]+ROTATIONS[d],high=180,low=-180)
        v=circvar(totals[d]*numpy.pi/180,high=numpy.pi,low=-numpy.pi)
        pylab.subplot(2,2,i+1,polar=True)
        orw,n,b,bc,ax = flypod.rose(totals[d]+ROTATIONS[d],360)
        #orw,n,b,bc,ax = flypod.rose(totals[d],360)
        pylab.hold('on')
        pylab.polar([0,numpy.pi/2-m*numpy.pi/180],[0,ax.get_rmax()*(1-v)])
        ax.set_rmax(.4)
        ax.set_rgrids([1],'')
        ax.set_thetagrids([0,90,180,270],['','','',''])
        pylab.title(d)
        ax.axesPatch.set_facecolor(COLORS[d])
        ax.axesPatch.set_alpha(0.2)
        worldTotals = numpy.concatenate((worldTotals,totals[d]+ROTATIONS[d]))

    M=circmean(worldTotals,high=180,low=-180)
    V=circvar(worldTotals*numpy.pi/180,high=numpy.pi,low=-numpy.pi)
    print M,V
