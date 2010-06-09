import pylab, time
pylab.ion()

def do_it(stopStartTimes,fly,sky=None):
    """
    example:
    stopTimes,startTimes=annotateStops.do_it(sky['stopStartTimes'],fly)
    fly['stopTimes']=stopTimes
    fly['startTimes']=startTimes
    
    outPklFile = open(os.path.join(dirName,fly['pklFileName']), 'wb')
    pickle.dump(fly, outPklFile)
    outPklFile.close()

    """
    if sky is not None:
        normM = sky['pixelMeans'] - min(sky['pixelMeans'])
        normM = normM/max(normM)
        normS = sky['pixelStds'] - min(sky['pixelStds'])
        normS = normS/max(normS)
    normO = fly['orientations'] - min(fly['orientations'])
    normO = normO/max(normO)
    pylab.ion()
    pylab.figure()
    if sky is not None:
        pylab.plot(sky['times'],normM)
        pylab.plot(sky['times'],normS)
    pylab.plot(fly['times'],normO)
    pylab.hold('on')
    stopTimes=[]
    startTimes=[]

    for stopTime in stopStartTimes:
        pylab.scatter(stopTime,0)
        pylab.xlim((stopTime-50,stopTime+50))
        pylab.draw()
        pt = pylab.ginput(1,0)
        stopTimes.append(pt[0][0])
        pylab.xlim((stopTime-50,stopTime+50))
        pt = pylab.ginput(1,0)
        startTimes.append(pt[0][0])
        
    return stopTimes, startTimes
    
