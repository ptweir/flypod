pylab.ion()

FPS = 291
SPEED = 0.5

COLORS = dict(N='b',E='y',S='r',W='g')
ROTATIONS = dict(N=0,E=90,S=180,W=270)

orientations = copy(fly['orientations'])
times = copy(fly['times'])
rotations = numpy.empty(fly['times'].shape)
rotations.fill(numpy.nan)
for i, cT in enumerate(sky['changeTimes'][:-1]):
    inds = (times > cT) & (times <= sky['changeTimes'][i+1])
    rotations[inds] = ROTATIONS[sky['directions'][i]]
    
worldOrientations = orientations+rotations

simVelX = numpy.sin(worldOrientations*numpy.pi/180)*SPEED/FPS
simVelY = numpy.cos(worldOrientations*numpy.pi/180)*SPEED/FPS
simTrajX = simVelX.cumsum()
simTrajY = simVelY.cumsum()

pylab.figure()
pylab.scatter(simTrajX[0],simTrajY[0])
hold('on')
for i, cT in enumerate(sky['changeTimes'][:-1]):
    inds = (times > cT) & (times <= sky['changeTimes'][i+1])
    pylab.plot(simTrajX[inds],simTrajY[inds],color=COLORS[sky['directions'][i]])
#pylab.plot(simTrajX,simTrajY)
#for i, d in enumerate(COLORS):
    #inds = rotations==ROTATIONS[d]
    #pylab.plot(simTrajX[inds],simTrajY[inds],color=COLORS[d])
ax=gca()
ax.set_aspect('equal')
title(fly['dirName'])
pylab.draw()
