MAX_NUM_STOPS = 4
MAX_TIME_STOPPED = 60
chunksPerTrial = 30
maxTimeStoppedChunk = 3
#chunksPerTrial = 5
#maxTimeStoppedChunk = 6

PVAL = .05
baseDirs = ['indoor/dark','circularPolarizer','noFilter','blueFilter','grayFilter']
    
#import statTestChunks
#reload(statTestChunks)
fig=pylab.figure()
ax = fig.add_subplot(111)
#testStatistics, p_values = statTestChunks.plotIt(ax,flies,skies,baseDirs,MAX_NUM_STOPS,MAX_TIME_STOPPED,CHANGEBUFFER,chunksPerTrial,maxTimeStoppedChunk)

import rayleighStatTestChunks
reload(rayleighStatTestChunks)
numChunks = 24
maxTimeStoppedChunk = 30
testStatistics, p_values = rayleighStatTestChunks.plotIt(ax,flies,skies,baseDirs,MAX_NUM_STOPS,MAX_TIME_STOPPED,CHANGEBUFFER,numChunks,maxTimeStoppedChunk)
for baseDirInd, baseDir in enumerate(baseDirs):
    print baseDir, np.sum(p_values[baseDirInd]<PVAL), " sig shift (p<", PVAL, ") out of N = ", len(testStatistics[baseDirInd])

