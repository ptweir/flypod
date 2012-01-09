import numpy as np
import pylab
pylab.ion()
from scipy.stats import ttest_ind
import matplotlib as mpl
mpl.rcParams['font.size'] = 8
mpl.rcParams['svg.embed_char_paths'] = False
"""
import pickle
inPklFile = open('allRotationFlies.pkl', 'rb')
flies = pickle.load(inPklFile)
inPklFile.close()
inPklFile = open('allRotationSkies.pkl', 'rb')
skies = pickle.load(inPklFile)
inPklFile.close()
"""

NUM_SUBPLOT_ROWS = 4
NUM_SUBPLOT_COLS = 6
CHANGEBUFFER = 9.5
MAX_NUM_STOPS = 4#8
MAX_TIME_STOPPED = 60

FIG_WIDTH_mm = 174 #188
fig0 = pylab.figure(figsize=(FIG_WIDTH_mm/25.4,5.75))
#fig0 = pylab.figure(figsize=(18,10))
fig0.set_facecolor('w')
ax0 = fig0.add_subplot(NUM_SUBPLOT_ROWS,NUM_SUBPLOT_COLS,1)
WIDTH_RIG_PICTURE = .35
ax0.set_position([-.05,.63,WIDTH_RIG_PICTURE+.1,.32])
arenaImg = pylab.imread('rotationArena.png')
ax0.imshow(arenaImg)
ax0.set_aspect('equal')
ax0.set_axis_off()

ax1 = fig0.add_subplot(NUM_SUBPLOT_ROWS,NUM_SUBPLOT_COLS,2)
ax1.set_position([.1+WIDTH_RIG_PICTURE,.81,.88-WIDTH_RIG_PICTURE,.16])
flyToPlot = flies['grayFilter'][2].copy()
skyToPlot = skies['grayFilter'][2].copy()
import plotIndividual
reload(plotIndividual)
plotIndividual.plotIt(ax1,flyToPlot,skyToPlot,CHANGEBUFFER)

ax2 = fig0.add_subplot(NUM_SUBPLOT_ROWS,NUM_SUBPLOT_COLS,3)
ax2.set_position([.1+WIDTH_RIG_PICTURE,.63,.88-WIDTH_RIG_PICTURE,.16])
import plotIndividualCompass
reload(plotIndividualCompass)
plotIndividualCompass.plotIt(ax2,flyToPlot,skyToPlot,CHANGEBUFFER)

import plotOrientationsAfterRotations
reload(plotOrientationsAfterRotations)
orientationChanges = plotOrientationsAfterRotations.plotIt(flies, skies, fig0, MAX_NUM_STOPS, MAX_TIME_STOPPED, CHANGEBUFFER)

axAll = fig0.add_subplot(NUM_SUBPLOT_ROWS,NUM_SUBPLOT_COLS,24)
axAll.set_position([.81,.08,.12,.21]) # boxplots
axIndv = fig0.add_subplot(NUM_SUBPLOT_ROWS,NUM_SUBPLOT_COLS,19)
axIndv.set_position([-.03,.0,.95,.28]) # trajectories
import newSimulatedTrajectories
reload(newSimulatedTrajectories)
allDistances = newSimulatedTrajectories.plotIt(flies,skies,axIndv,axAll, MAX_NUM_STOPS, MAX_TIME_STOPPED,CHANGEBUFFER)
for dist in allDistances:
    print "N = ", len(dist)
t, p_val = ttest_ind(allDistances[0],allDistances[1])
print p_val
t, p_val = ttest_ind(allDistances[2],allDistances[1])
print p_val 
t, p_val = ttest_ind(allDistances[3],allDistances[1])
print p_val
t, p_val = ttest_ind(allDistances[4],allDistances[1])
print p_val

#fig0.savefig('Fig4.svg', dpi=600)

