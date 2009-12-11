"""
Emulate an oscilloscope.  Requires the animation API introduced in
matplotlib 0.84.  See
http://www.scipy.org/wikis/topical_software/Animations for an
explanation.

This example uses gtk but does not depend on it intimately.  It just
uses the idle handler to trigger events.  You can plug this into a
different GUI that supports animation (GTKAgg, TkAgg, WXAgg) and use
your toolkits idle/timer functions.
"""
import gobject, gtk
import matplotlib
matplotlib.use('GTKAgg')
import numpy as np
from matplotlib.lines import Line2D
import time


class Scope:
	def __init__(self, ax, fd, maxt=100):
		self.ax = ax
		self.fd = fd
		self.canvas = ax.figure.canvas
		self.maxt = maxt
		self.Tstart = time.time()
		self.Tdata = [0]
		self.Pdata = [0]
		self.line = Line2D(self.Tdata, self.Pdata, animated=True)
		self.ax.add_line(self.line)
		self.background = None
		self.canvas.mpl_connect('draw_event', self.update_background)
		self.ax.set_ylim(0, 8000)
		self.ax.set_xlim(0, self.maxt)

	def update_background(self, event):
		self.background = self.canvas.copy_from_bbox(self.ax.bbox)

	def emitter(self, fd, Tstart):
		'return a random value with probability p, else 0'
		T=[]
		P=[]

		while 1:   # infinite loop
			line = fd.readline()
			if line:
				sline = line.split()
				T.append(float(sline[0])-Tstart)
				P.append(np.mod(int(sline[2]),8000))
				
			if not line:  # 'readline()' returns None at end of file.
				time.sleep(1)
				break
				
		return T, P

	def update(self, *args):
		if self.background is None: return True
		[T, P] = self.emitter(self.fd,self.Tstart)
		lastT = self.Tdata[-1]
		if lastT>self.Tdata[0]+self.maxt: # reset the arrays
			self.Tdata = [self.Tdata[-1]]
			self.Pdata = [self.Pdata[-1]]
			self.ax.set_xlim(self.Tdata[0], self.Tdata[0]+self.maxt)
			self.ax.figure.canvas.draw()

		self.canvas.restore_region(self.background)

		self.Tdata = self.Tdata + T
		self.Pdata = self.Pdata + P
		self.line.set_data(self.Tdata, self.Pdata)
		self.line.set_color('r')
		self.ax.draw_artist(self.line)

		self.canvas.blit(self.ax.bbox)
		return True


from pylab import figure, show
import sys
import os

fig = figure()

ax = fig.add_subplot(111)

ax.axesPatch.set_facecolor('k')
fig.figurePatch.set_facecolor('k')
ax.yaxis.set_ticks([0,1000,2000,3000,4000,5000,6000,7000,8000])
ax.yaxis.set_ticklabels([0, 45, 90, 135, 180, 225, 270, 315, 360],color='y')
ax.yaxis.grid(color='y', linestyle='-', linewidth=1)

fileNames = []
fileTimes = []
for fileName in os.listdir('./'):
    if fileName[12:16] == '2009':
    	fileNames.append(fileName)
        fileTimes.append(os.path.getmtime('./'+fileName))

sortedIndices = np.argsort(fileTimes)
newestFile = fileNames[sortedIndices[-1]]

fd = open(newestFile,'r')
scope = Scope(ax,fd)

gobject.idle_add(scope.update)

show()
fd.close()
