"""module to analyze magnetic tether outdoor experiments

PTW 12/10/2009
"""

from __future__ import division
import motmot.FlyMovieFormat.FlyMovieFormat as FMF
import sys
import matplotlib
matplotlib.use('tkagg')
import pylab
import numpy
import os

pylab.ion()

def get_centers(fileName,showFrames=0):
    """find center of fly in fmf, returns record array with x,y, and time
    
    arguments:      
    filename        .fmf movie
    showFrames      0 [default] or 1
    
    example:
    centers = get_centers('/home/cardini/2/movie20091202_165430.fmf',0)
    """ 
    fmf = FMF.FlyMovie(fileName)

    if showFrames:    
        ax = pylab.axes()

    top = 25
    bottom = 25
    left = 25
    right =  25

    nFrames = fmf.get_n_frames()

    centers = numpy.recarray((nFrames,), [('x',numpy.float),('y',numpy.float),('t',numpy.float)])

    for frameNumber in range(nFrames):
    #for frameNumber in range(0,10000,300):
        frame,timestamp = fmf.get_frame(frameNumber)
        
        ROIFrame = frame[top:-bottom,left:-right]
        #invertedFrame = 255 - ROIFrame
        #invertedFrame = invertedFrame - 120

        threshFrame = ROIFrame < numpy.mean(ROIFrame) - 2*numpy.std(ROIFrame)

        w,h = threshFrame.shape
        X = numpy.arange(w)
        Y = numpy.arange(h)
        cX = numpy.sum(numpy.sum(threshFrame,0)*X)/numpy.sum(threshFrame)
        cY = numpy.sum(numpy.sum(threshFrame,1)*Y)/numpy.sum(threshFrame)
        
        if showFrames and numpy.mod(frameNumber,100) == 0:
            #pylab.imshow(threshFrame)
            pylab.imshow(ROIFrame)
            pylab.hold(True)
            pylab.scatter(cX,cY)
            pylab.title('frame '+str(frameNumber)+' of '+str(nFrames))
            pylab.draw()
            pylab.hold(False)

        centers[frameNumber] = (cX,cY,timestamp)

    if showFrames:
        pylab.figure()
        pylab.scatter(centers['x'],centers['y'])
        pylab.draw()

    return centers

def circle_fit(dataX,dataY):
    """fit a circle to data, returns x,y,radius
    
    arguments:      
    dataX       numpy array containing x data
    dataY       numpy array containing y data (must be same size as dataX)
    
    example:
    cx, cy, r = analyze.circle_fit(x,y)
    """ 
    n = len(dataX)
    a = numpy.ones((n,3))
    a[:,0] = dataX
    a[:,1] = dataY
    b = -dataX**2 - dataY**2
    #out = numpy.linalg.solve(a,b)
    ai = numpy.linalg.pinv(a)
    out = numpy.dot(ai,b)
    circCenterX = -.5*out[0]
    circCenterY = -.5*out[1]
    circR  =  ((out[0]**2+out[1]**2)/4-out[2])**.5;
    return circCenterX, circCenterY, circR

def rose(data,wrapPoint=360,includeData=1,ax='current'):
    """makes polar histogram plot, returns wrapped data, n, bins, binCenters
    
    arguments:      
    data            numpy array containing data
    wrap point      scalar value at which data is wrapped [default 360]
    include data    TO ADD
    
    example:
    orw,n,b,bc = rose(orientations)
    """ 
    if ax == 'current':
        ax = pylab.gca()
        
    if not isinstance(ax,matplotlib.projections.polar.PolarAxes):
        ax = pylab.subplot(111,polar=True)

    wrappedData = numpy.mod(data,wrapPoint)*2*numpy.pi/wrapPoint
    n, bins, patches = pylab.hist(wrappedData)
    binCenters = bins[:-1] + (bins[1:] - bins[:-1])/2
    n = numpy.append(n,n[0])
    binCenters = numpy.append(binCenters,binCenters[0])
    pylab.cla()
    pylab.plot(binCenters,n)
    return wrappedData, n, bins, binCenters

def analyze_directory(dirName):
    """runs analyzes .fmf files in dirName and creates plots
    
    arguments:      
    dirName     directory path to analyze
    
    example:
    analyze_directory('/home/cardini/2/')
    """ 
    filenames = os.listdir(dirName)
    fig1 = pylab.figure()
    fig2 = pylab.figure()
    fig3 = pylab.figure()
    for f in range(len(filenames)):
        print f
        filename = dirName+filenames[f]
        centers = get_centers(filename,0)
        cx, cy, r = circle_fit(centers['x'],centers['y'])
        orientations = numpy.arctan2(centers['x']-cx,centers['y']-cy)*180/numpy.pi
        
        pylab.figure(fig1.number)
        pylab.subplot(221+f)
        pylab.plot(centers['x'],centers['y'])
        pylab.hold(True)
        th = range(0,700)
        th = numpy.array(th)/100.0
        pylab.plot(r*numpy.cos(th)+cx,r*numpy.sin(th)+cy)
        pylab.hold(False)
        pylab.draw()
        
        pylab.figure(fig2.number)
        pylab.subplot(221+f)
        pylab.hist(orientations)
        pylab.draw()
        
        pylab.figure(fig3.number)
        pylab.subplot(221+f,polar=True)
        orw,n,b,bc = rose(orientations,360,1)
        pylab.draw()  
    return


def check_orientations(fileName,orientations,cx,cy,frameStep=100):
    """plots input orientation from center (cx,cy) on images from .fmf movie file
    
    arguments:      
    filename        .fmf movie
    orientations    orientation in degrees clockwise from vertical for each frame
    cx              horizontal center of frame
    cy              vertical center of frame
    frame step      only plot 1 frame per frameStep [default 100]
    
    example:
    check_orientations('/home/cardini/2/movie20091202_165430.fmf',orientations,cx,cy)
    """ 
    fmf = FMF.FlyMovie(fileName)
 
    ax = pylab.axes()

    top = 25
    bottom = 25
    left = 25
    right =  25

    nFrames = fmf.get_n_frames()
    if nFrames != len(orientations):
        raise Exception('number of frames does not equal length of orientations')

    for frameNumber in range(0,nFrames,frameStep):
        frame,timestamp = fmf.get_frame(frameNumber)
        
        ROIFrame = frame[top:-bottom,left:-right]
        threshFrame = ROIFrame < numpy.mean(ROIFrame) - 2*numpy.std(ROIFrame)
        
        orx = (cx, numpy.sin(orientations[frameNumber]*numpy.pi/180)*cx + cx)
        ory = (cy, numpy.cos(orientations[frameNumber]*numpy.pi/180)*cy + cy)
    
        #pylab.imshow(threshFrame)
        pylab.imshow(ROIFrame)
        pylab.hold(True)
        pylab.plot(orx,ory)
        pylab.title('frame '+str(frameNumber)+' of '+str(nFrames))
        pylab.draw()
        pylab.hold(False)




