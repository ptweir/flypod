"""module to analyze magnetic tether outdoor experiments

PTW
"""

from __future__ import division
import motmot.FlyMovieFormat.FlyMovieFormat as FMF
from pygarrayimage.arrayimage import ArrayInterfaceImage
import motmot.imops.imops as imops
from pyglet import window
import sys
import matplotlib
#matplotlib.use('tkagg')
import pylab
import numpy
import os

pylab.ion()

def convert(frame,format):
    """ covert frame format """ 
    if format in ['RGB8','ARGB8','YUV411','YUV422']:
        frame = imops.to_rgb8(format,frame)
    elif format in ['MONO8','MONO16']:
        frame = imops.to_mono8(format,frame)
    elif (format.startswith('MONO8:') or
          format.startswith('MONO32f:')):
        # bayer
        frame = imops.to_rgb8(format,frame)
    return frame

def get_centers(filenames,showFrames=0,ROI=None):
    """find center of fly in fmf files, returns record array with x,y and time
    
    arguments:      
    filename        list of .fmf movie files or single .fmf movie file
    showFrames      0 [default] or 1
    ROI             4-tuple of region of interest (top, bottom, left, right) or None (whole frame)
    
    example:
    centers = get_centers('/home/cardini/data/fly07/flyN20100317_185940.fmf',0,(9,15,12,12))
    """ 
    THRESH = 1.5
    if isinstance(filenames,str):
        filenames = [filenames]
        
    if ROI is not None:
        top, bottom, left, right = ROI
    else:
        top, bottom, left, right = (1,0,0,1)
    
    COLOR = True
    ZOOM = 10
    FRAMESTEP = 100
    
    cents = []
    lenOutStr = 0
    for f,filename in enumerate(filenames):
        if filename[-3:] == 'fmf':
        
            fmf = FMF.FlyMovie(filename)

            nFrames = fmf.get_n_frames()
            
            frame,timestamp = fmf.get_frame(0)
            ROIFrame = frame[bottom:-top,left:-right]
            
            ## mask to eliminate points in corner -> 'ring' mask
            height,width = ROIFrame.shape
            R = numpy.floor(width/2)
            x0,y0 = R,R
            if numpy.mod(x0,2)==0:
                x0,y0 = x0-.5,y0-.5
            
            x = numpy.matrix(range(width))
            y = numpy.transpose(numpy.matrix(range(height)))

            xx = numpy.array((0*y+1)*x - y0)
            yy = numpy.array(y*(0*x+1) - x0)

            radius = numpy.sqrt(xx**2 + yy**2)
            ring = (radius<R) & (radius>R/4)
            ## end mask stuff
            
            if showFrames:
                
                dispFrame = convert(ROIFrame,fmf.format)
                if COLOR == True:
                    dispFrame = numpy.array([dispFrame,dispFrame,dispFrame])
                    dispFrame = numpy.swapaxes(dispFrame,0,2)
                    dispFrame = numpy.swapaxes(dispFrame,0,1)
                
                dispFrame = numpy.repeat(dispFrame,ZOOM,axis=0)
                dispFrame = numpy.repeat(dispFrame,ZOOM,axis=1)
                wnd = window.Window(visible=False, resizable=True)
                aii = ArrayInterfaceImage(dispFrame)
                img = aii.texture
                wnd.width = img.width
                wnd.height = img.height
                wnd.set_caption(filename[-23:])
                wnd.set_visible()

            for frameNumber in range(nFrames):
                frame,timestamp = fmf.get_frame(frameNumber)
                
                ROIFrame = frame[bottom:-top,left:-right]
                #invertedFrame = 255 - ROIFrame
                #invertedFrame = invertedFrame - 120

                threshFrame = (ROIFrame < numpy.mean(ROIFrame) - THRESH*numpy.std(ROIFrame)) & (ring) #mask stuff

                h,w = threshFrame.shape
                X = numpy.arange(w)
                Y = numpy.arange(h)
                cX = numpy.sum(numpy.sum(threshFrame,0)*X)/numpy.sum(threshFrame)
                cY = numpy.sum(numpy.sum(threshFrame,1)*Y)/numpy.sum(threshFrame)
                
                if numpy.mod(frameNumber,FRAMESTEP) == 0:
                    outStr = str(frameNumber)+' of '+ str(nFrames)
                    sys.stdout.write("%s%s\r" % (outStr, " "*lenOutStr ))
                    lenOutStr = len(outStr)
                    sys.stdout.flush()
                    if showFrames:
                        if wnd.has_exit:
                            break

                        wnd.dispatch_events()
                        dispFrame = convert(ROIFrame,fmf.format)
                        dispFrame[~ring]=dispFrame[~ring]/4
                        #dispFrame = threshFrame.astype(numpy.uint8)*155 +100
                        if ~numpy.isnan(cY) and ~numpy.isnan(cY):
                            dispFrame[round(cY),round(cX)] = 0
                        if COLOR == True:
                            dispFrame = numpy.array([dispFrame,dispFrame,dispFrame])
                            dispFrame = numpy.swapaxes(dispFrame,0,2)
                            dispFrame = numpy.swapaxes(dispFrame,0,1)
                            #dispFrame[ring,2]=255
                            if ~numpy.isnan(cY) and ~numpy.isnan(cY):
                                dispFrame[round(cY),round(cX),0] = 255

                        dispFrame = numpy.repeat(dispFrame,ZOOM,axis=0)
                        dispFrame = numpy.repeat(dispFrame,ZOOM,axis=1)
                        aii.view_new_array(dispFrame)
                        img.blit(0, 0, 0)
                        wnd.flip()

                cents.append((cX,cY,timestamp))
        if showFrames:
            wnd.close()

    centers = numpy.rec.fromarrays(numpy.transpose(cents), [('x',numpy.float),('y',numpy.float),('t',numpy.float)])
    if showFrames:
        pylab.figure()
        pylab.scatter(centers.x,centers.y,s=1)
        pylab.draw()

    return centers

def circle_fit(dataX,dataY):
    """fit a circle to data, returns x,y,radius
    
    arguments:      
    dataX       numpy array containing x data
    dataY       numpy array containing y data (must be same size as dataX)
    
    example:
    cx, cy, r = flypod.circle_fit(x,y)
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

def rose(data,wrapPoint=360,ax='current'):
    """makes polar histogram plot, returns wrapped data, n, bins, binCenters, axes
    
    arguments:      
    data            numpy array containing data
    wrap point      scalar value at which data is wrapped [default 360]
    ax              axes to plot in [default current]
    
    example:
    orw,n,b,bc,ax = rose(orientations)
    """ 
    NUMBINS = 10
    if ax == 'current':
        ax = pylab.gca()
        
    if not isinstance(ax,matplotlib.projections.polar.PolarAxes):
        ax = pylab.subplot(111,polar=True)

    wrappedData = numpy.mod(data,wrapPoint)*2*numpy.pi/wrapPoint
    n, bins, patches = pylab.hist(wrappedData,NUMBINS)
    binCenters = bins[:-1] + (bins[1:] - bins[:-1])/2
    n = numpy.append(n,n[0])
    binCenters = numpy.append(binCenters,binCenters[0])
    binCenters = -binCenters + numpy.pi/2 #this makes 0 straight up, with angles increasing clockwise
    pylab.cla()
    pylab.plot(binCenters,n)
    return wrappedData, n, bins, binCenters, ax

def analyze_directory(dirName):
    """runs analyzes .fmf files in dirName and creates plots
    
    arguments:      
    dirName     directory path to analyze
    
    example:
    analyze_directory('/home/cardini/data/fly07/')
    """ 
    filenames = os.listdir(dirName)
    flyFilenames = [f for f in filenames if f[:3] == 'fly']
    paramFileExists = False
    for f, filename in enumerate(filenames):
        if filename == 'params.txt':
            fd = open(os.path.join(dirName,filename),mode='r')
            line = fd.readline()
            sline = line.split()
            cx,cy,r = [float(i) for i in sline]
            line = fd.readline()
            sline = line.split()
            ROI = [int(i) for i in sline]
            fd.close()
            paramFileExists = True
    if not paramFileExists:
        ROI = input('please enter 4-tuple of (top,bottom,left,right): ')
        c = get_centers([os.path.join(dirName,f) for f in flyFilenames],0,ROI)
        cx, cy, r = circle_fit(c.x[~numpy.isnan(c.x)],c.y[~numpy.isnan(c.y)])
        fd = open(os.path.join(dirName,'params.txt'),mode='w')
        fd.write('%f %f %f\n'%(cx,cy,r))
        if ROI is not None:
            fd.write('%d %d %d %d'%ROI)
        else:
            fd.write('1 0 0 1')
        fd.flush()
        fd.close()
        

    spnum = -1
    fig1 = pylab.figure()
    fig2 = pylab.figure()
    fig3 = pylab.figure()
    for f, filename in enumerate(flyFilenames):
        csvFilename = filename[:-3]+'csv'
        csvFileExists = False
        if filename[-3:] == 'fmf':
            spnum = spnum+1
            if spnum == 4:
                spnum=0
                fig1 = pylab.figure()
                fig2 = pylab.figure()
                fig3 = pylab.figure()
            for fname in filenames:
                if fname == csvFilename:
                    csvFileExists = True
                    
            if csvFileExists:
                centers = pylab.csv2rec(os.path.join(dirName,csvFilename))
            else:
                centers = get_centers(os.path.join(dirName,filename),1,ROI)
                pylab.rec2csv(centers,os.path.join(dirName,csvFilename))
                
            #cx, cy, r = circle_fit(centers.x,centers.y)
            
            orientations = numpy.arctan2(centers.x[~numpy.isnan(centers.x)]-cx,centers.y[~numpy.isnan(centers.y)]-cy)*180/numpy.pi
            #these orientations are measured from 12 O'clock, increasing clockwise
            
            trueUpDirection = filename[3]
            if trueUpDirection == 'E':
                orientations = orientations+90
            elif trueUpDirection == 'S':
                orientations = orientations+180
            elif trueUpDirection == 'W':
                orientations = orientations+270
                
            orientations = numpy.mod(orientations+180,360)-180
            
            pylab.figure(fig1.number)
            pylab.subplot(221+spnum)
            pylab.scatter(centers.x,centers.y,s=1)
            pylab.hold(True)
            th = range(0,700)
            th = numpy.array(th)/100.0
            pylab.plot(r*numpy.cos(th)+cx,r*numpy.sin(th)+cy)
            pylab.plot(cx,cy,'+')
            pylab.hold(False)
            pylab.title(filename)
            pylab.draw()
            
            pylab.figure(fig2.number)
            pylab.subplot(221+spnum)
            pylab.hist(orientations)
            pylab.title(filename)
            pylab.draw()
            
            pylab.figure(fig3.number)
            pylab.subplot(221+spnum,polar=True)
            orw,n,b,bc,ax = rose(orientations,360)
            ax.set_rgrids([1],'')
            ax.set_thetagrids([0,90,180,270],['E','N','W','S'])
            pylab.title(filename)
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
    check_orientations('/home/cardini/12/flyS20100201_171506.fmf',orientations,cx,cy)
    """ 
    fmf = FMF.FlyMovie(fileName)
 
    ax = pylab.axes()
    pylab.gray()

    top = 20
    bottom = 5
    left = 15
    right =  10

    nFrames = fmf.get_n_frames()
    if nFrames != len(orientations):
        raise Exception('number of frames does not equal length of orientations')

    for frameNumber in range(0,nFrames,frameStep):
        frame,timestamp = fmf.get_frame(frameNumber)
        
        ROIFrame = frame[bottom:-top,left:-right]
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

#fileName = '/home/cardini/16/flyS20100215_174820.fmf'
#fileName = '/home/cardini/16/flyN20100215_173606.fmf'
#fileName = '/home/cardini/16/flyE20100215_174202.fmf'
#centers = get_centers(fileName,0)
#cx, cy, r = circle_fit(centers.x,centers.y)
#orientations = numpy.arctan2(centers.x[~numpy.isnan(centers.x)]-cx,centers.y[~numpy.isnan(centers.y)]-cy)*180/numpy.pi
#check_orientations(fileName,orientations,cx,cy)
