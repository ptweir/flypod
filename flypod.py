"""module to analyze magnetic tether outdoor experiments

PTW
"""

from __future__ import division
import pkg_resources
import motmot.FlyMovieFormat.FlyMovieFormat as FMF
from pygarrayimage.arrayimage import ArrayInterfaceImage
import motmot.imops.imops as imops
from pyglet import window
#matplotlib.use('tkagg')
import pylab, numpy, os, pickle, sys, matplotlib, time

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

def get_background(filename,FRAMESTEP=1000,ROI=None):
    """get a median background image
    
    arguments:      
    filename        .fmf movie file name
    FRAMESTEP       number of frames to skip in between ones used to calculate median (default 1000)
    ROI             4-tuple of region of interest (top, bottom, left, right) or None (whole frame)
    
    example:
    background = get_background('/home/cardini/data/fly07/flyN20100317_185940.fmf',(9,15,12,12))
    """ 
    if ROI is not None:
        top, bottom, left, right = ROI
    else:
        top, bottom, left, right = (1,0,0,1)
        
    if filename[-3:] == 'fmf':
        fmf = FMF.FlyMovie(filename)

        nFrames = fmf.get_n_frames()
        frameInds = range(0,nFrames,FRAMESTEP)
        #frameInds = range(1,round(nFrames*.55),FRAMESTEP)
        #frameInds = range(round(nFrames*.5),nFrames,FRAMESTEP)
        #frameInds = range(181000,257000,FRAMESTEP)
        
        frame,timestamp = fmf.get_frame(0)
        ROIFrame = frame[bottom:-top,left:-right]
        
        bg = numpy.empty((ROIFrame.shape[0],ROIFrame.shape[1],len(frameInds)))
        
        lenOutStr = 0
        for f,frameNumber in enumerate(frameInds):
            frame,timestamp = fmf.get_frame(frameNumber)
            ROIFrame = frame[bottom:-top,left:-right]
            
            if numpy.mod(frameNumber,FRAMESTEP) == 0:
                outStr = str(frameNumber)+' of '+ str(nFrames)
                sys.stdout.write("%s%s\r" % (outStr, " "*lenOutStr ))
                lenOutStr = len(outStr)
                sys.stdout.flush()

            bg[:,:,f] = ROIFrame

        background = numpy.median(bg,axis=2)

    return background

def get_centers(filenames,showFrames=0,ROI=None,THRESH=1.5,ringR=.25, background=None):
    """find center of fly in fmf files, returns record array with x,y and time
    
    arguments:      
    filename        list of .fmf movie files or single .fmf movie file
    showFrames      0 [default] or 1
    ROI             4-tuple of region of interest (top, bottom, left, right) or None (whole frame)
    THRESH          threshold for number of std for pixel to be counted as foreground (default 1.5)
    ringR           inner radius for mask (default .25)
    background      background image (default None)
    
    example:
    centers = get_centers('/home/cardini/data/fly07/flyN20100317_185940.fmf',0,(9,15,12,12))
    """ 

    if isinstance(filenames,str):
        filenames = [filenames]
        
    if ROI is not None:
        top, bottom, left, right = ROI
    else:
        top, bottom, left, right = (1,0,0,1)
    
    COLOR = True
    ZOOM = 10
    FRAMESTEP = 100
    
    ended_early = False
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
            ring = (radius<R) & (radius>R*ringR)
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
                if background is not None:
                    ROIFrame = ROIFrame - background

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
                            ended_early=True
                            break

                        wnd.dispatch_events()
                        #dispFrame = convert(ROIFrame,fmf.format)
                        #dispFrame = ROIFrame.astype(numpy.uint8)
                        #dispFrame[~ring]=dispFrame[~ring]/4
                        dispFrame = threshFrame.astype(numpy.uint8)*155 +100
                        if ~numpy.isnan(cY) and ~numpy.isnan(cY):
                            dispFrame[round(cY),round(cX)] = 0
                        if COLOR == True:
                            dispFrame = numpy.array([dispFrame,dispFrame,dispFrame])
                            dispFrame = numpy.swapaxes(dispFrame,0,2)
                            dispFrame = numpy.swapaxes(dispFrame,0,1)
                            dispFrame[ring,2]=255
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
        #pylab.scatter(centers.x,centers.y,s=1)
        pylab.plot(centers.x,centers.y)
        pylab.draw()
        
    if ended_early:
        centers = None
        
    return centers

def circle_fit(dataX,dataY):
    """fit a circle to data, returns x,y,radius
    
    arguments:      
    dataX       numpy array containing x data
    dataY       numpy array containing y data (must be same size as dataX)
    
    example:
    cx, cy, r = flypod2.circle_fit(x,y)
    """ 
    n = sum(~numpy.isnan(dataX))
    a = numpy.ones((n,3))
    a[:,0] = dataX[~numpy.isnan(dataX)]
    a[:,1] = dataY[~numpy.isnan(dataY)]
    b = -dataX[~numpy.isnan(dataX)]**2 - dataY[~numpy.isnan(dataY)]**2
    #out = numpy.linalg.solve(a,b)
    ai = numpy.linalg.pinv(a)
    out = numpy.dot(ai,b)
    circCenterX = -.5*out[0]
    circCenterY = -.5*out[1]
    circR  =  ((out[0]**2+out[1]**2)/4-out[2])**.5;
    return circCenterX, circCenterY, circR

def rose(data,wrapPoint=360,ax='current',plotArgs={}):
    """makes polar histogram plot, returns wrapped data, n, bins, binCenters, axes
    
    arguments:      
    data            numpy array containing data
    wrap point      scalar value at which data is wrapped [default 360]
    ax              axes to plot in [default current]
    plotArgs        dict of plot arguments (color, etc.)
    
    example:
    orw,n,b,bc,ax = rose(orientations)
    """ 
    NUMBINS = 36
    if ax == 'current':
        ax = pylab.gca()
        
    if not isinstance(ax,matplotlib.projections.polar.PolarAxes):
        ax = pylab.subplot(1,1,1,polar=True)

    wrappedData = numpy.mod(data,wrapPoint)*2*numpy.pi/wrapPoint
    n, bins, patches = pylab.hist(wrappedData,NUMBINS,range=(0,2*numpy.pi),visible=False)
    n = n/sum(n)
    binCenters = bins[:-1] + (bins[1:] - bins[:-1])/2
    n = numpy.append(n,n[0])
    binCenters = numpy.append(binCenters,binCenters[0])
    binCenters = -binCenters + numpy.pi/2 #this makes 0 straight up, with angles increasing clockwise
    pylab.plot(binCenters,n,**plotArgs)
    ax.set_rmax(.15)
    return wrappedData, n, bins, binCenters, ax

def compare_file_times(fn1, fn2):
    t1, t2 = int(fn1[4:12]+fn1[13:19]),int(fn2[4:12]+fn2[13:19])
    return cmp(t1,t2)

def analyze_directory(dirName):
    """runs analyzes .fmf files in dirName and creates plots
    
    arguments:      
    dirName     directory path to analyze
    
    example:
    analyze_directory('/home/cardini/data/fly07/')
    """ 
    filenames = os.listdir(dirName)
    flyFilenames = [f for f in filenames if f[:3] == 'fly' and f[-3:] == 'fmf']
    flyFilenames.sort(compare_file_times)
    pklFilenames = [f for f in filenames if f[:4] == 'aFly' and f[-3:] == 'pkl']
    pklFilenames.sort(compare_file_times)
    if len(pklFilenames) > 0:
        pklFilename = pklFilenames[-1]
    else:
        date_time = time.strftime("%Y%m%d_%H%M%S")
        pklFilename = 'aFly'+ date_time +'.pkl'

    for f, filename in enumerate(filenames):
        if filename == pklFilename:
            inPklFile = open(os.path.join(dirName,filename), 'rb')
            fly = pickle.load(inPklFile)
            inPklFile.close()
            break
    else:
        fly = {}
        ROI = input('please enter 4-tuple of integers (top,bottom,left,right): ')
        THRESH = input('please enter scalar THRESH: ')
        ringR = input('please enter scalar ringR: ')
        useBackground = input('use background subtraction (0 or 1)?: ')
        
        fullFilenames = [os.path.join(dirName,filename) for filename in flyFilenames]
        background = None
        if useBackground:
            background = get_background(fullFilenames[-1],FRAMESTEP=1000,ROI=ROI)
        centers = get_centers(fullFilenames,1,ROI,THRESH,ringR,background)
        if centers is None:
            return centers
        cx, cy, r = circle_fit(centers.x,centers.y)

        orientations = numpy.arctan2(centers.x-cx,centers.y-cy)*180/numpy.pi
        
        #these orientations are measured from 12 O'clock, increasing clockwise
            
        fly['dirName'], fly['fileName'] = dirName, flyFilenames[0]
        fly['pklFileName'] = pklFilename
        fly['ROI'], fly['THRESH'], fly['ringR'] = ROI, THRESH, ringR
        fly['x'] = centers.x
        fly['y'] = centers.y
        fly['times'] = centers.t
        fly['orientations'] = orientations
        fly['background'] = background
        
        pylab.figure()
        pylab.plot((fly['times']-fly['times'][0])/60,fly['orientations'])
        pylab.title('total pts missed tracking = ' + str(sum(numpy.isnan(fly['x']))))

    outPklFile = open(os.path.join(dirName,pklFilename), 'wb')
    pickle.dump(fly, outPklFile)
    outPklFile.close()
    return fly

def check_orientations(fileName,orientations,cx,cy,ROI=None,frameStep=100):
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

    fig = pylab.figure()
    ax = pylab.axes()
    pylab.gray()

    if ROI is not None:
        top, bottom, left, right = ROI
    else:
        top, bottom, left, right = (1,0,0,1)

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

