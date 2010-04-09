"""sky_times

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

import tools

pylab.ion()

skies = {}
class Sky:
    def __init__(self, dirName, filename):
        self.dirName = dirName
        self.fileName = filename

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

def get_pixel_stats(filenames,showFrames=0,ROI=None):
    """get pixel stats
    
    arguments:      
    filename        list of .fmf movie files or single .fmf movie file
    showFrames      0 [default] or 1
    ROI             4-tuple of region of interest (top, bottom, left, right) or None (whole frame)
    """ 

    if isinstance(filenames,str):
        filenames = [filenames]
        
    if ROI is not None:
        top, bottom, left, right = ROI
    else:
        top, bottom, left, right = (1,0,120,600)
        #top, bottom, left, right = (1,0,0,1)
    
    FRAMESTEP = 1000
    
    pixStats = []
    lenOutStr = 0
    for f,filename in enumerate(filenames):
        if filename[-3:] == 'fmf':
        
            fmf = FMF.FlyMovie(filename)

            nFrames = fmf.get_n_frames()
            
            frame,timestamp = fmf.get_frame(0)
            ROIFrame = frame[bottom:,left:right]
            
            ## mask to eliminate points in corner -> 'ring' mask
            height,width = ROIFrame.shape
            R = numpy.floor(height/2)
            x0,y0 = R,R
            if numpy.mod(x0,2)==0:
                x0,y0 = x0-.5,y0-.5
            
            x = numpy.matrix(range(width))
            y = numpy.transpose(numpy.matrix(range(height)))

            xx = numpy.array((0*y+1)*x - y0)
            yy = numpy.array(y*(0*x+1) - x0)

            radius = numpy.sqrt(xx**2 + yy**2)
            ring = (radius<R)& (radius>R/5)
            ## end mask stuff
            
            if showFrames:
                
                dispFrame = convert(ROIFrame,fmf.format)

                wnd = window.Window(visible=False, resizable=True)
                aii = ArrayInterfaceImage(dispFrame)
                img = aii.texture
                wnd.width = img.width
                wnd.height = img.height
                wnd.set_caption(filename[-23:])
                wnd.set_visible()

            for frameNumber in range(nFrames):
                frame,timestamp = fmf.get_frame(frameNumber)
                
                ROIFrame = frame[bottom:,left:right]
                #invertedFrame = 255 - ROIFrame
                #invertedFrame = invertedFrame - 120

                pixM = numpy.mean(ROIFrame[ring])
                pixS = numpy.std(ROIFrame[ring])
                
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
                        dispFrame[~ring]=255

                        aii.view_new_array(dispFrame)
                        img.blit(0, 0, 0)
                        wnd.flip()

                pixStats.append((pixM,pixS,timestamp))
        if showFrames:
            wnd.close()

    pixelStats = numpy.rec.fromarrays(numpy.transpose(pixStats), [('m',numpy.float),('s',numpy.float),('t',numpy.float)])
    if showFrames:
        pylab.figure()
        pylab.plot(pixelStats.t,pixelStats.m)
        pylab.hold('on')
        pylab.plot(pixelStats.t,pixelStats.s)
        pylab.draw()

    return pixelStats


def compare_file_times(fn1, fn2):
    t1, t2 = int(fn1[4:12]+fn1[13:19]),int(fn2[4:12]+fn2[13:19])
    return cmp(t1,t2)

def analyze_directory(dirName):
    """pixel mean analysis
    
    arguments:      
    dirName     directory path to analyze
    
    example:
    analyze_directory('/home/cardini/data/fly07/')
    """ 
    filenames = os.listdir(dirName)
    skyFilenames = [f for f in filenames if f[:3] == 'sky' and f[-3:] == 'fmf']    
    skyFilenames.sort(compare_file_times)
    for f, filename in enumerate(skyFilenames):
        csvFilename = filename[:-3]+'csv'
        csvFileExists = False
        if filename[-3:] == 'fmf':
            for fname in filenames:
                if fname == csvFilename:
                    csvFileExists = True
                    
            if csvFileExists:
                pixelStats = pylab.csv2rec(os.path.join(dirName,csvFilename))
            else:
                pixelStats = get_pixel_stats(os.path.join(dirName,filename),1)
                pylab.rec2csv(pixelStats,os.path.join(dirName,csvFilename))
            
            trueUpDirection = filename[3]
            
            pylab.figure()
            pylab.plot(pixelStats.t,tools.smooth(pixelStats.m,200))
            pylab.draw()
            pylab.title('click two points for threshold')
            pylab.hold('on')
            pts = pylab.ginput(2)
            x = [pt[0] for pt in pts]
            y = [pt[1] for pt in pts]
            slope = (y[0] - y[1])/(x[0] - x[1])
            threshold = slope * (pixelStats.t - x[0]) + y[0]
            covered = tools.smooth(pixelStats.m,200) < threshold
            pylab.plot(pixelStats.t,threshold)
            pylab.plot(pixelStats.t,covered)
            
            skies[f] = Sky(dirName,filename)
            skies[f].pixelMeans = pixelStats.m
            skies[f].pixelStds = pixelStats.s
            skies[f].times = pixelStats.t
            skies[f].covered = covered
    
    return skies
