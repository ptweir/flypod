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

def polarimetry(filenames,showFrames=0):
    """get polarimetry data from fmf movies filenames
    
    arguments:      
    filename        list of .fmf movie files or single .fmf movie file
    showFrames      0 [default] or 1
    
    example:
    pol, time = polarimetry('/home/cardini/2/movie20091202_165430.fmf',0)
    """ 
    if isinstance(filenames,str):
        filenames = [filenames]
        
    top = 250
    bottom = 100
    left = 350
    right = 300
    
    ZOOM = 1
    FRAMESTEP = 1
    
    #pix = []
    oldFrameNumber = 0
    for f,filename in enumerate(filenames):
        if filename[-3:] == 'fmf':
        
            fmf = FMF.FlyMovie(filename)

            nFrames = fmf.get_n_frames()
            
            #pix = []
            frame,timestamp = fmf.get_frame(0)
            ROIFrame = frame[top:-bottom,left:-right]
            p = frame[::10,:10*numpy.floor(frame.shape[-1]/10.0):10]
            
            if showFrames:
                DnSmpFrame = p
                #dispFrame = convert(ROIFrame,fmf.format)
                dispFrame = convert(DnSmpFrame,fmf.format)
                
                dispFrame = numpy.repeat(dispFrame,ZOOM,axis=0)
                dispFrame = numpy.repeat(dispFrame,ZOOM,axis=1)
                wnd = window.Window(visible=False, resizable=True)
                aii = ArrayInterfaceImage(dispFrame)
                img = aii.texture
                wnd.width = img.width
                wnd.height = img.height
                wnd.set_caption(filename)
                wnd.set_visible()
                
            pix = numpy.empty([p.shape[-2],p.shape[-1],nFrames])
            pix.fill(numpy.nan)
            time = numpy.empty(nFrames)
            
            for frameNumber in range(nFrames):
                frame,timestamp = fmf.get_frame(frameNumber)
                
                ROIFrame = frame[top:-bottom,left:-right]
                
                DnSmpFrame = numpy.empty([p.shape[-2],p.shape[-1],10])
                for d in range(10):
                    DnSmpFrame[:,:,d] = frame[d::10,d:10*numpy.floor(frame.shape[-1]/10.0):10]
                    
                p = numpy.mean(DnSmpFrame,axis=2)
                
                if numpy.mod(frameNumber,FRAMESTEP) == 0:
                    sys.stdout.write('\b'*(len(str(oldFrameNumber))+4+len(str(nFrames)))+str(frameNumber)+' of '+ str(nFrames))
                    oldFrameNumber = frameNumber
                    sys.stdout.flush()
                    if showFrames:
                        if wnd.has_exit:
                            break

                        wnd.dispatch_events()
                        #dispFrame = convert(ROIFrame,fmf.format)
                        dispFrame = p.astype(numpy.uint8)
                        #dispFrame = convert(DnSmpFrame,fmf.format)

                        dispFrame = numpy.repeat(dispFrame,ZOOM,axis=0)
                        dispFrame = numpy.repeat(dispFrame,ZOOM,axis=1)
                        aii.view_new_array(dispFrame)
                        img.blit(0, 0, 0)
                        wnd.flip()

                #pix.append((p,timestamp))
                pix[:,:,frameNumber] = p
                time[frameNumber] = timestamp
                
        if showFrames:
            wnd.close()
    pixel = pix

    return pixel, time

