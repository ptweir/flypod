from __future__ import division
import motmot.FlyMovieFormat.FlyMovieFormat as FMF
from pygarrayimage.arrayimage import ArrayInterfaceImage
import motmot.imops.imops as imops
from pyglet import window
import sys
import pylab
import numpy as np
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

def get_pixels(filenames,showFrames=0):
    """get pixel value data from fmf movies filenames
    
    arguments:      
    filename        list of .fmf movie files or single .fmf movie file
    showFrames      0 [default] or 1
    
    example:
    pol, time = get_pixels('/home/cardini/2/movie20091202_165430.fmf',0)
    """ 
    if isinstance(filenames,str):
        filenames = [filenames]
        
    top = 250
    bottom = 100
    left = 350
    right = 300
    
    ZOOM = 10
    FRAMESTEP = 1
    dnSmpRate = 5.0
    
    #pix = []
    oldFrameNumber = 0
    for f,filename in enumerate(filenames):
        if filename[-3:] == 'fmf':
        
            fmf = FMF.FlyMovie(filename)

            nFrames = fmf.get_n_frames()
            
            #pix = []
            frame,timestamp = fmf.get_frame(0)
            ROIFrame = frame[top:-bottom,left:-right]
            p = frame[::dnSmpRate,:dnSmpRate*np.floor(frame.shape[-1]/dnSmpRate):dnSmpRate]
            
            if showFrames:
                DnSmpFrame = p
                #dispFrame = convert(ROIFrame,fmf.format)
                dispFrame = convert(DnSmpFrame,fmf.format)
                
                dispFrame = np.repeat(dispFrame,ZOOM,axis=0)
                dispFrame = np.repeat(dispFrame,ZOOM,axis=1)
                wnd = window.Window(visible=False, resizable=True)
                aii = ArrayInterfaceImage(dispFrame)
                img = aii.texture
                wnd.width = img.width
                wnd.height = img.height
                wnd.set_caption(filename)
                wnd.set_visible()
                
            pix = np.empty([p.shape[-2],p.shape[-1],nFrames])
            pix.fill(np.nan)
            time = np.empty(nFrames)
            
            for frameNumber in range(nFrames):
                frame,timestamp = fmf.get_frame(frameNumber)
                
                ROIFrame = frame[top:-bottom,left:-right]
                
                DnSmpFrame = np.empty([p.shape[-2],p.shape[-1],dnSmpRate])
                for d in range(dnSmpRate):
                    DnSmpFrame[:,:,d] = frame[d::dnSmpRate,d:dnSmpRate*np.floor(frame.shape[-1]/dnSmpRate):dnSmpRate]
                    
                p = np.mean(DnSmpFrame,axis=2)
                
                if np.mod(frameNumber,FRAMESTEP) == 0:
                    sys.stdout.write('\b'*(len(str(oldFrameNumber))+4+len(str(nFrames)))+str(frameNumber)+' of '+ str(nFrames))
                    oldFrameNumber = frameNumber
                    sys.stdout.flush()
                    if showFrames:
                        if wnd.has_exit:
                            break

                        wnd.dispatch_events()
                        #dispFrame = convert(ROIFrame,fmf.format)
                        dispFrame = p.astype(np.uint8)
                        #dispFrame = convert(DnSmpFrame,fmf.format)

                        dispFrame = np.repeat(dispFrame,ZOOM,axis=0)
                        dispFrame = np.repeat(dispFrame,ZOOM,axis=1)
                        aii.view_new_array(dispFrame)
                        img.blit(0, 0, 0)
                        wnd.flip()

                #pix.append((p,timestamp))
                pix[:,:,frameNumber] = p
                time[frameNumber] = timestamp
                
        if showFrames:
            wnd.close()
    pixels = pix

    return pixels, time

def do_polarimetry(pixels,time=None):

    data = pixels
    if time is None:
        t = arange(data.shape[-1])
    else:
        t = time
    d = np.median(np.diff(time))

    N=None
    sp = np.fft.fft(data,N)
    rsp = np.fft.rfft(data,N)
    fp = np.fft.fftfreq(sp.shape[-1],d)
    #fp = np.fft.fftfreq(data.shape[-1],d)
    fpp = fp[fp>=0]
    freq = np.empty(rsp.shape[-1])
    freq[:fpp.shape[-1]] = fpp
    if fp.shape[-1] != rsp.shape[-1]:
        freq[-1] = -np.min(fp)


    amp = np.abs(rsp)
    pwr = amp**2
    phs = np.angle(rsp) # something is wrong here
    #phs = np.arctan2(rsp.real,rsp.imag)

    ind = np.argmax(pwr*(freq!=0),axis=2)
    #m = pwr[ind]

    pw = np.empty(ind.shape)
    ph = np.empty(ind.shape)

    for r in range(pwr.shape[-3]):
        for c in range(pwr.shape[-2]):
            #plot(freq, pwr[r,c,:])#,[fp[ind], fp[ind]],[0,np.sqrt(m[r,c])])
            pw[r,c] = pwr[r,c,ind[r,c]]
            ph[r,c] = phs[r,c,ind[r,c]]

    s = np.sum(pwr,axis=2)
    i = np.median(ind)
    pylab.figure()
    pylab.imshow(pwr[:,:,i]/s)
    #pylab.imshow(p/s)
    pylab.show()

    pylab.figure()
    #pylab.imshow(ph)
    pylab.imshow(phs[:,:,i])
    pylab.show()

