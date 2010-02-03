import numpy as np
import motmot.FlyMovieFormat.FlyMovieFormat as FMF
import colormapTools as cmt
import pylab
import sys
import time
import tables

def do_fft(pixels,time=None):

    data = pixels
    if time is None:
        t = arange(data.shape[-1])
    else:
        t = time
    d = np.median(np.diff(time))

    rsp = np.fft.rfft(data)
    fp = np.fft.fftfreq(data.shape[-1],d)
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
    i = np.median(ind)
    
    s = np.sum(pwr,axis=2)
    
    power = pwr[:,:,i]/s
    phase = phs[:,:,i]
    
    return power, phase

def do_polarimetry(filename,nFrames=100):

    fmf = FMF.FlyMovie(filename)

    if fmf.get_n_frames() < nFrames:
        nFrames = fmf.get_n_frames()
        print filename + " only has " + str(nFrames) + " frames"

    frame,timestamp = fmf.get_frame(0)

    N = 3
    Nx = N
    Ny = N
    LEFT = 100
    RIGHT = 600
    TOP = 0
    BOTTOM = frame.shape[-2]

    #X = np.round(np.linspace(0,frame.shape[-1],Nx+1))
    X = np.round(np.linspace(LEFT,RIGHT,Nx+1))
    Y = np.round(np.linspace(TOP,BOTTOM,Ny+1))

    power = np.empty(frame.shape)
    power.fill(np.nan)
    phase = np.empty(frame.shape)
    phase.fill(np.nan)
    intensity = np.empty(frame.shape)
    intensity.fill(np.nan)

    for i,x in enumerate(X[:-1]):
    #for i,x in enumerate(X[:2]):
        sys.stdout.write('i='+str(i)+'\n')
        sys.stdout.flush()
        for j,y in enumerate(Y[:-1]):
        #for j,y in enumerate(Y[:2]):
            sys.stdout.write('j='+str(j)+'\n')
            sys.stdout.flush()
            ROIFrames = np.empty([Y[j+1]-y,X[i+1]-x,nFrames])
            timestamps = np.empty(nFrames)
            for frameNumber in range(nFrames):
                frame,timestamps[frameNumber] = fmf.get_frame(frameNumber)
                ROIFrames[:,:,frameNumber] = frame[y:Y[j+1],x:X[i+1]]
                
            power[y:Y[j+1],x:X[i+1]], phase[y:Y[j+1],x:X[i+1]] = do_fft(ROIFrames,timestamps)
            intensity[y:Y[j+1],x:X[i+1]] = np.mean(ROIFrames,axis = 2)
            
    #polarimetry = np.rec.fromarrays([power,phase,intensity],names='power,phase,intensity')
    return power, phase, intensity

def analyze_directory(dirName):
    filenames = os.listdir(dirName)
    skyFilenames = [f for f in filenames if f[:3] == 'sky']

    spnum = -1
    fig1 = pylab.figure()
    fig2 = pylab.figure()
    fig3 = pylab.figure()
    for f, filename in enumerate(skyFilenames):
        h5Filename = filename[:-3]+'h5'
        h5FileExists = False
        if filename[-3:] == 'fmf':
            spnum = spnum+1
            for fname in filenames:
                if fname == h5Filename:
                    h5FileExists = True
                    
            if h5FileExists:
                "reading polarimetry from "+os.path.join(dirName,h5Filename)
                h5File = tables.openFile(os.path.join(dirName,h5Filename),'r')
                power = h5File.getNode("/power").read()
                phase = h5File.getNode("/phase").read()
                intensity = h5File.getNode("/intensity").read()
            else:
                print "doing polarimetry on "+os.path.join(dirName,filename)
                power, phase, intensity = do_polarimetry(os.path.join(dirName,filename))
                h5File = tables.openFile(os.path.join(dirName,h5Filename), mode = "w", title = h5Filename + "pol data")
                h5File.createArray(h5File.root,"power",power)
                h5File.createArray(h5File.root,"phase",phase)
                h5File.createArray(h5File.root,"intensity",intensity)
                h5File.close()
                
            phase = phase - np.mean(phase[222:261,338:382])
            phase = cmt.add_colordisc(phase)
            
            #trueUpDirection = filename[3]
            #if trueUpDirection == 'E':
            #    phase = phase+np.pi/2
            #elif trueUpDirection == 'S':
            #    phase = phase+np.pi
            #elif trueUpDirection == 'W':
            #    phase = phase+3*np.pi/2
            
            phase = np.mod(phase,2*np.pi)
            
            pylab.figure(fig1.number)
            ax = pylab.subplot(221+spnum)
            pylab.imshow(intensity,cmap='gray')
            pylab.title(filename)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            pylab.show()
            
            pylab.figure(fig2.number)
            ax = pylab.subplot(221+spnum)
            pylab.imshow(power,cmap='jet')
            pylab.colorbar()
            pylab.title(filename)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            pylab.show()
            
            gb180 = cmt.get_cmap('gb180')
            pylab.figure(fig3.number)
            ax = pylab.subplot(221+spnum)
            pylab.imshow(phase,cmap=gb180)
            pylab.title(filename)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            pylab.show()  
    return

dirName = '/home/cardini/12/'
analyze_directory(dirName)

#filename = '/home/cardini/12/skyN20100201_173114.fmf'

#power, phase, intensity = do_polarimetry(filename)
#phase = phase - np.mean(phase[222:261,338:382])
#phase = cmt.add_colordisc(phase)
#phase = np.mod(phase,2*np.pi)

#pylab.figure()
#pylab.imshow(intensity,cmap='gray')
#pylab.show()

#pylab.figure()
#pylab.imshow(power,cmap='jet')
#pylab.colorbar()
#pylab.show()

#gb180 = cmt.get_cmap('gb180')

#pylab.figure()
#pylab.imshow(phase,cmap=gb180)
#pylab.show()
