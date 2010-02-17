import numpy as np
import motmot.FlyMovieFormat.FlyMovieFormat as FMF
import colormapTools as cmt
import pylab
import sys
import time
import tables
from scipy.stats.morestats import circmean
from scipy.signal import sepfir2d, gaussian #, convolve2d

def show_angle(angle,power):
    ARROW_STEP = 40
    kernel = np.ones((ARROW_STEP,ARROW_STEP))
    rowFilter = gaussian(ARROW_STEP,ARROW_STEP/5)
    colFilter = rowFilter
    
    gb180 = cmt.get_cmap('gb180')
    
    #X = np.arange(0,angle.shape(-1),ARROW_STEP)
    #Y = np.arange(0,angle.shape(-2),ARROW_STEP)

    x = np.matrix(np.arange(ARROW_STEP/2,angle.shape[-1],ARROW_STEP))
    y = np.transpose(np.matrix(np.arange(ARROW_STEP/2,angle.shape[-2],ARROW_STEP)))

    X = np.array((0*y+1)*x)
    Y = np.array(y*(0*x+1))
    
    #u = convolve2d(sin(angle),kernel,mode='same')
    #v = convolve2d(cos(angle),kernel,mode='same')
    #p = convolve2d(power,kernel,mode='same')
    
    u = sepfir2d(sin(angle),rowFilter,colFilter)
    v = sepfir2d(cos(angle),rowFilter,colFilter)
    p = sepfir2d(power,rowFilter,colFilter)
    
    U = u[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP]*p[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP]
    V = v[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP]*p[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP]
    
    #U = sin(angle[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP])*(power[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP])
    #V = cos(angle[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP])*(power[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP])

    #X = X[(power[::ARROW_STEP,::ARROW_STEP]>.016)]
    #Y = Y[(power[::ARROW_STEP,::ARROW_STEP]>.016)]
    #U = U[(power[::ARROW_STEP,::ARROW_STEP]>.016)]
    #V = V[(power[::ARROW_STEP,::ARROW_STEP]>.016)]
    
    ua = ARROW_STEP/1.5*np.nansum(sin(angle)*power)/np.nansum(power)
    va = -ARROW_STEP/1.5*np.nansum(cos(angle)*power)/np.nansum(power)
    xc, yc = angle.shape[-1]/2, angle.shape[-2]/2
    
    pylab.imshow(angle,cmap=gb180)
    #pylab.imshow(angle,cmap='hsv')
    pylab.quiver(X,Y,U,V,pivot='middle',color='w',headwidth=1,headlength=0)
    pylab.arrow(xc,yc,ua,va,color='w',linewidth=2)
    pylab.arrow(xc,yc,-ua,-va,color='w',linewidth=2)
    pylab.show()

    A = np.arctan2(-ua,va)
    return A

def do_fft(pixels,time=None):

    data = pixels
    if time is None:
        t = arange(data.shape[-1])
    else:
        t = time
    d = np.median(np.diff(time))

    rsp = np.fft.rfft(data)
    fp = np.fft.fftfreq(data.shape[-1],d)
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
    print freq[i]
    
    s = np.sum(pwr,axis=2)
    
    power = pwr[:,:,i]/s
    phase = phs[:,:,i]
    
    return power, phase

def do_polarimetry(filename,nFrames=500):
    
    FRAMES_TO_SKIP = 10 #number of rames at the beginning of the movie to skip
    PLOTPIX = True
    if PLOTPIX is True:
        fig = pylab.figure()
        
    fmf = FMF.FlyMovie(filename)

    if fmf.get_n_frames() < nFrames:
        nFrames = fmf.get_n_frames()
        print filename + " only has " + str(nFrames) + " frames"

    frame,timestamp = fmf.get_frame(0)

    N = 3
    Nx = N
    Ny = N
    LEFT = 120
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
                frame,timestamps[frameNumber] = fmf.get_frame(frameNumber+FRAMES_TO_SKIP) # skip first FRAMES_TO_SKIP frames
                ROIFrames[:,:,frameNumber] = frame[y:Y[j+1],x:X[i+1]]
                
            power[y:Y[j+1],x:X[i+1]], phase[y:Y[j+1],x:X[i+1]] = do_fft(ROIFrames,timestamps)
            intensity[y:Y[j+1],x:X[i+1]] = np.mean(ROIFrames,axis = 2)
            if PLOTPIX is True:
                pylab.figure(fig.number)
                pylab.plot(timestamps,ROIFrames[0,0,:], label='p=' + str(power[y,x])[:4] + ' a=' + str(phase[y,x]*180/np.pi)[:4])
                pylab.legend()
                pylab.show()
                pylab.draw()
            
    power = power[Y[0]:Y[-1],X[0]:X[-1]] # not checked
    phase = phase[Y[0]:Y[-1],X[0]:X[-1]]
    intensity = intensity[Y[0]:Y[-1],X[0]:X[-1]]
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
                print "reading polarimetry from "+os.path.join(dirName,h5Filename)
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
                            
            phase = phase - circmean(np.ravel(phase[222:261,222:272]),high=np.pi,low=-np.pi)
            
            angle = phase/2.0 # because phase offset of intensity values is twice angle between overlapping polarizers
            
            trueUpDirection = filename[3]
            if trueUpDirection == 'E':
                power = np.rot90(power,1)
                angle = np.rot90(angle,1)
                intensity = np.rot90(intensity,1)
                angle = angle + np.pi/2
            elif trueUpDirection == 'S':
                power = np.rot90(power,2)
                angle = np.rot90(angle,2)
                intensity = np.rot90(intensity,2)
                angle = angle + np.pi
            elif trueUpDirection == 'W':
                power = np.rot90(power,3)
                angle = np.rot90(angle,3)
                intensity = np.rot90(intensity,3)
                angle = angle + 3*np.pi/2

            mask = intensity>(np.mean(intensity)-.45*np.std(intensity)) #hack
            #mask[100:300,100:300] = True not sure if central dot (polarizer) should be in or not... if so - threshold should be ~1 std below mean intensity
            power[~mask] = nan
            angle[~mask] = nan
            
            angle = cmt.add_colordisc(angle,width=71)
            angle = np.mod(angle+np.pi,2*np.pi)-np.pi
            
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
            
            pylab.figure(fig3.number)
            ax = pylab.subplot(221+spnum)
            A = show_angle(angle,power)
            pylab.title(filename)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            pylab.show()  
    return power, angle, intensity

dirName = '/home/cardini/15/'
pwr, ang, ints = analyze_directory(dirName)

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
