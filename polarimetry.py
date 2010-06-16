import numpy as np
import motmot.FlyMovieFormat.FlyMovieFormat as FMF
import colormapTools as cmt
import pylab
import sys, os, time
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
    
    u = sepfir2d(np.sin(angle),rowFilter,colFilter)
    v = sepfir2d(np.cos(angle),rowFilter,colFilter)
    p = sepfir2d(power,rowFilter,colFilter)
    
    U = u[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP]*p[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP]
    V = v[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP]*p[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP]
    
    #U = sin(angle[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP])*(power[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP])
    #V = cos(angle[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP])*(power[ARROW_STEP/2::ARROW_STEP,ARROW_STEP/2::ARROW_STEP])

    #X = X[(power[::ARROW_STEP,::ARROW_STEP]>.016)]
    #Y = Y[(power[::ARROW_STEP,::ARROW_STEP]>.016)]
    #U = U[(power[::ARROW_STEP,::ARROW_STEP]>.016)]
    #V = V[(power[::ARROW_STEP,::ARROW_STEP]>.016)]
    
    ua = ARROW_STEP/1.5*np.nansum(np.sin(angle)*power)/np.nansum(power)
    va = -ARROW_STEP/1.5*np.nansum(np.cos(angle)*power)/np.nansum(power)
    xc, yc = angle.shape[-1]/2, angle.shape[-2]/2
    
    pylab.imshow(angle,cmap=gb180)
    #pylab.imshow(angle,cmap='hsv')
    pylab.quiver(X,Y,U,V,pivot='middle',color='w',headwidth=1,headlength=0)
    pylab.arrow(xc,yc,ua,va,color='w',linewidth=2)
    pylab.arrow(xc,yc,-ua,-va,color='w',linewidth=2)
    ax=pylab.gca()
    ax.set_axis_off()
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
    #i=9
    print i, freq[i]
    
    s = np.sum(pwr,axis=2)
    
    power = pwr[:,:,i]/s
    phase = phs[:,:,i]
    
    return power, phase

def do_polarimetry(fmf,firstFrame=0,nFrames=500):
    
    PLOTPIX = False
    if PLOTPIX is True:
        fig = pylab.figure()
        fig.hold('on')
        
    if fmf.get_n_frames() < nFrames+firstFrame:
        nFrames = fmf.get_n_frames()-firstFrame
        print "fmf only has " + str(fmf.get_n_frames()) + " frames"

    frame,timestamp = fmf.get_frame(firstFrame)

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
        for j,y in enumerate(Y[:-1]):
            ROIFrames = np.empty([Y[j+1]-y,X[i+1]-x,nFrames])
            timestamps = np.empty(nFrames)
            for frameNumber in range(nFrames):
                frame,timestamps[frameNumber] = fmf.get_frame(frameNumber+firstFrame) # start at firstFrame
                ROIFrames[:,:,frameNumber] = frame[y:Y[j+1],x:X[i+1]]
                
            power[y:Y[j+1],x:X[i+1]], phase[y:Y[j+1],x:X[i+1]] = do_fft(ROIFrames,timestamps)
            intensity[y:Y[j+1],x:X[i+1]] = np.mean(ROIFrames,axis = 2)
            if PLOTPIX is True:
                pylab.figure(fig.number)
                fig.hold('on')
                pylab.plot(timestamps,ROIFrames[0,0,:], label='p=' + str(power[y,x])[:4] + ' a=' + str(phase[y,x]*180/np.pi)[:4])
                pylab.legend()
                pylab.show()
                pylab.draw()
            
    power = power[Y[0]:Y[-1],X[0]:X[-1]] # not checked
    phase = phase[Y[0]:Y[-1],X[0]:X[-1]]
    intensity = intensity[Y[0]:Y[-1],X[0]:X[-1]]
    return power, phase, intensity

def compare_file_times(fn1, fn2):
    t1, t2 = int(fn1[4:12]+fn1[13:19]),int(fn2[4:12]+fn2[13:19])
    return cmp(t1,t2)
    
def analyze_file(sky,fname=None):
    """
    example:
    power,angle,intensity=polarimetry.analyze_file(sky)
    """
    WAIT_TIME = 20 #seconds after changeTimes to start polarimetry
    ROT180 = True #because camera returns rotated image
    
    if fname is None:
        fileName = sky['fileName']
    else:
        fileName = fname
    dirName = sky['dirName']
    
    N = len(sky['changeTimes'][:-1])
    
    fmf = FMF.FlyMovie(os.path.join(dirName,fileName))
    frame,timestamp = fmf.get_frame(0)
    
    if not sky.has_key('times'):
        timestamps = fmf.get_all_timestamps()
    else:
        timestamps = sky['times']
    
    for i, startTime in enumerate(sky['changeTimes'][:-1]):
        startInd = np.argmin(abs(startTime + WAIT_TIME - timestamps))
        sys.stdout.write(time.ctime(startTime + WAIT_TIME)+'\n')
        #sys.stdout.write("%s\n" % (str(i)))
        sys.stdout.flush()
        pwr, phs, ints = do_polarimetry(fmf,firstFrame=startInd,nFrames=500)
                        
        phs = phs - circmean(np.ravel(phs[222:261,222:272]),high=np.pi,low=-np.pi)
        
        ang = phs/2.0 # because phase offset of intensity values is twice angle between overlapping polarizers
        
        if ROT180:
            pwr = np.rot90(pwr,2)
            ang = np.rot90(ang,2)
            ints = np.rot90(ints,2)
        """
        trueUpDirection = filename[3]
        if trueUpDirection == 'E':
            pwr = np.rot90(pwr,1)
            ang = np.rot90(ang,1)
            ints = np.rot90(ints,1)
            ang = ang + np.pi/2
        elif trueUpDirection == 'S':
            pwr = np.rot90(pwr,2)
            ang = np.rot90(ang,2)
            ints = np.rot90(ints,2)
            ang = ang + np.pi
        elif trueUpDirection == 'W':
            pwr = np.rot90(pwr,3)
            ang = np.rot90(ang,3)
            ints = np.rot90(ints,3)
            ang = ang + 3*np.pi/2
        """
        
        mask = ints>(np.mean(ints)-.45*np.std(ints)) #hack
        mask[100:300,100:300] = True #not sure if central dot (polarizer) should be in or not... if so - threshold should be ~1 std below mean intensity
        pwr[~mask] = np.nan
        ang[~mask] = np.nan
        
        ang = np.mod(ang+np.pi/2,2*np.pi)-np.pi/2
        ang = cmt.add_colordisc(ang,width=71)
        #ang = np.mod(ang+np.pi,2*np.pi)-np.pi
        ang = np.mod(ang,2*np.pi)
        
        if i==0:
            w,h=ang.shape
            power = np.empty([w,h,N])
            power.fill(np.nan)
            angle = np.empty([w,h,N])
            angle.fill(np.nan)
            intensity = np.empty([w,h,N])
            intensity.fill(np.nan)
            
        power[:,:,i],angle[:,:,i],intensity[:,:,i]=pwr,ang,ints
    return power, angle, intensity


