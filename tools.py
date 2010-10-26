import numpy
import motmot.FlyMovieFormat.FlyMovieFormat as FMF
import pylab

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

def get_circle(fileName):
    fmf = FMF.FlyMovie(fileName)
    
    frame,timestamp = fmf.get_frame(0)
    pylab.imshow(frame)
    pylab.title('click on circle 8 times')
    pts = pylab.ginput(8,100)
    x = [pt[0] for pt in pts]
    y = [pt[1] for pt in pts]
    cx, cy, r = circle_fit(numpy.asarray(x),numpy.asarray(y))
    
    th = range(0,700)
    th = numpy.array(th)/100.0
    pylab.plot(r*numpy.cos(th)+cx,r*numpy.sin(th)+cy)
    pylab.plot(cx,cy,'+')
    pylab.draw()
    return cx,cy,r

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """
    
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
        
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
        
    if window_len<3:
        return x
        
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        
    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')
        
    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]
    
def zigzag(x):
    """compute a vector of difference between reversal points of the input
    example:
    x1 = randn(200)
    x=x1.cumsum()
    y=tools.zigzag(x)
    f=figure()
    a=f.add_subplot(111)
    a.plot(x)
    a.plot(y.cumsum())
    """
    x = numpy.array(x)
    s = numpy.diff(numpy.sign(numpy.diff(x)))
    inds = (s != 0)
    
    out = numpy.zeros(len(x))
    lastx = 0
    for i, ind in enumerate(inds[1:]):
        if ind:
            out[i+2] = x[i+2]-lastx
            lastx = x[i+2]
    
    return out

def local_minima(x,window=1):
    N = len(x)
    minima_inds = numpy.zeros(N)
    
    for i in range(window,N-window):
        for w in range(window):
            if not x[i] < x[i+w+1] or not x[i] < x[i-w-1]:
                break
        else:
            minima_inds[i] = 1

    minima_inds = (minima_inds==1)
    return minima_inds
