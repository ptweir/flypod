import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib._cm as _cm


def loop_cmap(cmap_name, num_loops=2, reflect=False, plot=False):
    '''
    tool to loop colormap, returns colormap
    
    arguments:
    cmap_name       colormap to start with
    num_loops       number of times to repeat colormap (integer, default 2)
    reflect         reflect every other iteration of colormap (boolean, default False)
    plot            plot RGB data of new colormap (boolean, default False)
    
    example:
    my_hot = loop_cmap('hot',num_loops=4, reflect=True, plot=True)
    
    Based on code by Gary Ruben found at http://www.mail-archive.com/matplotlib-users@lists.sourceforge.net/msg14237.html
    PTW 2010
    '''
    if type(cmap_name) is str:
        cmap = eval('_cm._%s_data' % cmap_name)
        cmap_str = cmap_name
    else:
        try:    
            cmap = eval('_cm._%s_data' % cmap_name.name)
        except AttributeError:
            cmap = cmap_name._segmentdata
        cmap_str = cmap_name.name
        
    num_loops = round(num_loops)
    LUTSIZE = plt.rcParams['image.lut']
    r = np.array(cmap['red'])
    g = np.array(cmap['green'])
    b = np.array(cmap['blue'])
    
    Nr = r[:,0].shape[-1]
    Ng = g[:,0].shape[-1]
    Nb = b[:,0].shape[-1]
    R = np.empty((Nr*num_loops,3))
    G = np.empty((Ng*num_loops,3))
    B = np.empty((Nb*num_loops,3))
    for n in range(num_loops):
        if reflect and np.mod(n,2)==1:
            R[n*Nr:(n+1)*Nr,0] = (n + 1 - r[::-1,0])/num_loops
            R[n*Nr:(n+1)*Nr,1:] = r[::-1,1:]
            G[n*Ng:(n+1)*Ng,0] = (n + 1 - g[::-1,0])/num_loops
            G[n*Ng:(n+1)*Ng,1:] = g[::-1,1:]
            B[n*Nb:(n+1)*Nb,0] = (n + 1 - b[::-1,0])/num_loops
            B[n*Nb:(n+1)*Nb,1:] = b[::-1,1:]
        else:
            R[n*Nr:(n+1)*Nr,0] = r[:,0]/num_loops + n/num_loops
            R[n*Nr:(n+1)*Nr,1:] = r[:,1:]
            G[n*Ng:(n+1)*Ng,0] = g[:,0]/num_loops + n/num_loops
            G[n*Ng:(n+1)*Ng,1:] = g[:,1:]
            B[n*Nb:(n+1)*Nb,0] = b[:,0]/num_loops + n/num_loops
            B[n*Nb:(n+1)*Nb,1:] = b[:,1:]
        
    _my_data = {'red':   tuple(map(tuple,R)),
                'green': tuple(map(tuple,G)),
                'blue':  tuple(map(tuple,B))
               }
    my_cmap = colors.LinearSegmentedColormap('my_'+cmap_str, _my_data,  LUTSIZE)

    if plot:
        plt.figure()
        plt.plot(R[:,0], R[:,1], 'r', G[:,0], G[:,1], 'g', B[:,0],  B[:,1], 'b', lw=3)
        plt.axis(ymin=-0.2, ymax=1.2)

    return my_cmap
def reflect_cmap(cmap_name, plot=False):
    '''
    reflect a colormap once
    DEPRECATED: use loop_cmap
    
    Based on code by Gary Ruben found at http://www.mail-archive.com/matplotlib-users@lists.sourceforge.net/msg14237.html
    PTW 2010
    '''
    if type(cmap_name) is str:
        cmap = eval('_cm._%s_data' % cmap_name)
        cmap_str = cmap_name
    else:
        try:    
            cmap = eval('_cm._%s_data' % cmap_name.name)
        except AttributeError:
            cmap = cmap_name._segmentdata
        cmap_str = cmap_name.name
        
    LUTSIZE = plt.rcParams['image.lut']
    r = np.array(cmap['red'])
    g = np.array(cmap['green'])
    b = np.array(cmap['blue'])
    
    print g
    Nr = r[:,0].shape[-1]
    Ng = g[:,0].shape[-1]
    Nb = b[:,0].shape[-1]
    R = np.empty((Nr*2,3))
    G = np.empty((Ng*2,3))
    B = np.empty((Nb*2,3))

    R[:Nr,0] = r[:,0]/2
    R[:Nr,1:] = r[:,1:]
    G[:Ng,0] = g[:,0]/2
    G[:Ng,1:] = g[:,1:]
    B[:Nb,0] = b[:,0]/2
    B[:Nb,1:] = b[:,1:]
    
    R[Nr:,0] = 1 - r[::-1,0]/2
    R[Nr:,1:] = r[::-1,1:]
    G[Ng:,0] = 1 - g[::-1,0]/2
    G[Ng:,1:] = g[::-1,1:]
    B[Nb:,0] = 1 - b[::-1,0]/2
    B[Nb:,1:] = b[::-1,1:]
        
    _my_data = {'red':   tuple(map(tuple,R)),
                'green': tuple(map(tuple,G)),
                'blue':  tuple(map(tuple,B))
               }
    my_cmap = colors.LinearSegmentedColormap('my_'+cmap_str, _my_data,  LUTSIZE)

    if plot:
        plt.figure()
        plt.plot(R[:,0], R[:,1], 'r', G[:,0], G[:,1], 'g', B[:,0],  B[:,1], 'b', lw=3)
        plt.axis(ymin=-0.2, ymax=1.2)

    return my_cmap
    
def rescale_cmap(cmap_name, rLow=0.0, rHigh=1.0, gLow=0.0, gHigh=1.0, bLow=0.0, bHigh=1.0, plot=False):
    '''
    Example 1:
        my_hsv = rescale_cmap('hsv', low = 0.3)     # equivalent scaling  
                    to cplot_like(blah, l_bias=0.33, int_exponent=0.0)
    Example 2:
        my_hsv = rescale_cmap(cm.hsv, low = 0.3)
    based on code by Gary Ruben found at http://www.mail-archive.com/matplotlib-users@lists.sourceforge.net/msg14237.html
    '''
    if type(cmap_name) is str:
        cmap = eval('_cm._%s_data' % cmap_name)
        cmap_str = cmap_name
    else:
        try:    
            cmap = eval('_cm._%s_data' % cmap_name.name)
        except AttributeError:
            cmap = cmap_name._segmentdata
        cmap_str = cmap_name.name
    LUTSIZE = plt.rcParams['image.lut']
    r = np.array(cmap['red'])
    g = np.array(cmap['green'])
    b = np.array(cmap['blue'])
    rRange = rHigh - rLow
    gRange = gHigh - gLow
    bRange = bHigh - bLow
    r[:,1:] = r[:,1:]*rRange+rLow
    g[:,1:] = g[:,1:]*gRange+gLow
    b[:,1:] = b[:,1:]*bRange+bLow
    _my_data = {'red':   tuple(map(tuple,r)),
                'green': tuple(map(tuple,g)),
                'blue':  tuple(map(tuple,b))
               }
    my_cmap = colors.LinearSegmentedColormap('my_'+cmap_str, _my_data,  LUTSIZE)

    if plot:
        plt.figure()
        plt.plot(r[:,0], r[:,1], 'r', g[:,0], g[:,1], 'g', b[:,0],  b[:,1], 'b', lw=3)
        plt.axis(ymin=-0.2, ymax=1.2)

    return my_cmap

def get_cmap(cmap_name):
    if cmap_name is 'gb180':
        #180 degree rotationally symmetric colorblind-frienly colormap
        mygbdata = {
        'red' : ((0., 0., 0.), (1., 0., 0.)),
        'green': ((0., 0., 0.), (0.25, 1., 1.), (.5, 1., 1.), (.75, 0., 0.), (1., 0., 0.)),
        'blue' : ((0., 0., 0.), (0.25, 0., 0.), (.5, 1., 1.), (.75, 1., 1.), (1., 0., 0.))
        }
        mygb = colors.LinearSegmentedColormap('mygb', mygbdata)

        cmap = loop_cmap(mygb)
    elif cmap_name is 'hsv180':
        #180 degree rotationally symmetric colorblind-frienly colormap
        cmap = loop_cmap('hsv')

    #my_hot = cmt.loop_cmap('hot',num_loops=2, reflect=True)
    #my_gray = cmt.loop_cmap('gray',num_loops=2, reflect=False)
    
    return cmap
    
def add_colordisc(image,width=101):

    height = width
    R = np.floor(width/2)
    y0=R
    x0=R
    
    x = np.matrix(range(width))
    y = np.transpose(np.matrix(range(height)))

    X = np.array((0*y+1)*x - y0)
    Y = np.array(y*(0*x+1) - x0)

    angle = np.arctan2(Y,X) + np.pi/2
    radius = np.sqrt(X**2 + Y**2)
    
    ring = (radius<R)
    radius[~ring]=np.nan
    radius[ring]=1
    
    image[:height,:width] = angle*radius
    
    return image
