# Select a few points on your curve and have data to be scaled and offset to
# those points. Done in Igor frequently.
# Start with one point which you'll offset multiple curves to
# Choose one x-value, call that y-value zero
# Choose one x-value, call that y-value one



"""A module of graphing routines which make use of hawpy."""

import hawpy
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize as opt

from matplotlib import cm

mpl.rcParams['lines.linewidth'] = 0.5
hawpy.__verbose__ = False


def lorentzian(x, x0, gamma):
    """Returns the y value of the Lorentzian distribution for the given x."""
    y = (1.0 / (math.pi*gamma)) * (1.0 / (1 + ((x-x0) / float(gamma))**2))
    return y

def plot_scan_range(filename, start, end, 
                    xcol='TwoTheta', ycol='ChT_REIXS', mcol='I0_BD3'):
    """A routine for plotting a range of consecutive scans."""
    scanrange = [x for x in range(start, end+1)]
    try:
        plot_multiple_traces(filename, scanrange, xcol, ycol, mcol)
    except ValueError:
        print 'Please esnsure that all scans in the range are of the same type.'
        
    
def plot_multiple_traces(filename, scanlist,
                         xcol='TwoTheta', ycol='ChT_REIXS', mcol='I0_BD3'):
    """A routine for plotting a range of scans, coloured according to order."""
    specfile = hawpy.SpecDataFile(filename)

    i = 0

    for scan_no in scanlist:
        scan = specfile[scan_no]
        xdata = getattr(scan, xcol)
        ydata = getattr(scan, ycol)

        if mcol is not None:
            mdata = getattr(scan, mcol)
            ydata = ydata / mdata

            
        color_map = cm.get_cmap(name='gist_rainbow')

        plt.plot(xdata, ydata, label='S{}'.format(scan_no),
                 color=color_map(1 - float(i)/len(scanlist)))
        i += 1

    plt.title('Multiple scans from {}'.format(filename))
    plt.xlabel(xcol)
    plt.legend()

    if mcol is not None:
        plt.ylabel('{} / {}'.format(ycol, mcol))
    else:
        plt.ylabel('{}'.format(ycol))   
    
    plt.show()
    
def plot_lorentzian_fit(filename, scan_no, xcol='TwoTheta', 
                        ycol='ChT_REIXS', mcol='I0_BD3'):
    """Perform and plot a Lorentzian fit of the given xy plot."""
    specfile = hawpy.SpecDataFile(filename)
    
    scan = specfile[scan_no]
    xdata = getattr(scan, xcol)
    ydata = getattr(scan, ycol)
    
    if mcol is not None:
        mdata = getattr(scan, mcol)
        ydata /= mdata

    ymax = max(ydata)
    
    ydata /= ymax
        
    popt, pcov = opt.curve_fit(lorentzian, xdata, ydata)

    fit_ydata = lorentzian(xdata, *popt)
    
    ydata *= ymax
    fit_ydata *= ymax
    
    plt.title('Lorentzian fit of S{}'.format(scan_no))
    plt.xlabel('{}'.format(xcol))
    
    if mcol is not None:
        plt.ylabel('{} / {}'.format(ycol, mcol))
    else:
        plt.ylabel('{}'.format(ycol))    
    
    plt.plot(xdata, fit_ydata, label='Lorentzian fit.\nx0={}\n$\gamma={}$'.format(popt[0], popt[1]))
    plt.plot(xdata, ydata, label='S{}'.format(scan_no))
    
    plt.legend()
    plt.show()
    

# Function aliases.
plotR = plot_scan_range
plotMT = plot_multiple_traces
plotLF = plot_lorentzian_fit
    
    
    
if __name__ == '__main__':
    SCANLIST = [173, 169, 165, 161, 155, 151, 147, 143, 139,
                135, 103, 131, 127, 125, 121, 117, 109, 113]

    plotMT('LNSCO', SCANLIST, xcol='TwoTheta',
                         ycol='ChT_REIXS', mcol='I0_BD3')

    plotLF('LNSCO', 298, xcol='TwoTheta', ycol='ChT_REIXS', mcol='I0_BD3')

    plotR('LNSCO', 6, 9, xcol='MonoEngy', ycol='MCP_REIXS', mcol=None)