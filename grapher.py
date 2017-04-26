# hawpy.py (c) Silas Vriend 2017
#
# Silas Vriend. Hawthorn Research Group. Winter 2017.
# University of Waterloo Department of Physics and Astronomy.

"""A module of graphing routines which make use of hawpy."""

import hawpy
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize as opt

from matplotlib import cm


def lorentzian(x, x0, gamma):
    """Returns the y value of the Lorentzian distribution for the given x."""
    y = (1.0 / (math.pi*gamma)) * (1.0 / (1 + ((x-x0) / float(gamma))**2))
    return y

def plot_scan_range(filename, start, end,
                    xcol='TwoTheta', ycol='ChT_REIXS', mcol='I0_BD3', **kwargs):
    """A routine for plotting a range of consecutive scans."""
    scanrange = [x for x in range(start, end+1)]
    try:
        fig, ax, lines = plot_multiple_traces(filename, scanrange, xcol, 
                                              ycol, mcol, **kwargs)
        return fig, ax, lines
    except ValueError:
        print 'Please esnsure that all scans in the range are of the same type.'
    

def plot_multiple_traces(filename, scanlist, xcol='TwoTheta',
                         ycol='ChT_REIXS', mcol='I0_BD3', **kwargs):
    """A routine for plotting a range of scans, coloured according to order.
    
    Any keyword arguments passed through the **kwargs paramter must be
    keyword arguments for the pyplot.plot() function. See the Matplotlib API 
    docs for more details.
    """
    specfile = hawpy.SpecDataFile(filename)
    
    fig, ax = plt.subplots()
    
    i = 0
    
    lines = tuple()
    
    for scan_no in scanlist:
        scan = specfile[scan_no]
        xdata = getattr(scan, xcol)
        ydata = getattr(scan, ycol)

        if mcol is not None:
            mdata = getattr(scan, mcol)
            ydata = ydata / mdata

        color_map = cm.get_cmap(name='gist_rainbow')

        line, = plt.plot(xdata, ydata, label='S{}'.format(scan_no),
                         color=color_map(1 - float(i)/len(scanlist)), **kwargs)
        lines += line,
        i += 1

    ax.set_title('Multiple scans from {}'.format(filename))
    ax.set_xlabel(xcol)
    ax.legend()

    if mcol is not None:
        ax.set_ylabel('{} / {}'.format(ycol, mcol))
    else:
        ax.set_ylabel('{}'.format(ycol))
    
    return fig, ax, lines
    
def plot_lorentzian_fit(filename, scan_no, xcol='TwoTheta',
                        ycol='ChT_REIXS', mcol='I0_BD3'):
    """Perform and plot a Lorentzian fit of the given xy plot."""
    specfile = hawpy.SpecDataFile(filename)

    scan = specfile[scan_no]
    xdata = getattr(scan, xcol)
    ydata = getattr(scan, ycol)
    
    fig, ax = plt.subplots()
    
    if mcol is not None:
        mdata = getattr(scan, mcol)
        ydata /= mdata

    ymax = max(ydata)

    ydata /= ymax

    popt, pcov = opt.curve_fit(lorentzian, xdata, ydata)

    fit_ydata = lorentzian(xdata, *popt)

    ydata *= ymax
    fit_ydata *= ymax

    ax.set_title('Lorentzian fit of S{}'.format(scan_no))
    ax.set_xlabel('{}'.format(xcol))

    if mcol is not None:
        ax.set_ylabel('{} / {}'.format(ycol, mcol))
    else:
        ax.set_ylabel('{}'.format(ycol))

    line1, = plt.plot(xdata, fit_ydata, label='Lorentzian fit.\nx0={}\n$\gamma={}$'.format(popt[0], popt[1]))
    line2, = plt.plot(xdata, ydata, label='S{}'.format(scan_no))

    lines = line1, line2
        
    ax.legend()
    
    return fig, ax, lines

def scale_offset(plot, region1, region2):
    """This function performs a scale and offset on a given plot.
    
    The basic algorithm for a scale and offset is this:
    
        Select region 1. Region 1 is the region to offset relative to. For
        simplicity, this will be a single x-value.
        
        Select region 2. Region 2 is the region to scale relative to. For
        simpliticty, this will also be a single x-value.
        
        For each trace on the plot: determine the y-value of the trace at the
        x-value in region 1. Subtract that y-value from the entire trace.
        
        For each trace on the plot: determine the y-value of the trace at the
        x-value in region 2. Divide the entire trace by that y-value, unless
        the y-value is zero.
    """
    
    pass
    

# Function aliases.
plotR = plot_scan_range
plotMT = plot_multiple_traces
plotLF = plot_lorentzian_fit
zero2one = scale_offset


if __name__ == '__main__':
    hawpy.__verbose__ = False
    
    SCANLIST = [173, 169, 165, 161, 155, 151, 147, 143, 139,
                135, 103, 131, 127, 125, 121, 117, 109, 113]
    
    plotMT('LNSCO', SCANLIST, xcol='TwoTheta',
                         ycol='ChT_REIXS', mcol='I0_BD3', linewidth=0.5)

    plotLF('LNSCO', 298, xcol='TwoTheta', ycol='ChT_REIXS', mcol='I0_BD3')

    plotR('LNSCO', 6, 9, xcol='MonoEngy', ycol='MCP_REIXS', mcol=None)
    
    plt.show()