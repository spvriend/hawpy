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
from bisect import bisect_left


def plot_scan_range(filename, start, end,
                    xcol='TwoTheta', ycol='ChT_REIXS', mcol='I0_BD3', **kwargs):
    """A routine for plotting a range of consecutive scans."""
    scanrange = [x for x in range(start, end+1)]
    try:
        plot = plot_multiple_traces(filename, scanrange, 
                                    xcol, ycol, mcol, **kwargs)
        return plot
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
    
    plot = (fig, ax, lines)
    
    return plot
    

# Function aliases.
plotR = plot_scan_range
plotMT = plot_multiple_traces
plotLF = plot_lorentzian_fit
zero2one = scale_offset


if __name__ == '__main__':
    hawpy.__verbose__ = False
    
    SCANLIST = [173, 169, 165, 161, 155, 151, 147, 143, 139,
                135, 103, 131, 127, 125, 121, 117, 109, 113]
    
    PLOT1 = plotMT('LNSCO', SCANLIST, xcol='TwoTheta',
                         ycol='ChT_REIXS', mcol='I0_BD3', linewidth=0.5)
                         
    PLOT2 = plotMT('LNSCO', SCANLIST, xcol='TwoTheta',
                         ycol='ChT_REIXS', mcol='I0_BD3', linewidth=0.5)

    PLOT2 = plotLF('LNSCO', 298, xcol='TwoTheta', ycol='ChT_REIXS', mcol='I0_BD3')

    PLOT3 = plotR('LNSCO', 6, 9, xcol='MonoEngy', ycol='MCP_REIXS', mcol=None)
    
    scale_offset(PLOT2, 58, 60)
    
    plt.show()