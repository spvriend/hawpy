"""A module of graphing routines which make use of hawpy."""

import hawpy
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize as opt

from matplotlib import cm

hawpy.__verbose__ = False


def lorentzian(x, x0, gamma):
    """Returns the y value of the Lorentzian distribution for the given x."""
    y = (1.0 / (math.pi*gamma)) * (1.0 / (1 + ((x-x0) / float(gamma))**2))
    return y
    
def plot_multiple_traces(filename, scanlist,
                         xcol='TwoTheta', ycol='ChT_REIXS', mcol=None):
    """A routine for plotting a range of scans, coloured according to order."""
    specfile = hawpy.SpecDataFile(filename)

    mpl.rcParams['lines.linewidth'] = 0.5

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

    plt.xlim((58, 62))
    plt.title('Cu edge 001\nWednesday (From ATS CU)')
    plt.xlabel('Two Theta (degrees)')
    plt.ylabel('Channeltron (REIXS) / I0 (arb. units)')

    plt.legend()
    plt.show()
    
def plot_lorentzian_fit(filename, scan_no, xcol='TwoTheta', 
                        ycol='ChT_REIXS', mcol=None):
    """Perform and plot a Lorentzian fit of the given x-y plot."""
    specfile = hawpy.SpecDataFile(filename)
    
    mpl.rcParams['lines.linewidth'] = 0.5
    
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
    

if __name__ == '__main__':
    # SCANLIST = [173, 169, 165, 161, 155, 151, 147, 143, 139,
                # 135, 103, 131, 127, 125, 121, 117, 109, 113]

    # plot_multiple_traces('LNSCO', SCANLIST, xcol='TwoTheta',
                         # ycol='ChT_REIXS', mcol='I0_BD3')
    
    plot_lorentzian_fit('LNSCO', 298, xcol='TwoTheta', ycol='ChT_REIXS', mcol='I0_BD3')
    
