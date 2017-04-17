"""A module of graphing routines which make use of hawpy."""

import hawpy

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

hawpy.__verbose__ = False

def plot_multiple_scans(filename, scanlist, 
                        xcol='TwoTheta', ycol='ChT_REIXS', mcol=None):
    """A routine for plotting a range of scans, coloured by temperature."""
    specfile = hawpy.SpecDataFile(filename)
    
    mpl.rcParams['lines.linewidth'] = 0.5
    
    i = 0
    
    for scan_no in scanlist:
        scan = specfile[scan_no]
        xdata = getattr(scan, xcol)
        ydata = getattr(scan, ycol)
        
        if mcol is not None:
            ydata = ydata / getattr(scan, mcol)
        
        plt.plot(xdata, ydata, label='S{}'.format(scan_no), color=plt.cm.gist_rainbow((len(scanlist)-i)/float(len(scanlist))))
        i += 1
    
    plt.xlim((58,62))
    plt.title('Cu edge 001\nWednesday (From ATS CU)')
    plt.xlabel('Two Theta (degrees)')
    plt.ylabel('Channeltron (REIXS) / I0 (arb. units)')

    plt.legend()
    plt.show()


if __name__ == '__main__':
    scanlist = [173, 169, 165, 161, 155, 151, 147, 143, 139, 
                135, 103, 131, 127, 125, 121, 117, 109, 113]

    plot_multiple_scans('LNSCO', scanlist, xcol='TwoTheta', 
                        ycol='ChT_REIXS', mcol='I0_BD3')
