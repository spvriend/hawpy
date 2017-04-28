"""A test script to demonstrate the features of the hawpy module.

Silas Vriend. Hawthorn Research Group. Winter 2017.
University of Waterloo. Department of Physics of Astronomy.

This script makes use of the function aliases defined within the SpecPlot class.
The aliases are as follows:

    plR = plot_scan_range
    plMT = plot_multiple_traces
    Lfit = lorentz_fit
    zero2one = offset_and_scale
    
    setTL = set_title
    setXL = set_xlabel
    setYL = set_ylabel
    setCBL = set_clb_label
    
    appTL = title_append
    appXL = xlabel_append
    appYL = ylabel_append

Note that any changes to the labels must occur after any of the first four
functions listed above are called. Otherwise, the changes to the labels will
be undone.

"""

import hawpy as hp
import matplotlib.pyplot as plt

__verbose__ = False


# Instantiate SpecDataFile instances from the spec files you are dealing with.
LNSCO = hp.SpecDataFile('LNSCO')
YBCO = hp.SpecDataFile('YBCO_XAS')

# Create SpecScan objects by passing SpecDataFile objects the scan number.
SCAN1 = LNSCO[298]
SCAN2 = YBCO[11]
SCAN3 = YBCO[9]


# Standard plot. 
# 
# The user must provide the y-data label using the exact characters from 
#   the original spec data file.
SCAN1.pl(ycol='ChT_REIXS')


# Standard plot without normalization.
#
# To suppress data normalization, the keyword argument mcol must be set equal
#   to None. The default value of mcol is 'I0_BD3'.
SCAN1.pl(ycol='ChT_REIXS', mcol=None)


# Standard plot using object-oriented methods.
#
# Using the object-oriented methods demonstrated here allows for the editing of
#   plot labels. The labels are generated automatically otherwise.
PLOT3 = SCAN1.pl(ycol='ChT_REIXS')
PLOT3.setTL('This is a standard plot with a custom title.')
PLOT3.appTL(' Neat!')
PLOT3.setXL('Two Theta (degrees)')
PLOT3.appYL(' (arb.units)')


# Standard plot with Lorentzian fit.
#
# The lorentz_fit() method loops through each trace on the plot and tries to
#   fit the curve to the Lorentz distribution. The fit parameters are displayed
#   in the legend.
#
PLOT4 = SCAN1.pl(ycol='ChT_REIXS')
PLOT4.Lfit()


# Standard plot with multiple scans.
#
# Ordering the scans by temperature is up to the user. The colour gradient
#   will flow from purple at the beginning of the scan list to red at the end
#   of the scan list.
#
SCANLIST = [173, 169, 165, 161, 155, 151, 147, 143, 139,
            135, 103, 131, 127, 125, 121, 117, 109, 113]

PLOT6 = hp.SpecPlot()
PLOT6.plMT(LNSCO, SCANLIST)
PLOT6.setTL('Multiple scans, coloured in order.')


# Standard plot with multiple scans, then zero-to-one offset-and-scale.
SCANLIST = [173, 169, 165, 161, 155, 151, 147, 143, 139,
            135, 103, 131, 127, 125, 121, 117, 109, 113]
PLOT5 = hp.SpecPlot()
PLOT5.plMT(LNSCO, SCANLIST)
PLOT5.zero2one(58, 60.5)
PLOT5.setTL('This is a zero-to-one offset-and-scale test.')


# Standard plot with a range of scans.
PLOT7 = hp.SpecPlot()
PLOT7.plR(LNSCO, 6, 9, ycol='MCP_REIXS', mcol=None)
PLOT7.setTL('This is a range of scans.')


# Mesh plot.
SCAN2.pl(ycol='TEY_REIXS')


# Mesh plot without normalization.
SCAN2.pl(ycol='TEY_REIXS', mcol=None)


# Mesh plot using object oriented methods.
PLOT10 = SCAN3.pl(ycol='TEY_REIXS')
PLOT10.setTL('This is a mesh plot with a custom title.')
PLOT10.appTL(' Wow!')
PLOT10.setCBL('TEY_REIXS / I0_BD3 (arb. units)')
PLOT10.appXL(' (mm)')
PLOT10.appXL(' (mm)')


# Calling plt.show() at the very end shows all of the plots at once.
plt.show()
