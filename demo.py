# The following is a test script to demonstrate the features of this module.
# This code is not executed if this module is imported in a different
#   script or module.
__verbose__ = False

LNSCO = SpecDataFile('LNSCO')
YBCO = SpecDataFile('YBCO_XAS')

SCAN1 = LNSCO[298]
SCAN2 = YBCO[11]
SCAN3 = YBCO[9]


# STANDARD PLOT TEST.
SCAN1.do_plot(ycol='ChT_REIXS')


# STANDARD PLOT TEST without normalization.
SCAN1.do_plot(ycol='ChT_REIXS', mcol=None)


# STANDARD PLOT TEST using object oriented methods.
#
# This allows for editing of the plot title and the axis labels.
#
# The SpecPlot class sets the title and x- and y-axis labels automatically,
#   but they can be modified with the getter and setter methods shown below.
#
# Function aliasing could be used to make the commands shorter.
#   For example, one could add the following lines of code to the body of
#   the SpecPlot class:
#
#       setXL = set_xlabel
#       setYL = set_ylabel
#       appXL = xlabel_append
#       appYL = ylabel_append
#
PLOT3 = SCAN1.do_plot(ycol='ChT_REIXS')
PLOT3.set_title('This is a standard plot with a custom title.')
PLOT3.title_append(' Neat!')
PLOT3.set_xlabel('Two Theta (degrees)')
PLOT3.ylabel_append(' (arb.units)')


#STANDARD PLOT TEST with Lorentzian fit.
PLOT4 = SCAN1.do_plot(ycol='ChT_REIXS')
PLOT4.lorentz_fit()


# STANDARD PLOT TEST with Lorentzian fit and zero-to-one offset/scale.
PLOT5 = SCAN1.do_plot(ycol='ChT_REIXS')
PLOT5.lorentz_fit()
PLOT5.scale_offset(124, 125.5)


# STANDARD PLOT TEST with multiple scans.
SCANLIST = [173, 169, 165, 161, 155, 151, 147, 143, 139,
            135, 103, 131, 127, 125, 121, 117, 109, 113]

PLOT6 = SpecPlot()
PLOT6.plotMT(LNSCO, SCANLIST)
PLOT6.set_title('This is multiple scans, coloured in order.')


# STANDARD PLOT TEST with a range of scans.
PLOT7 = SpecPlot()
PLOT7.plotR(LNSCO, 6, 9, ycol='MCP_REIXS', mcol=None)
PLOT7.set_title('This is a range of scans.')


# MESH PLOT TEST.
SCAN2.do_plot(ycol='TEY_REIXS')


# MESH PLOT TEST without normalization.
SCAN2.do_plot(ycol='TEY_REIXS', mcol=None)


# MESH PLOT TEST using object oriented methods.
#
# This allows for editing of the title, colorbar label and axis labels.
#
# Note that there is not a get_clb_label or clb_label_append method for the
#   SpecPlot class. This is because the current implementation of the
#   matplotlib.colorbar.Colorbar class does not have a get_label() method.
#
PLOT10 = SCAN3.do_plot(ycol='TEY_REIXS')
PLOT10.set_title('This is a mesh plot with a custom title.')
PLOT10.title_append(' Wow!')
PLOT10.set_clb_label('TEY_REIXS / I0_BD3 (arb. units)')
PLOT10.xlabel_append(' (mm)')
PLOT10.ylabel_append(' (mm)')


# Calling plt.show() at the very end shows all of the plots at once.
#   Calling it earlier would mean that the plots created after the call to
#   plt.show() would not be displayed.
plt.show()
