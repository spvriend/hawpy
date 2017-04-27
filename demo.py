import hawpy as hp
import matplotlib.pyplot as plt

hp.__verbose__ = False

LNSCO = hp.SpecDataFile('LNSCO')
YBCO = hp.SpecDataFile('YBCO_XAS')

SCAN1 = LNSCO[298]
SCAN2 = YBCO[11]
SCAN3 = YBCO[9]

# STANDARD PLOT.
SCAN1.do_plot(ycol='ChT_REIXS')

# STANDARD PLOT without normalization.
SCAN1.do_plot(ycol='ChT_REIXS', mcol=None)

# STANDARD PLOT using object oriented methods.
# The SpecPlot class sets the title and x- and y-axis labels automatically,
#   but they can be modified with the getter and setter methods shown below.
PLOT3 = SCAN1.do_plot(ycol='ChT_REIXS')
PLOT3.set_title('This is a standard plot with a custom title.')
PLOT3.title_append(' Neat!')
PLOT3.set_xlabel('Two Theta (degrees)')
PLOT3.ylabel_append(' (arb.units)')

# MESH PLOT.
SCAN2.do_plot(ycol='TEY_REIXS')

# MESH PLOT without normalization.
SCAN2.do_plot(ycol='TEY_REIXS', mcol=None)

# MESH PLOT using object oriented methods.
# This allows for editing of the title, colorbar label and axis labels.
PLOT6 = SCAN3.do_plot(ycol='TEY_REIXS')
PLOT6.set_title('This is a mesh plot with a custom title.')
PLOT6.title_append(' Wow!')
PLOT6.set_clb_label('TEY_REIXS / I0_BD3 (arb. units)')
PLOT6.xlabel_append(' (mm)')
PLOT6.ylabel_append(' (mm)')

# Calling plt.show() at the very end shows all of the plots at once.
#   Calling it earlier would mean that the plots created after the call to
#   plt.show() would not be displayed.
plt.show()
