# hawpy
Silas Vriend. Winter 2017.

Hawthorn Research Group.

Department of Physics and Astronomy.

University of Waterloo.

A Python module for reading, storing, and plotting spec data files.

## Project Overview
 
The objective is to write a Python package which can read in data
from a spec file, store it as a Python object, and then carry out plotting
routines according to user-specified parameters.

Ideally, this Python module would mimic the functionality of the 
SPEC_REIXSloader and ATools Igor procedures which were written by 
previous members of the Hawthorn Research Group.

The package will consist of a basic, well-documented set of modules
which may be improved upon by future members of the group.

For the first iteration of the package, the goal is to have a command
line interface which would be able to perform some of the essential
operations which the group needs to analyze data.

Code from the [pyspec package](https://github.com/stuwilkins/pyspec) written 
by Stuart Wilkins at Brookhaven National Laboratory could be modified to suit 
the needs of the group.
    
    
## Technical Specifications
    
The package will be written in Python 2.7.13.  The package will
take advantage of several modules from the SciPy stack, via the
[Anaconda distribution](https://www.continuum.io/downloads).
    
    
## spec Data File Structure

spec data files are ASCII files containing two types of lines:
header lines, and raw numerical data lines. Header lines are preceded 
by a # symbol. There are several different types of headers:
    
    #F      filename
    #E      Unix epoch time
    #D      Gregorian calendar date
    #C      comment line
    
    #O_     motor names
    #o_     motor name mnemonics
    
    #J_     counter names
    #j_     counter name mnemonics
    
    #S      scan header info
    #T      time interval per data point
    #G_     initial scan geometry information
    #Q      initial reciprocal lattice coordinates
    #P_     initial motor settings (given in same order as in #O_ lines)
    #N      number of data columns in the scan
    #L      data column labels
    
The remaining lines of the spec data file contain raw numerical data
for all motors and counters at each time interval during the scan. The 
values are separated by spaces and organized into columns.

### Scan types: 

- ascan is a single-motor absolute-position scan. 
- a2scan is a two-motor absolute-position scan.
- Escan is an energy scan.
- mesh is a nested two-motor absolute-position scan.
- rscan is a one-motor relative-position scan.
- timescan is an indefinite time-lapse scan.
    
### Scan syntaxes:

Note the spacing between select words.

    ascan  motor start finish  intervals time
    a2scan  m1 s1 f1  m2 s1 f2  intervals time
    Escan start end intervals time
    mesh  m1 s1 f1 intervals1  m2 s2 f2 intervals2  time
    rscan motor s1 f1 intervals1 f2 intervals2 f3 intervals3 [...] time
    timescan count_time sleep_time

## Desired Functionality

"*" indicates functionality which could be implemented in future iterations.

### From SPEC_REIXSloader:

- Command to select a folder which contains the spec data file.
- Command to load a spec data file.
- Progress indicator is shown as data file is loaded.
- A list of scan headers from the file is generated and displayed.
- Several scan parameters are read and displayed from the
    currently highlighted scan header.
- Loader automatically determines the appropriate x-axis based on
    the scan header.
- A list of possible y-axes is presented. 
- The user can opt to generate multiple plots by selecting multiple
    y-axis options and using the "Display" command.
- The user can opt to append an additional trace to an existing plot
    by selecting a new y-axis and using the "Append" command.
        
### From ATools:

- Command to merge arbitrarily many graphs with identical axes 
    onto the same plot. The user can select the graph from which
    the new merged graph should inherit features.
- *Command to plot T-scan.
- *Command to fix theta scan relative to H, K, or L. (Switch x-axis)
- *Command to order scans by T. (For vertical scan stack.)
- *Commands to format graphs in predetermined ways.
- *Command to subtract background noise from graph.
- *Command to subtract hidden polynomial noise from graph.
- *Command to perform Lorentzian fit to graph.
- *Command to return the fit parameters.
- *Command to export visible graphs as PDFs.
