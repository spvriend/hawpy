# hawpy.py (c) Silas Vriend 2017
# Modified from the pyspec package (c) Stuart B. Wilkins 2008
#
# Silas Vriend. Hawthorn Research Group. Winter 2017.
# University of Waterloo Department of Physics and Astronomy.

"""A module for reading, storing, and plotting spec data files.

This module defines classes for reading spec data files into Python,
extracting data from specific scans, and creating standard x- vs. y-axis
plots or two-dimensional mesh images.

The __verbose__ parameter controls whether the module prints output to the
terminal when a script is run. To turn output off, add the following line to
your script:
    
    hawpy.__verbose__ = False

Classes:

    SpecDataFile :
        Acts as a spec-specific file object, with scan indexing built in.

    SpecScan :
        A class which represents one or more scans, which may be plotted.

    SpecScanHeader :
        Contains information from the scan header for a SpecScan object.

    SpecScanData :
        Contains the raw numerical data for a SpecScan object.

    SpecPlot :
        Contains methods to plot standard x- vs. y-axis data and to plot
        two-dimensional mesh images.

"""

import time
import copy

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import griddata

__verbose__ = True


def get_mesh_dims(scan):
    """Return a tuple containing the dimensions of the mesh grid.
    
    Parameters
    ----------
    scan : SpecScan
        `scan` is the scan to be checked.
        
    Returns
    -------
    mesh_dims : tuple of int
        The x and y dimensions of the mesh grid.
    
    """
    
    mesh_cmd = scan.header.scan_cmd.split()
    xint = int(mesh_cmd[4])
    yint = int(mesh_cmd[8])
    mesh_dims = (xint, yint)
    return mesh_dims

def istwod(scan):
    """Return True if the scan is a mesh scan.
    
    Parameters
    ----------
    scan : SpecScan
        The SpecScan object to check.
        
    Returns
    -------
    bool
        A bool indicating whether the scan is a mesh scan.
        
    """
    
    return scan.header.scan_type.strip() == 'mesh'


class SpecDataFile(object):
    """Acts as a spec-specific file object.

    To create a SpecScan object, call SpecDataFile[num], where `num` is
    the scan number as an integer, or a list of integer scan numbers.

    To print a string of statistics on the spec file, pass the SpecDataFile
    instance in a print statement.

    Parameters
    ----------
    fn : str
        The name of the spec data file from which to read data.
    
    Attributes
    ----------
    file : file object
        A temporary file object, created whenever the file is read.    
        
    filename : str
        The name of the spec data file.
        
    scan_index : dict
        A dict mapping scan numbers to scan locations within the file.

    motors : list
        A list of the motors found in #O headers within the spec file.

    motormap : dict
        A dict which maps motor mnemonics to motor names.
        
    scan_objects : dict
        A dict mapping scan numbers to SpecScan objects.

    """

    def __init__(self, fn):
        """Initialize an instance of the SpecDataFile class."""
        self.filename = fn
        self.file = open(self.filename, 'rb')
        self.scan_index = {}
        self.motors = []
        self.motormap = {}
        self.scan_objects = {}

        self._load_spec_file()
        return

    def __getitem__(self, item):
        """Call get_scan(item) when self[item] is called.
        
        Parameters
        ----------
        item : int
            The scan to be returned.
            
        Returns
        -------
        SpecScan
            A scan object corresponding to the scan number in item.
            
        """

        return self.get_scan(item, set_labels=True)

    def __str__(self):
        """Call show() when self is passed in a print statement."""
        return self.show()

    def _load_spec_file(self):
        # Call methods to index the file, read the header, and get stats.
        if __verbose__:
            print "**** Opening spec data file {}.".format(self.filename)

        if self.file.closed:
            self.file = open(self.filename, 'rb')

        self.index()
        self.read_header()
        
        self.file.close()
        return

    def _moveto(self, item):
        # Move to the location of the start of the scan in the data file.
        if item in self.scan_index:
            self.file.seek(self.scan_index[item])
        else:
            # Try re-indexing the file here.
            if __verbose__:
                print '**** Re-indexing the spec data file.\n'
            self.index()
            if item in self.scan_index:
                self.file.seek(self.scan_index[item])
            else:
                raise Exception('Scan {} is not in the file.'.format(item))

    def get_line(self):
        """Return the next line from the data file."""
        line = self.file.readline()
        return line

    def get_scan(self, item, set_labels=True, reread=False):
        """Create a single SpecScan object for the desired scan."""
        if self.file.closed:
            self.file = open(self.filename, 'rb')
            
        if __verbose__:
            print '**** Reading scan {}.'.format(item)
            
        if (item not in self.scan_objects) or (reread is True):
            self._moveto(item)
            self.scan_objects[item] = SpecScan(self, set_labels)

        self.file.close()
        return self.scan_objects[item]

    def index(self):
        """Read the file to generate the scan_index dict."""
        if __verbose__:
            print '---- Indexing scan :       ',

        self.file.seek(0, 0)

        pos = self.file.tell()
        line = self.file.readline()

        while line != '':
            if line[0:3] == '#S ':
                temp = line.split()
                num = int(temp[1])
                if (num % 5) == 0 and __verbose__:
                    print '\b\b\b\b\b\b\b{:5d} '.format(num),
                self.scan_index[num] = pos

            pos = self.file.tell()
            line = self.file.readline()

        if __verbose__:
            print '\b\b\b\b\b\b\bDONE '

        return

    def read_header(self):
        """Set attributes based on information in the scan header.

        Currently, this only reads in the names of the motors as contained
        on #O control lines.

        """

        if __verbose__:
            print '---- Reading file header.'

        self.file.seek(0, 0)
        line = self.file.readline()
        mnems = []
        names = []

        while (line[0:2] != '#S') and (line != ''):
            if line[0:2] == '#o':
                mnems.extend(line.split()[1:])
            elif line[0:2] == '#O':
                names.extend(line.split()[1:])
                self.motors.extend(line.split()[1:])
            line = self.file.readline()

        for mnem, name in zip(mnems, names):
            self.motormap[mnem] = name

        return

    def reset(self):
        """Reset the scan_objects dict as though no scans had been read."""
        self.scan_objects.clear()
        if __verbose__:
            print '**** SpecDataFile reset (all scan objects removed).'

    def reread(self):
        """Reload the data file into the SpecDataFile object."""
        if __verbose__:
            print '**** Reloading the spec file.'
            self._load_spec_file()

    def show(self, head='---- '):
        """Return a string of statistics on the data file.
        
        Parameters
        ----------
        head : str, optional
            The sring of characters which precedes each line of information.
            
        Returns
        -------
        statistics : str
            A string giving the number of scans, and the first and last scans.
            
        """
        
        file_length = len(self.scan_index)
        start_scan = min(self.scan_index.keys())
        end_scan = max(self.scan_index.keys())

        statistics = ''
        statistics += head + 'Spec file contains {} scans.\n'.format(file_length)
        statistics += head + 'Start scan = {}\n'.format(start_scan)
        statistics += head + 'End scan = {}\n'.format(end_scan)

        return statistics

class SpecScan(object):
    """This class represents a single scan from the data file.

    Currently, #G control lines are not read in, as more information is
    needed about the information that they contain. Information about
    #G geometry control lines is contained in files within the macros
    folder for the version of spec at the REIXS beamline.

    macros/fourc.src and macros/ub.mac should have the
    identifiers for #G0 and #G1. The interpretation of the remaining #G
    control lines could come straight from Stuart Wilkins' pyspec package.

    Attributes:

        set_labels : bool
            Decides whether motor labels are passed from data file to scan.

        specfile : SpecDataFile
            The spec data file from which the scan is to be read.

        header : SpecScanHeader
            An object containing information from the header.

        data : SpecScanData
            An object containing the scan data.

    """

    def __init__(self, specfile, set_labels=True, mask=None):
        """Initialize an instance of the SpecScan class.

        Parameters
        ----------
        specfile : SpecDataFile         Instance of a spec data file.
        item : int or list              The scan number(s) to read in.
        set_labels : bool               If true, set motor labels as keys in
                                            the dict of class variables.

        """

        self.specfile = specfile
        self.header = SpecScanHeader()
        self.data = SpecScanData()
        self.set_labels = set_labels

        line = self.specfile.get_line()

        sline = line.strip().split()

        self.header.scan_no = int(sline[1])
        self.header.scan_type = sline[2]
        self.header.scan_cmd = ' '.join(sline[2:])

        self.header.fullheader = line

        line = self.specfile.get_line()
        self.header.fullheader += line

        while (line[0:2] != '#L') and (line != ''):
            self.header.parse_header_line(self, line)
            line = self.specfile.get_line()

        if line[0:2] == '#L':
            self.header.labels = line[3:].split()

        line = self.specfile.get_line()

        while (line[0:2] != '#S') and (line != ''):
            if line[0] != '#':
                self.data.parse_data_line(line)
            elif line[0:2] == '#C':
                self.header.parse_comment_line(line)
            else:
                self.header.fullheader += line

            line = self.specfile.get_line()

        if isinstance(self.data.raw, list):
            self.data.raw = np.asarray(self.data)

        if mask != None:
            if __verbose__:
                print '---- Removing rows {} from data.'.format(str(mask))
            self.data = np.delete(self.data, mask, axis=0)


        self._setcols()

        if __verbose__:
            print '---- Data is {} rows x {} cols.'.format(self.data.raw.shape[0],
                                                           self.data.raw.shape[1])
        return None

    def __str__(self):
        """Call show() when the class is passed in a print statement."""
        return self.show()

    def _setcols(self):
        """Map column labels to column data arrays.

        This is done both within the SpecScan object and within its
        corresponding SpecScanData object. Creates SpecScan attributes
        from all of the key-valuepairs.

        """

        # If there is at least one row of data,
        if self.data.raw.shape[0] > 0:
            # go through each column label
            for i in range(len(self.header.labels)):
                # and if there is more than one row of data
                if len(self.data.raw.shape) == 2:
                    # map the column label to that column of data
                    motor_name = self.header.labels[i]
                    data_column = self.data.raw[:, i]
                    self.data.cols[motor_name] = data_column
                else:
                    # Otherwise, map that column label to the one data point.
                    motor_name = self.header.labels[i]
                    data_value = np.array([self.data[i]])
                    self.data.cols[motor_name] = data_value

            # Then update the local copy of the values dict.
            self.values = self.data.cols.values

            # If set_labels is true, create SpecScan object attributes from
            # all of the key-value pairs contained in the SpecScanData
            # values dict.
            if self.set_labels:
                for i in self.data.cols:
                    new_attr = i
                    if __verbose__:
                        print 'oooo Setting attribute {}'.format(i)
                    setattr(self, new_attr, self.data.cols[i])

    def concat(self, newscan):
        """Concatenates self with the SpecScan passed to the method."""
        if self.header.labels != newscan.header.labels:
            raise Exception('Scan column headers are not the same.')

        self.header += newscan.header

        self.data = np.concatenate((self.data, newscan.data))
        self.data.scan_nums = np.concatenate((self.data.scan_nums, newscan.data.scan_nums))
        self.data.row_nums = np.concatenate((self.data.row_nums, newscan.data.row_nums))

        self._setcols()

    def plot(self, **kwargs):
        """Return a SpecPlot object according to the arguments."""
        plot = SpecPlot(self, **kwargs)
        plot.doplot()

    def show(self, prefix='', nperline=6):
        """Return a string of statistics about this SpecScan instance."""
        info = ''
        info += 'Scan:\n\n'
        info += '\t{}\n\n'.format(self.header.scan_no)
        info += 'Datafile:\n\n'
        info += '\t{}\n\n'.format(self.specfile.filename)
        info += 'Scan Command:\n\n'
        info += '\t{}\n\n'.format(self.header.scan_cmd)
        info += 'Scan Constants:\n\n'

        j = nperline
        types_to_print = [float, str, np.ndarray, int, np.float]

        for key in self.__dict__:
            if not self.data.cols.has_key(key):
                if types_to_print.count(type(getattr(self, key))):
                    info += '{:19} '.format(key)
                    j -= 1
                    if j == 0:
                        info += '\n'
                        j = nperline

        info += '\n\n'
        info += self.data.show(prefix, nperline)
        return info

class SpecScanHeader(object):
    """This class contains the data from the given scan header.

    Attributes:

        fullheader : string
            A string containing the entire scan header.

        scan_no : int
            The scan number as an integer.

        scan_cmd : string
            The scan command as a string.

        scan_type : string
            The scan type as a string.

        date : string
            A string containing the time as it appears in the scan file.

        labels : list
            A list of all the column labels for the scan data.

        comments : string
            A string containing all of the comment lines from the scan.
    """

    def __init__(self):
        # Declare all of the attributes ahead of time. Default them to None.
        self.fullheader = ''
        self.scan_no = None
        self.scan_cmd = ''
        self.scan_type = ''
        self.date = ''
        self.labels = ''
        self.comments = ''
        return

    def parse_header_line(self, specscan, line):
        """Read the line and set attributes accordingly."""
        if line[0:2] == '#P':
            self.parse_motor_line(specscan, line)

        elif line[0:2] == '#C':
            self.parse_comment_line(line)

        elif line[0:2] == '#D':
            self.parse_date_line(line)

    def parse_motor_line(self, specscan, line):
        """Read #P lines and set motor positions."""
        self.fullheader += line
        counter = 0
        positions = line.strip().split()
        for i in range(1, len(positions)):
            if counter < len(specscan.specfile.motors):
                key = specscan.specfile.motors[counter]
                value = np.array([float(positions[i])])
                specscan.data.set_value(key, value)
            else:
                print '**** WARNING: More motors in scan than'
                print '****     defined in the scan header. Someone'
                print '****     changed the motor order without'
                print '****     updating the file. CHECK that the'
                print '****     assignments are correct.'
            counter += 1

    def parse_comment_line(self, line):
        """Read in the #C comment line."""
        self.fullheader += line
        self.comments += line

    def parse_date_line(self, line):
        """Read in the #D date line."""
        self.fullheader += line
        date_obj = time.strptime(line[2:].strip())
        self.date = time.asctime(date_obj)

class SpecScanData(object):
    """This class defines the data contained within a scan.

    Attributes:

        values : dict
            A dict mapping motor names to their initial positions, and column
            labels to their corresponding columns (arrays) of data.

        raw : NumPy array
            A NumPy array containing the scan data as it appears in the file.

        cols : dict
            A dict which maps motor names to their initial positions, and
            column labels to their corresponding column arrays of data.

        row_nums : NumPy array
            An nx1 array where the nth entry is the integer n, starting at 0.
            (n is the number of rows of data in the scan.)

        scan_nums : NumPy array
            An nx1 array where each entry is the scan number as an int.
            (n is the number of rows of data in the scan.)

    """

    def __init__(self):
        """Initializes an instance of the SpecScanData class."""
        self.values = {}
        self.raw = np.array([])
        self.cols = {}
        self.scan_nums = np.array([])
        self.row_nums = np.array([])

    def __str__(self):
        """Call show() when the class is passed in a print statement."""
        return self.show()

    def parse_data_line(self, line):
        """Read in a line of data and add it to the full array."""
        dataline = map(float, line.strip().split())
        newrow = np.asarray(dataline)
        if len(dataline) != 0:
            if len(self.raw) == 0:
                # If the SpecScan object has no data already:
                self.raw = newrow
            else:
                # Otherwise, add the new row below the old data.
                # Note that np.vstack is deprecated in favor of
                # np.stack and np.concatenate.
                self.raw = np.vstack((self.raw, newrow))

    def get(self, key):
        """Return the value for the given key, if possible."""
        if self.values.has_key(key):
            return self.values[key]
        else:
            return None

    def set_value(self, key, data):
        """Add the given key-value pair to the values dict."""
        self.values[key] = data
        if __verbose__:
            print 'oooo Setting key {}'.format(key)

    def show(self, prefix='', nperline=6):
        """Return a string which reads out all motors and scan variables."""
        j = nperline
        info = ''
        info += prefix + 'Motors:\n\n'
        info += prefix
        for i in self.values:
            if self.values[i].size == 1:
                info += '{:19}'.format(i)
                j -= 1
                if j == 0:
                    info += '\n' + prefix
                    j = nperline

        if j != nperline:
            info += '\n' + prefix

        info += '\n'

        info += prefix + '\n'
        info += prefix + 'Scan Variables:\n'
        info += prefix + '\n'
        j = nperline
        for i in self.values:
            if self.values[i].size > 1:
                info += '{:19}'.format(i)
                j -= 1
                if j == 0:
                    info += '\n' + prefix
                    j = nperline

        info += '\n'
        return info


class SpecPlot(object):
    """This class represents a standard x-axis vs. y-axis plot.

    Attributes:

        plot : None or pyplot plot
            The matplotlib plot associated with the object.

        plotted : bool
            A boolean indicating whether the scan has been plotted.

        plotx : NumPy array
            The x-data for the plot.

        ploty : NumPy array
            The y-data for the plot. May be normalized.

        scan : SpecScan object
            The SpecScan object from which to graph data.

    """

    def __init__(self, specscan):
        """Initialize an instance of the SpecPlot class."""
        self.scan = specscan
        self.plt = None
        self.plotted = False

        self.xcol = None
        self.x2col = None

    def doplot(self, norm=False, ycol='ChT_REIXS', mcol='I0_BD3', fmt=''):
        """Generates a plot according to the provided kwargs."""
        twod = istwod(self.scan)
        
        if twod:
            ycol = 'TEY_REIXS'
        
        self.set_x_axes()
        self.check_x_type()
        self.labels_to_indices()

        if isinstance(ycol, str):
            ycol = self.scan.header.labels.index(ycol)
        if isinstance(mcol, str):
            mcol = self.scan.header.labels.index(mcol)

        plotx = self.scan.data.raw[:, self.xcol]
        if self.x2col != None:
            plotx = np.vstack((plotx, self.scan.data.raw[:, self.x2col]))
            plotx = plotx.transpose()

        if __verbose__:
            print '**** Plotting scan {} ({})'.format(self.scan.header.scan_no,
                                                      self.scan.header.scan_cmd)
            print '---- x = {}'.format(self.scan.header.labels[self.xcol])
            if self.x2col != None:
                print '---- x2 = {}'.format(self.scan.header.labels[self.x2col])

        if norm:
            ploty = self.scan.data.raw[:, ycol]/self.scan.data.raw[:, mcol]
            y_label = '{} / {}'.format(self.scan.cols[ycol],
                                       self.scan.cols[mcol])
            if __verbose__:
                print '---- y = {} / {}'.format(self.scan.cols[ycol],
                                                self.scan.cols[mcol])
        else:
            ploty = self.scan.data.raw[:, ycol]
            y_label = self.scan.header.labels[ycol]
            if __verbose__:
                print '---- y = {}'.format(self.scan.header.labels[ycol])

        plt.figure()

        if twod:
            self.do_mesh_plot(plotx, ploty)

        else:
            self.do_std_plot(plotx, ploty, fmt, y_label)


        plt.title('{} {} {} {}\n{}'.format(self.scan.specfile.filename,
                                           self.scan.header.scan_no, ' - ',
                                           self.scan.header.date,
                                           self.scan.header.scan_cmd))
        self.plotted = True
        plt.show()

    def set_x_axes(self):
        """Set values for the x axes according to the scan type."""
        if self.xcol is None:
            if self.scan.header.scan_type.strip() == 'mesh':
                self.xcol = 0
                self.x2col = 1
            else:
                self.xcol = 0
                self.x2col = None

    def check_x_type(self):
        """Check that the x-axes have the right type."""
        if isinstance(self.xcol, list) or isinstance(self.xcol, np.ndarray):
            self.xcol = self.xcol[0]
            self.x2col = self.xcol[1]
        else:
            if not isinstance(self.xcol, int) and not isinstance(self.xcol, str):
                raise Exception("Illegal {} : xcol can only be 'int', 'list'\
                                 or 'np.ndarray'.".format(type(self.xcol)))

    def labels_to_indices(self):
        """Converts the axis labels to indices."""
        if isinstance(self.xcol, str):
            self.xcol = self.scan.cols.index(self.xcol)
        if isinstance(self.x2col, str):
            self.x2col = self.scan.cols.index(self.x2col)

    def do_mesh_plot(self, plotx, ploty):
        """Plot the mesh scan data."""
        j = 1j
        xint, yint = get_mesh_dims(self.scan)
        grid_x, grid_y = np.mgrid[min(plotx[:, 0]):
                                  max(plotx[:, 0]):xint*j,
                                  min(plotx[:, 1]):
                                  max(plotx[:, 1]):yint*j]
        if __verbose__:
            print '---- xint = {}'.format(xint)
            print '---- yint = {}'.format(yint)

        grid_z = griddata(plotx, ploty,
                          (grid_x, grid_y), method='linear')

        plt.pcolormesh(grid_x, grid_y, grid_z, cmap='inferno')
        plt.colorbar()

        plt.xlim(min(plotx[:, 0]), max(plotx[:, 0]))
        plt.ylim(min(plotx[:, 1]), max(plotx[:, 1]))

        plt.xlabel(self.scan.header.labels[self.xcol])
        plt.xlabel(self.scan.header.labels[self.x2col])

    def do_std_plot(self, plotx, ploty, fmt, y_label):
        """Plot a standard x- vs. y-axis plot."""
        self.plt = plt.plot(plotx, ploty, fmt)

        plt.xlabel(self.scan.header.labels[self.xcol])
        plt.ylabel(y_label)

        plt.xlim(min(plotx), max(plotx))

        
if __name__ == '__main__': 
    YBCO = SpecDataFile('YBCO_XAS')
    LNSCO = SpecDataFile('LNSCO')

    # Test standard plot.
    SCAN1 = LNSCO[298]
    SCAN1.plot()

    # Test mesh plot.
    SCAN2 = YBCO[11]
    SCAN2.plot()
