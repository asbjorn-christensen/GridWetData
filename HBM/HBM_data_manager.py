#!/usr/bin/env python
# -*- coding: utf-8 -*-  
##########################################################################################
#        @package data_manager
#        Light weight data manager
#
#        Manage an inventory list of a provided data folder with arbitrary content and
#        service availability requests. DataManager does not know about specific
#        grid/griddata classes, but parses folder content based on structured file names.
#        File names for needed data can be retrieved from the a DataManager instance
#
##########################################################################################
import sys
import os
import exceptions
from numpy    import *  # searchsorted 
from datetime import *
import re

import copy                     # for grid copy, notice numpy also provides a copy which is shadowed by this import

_thisdir         = os.path.dirname(__file__)  # allow remote import
_verbose         = False # True  # log info
_toffset         = datetime(2000,1,1)

# ----------------------------------------------------------------------------------
# In a generalization, these regexes could be provided at DataManager instantiation
# ----------------------------------------------------------------------------------
#
#            PROPTAG_GRIDTAG_YEAR_MONTH_DAY_HOUR_MINUTE_SECOND.nc
#
data_file_template  = re.compile(r"""
(?P<proptag> \S+    ) _
(?P<gridtag> \w+    ) _
(?P<year>    \d{4,4}) _
(?P<month>   \d{2,2}) _
(?P<day>     \d{2,2}) _
(?P<hour>    \d{2,2}) _
(?P<minute>  \d{2,2}) _
(?P<second>  \d{2,2}).nc
 """, re.VERBOSE)
#
#            GRIDTAG_grid.nc
#
grid_file_template  = re.compile(r"""
(?P<gridtag> \w+    ) _
grid.nc
 """, re.VERBOSE)

## create a time hash of a datetime instance for easier indexing/searching
#  @param t  datetime object
def _create_time_hash(t):
    tdiff = t - _toffset
    return tdiff.total_seconds()

#  ==================================================================================
## Manage an inventory list of a provided data folder with arbitrary content and service availability requests
#
# DataManager does not know about grid/griddata classes desired
# for representation, but does just represent source files for grid/griddata
# grid/griddata are associated with matching GRIDTAGs. GRIDTAGs represent 
# a combination of spatial domain and grid type
#
class DataManager:
    #  -----------------------------------------------------
    ## constructor
    #  scan for grid/griddata files using regular expressions defined in module for grid/data files, respectively 
    #  consider all files in path, but discard files not mathcing satisfactory to grid/data regular expressions pattern
    #  For easy access, create a hast table of time since _toffset
    #  @param self    The object pointer
    #  @param path    path to directory where data is
    #  @param gclass  grid class through which grids should be viewed
    #
    def __init__(self, path, gclass=str):
        ##  path to directory
        self.path     = path # @var path to directory
        ##  datasets [PROPTAG] [GRIDTAG] = (times_of_data, files, datetime) where time_of_data is a time hash corresponding to datetime
        #            In the plain lists (times_of_data, files), entries correspond pairwise           
        self.datasets = {}
        ## grids [GRIDTAG] = grid descriptors corresponding to GRIDTAG, instantiated by provided instantiator gclass
        self.grids    = {}
        #
        # -------- scan for grids     --------
        #
        for fname in os.listdir( path ):
            amatch = grid_file_template.match(fname)
            #
            # --- extract set info from filename
            #
            if amatch is not None:
                grid   = amatch.group("gridtag")
                self.grids[grid] = gclass(os.path.join(path, fname)) # archive instance
            else:
                continue # invalid name, silently discard this entry
        #   
        # -------- scan for grid data --------
        #
        for fname in os.listdir( path ): # os.listdir skips head of path
            amatch = data_file_template.match(fname)
            #
            # --- extract set info from filename
            #
            if amatch is not None:
                prop   = amatch.group("proptag")
                grid   = amatch.group("gridtag")
                year   = int(amatch.group("year"))
                month  = int(amatch.group("month"))
                day    = int(amatch.group("day"))
                hour   = int(amatch.group("hour"))
                minute = int(amatch.group("minute"))
                second = int(amatch.group("second"))
            else:
                continue # invalid name, silently discard this entry
            #
            # --- archive reference
            #
            if not self.datasets.has_key(prop):
                self.datasets[prop] = {}
            if not self.datasets[prop].has_key(grid):
                self.datasets[prop][grid] = ([],[],[])
            try:
                thash = _create_time_hash( datetime(year, month, day, hour, minute, second) )
            except exceptions.ValueError:
                continue # silently discard entry
            #
            # --- insert entry for this (prop,grid) at right place ---
            #
            ix = searchsorted( self.datasets[prop][grid][0], thash )
            self.datasets[prop][grid][0].insert(ix, thash)
            self.datasets[prop][grid][1].insert(ix, os.path.join(path, fname))  # store full path
            self.datasets[prop][grid][2].insert(ix, datetime(year, month, day, hour, minute, second))
        #

    # ------------------------------------------------------------------------------------
    ## Resolve filenames suitable for (linear) time interpolation at attime of prop in grid
    # Return (file0, zfile0, w0), (file1, zfile1, w1) ordered according to time
    # or raise LookupError where (w0, w1) are relative weights to be applied for linear time
    # (zfile0, zfile1) are corresponding sea level elevation that
    # is needed to update grids corresponding to a given time.
    # sea level elevation is given by property z
    #  @param self    The object pointer
    #  @param prop    property tag of solicited data
    #  @param grid    grid tag (spatial grid + resolution) of solicited data
    #  @param attime  datetime instance corresponding to interpolation time
    #
    def get_time_interpolation_files(self, prop, grid, attime):
        
        (file0,  tprop0), (file1,  tprop1) = self.get_interpolation_bracket(prop, grid, attime)
        (zfile0, tz0),    (zfile1, tz1)    = self.get_interpolation_bracket("z",  grid, attime) # sea level elevation
        if (tprop0 != tz0) or (tprop1 != tz1):
            raise exceptions.LookupError("(prop,z) frames does not match time wise: %d vs %d ; %d vs %d" % (tz0,tprop0,tz1,tprop1))
        # --- compute relative weigths for left/right bracket ---
        dt    = tprop1 - tprop0
        tseek = _create_time_hash(attime)
        w0     = (tprop1-tseek)/dt
        w1     = 1.0 - w0 
        return (file0, zfile0, w0),  (file1, zfile1, w1)
        
        
    # ---------------------------------------------------------------
    ## Resolve closest pair of filenames suitable for (linear) time interpolation
    # of prop in grid. Return (file0, thash0), (file1, thash1) ordered according to time
    # or raise LookupError
    #  @param self    The object pointer
    #  @param prop    property tag of solicited data
    #  @param grid    grid tag (spatial grid + resolution) of solicited data
    #  @param attime  datetime instance corresponding to interpolation time
    #
    def get_interpolation_bracket(self, prop, grid, attime):
        tseek = _create_time_hash(attime)
        try:
            i1  = searchsorted( self.datasets[prop][grid][0], tseek )
        except:
            print tseek
            print self.datasets[prop][grid][0]
            raise exceptions.LookupError("searchsorted failed for %s in %s for %s" % (prop, grid, `attime`))
        nt   = len(self.datasets[prop][grid][1])
        if nt<2:
            raise exceptions.LookupError("too few fields for %s in %s to time interpolate" % (prop, grid))
        if i1 >= nt:
            raise exceptions.ValueError("i1 >= nt", attime)
        i0 = i1-1
        if i0<0:
            raise exceptions.ValueError("i0<0", attime)
        # --- we've identified the proper time bracket ---
        ths    = self.datasets[prop][grid][0]
        fnames = self.datasets[prop][grid][1]
        return (fnames[i0], ths[i0]),  (fnames[i1], ths[i1])
        #

    # ---------------------------------------------------------------
    ## Generator function for looping over data sets in inventory
    #  Staggering not supported (ignored)
    #  Return (time, griddata) pair where time is a date_time instance and griddata is a griddata_class instance
    #  @param self              The object pointer
    #  @param griddata_class    spectator class to be applied for grid data - must be instantiated as a GridData_3D
    #  @param proptag           tag for data property to be looped over
    #  @param gridtag           tag for grid to be applied 
    def loop_over(self, griddata_class, proptag, gridtag):
        if _verbose: print "load_frame: %s" % fname_data
        for (thash, fname, dtime) in zip(*self.datasets[proptag][gridtag]):    
            izf     = self.datasets["z"][gridtag][2].index(dtime) # pick corresponding sea level elevation
            zfname  = self.datasets["z"][gridtag][1][izf]         # 
            newgrid = copy.copy(self.grids[gridtag])              # shallow copy of ancestor grid
            newgrid.update_sea_level(zfname)                
            yield dtime, griddata_class(newgrid, fname)
        
        
    # ---------------------------------------------------------------
    ## Print inventory lists of grids/griddata sets available to stdout
    #  mainly for debugging
    def print_inventory(self):
        print "Inventory of %s :" % self.path
        print
        print "Grids:"
        print "%5s   %s" % ("grid", "string rep")
        for grid in self.grids.keys():
            print "%5s   %s" % (grid, str(self.grids[grid]))
        print
        print "Data sets:"
        print "%10s %5s   %14s     %s" % ("property", "grid", "time hash", "file name")
        for prop in self.datasets.keys():
            for grid in self.datasets[prop].keys():
                for (thash, fname) in zip(self.datasets[prop][grid][0], self.datasets[prop][grid][1]):
                    print "%10s %5s   %14.2f     %s" % (prop, grid, thash, fname)

        

#################### self test #########################
if __name__ == "__main__":
    _verbose = True
    dmg = DataManager("../DMI_data")
    dmg.print_inventory()
    print
    now = datetime(2015, 05, 20, 00, 20, 00) 
    print "file mapping and interpolation weights for ", now
    intp_set = dmg.get_time_interpolation_files("t", "ns", datetime(2015, 05, 20, 00, 20, 00))
    print "left :", intp_set[0]
    print "rigth:", intp_set[1]
