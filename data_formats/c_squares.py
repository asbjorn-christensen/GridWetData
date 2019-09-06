#! /usr/bin/python
# ---------------------------------------------------------------
#  @package c_squares
#
#  Standalone implementation of c-square specification
#  at
#     http://www.marine.csiro.au/csquares/ 
#     http://www.cmar.csiro.au/csquares/spec1-1.htm
#  or 
#     Rees, Tony, 2003.
#     "C-squares", a new spatial indexing system and its applicability to the description of oceanographic datasets.
#     Oceanography, vol.16(1): 11-19.
#
#  This module provides:
#
#     1) raw transformations:
#           tag2bbox
#           tag_of_enclosing_csquare
#           tag2centerpoint
# 
#     2) classes for structured handling of CSquares
#           CSquare
#           CSquareList (can be instantiated with a list, wildcards or multiple c-square tags)
#                  
#
#  Longitudes are decimal degree specification in range [-180:180]
#  Latitudes are decimal degree specification in range [-90:90]
#  A bounding box is a lon-lat specification  (W,S,E,N) = (lonmin, latmin, lonmax, latmax) 
#  where (lonmin < lonmax) and (latmin < latmax); this corresponds to (SW corner, NE corner) of
#  the bounding box. Provide some minimum level of input error checking
#
#  Examples:
#
#     (W,S,E,N), resolution = tag2bbox("5704:143:373:225")         # probe raw bounding box corresponding to c-square tag
#     tag_of_enclosing_csquare(0, 55, 0.1)                         # generate raw tag of c-square surrounding (lon=0,lat=55) with resolution 0.1 deg
#
#     csqr = CSquare(0, 55, 0.1)                                   # create csquare surrounding (lon=0,lat=55) with resolution 0.1 deg
#     csqr = CSquare("7500:100:100")                               # create csquare corresponding to tag
#     print csqr.area()                                            # assess area of a csquare
#     print csqr.bbx                                               # assess bounding box (W,S,E,N) of a csquare
#     csqr.my_stuff = 17                                           # attach some custom data to a csquare
#
#     csqrlist = CSquareList("5704:143:373:***|7016:466:469:131")  # create a list of csquares by parsing input string
#
# ---------------------------------------------------------------

import exceptions
from math import sin, pi     # for CSquare.area
earth_radius  = 6371.0       # Earth average radius in km
deg2rad       = pi/180.

def nint(x):
    ## Round float to nearest integer
    #  @param x : float
    #  @return  : float rounded to nearest integer
    return int(round(x))


def tag2bbox(csquare):
    ## Identify bounding box of primitive csquare tag at any resolution
    #  @param csquare : c-square tag, e.g. "5704:143:373:225"
    #  @return  : bounding box (W,S,E,N), resolution[deg] 
    # ===========================================================================
    # csquare = WMOsqr[:subsqr[:subsqr[...]]]
    #
    # return bounding box of csquare and resolution [deg]
    # (resolution is implied by bounding box, return it as a convenience)
    # 
    # A primitive csquare tag is a designation identifying a single csquare
    # enforce that intermediate sub squares (after WMOsqr) are a full (3 digits)
    # only last sub square may be half (one digit, 1-4)
    # ignore leading whitespaces (if any) on WMOP tag
    # ===========================================================================
    #
    # --- check tag validity: deny tags containing "*" or "|"
    #
    if csquare.find("*") > -1: # not found => -1
        raise exceptions.ValueError("%s: wildcard designation not supported here" % csquare)
    if csquare.find("|") > -1: # not found => -1
        raise exceptions.ValueError("%s: multiple csquare parsing not supported here" % csquare)
    tags = csquare.split(":")
    #
    # --- identify bounding box corner (lon0,lat0) closest to (lon,lat) = (0,0)
    #
    #     1) parse the WMO tag (10 deg resolution)
    #
    WMOtag = tags[0].lstrip() # ignore leading whitespaces (if any) on WMOP tag
    if len(WMOtag) != 4:
        raise exceptions.ValueError("%s: invalid WMO tag %s encountered" % (csquare, WMOtag))
    if not WMOtag[0] in ("1", "3", "5", "7"):
        raise exceptions.ValueError("%s: invalid WMO tag %s encountered" % (csquare, WMOtag))
    # identify hemisphere
    if WMOtag[0] in ("1", "3"):
        lonsign =  1  # Eastern hemisphere
    else:
        lonsign = -1  # Western hemisphere
    if WMOtag[0] in ("1", "7"):
        latsign =  1  # Northern hemisphere
    else:
        latsign = -1  # Southern hemisphere
    # assert valid range of lon/lat fields
    if int(WMOtag[1])>8:
        raise exceptions.ValueError("%s: WMO tag %s exceed latmax" % (csquare, WMOtag))
    if int(WMOtag[2])>17:
        raise exceptions.ValueError("%s: WMO tag %s exceed lonmax" % (csquare, WMOtag))
    #
    lat0  = float(WMOtag[1])   * 10*latsign
    lon0  = float(WMOtag[2:4]) * 10*lonsign
    resol = 10. # degrees
    #
    #     2) parse sub square tags (if any)
    #
    for (itg,subsqr) in enumerate(tags[1:]): # loop possibly void
        if len(subsqr) == 1: # half division
            if (itg < len(tags)-2): # WMO tag already parsed
                raise exceptions.ValueError("%s: only last subsquare can be half-divided" % csquare)
            if not subsqr in ("1", "2", "3", "4"):
                raise exceptions.ValueError("%s: invalid half square %s encountered" % (csquare,subsqr))
            resol *= 0.5
            if subsqr in ("2", "4"):
                lon0 += resol*lonsign
            if subsqr in ("3", "4"):
                lat0 += resol*latsign
        elif len(subsqr) == 3: # deca division (ignore subsqr[0] )
            resol *= 0.1
            lat0  += resol*float(subsqr[1])*latsign
            lon0  += resol*float(subsqr[2])*lonsign
        else: # invalid subsquare length
            raise exceptions.ValueError("%s: subsquare %s must have len 1|3" % (csquare,subsqr))
    #
    #    3) construct and normalize bounding box
    #
    lon1 = lon0 + resol*lonsign
    lat1 = lat0 + resol*latsign
    bbx  = min(lon0,lon1), min(lat0,lat1), max(lon0,lon1), max(lat0,lat1) # standard order
    return bbx, resol       


def tag2centerpoint(csquare):
    ## Identify center of primitive csquare tag at any resolution
    #  @param csquare : c-square tag, e.g. "5704:143:373:225"
    #  @return  :       center of c-square (lonc[deg],latc[deg]), resolution[deg] 
    # ===========================================================================
    # csquare = WMOsqr[:subsqr[:subsqr[...]]]
    # ===========================================================================
    bbx, resol = tag2bbox(csquare)
    lonc       = 0.5*(bbx[0]+bbx[2])
    latc       = 0.5*(bbx[1]+bbx[3])
    return (lonc,latc), resol
    

def tag_of_enclosing_csquare(lon, lat, resol):
    ## Identify enclosing csquare tag for point (lon, lat) at resolution resol degrees
    #  @param lon   : longitude[deg] of target point
    #  @param lat   : latitude[deg] of target point
    #  @param resol : resolution[deg] desired c-square
    #  @return      : c-square tag corresponding to input
    # ===============================================================================
    # resol must be a standard step: 10, 5, 1, 0.5, 0.1, ....
    # handle boundary cases
    # ===============================================================================
    if not -180 < lon < 180:
        raise exceptions.ValueError("range error: lon = %f " % lon)
    if not -90  < lat < 90:
        raise exceptions.ValueError("range error: lat = %f " % lat)
    #
    # identify global WMO quadrant and tag
    #
    if lon>0:
        if lat>0:
            wmo_quadrant = 1
        else:
            wmo_quadrant = 3
    else:
        if lat>0:
            wmo_quadrant = 7
        else:
            wmo_quadrant = 5
    # --- fold to NE quadrant
    lon  = abs(lon)/10
    lat  = abs(lat)/10
    ilon = int(lon)
    ilat = int(lat)
    this_resol = 10.0 # avoid integer division below
    tag  = "%1d%1d%02d" % (wmo_quadrant, ilat, ilon) # include leading zeros from lon designator
    #
    # csquare quadrant labelling is different from WMO
    #
    while this_resol/resol > 1.001: 
        lon = 10*(lon-ilon) # reap next lon digit
        lat = 10*(lat-ilat) # reap next lat digit
        csq_quadrant = 1
        if lon>5:
                csq_quadrant += 1
        if lat>5:
                csq_quadrant += 2
        ilon = int(lon)
        ilat = int(lat)
        tag  = tag + ":%1d%1d%1d" % (csq_quadrant, ilat, ilon)
        this_resol /= 10. # corresponding to current tag
    
    if this_resol/resol < 0.201: # last step was half division
        tag = tag[:-2]   # drop last two digits to back up half step
    return tag



# ===========================================================================
#
#                       Represent a CSquare
#
#   attributes:
#      id         : CSquare tag
#      bbx        : (W,S,E,N) of C-Square box  (strictly implied by id)
#      resolution : resolution of this C-Square (strictly implied by bbx)
#
# ===========================================================================
class CSquare:
    ## Represent a CSquare
    def __init__(self, *args):
        ## constructor for class CSquare:
        #  @param args   : eiter a primitive C-Square tag (precedence) or a (lon[deg],lat[deg],resolution[deg]) of an enclosed point
        # ===========================================================================
        # args is eiter a primitive C-Square tag (precedence) or a (lon,lat,resolution) tuple
        # corresponding to a point enclosed by the desired C-Square
        # ===========================================================================
        if isinstance(args[0], basestring):
            if len(args)>1:
                raise exceptions.ValueError("CSquare instantiation: expect no args beyond tag" )
            self.bbx, self.resolution = tag2bbox(args[0])
            self.id = args[0]
        else:
            if len(args) != 3:
                raise exceptions.ValueError("CSquare instantiation: expect args (lon,lat,resolution)")
            self.id = tag_of_enclosing_csquare(*args)
            self.bbx, self.resolution = tag2bbox(self.id)
            
    def center(self):
        ## return center position of C-Square
        #  @return (lon[deg],lat[deg]) of center position of C-Square
        return 0.5*(self.bbx[0]+self.bbx[2]), 0.5*(self.bbx[1]+self.bbx[3])

    def area(self):
        ## compute C-Square area in km^2
        # ===========================================================================
        # When bounding box (W,S,E,N) is in radians
        #    area = Integrate(S < phi < N, W < lamda < E, R^2*cos(phi))
        #         = R^2*(E-W)*(sin(N)-sin(S))
        #         < 4 pi * R^2
        # ===========================================================================
        (W,S,E,N) = self.bbx
        return  earth_radius**2 * deg2rad*(E-W) * (sin(N*deg2rad) - sin(S*deg2rad))
    
    def __eq__(self, other):
        ## compare CSquare instances (for list operations in CSquareList)
        #  @param other : a CSquare instance
        # A CSquare is uniquely determined by its primitive c-square tag
        return self.id == other.id

# ===========================================================================
#
#             Represent ordered list of CSquare instances
#   
# ===========================================================================

 
class CSquareList(list):
    ## Represent ordered list of CSquare instances
    def __init__(self, input):
        ## constructor for class CSquareList
        #  @param input : designation[ | designation [ | designation [ ... ]]] (see below)
        # ===========================================================================
        # input is either a parsable string or list of CSquares
        # If string, instantiate list of CSquare corresponding to string pattern
        # input string may contain wildcards ("*" or "***") and reference multiple 
        # c-squares (separated by "|") as described in format specification
        #
        # Parse input string as follows: 
        #   input       = designation[ | designation [ | designation [ ... ]]]
        #   designation = tag | tag\* | tag\*\*\*
        #
        # List order is left-to-right, with wildcards expanded in reproducable order
        # set by loops below. Do not check for duplicates nor resolution consistency
        # ===========================================================================
        if isinstance(input, basestring):
            for designation in input.split("|"):
                if   designation[-4:] == ":***": # all subsquares of level above
                    basetag = designation[:-4]
                    for ix in range(10):
                        for iy in range(10):
                            csq_quadrant = 1
                            if ix>4:
                                csq_quadrant += 1
                            if iy>4:
                                csq_quadrant += 2
                            tag = basetag + (":%d%d%d" % (csq_quadrant,iy,ix))
                            self.append(CSquare(tag))
                    
                elif designation[-2:] == ":*":   # append corresponding 4 half-squares
                    basetag = designation[:-2]
                    for csq_quadrant in range(1,5):
                        tag = basetag + (":%d" % csq_quadrant)
                        self.append(CSquare(tag))
                        
                else: # assume it is a primitive tag
                    self.append(CSquare(designation))
        else:
            if not hasattr(input, "__len__"):
                raise exceptions.ValueError("CSquareList instantiation: input not a sequence")
            
            for i,member in enumerate(input): # validate all members in list
                if not isinstance(member, CSquare):
                    raise exceptions.ValueError("CSquareList instantiation: element %d is not CSquare instance" % i)
                self.append(member)
                
    def bounding_box(self):
        ## Identify Overall minimal bounding box for collection of CSquares
        #  @return : bounding box as (W,S,E,N) in degrees
        # ===========================================================================
        # Overall minimal bounding box for collection of CSquares
        # (possibly at different resolution)
        # ===========================================================================
        if len(self) == 0:
            return None
        else:
            (W,S,E,N) = self[0].bbx
            for member in self[1:]: # possibly void
                 W = min(W, member.bbx[0])
                 S = min(S, member.bbx[1])
                 E = max(E, member.bbx[2])
                 N = max(N, member.bbx[3])
            return (W,S,E,N)

        
    def minimal_grid(self):
        ## Identify minimal common regular grid of a collection of CSquares
        #  @return : nx, ny, lonmin, latmin, resolution (see below)
        # ===========================================================================
        # Explore, if collection of CSquares can be cast onto a common regular grid
        # requires all members are on same resolution
        #
        # return (nx, ny, lonmin, latmin, resolution) if possible
        # (lonmin, latmin) are center of the SW corner cell;
        # (nx, ny) are grid points along longitude and latitude axes, respectively
        # resolution is the common longitude/latitude grid step
        #
        # return None, if list is empty or members have different resolution
        # ===========================================================================
        if len(self) == 0:
            return None
        else:
            res = self[0].resolution
            lonmin,latmin = self[0].center()
            lonmax,latmax = lonmin,latmin
            for member in self[1:]: # possibly void
                if res != member.resolution:
                    return None 
                else:    
                    lon,lat = member.center()
                    lonmin  = min(lon,lonmin)
                    latmin  = min(lat,latmin)
                    lonmax  = max(lon,lonmax)
                    latmax  = max(lat,latmax)
            #
            nx = nint(1 + (lonmax-lonmin)/res)
            ny = nint(1 + (latmax-latmin)/res)
            return nx, ny, lonmin, latmin, res
            
    def summary(self):
        ## generate a printable summary of collection of CSquares
        #  @return : printable summary
        # ===========================================================================
        # Mainly for debugging
        # ===========================================================================
        str = "CSquareList instance with %d members:\n" % len(self)
        for (i,member) in enumerate(self):
            str = str + ("CSquare #%d = %s\n" % (i, member.id))
        return str
        
###################################################################################33
#   simple module self test
###################################################################################
if __name__=="__main__":
    #
    print 20*">", "testing c_squares.py ", 20*"<"
    
    #
    # SE global quadrant test case from http://www.cmar.csiro.au/csquares/c-squares-encoder.xls
    #
    #tag                 N limit	S limit	W limit	E limit	 Centre lat long
    #3414	            -40	    -50	    140	    150	    -45	    145
    #3414:2	            -40	    -45	    145	    150	    -42.5	147.5
    #3414:227	        -42	    -43	    147	    148	    -42.5	147.5
    #3414:227:3	        -42.5	-43	    147	    147.5	-42.75	147.25
    #3414:227:382	    -42.8	-42.9	147.2	147.3	-42.85	147.25
    #3414:227:382:4	    -42.85	-42.9	147.25	147.3	-42.875	147.275
    #3414:227:382:458	-42.85	-42.86	147.28	147.29	-42.855	147.285
    
    print "\n---- %s ----" % "testing tag2bbox(1):"
    for tag in ("3414", "3414:2", "3414:227", "3414:227:3", "3414:227:382", "3414:227:382:4", "3414:227:382:458"):
        bbx,resol = tag2bbox(tag)
        print "%20s : " % tag, (4*"%9.3f") % bbx, "resol = %f deg" % resol
        
    #
    # test cases from http://www.marine.csiro.au/marq/csq_builder.init
    #
    # NE global quadrant: 1312:227:382:458  -> (W,S,E,N) = 127.28  32.85  127.29  32.86
    # SW global quadrant: 5704:143:373:225  -> (W,S,E,N) = -43.36 -74.73  -43.35 -74.72
    # NW global quadrant: 7016:466:469:131  -> (W,S,E,N) =-166.92   6.63 -166.91   6.64
    #
    
    print "\n---- %s ----" % "testing tag2bbox(2):"
    for tag in ("1312:227:382:458", "5704:143:373:225", "7016:466:469:131"):
        bbx,resol = tag2bbox(tag)
        print "%20s : " % tag, (4*"%9.3f") % bbx, "resol = %f deg" % resol

    #
    # test cases from http://www.marine.csiro.au/marq/csq_builder.init
    #
    # SE global quadrant: ( 147.285 -42.855, resol=0.01) -> 3414:227:382:458
    # NE global quadrant: ( 127.285  32.855, resol=0.01) -> 1312:227:382:458  
    # SW global quadrant: ( -43.355 -74.725, resol=0.01) -> 5704:143:373:225 
    # NW global quadrant: (-166.915   6.635, resol=0.01) -> 7016:466:469:131
    #
    print "\n---- %s ----" % "testing tag_of_enclosing_csquare:"
    for (lon, lat) in ((147.285, - 42.855),
                       (127.285,   32.855),
                       ( -43.355, -74.725),
                       (-166.915,   6.635)):
        for resol in (10, 5, 1, 0.5, 0.1, 0.05, 0.01):
            print "lon=%8.3f lat=%8.3f resol=%6.2f - tag=%s" % (lon, lat, resol, tag_of_enclosing_csquare(lon, lat, resol))
    
    print "\n---- %s ----" % "testing : CSquare.area"
    csqr = CSquare(0, 55, 0.1)
    print "area of %s = %f km2" % (csqr.id, csqr.area())
    csqr = CSquare("7500:100:100") # alternative instantiation
    print "area of %s = %f km2" % (csqr.id, csqr.area())
    
    print "\n---- %s ----" % "testing : CSquareList"
    print CSquareList("5704:143:373:225").summary()
    print CSquareList("5704:143:373:225|7016:466:469:131").summary()
    print CSquareList("5704:143:373:*|7016:466:469:131").summary()
    print CSquareList("5704:143:373:***|7016:466:469:131").summary()
