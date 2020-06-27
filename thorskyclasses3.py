#!/usr/bin/env python
"""thorskyclasses.py -- classes used for observational
   circumstances.  
"""

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, EarthLocation, CartesianRepresentation
from astropy.coordinates import FK5, PrecessedGeocentric
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric
from astropy.coordinates import get_body, get_moon, get_sun
from astropy.time import Time, TimeDelta
from datetime import datetime
from pytz import timezone
import pytz

import time as ttime  
import dateutil.parser

from thorsky.thorskyutil import altazparang,lpsidereal,lpmoon,lpsun,accumoon,min_max_alt,phase_descr
from thorsky.thorskyutil import ha_alt, jd_sun_alt, jd_moon_alt, local_midnight_Time, hrs_up
from thorsky.thorskyutil import ztwilight, lunskybright, true_airmass
from thorsky.thorskyutil import currentgeocentframe, currentFK5frame, thorconsts
from thorsky.thorskyutil import precessmatrix, cel2localmatrix, cel2xyz
from thorsky.thorskyutil import getbrightest, skyproject, angle_between
from thorsky.thorskyutil import time_rounded_to_minute

# import matplotlib.pyplot as plt   # for tests
import sys  # for tests

import pkgutil   # to get fixed data

# import warnings  # to filter ERFA warnings about dubious years 
# Don't have this working.

# this is slightly modified cribbed code.
   
#def fxn():
#    warnings.warn("dubious", ErfaWarning)

#with warnings.catch_warnings():
#    warnings.simplefilter("ignore")
#    fxn()

class obs_site :  # defaults to MDM.  Arguments are all keyworded.
    """
    Site from which observations are taken.

    This initializes by default to Kitt Peak.  In practice it is 
    usually specified by using an abbreviation from the list of observatories
    distributed with the package, 'observatories_rev.dat'.  
    Other arguments are keyworded.

    Notes
    -----

    The core of this is an EarthLocation, but it also includes time zone 
    information and a guess at the height of the observatory above surrounding
    terrain for adjusting rise and set times for the sun and moon.

    """
    
    def __init__(self, name = "MDM Observatory [Kitt Peak]", loc = (-1996199., -5037542., 3356753.),
                 tzstr = "America/Phoenix", tzstdabbr = "MST", tzdayabbr = "MDT",
                 westheight = 700., eastheight = 0.) :
        """Initialize an obs_site. 

        Parameters
        ----------
        name : str, optional
             The full name of the observatory.  Defaults to "MDM Observatory [Kitt Peak]"

        loc : tuple or list 
             Geocentric cartesian location, in meters.  Default (-1996199., -5037542., 3356753.)
             If you need to compute these for a new site, see the documentation for 
             astropy.coordinates.EarthLocation, in particular the "from_geodetic" method.

        tzstr : str
             time zone name in format used by unix programs; default "America/Phoenix".

        tzstdabbr, tzdayabbr : str
             abbreviations for standard and daylight time in this zone; defaults "MST","MDT"

        westheight, eastheight :
             rough height of observatory above surrounding terrain, in meters: defaults 700. and 0.
        
        Attributes: 
        
        name : str
             Name of observatory.

        location : EarthLocation
             EarthLocation of the observation site.

        localtz : timezone
             Local time zone.

        localtzabbrev, localtzdayabbrev :
             Abbreviations for standard and daylight time.

        risealt, setalt : Angle
             elevations at which the sun/moon rises or sets.

        axisdist : Quantity in distance units
             distance of site from earth's rotation axis

        diurnalspeed : Quantity in velocity units
             speed of diurnal rotation at this site, for use in velocity corrections.
              
        """

        self.name = name

        if isinstance(loc, EarthLocation) : 
            self.location = loc
        elif isinstance(loc,tuple) or isinstance(loc, list) :  # geocentric xyz
            self.location = EarthLocation.from_geocentric(loc[0], loc[1], loc[2], 
              unit = u.m)  # MDM

        self.tzstr = tzstr
        self.localtz = timezone(tzstr)
        self.localtzabbrev = tzstdabbr   # abbrevn for standard time zone
        self.localtzdayabbrev = tzdayabbr  # abbrevn for daylight time zone

        self.height_above_west = westheight * u.m    # height of obs over its west horizon
        self.height_above_east = eastheight * u.m    # height of obs over its east horizon

        # compute the altitude at which sunrise and sunset occur, defined as the
        # top limb of the sun or moon coinciding with the horizon.  In the
        # absence of terrain effects, this is conventionally at zenith dist
        # 90 degrees 50 minutes.  In real life higher precision is not useful
        # because atmospheric refraction varies a great deal at such low
        # altitudes.  Include in the calculation the possibility that the 
        # observatory is above or below the terrain that forms its horizon,
        # using a simple approximation for the depression of the horizon.

        # These altitudes also serve for the moon.  The angular diameter
        # variation is insignificant for this purpose, again because of 
        # atmospheric refraction.

        self.risealt = None
        self.setalt = None

        west_depression = Angle(np.sqrt(2. * ( abs(self.height_above_west)
                            / (thorconsts.EQUAT_RAD))) * u.rad)
#        print "west_depression",west_depression

        if self.height_above_west > 0. :
           self.setalt = Angle(-0.833,unit=u.deg) - west_depression
                              # zd = 90 deg 50 arcmin
        elif self.height_above_west <= 0. :
           self.setalt = Angle(-0.833,unit=u.deg) + west_depression

        east_depression = Angle(np.sqrt(2. * (abs(self.height_above_east)
                            / (thorconsts.EQUAT_RAD))) * u.rad)
#        print "east_depression",east_depression

        if self.height_above_east > 0. :
           self.risealt = Angle(-0.833 * u.deg) - east_depression
                              # zd = 90 deg 50 arcmin
        elif self.height_above_east <= 0. :
           self.risealt = Angle(-0.833 * u.deg) + east_depression
#        print "self.risealt = ",self.risealt

# compute and store the speed of diurnal rotation at this site
# for later use in barycentric velocity correction.

        axisdist = np.sqrt(self.location.x ** 2 + self.location.y ** 2)
        axisdist = axisdist.to(u.km)
        # one sidereal day is about 86164.09 seconds.
        self.diurnalspeed = (2. * np.pi * axisdist) / (86164.09 * u.s) 

def get_sites() :  # reads from a file of site names.
    """Creates a dictionary of obs_site instances from the file of site parameters
    provided in the package.  Runs at import time and will probably never need to be
    called by a user

    """
#    inf = open("observatories_rev.dat")
    # use pkgutil since I'm packaging this ... 
    #  https://stackoverflow.com/questions/6028000/how-to-read-a-static-file-from-inside-a-python-package
    # print('package name:',__name__)   # was __package__, which broke in Python 3.
    bytestr = pkgutil.get_data(__name__, 'observatories_rev.dat') 
    # this was just a file bu tis now a byte string containing the raw contents of the 
    # file.  Turn it into ASCII and cut into lines, and it'll read as if it were a file.
    inf = bytestr.decode('ASCII').splitlines()
    sitedict = {}
    for l in inf :
        # print("l = ",l)  # diagnostic -- this took some doing!
        if l[0] != '#' :  # allow commenting out of lines
            try :
                x = l.split('\t')
                y = x[2].split(",")  # xyz coordinates in meters as characters
                loc = (float(y[0]),float(y[1]),float(y[2]))
                sitedict[x[0]] = obs_site(name = x[1], loc = loc, tzstr = x[3], tzstdabbr = x[4],
                    tzdayabbr = x[5], eastheight = float(x[6]), westheight = float(x[7]))
            except :
                print("Bad input line in get_sites, may not be tab-separated:")
                print(l)
    return sitedict

sitedict = get_sites()   # reads site dictionary from package data file 'observatories_rev.dat'

class Observation :
    """Circumstances of an observation at a given SkyCoord, obs_site, and time.

    This is the main class, used to compute time-and-the-sky information.

    A minimal use example ::
       
       import astropy.units as u
       import thorsky.thorskyutil as tsu
       import thorsky.thorskyclasses3 as tsc3

       o = tsc3.Observation()      # defaults to zenith at default site at present moment
  
       o.computesky()              # hour angle, airmass, and more
       o.computesunmoon()          # note that the order in which these are called is important
       o.computequickbary()
       o.setnightevents()
       o.computeplanets()

       o.printnow()                # see what you have
       o.printnight() 


    This class has a huge number of attributes.  To avoid unneeded computation,
    these are mostly left blank unless computed explicitly using the methods 
    provided.

    The order in which you call the methods that compute results is important
    since results are cumulative:  In sum:

    * computesky needs to be first
    * computesunmoon, if needed, is second.
    * computequickbary depends on computesky but nothing else
    * setnightevents depends on both computesky and computesunmoon
    * hours_up depends on computesky, computesunmoon, and setnightevents
    * computeplanets depends on computsky and computesunmoon.

    These have been split out to save compute time.  For example, if you only 
    need airmass, all you need to run is computesky, but if you want to know 
    if the moon is up you need to run computesunmoon as well.
    
    """

    def __init__(self,celest = None, t = None, site = None, default_site = 'mdm', 
          use_local_time = True, ra_unit = u.hourangle, dec_unit = u.deg) :
        """
        
        Parameters :

        celest : SkyCoord, or string, or list, or tuple; optional
            Celestial location of the observation.  Defaults to the zenith for the site and time.
            A wide variety of input formats are permissible; see 'setcelest' below.

        t : astropy Time, or string, or tuple, or float; optional
            Instant at which observation occurs.  Default is time at which code is executed.
            Also accepts a wide variety of inputs.

        site : obs_site, or str; optional
            Site on earth from which observations are taken.  Usually passed as a key from
            the builtin dictionary of sites, but can be a obs_site instance.
            Default is 'mdm'.

        use_local_time : boolean; optional
            Report results in local time; alternative is UT.  Default is True.

        ra_unit : astropy.unit; optional 
            unit for ra input; default is u.hourangle, which I wish they'd named something 
            different since it's a term of art.   

        dec_unit : astropy.unit; optional
            unit for dec input; default is u.deg

        Attributes :
        
        site : obs_site 
            The observer's site.

        celest : astropy SkyCoord 
            celestial position observed.

        t : astropy Time 
            instant in time of the observation.

        julyear : float 
            julian year corresponding to t.

        lst : astropy Angle 
            local sidereal time, computed to about 1 sec bypassing slow astropy routines.

        nowfr : astropy coordinate frame 
            FK5 frame for equinox of obsevation ('of date')

        hanow : astropy Angle 
            hour angle (from meridian)
     
        airmass : float 
            true airmass (not secant-z)

        cel_J2000 : astropy SkyCoord  
            celestial location for J2000.

        icrs2now : numpy array 
            rotation matrix to transform cartesian icrs to present equinox
        
        current2topoxyz : numpy array 
            rotation matrix to transform current equinox to topocentric XYZ

        icrs2topoxyz : numpy array
            matrix product of icrs2now and current2topoxyz; tranforms icrs to topocentric XYZ

        constel : string
            IAU abbreviation of the constellation, e.g. 'UMa'

        previouslstmid : astropy Angle 
            sidereal time of previous night's midnight 

        moonpos : SkyCoord 
            topocentric moon celestial position to about 10 arcsec precision in current equinox

        moonobjang : astropy Angle 
            angle subtended by Observation's celest and the moon

        moonphasedescr : string 
            human-readable description of the moon phase

        moonillumfrac : float 
            fraction of moon's face that is illuminated (0. to 1.)

        moonaltit : astropy Angle 
            elevation of moon above horizon

        lunskybrightness : float 
            predicted sky brightness added by the moon in mag per square arcsec

        sunpos :  astropy SkyCoord 
            sun position to about 0.01 degrees

        sunobjang : astropy Angle 
            angle subtended by sun and moon

        sunaltit : astropy Angle 
            elevation of sun center above horizon (no refraction)

        sunaz : astropy Angle 
            azimuth of sun

        twi : float 
            rough zenight twilight brightness in blue mag per square arcsec

        lstmid : astropy Angle 
            local sidereal time at midnight

        tsunset, tsunrise : astropy Time 
            times of sunset and sunrise on present night

        tevetwi, tmorntwi : astropy Time 
            end and start of astronomical (18-degr) twilight 
        
        teve12,tmorn12 : astropy Time 
            end and start of nautical (12-degr) twilight

        tnightcenter : astropy Time 
            time of lower culmination of the sun

        uptime30, uptime20, uptime15 : astropy TimeDelta 
            How long object has <3, <2, and <1.5 airmasses during night defined by 12-deg twilight.

        barytcorr : astropy TimeDelta 
            Amount to add to Observation.t to get arrival time at barycenter

        baryvcorr : astropy Quantity 
            velocity to add to shift velocity to solar system barycenter frame

        planetdict : Dictionary of astropy SkyCoords 
            sky positions of planets, keyed by lowercase name

        planetmags : Dictionary of floats  
            apparent magnitudes of planets

        """
            # set up the site.

        self.setsite(site, default_site = default_site) 

        # set up the time.

        if t == None :   # set to system clock if no time specified.
            self.t = Time(ttime.time(),format='unix')

        else : 
            self.settime(t,use_local_time = use_local_time)

        # set up the celestial coordinates.

        if celest != None and celest != 'ZENITH' :  
            self.setcelest(celest)    # many legal formats -- see below. 
           
        else : # default to zenith
            self.lst = lpsidereal(self.t, self.site.location)    # sidereal
            self.nowfr = currentFK5frame(self.t)
            self.celest = SkyCoord(self.lst,self.site.location.lat,frame=self.nowfr)  # zenith

        self.cel_J2000 = self.celest.transform_to('icrs')
        # print("after xform, self.celest is",self.celest)
    
        self.constel = self.celest.get_constellation(short_name = True)

        self.previouslstmid = None

        self.moonpos = None   # SkyCoord of moon
        self.moonobjang = None
        self.moonphasedescr = None
        self.moonillumfrac = None
        self.moonaltit = None
        self.moonaz = None
        self.lunskybright = None  # approximate lunar contrib to sky 
                                  # brightness at this location
        self.sunpos =  None   # SkyCoord of sun
        self.sunobjang = None
        self.sunaltit = None
        self.sunaz = None
        self.twi = None  # rough magnitude difference between dark
                # night sky and present blue twilight contrib in zenith.

        self.lstmid = None
        self.tsunset = None
        self.tevetwi = None
        self.tnightcenter = None
        self.tmorntwi = None
        self.tsunrise = None

#        self.moonrise = None
#        self.sunrise = None

        self.barytcorr = None
        self.baryvcorr = None

        self.planetdict = {}  # dictionary of planet SkyCoords by name.
        self.planetmags = {}  # dictionary of planet visual mags by name.

    def setsite(self, site, default_site = 'mdm') :
        """Set the site on earth of the observation.

        Parameters
        ----------
        site : obs_site instance or string
            if string, it is the key used for the site in the site dictionary.

        """
        if isinstance(site, obs_site) : 
            self.site = site

        else : 
#            sitedict = get_sites()
            # print "got sites; ",sitedict.keys()
            if site == None :
                if default_site != None :
                     self.site = sitedict[default_site]
                    # 'backup default' is mdm
                else : self.site = sitedict['mdm']
            else :
                self.site = sitedict[site]
    

    def settime(self, t, use_local_time = True) :
        """set Observation time t using the input 't' 

        Note the instance variable 't' is always stored as 
        in utc.  If use_local_time, then the input is assumed to be local
        zone time (though if t is already a Time instance this is ignored).
        Input t can be an ISO string, or a tuple of at least yr, month,
        day, hour, and minute.

        Note also that if use_local_time is True (default) and the site 
        switches back and forth from standard to daylight saving time, then
        one hour of clock time gets repeated in the autumn, and an hour gets
        skipped in the spring, which is tricky if localtime is input.  For
        example, specifying 2:30 AM Eastern time on the 'spring forward' night 
        sets to 3:30 AM Eastern Daylight time; 2:30 AM Eastern Standard does
        not happen on that night.  When Daylight Saving ends, a time in the 
        ambiguous hour reverts to standard, so there's an hour of time that 
        you just can't specifiy using local time input.
  
        Parameters
        ----------
        t : a Time, or a string, or a tuple
            time to which to set Observation attribute t.  If not already a Time,
            
            * if a string, can be anything the dateutil parser accepts,e.g."2023-04-17 19:08:23"
            * if a tuple, must be (year, month, day, hour, minute, [sec]) as ints; seconds are optional.

        use_local_time : boolean 
            whether input time is interpreted as local or not; default True.

        """

        if isinstance(t, Time) : 
            self.t = t

        elif isinstance(t,str) :  # generally an ISO string like '2018-07-22 14:13:22'
            if use_local_time : 
               localt = self.site.localtz.localize(dateutil.parser.parse((t)))
               self.t = Time(localt)
            else : 
               self.t = Time(dateutil.parser.parse(t))

        elif isinstance(t,tuple) or isinstance(t,list) :
            if len(t) == 5 : 
                dtin = datetime(t[0],t[1],t[2],t[3],t[4])
            else :
                dtin = datetime(t[0],t[1],t[2],t[3],t[4],t[5])
            if use_local_time : 
               # print "using local"
               localt = self.site.localtz.localize(dtin)
               self.t = Time(localt)
            else : 
               # print "not using local"
               self.t = Time(dtin)

        elif isinstance(t, float) : 
            self.t = Time(t,format = 'jd') 

        # Keep a 'pure number' version of the julian epoch.
        self.julyear = 2000. + (self.t - thorconsts.J2000_Time).jd / 365.25 

    def advancetime(self, delta, forward = True) :
        """Advance Observation.t by the specified amount.

        Parameters
        ----------
        delta : Amount to advance the time, as an astropy TimeDelta, or a string 
               such as "143 s" or "5 d" (see below), or a float in sec.
        forward : boolean : if True, go forward; if False, go back.
        
        When specifying the interval with a string, the codes for units are

        * s - second
        * m - minute
        * h - hour
        * t - 1 sidereal day (transit-to-transit)
        * d - 1 solar day
        * w - one week of solar days
        * l - one lunation
        * y - one Julian year of 365.25 days

        """
            # non-obvious: 't' for 'transit' is 1 sidereal day, 
            # 'w' is a week, 'l' is exactly 30 days -- a very rough 'lunation' that 
            # will land on the same time of night about one month later -- and
            # 'y' is a 365-day 'year', again because it's more likely you want to
            # land on the same time of night next year than keep track over many years.
 
        if isinstance(delta, TimeDelta) :
            if forward :
                self.t = self.t + delta
            else :
                self.t = self.t - delta

        elif isinstance(delta, str) :   # e.g., "123. or 2. d"

            # codes for time delta intervals, and their values in seconds.
            # non-obvious: 't' for 'transit' is 1 sidereal day, 
            # 'w' is a week, 'l' is exactly 30 days -- a very rough 'lunation' that 
            # will land on the same time of night about one month later -- and
            # 'y' is a 365-day 'year', again because it's more likely you want to
            # land on the same time of night next year than keep track over many years.
            t_unit_dict = {'s' : 1., 'm' : 60., 'h' : 3600., 't' : 86164.0905352, 
               'd' : 86400, 'w' : 604800., 'l' : 2592000., 'y' : 31536000.}
               
            
            x = delta.split()
            try : 
                deltafloat = float(x[0])  
                if len(x) > 1 : 
                    deltafloat = deltafloat * t_unit_dict[x[1][0]]
                    dt = TimeDelta(deltafloat, format = 'sec')
                if forward :
                    self.t = self.t + dt
                else : 
                    self.t = self.t - dt
            except :
                print("Bad time step string.")

        elif isinstance(delta, float) :  # float defaults to seconds.
            dt = TimeDelta(delta, format = 'sec')  
            if forward : 
                self.t = self.t + dt
            else : 
                self.t = self.t - dt

        else : 
            print("Time step must be a float, an astropy TimeDelta or a string.")

        self.julyear = 2000. + (self.t - thorconsts.J2000_Time).jd / 365.25 
 
    def setcelest(self,celestin,ra_unit = u.hourangle, dec_unit = u.deg) :  
        """Sets the celestial coords.  

        Input can be a SkyCoord instance, a list or tuple of (ra,dec) or
        (ra,dec,equinox), or a character string.  

        Parameters 
        ----------

        celestin : Skycoord, or list/tuple, or string. Celestial coordinates to set.  See below.
        ra_unit : astropy unit; input unit for RA, defaults to u.hourangle  (i.e, hrs min sec)
        dec_unit : astropy unit; input unit for dec, defaults to u.deg 

        If celestin is a string, some permissible formats are 

        * '18:22:22.3  -0:18:33'   (defaults to ICRS or basically J2000)
        * '18:22:22.3  -0:18:33 2015'
        * '18 22 22.3  -0 18 33'
        * '18 22 22.3  -0 18 33 2015'

        Input units only apply if celestin is not already a SkyCoord.

        """
        if isinstance(celestin, SkyCoord) :  # if it's a SkyCoord, just copy it.
            self.celest = celestin
        elif isinstance(celestin, tuple) or isinstance(celestin, list) :   # 
            if len(celestin) == 2 :
                self.celest = SkyCoord(celestin[0],celestin[1], unit=(ra_unit, dec_unit))   # parse a tuple
            elif len(celestin) == 3 :
               # print("celestin[2] = ",celestin[2])
               eq = float(celestin[2])
               if eq == 2000. : 
                   self.celest = SkyCoord(celestin[0],celestin[1],unit = (ra_unit, dec_unit))
               else : 
                   eq = "J%7.2f" % (eq)
                   self.celest = SkyCoord(celestin[0],celestin[1],unit = (ra_unit, dec_unit),frame = FK5(equinox = eq))
        elif isinstance(celestin, str) :   # str includes unicode in python3.
            pieces = celestin.split() 
            if len(pieces) >= 6 :  # space-separated fields - glue back together.
                rastr = pieces[0] + ":" + pieces[1] + ":" + pieces[2]
                decstr = pieces[3] + ":" + pieces[4] + ":" + pieces[5]
                if len(pieces) == 7 :   # if there's an equinox ...
                   eqstr = pieces[6]
                else : eqstr = '2000.'
            else :    # e.g. 21:29:36.2   so split gets ra and dec separately
                rastr = pieces[0]  
                decstr = pieces[1]
                if len(pieces) > 2 : eqstr = pieces[2]  # if there's an equinox ...
                else : eqstr = '2000.'
            # print "rastr decstr eqstr",rastr,decstr,eqstr
            if float(eqstr) == 2000. :  # default to ICRS if 2000.
                self.celest = SkyCoord(rastr,decstr,unit = (ra_unit, dec_unit))
            else :   # or set in FK5 frame of date if not 2000.
                eq = "J"+eqstr
                self.celest = SkyCoord(rastr,decstr,unit = (ra_unit, dec_unit),frame = FK5(equinox = eq))

#        print(" ************ IN SETCELEST ************ ")
#        print( "self.celest:",self.celest)
#        print( "frame.name:",self.celest.frame.name)
                
        if self.celest.frame.name == 'icrs' or self.celest.frame.name == 'gcrs' : 
            self.cel_J2000 = self.celest
        elif self.celest.frame.name == 'precessedgeocentric' : 
            self.cel_J2000 = self.celest.transform_to('gcrs')
        else :
            self.cel_J2000 = self.celest.transform_to('icrs')

#        print "self.cel_J2000:",self.cel_J2000

        
    def computesky(self, redo_coords = True) :  
        """Compute basic observational circumstances.

        Quantities computed include lst, current-equinox coords, 
        hour angle, alt, az, parallactic angle, midhight, LST mid,.  Precession to
        current equinox can be suppressed for speed.

        Parameters
        ----------

        redo_coords : boolean, defaults to true. If False, turns off precession.

        """

        # current equinox etc., useful for repeat calls tracing out a single
        # night

        # Also, if redo_coords is on, computes rotation matrices to 
        # take celestial cartesian coords in ICRS directly into the 
        # observer's topocentric frame.  This speeds up such things
        # as the star display by a large factor.

#        print("entering computesksy: self.celest = ", self.celest)

        self.lst = lpsidereal(self.t, self.site.location)    # sidereal

        # use these to test for 'this is the same night.'
        self.midnight = local_midnight_Time(self.t,self.site.localtz) 
        self.lstmid = lpsidereal(self.midnight,self.site.location)

        if redo_coords : 
            if self.celest.frame.name == 'gcrs' or self.celest.frame.name == 'precessedgeocentric' :
               self.nowfr = currentgeocentframe(self.t)
 
            else : self.nowfr = currentFK5frame(self.t)

#        print("self.nowfr = ",self.nowfr)
        
            self.celnow = self.celest.transform_to(self.nowfr)

        self.hanow = (self.lst - self.celnow.ra).wrap_at(12.*u.hour)     # hour angle
        (self.altit, self.az, self.parang) = altazparang(self.celnow.dec, self.hanow, 
               self.site.location.lat)  
        self.airmass = true_airmass(self.altit)  # polynomial expansion.

        if redo_coords : 
            self.constel = self.celest.get_constellation(short_name = True)

            # compute some convenient rotation matrices:    
    
            # compute matrix for precession from J2000 (i.e., icrs almost exactly) 
            # to now, for later use.
            self.icrs2now = precessmatrix(thorconsts.J2000_Time,self.t)
            # and matrix to rotate a current-equinox celestial xyz into topocentric
            self.current2topoxyz = cel2localmatrix(self.lst, self.site.location.lat)
            # and matrix product to rotate icrs into topocentric xyz.
            self.icrs2topoxyz = self.current2topoxyz.dot(self.icrs2now) # matrix product 

            # If we're in the south, rotate by 180 degrees around vertical axis
            # to invert display.  Simply negate top two rows of the matrix.
            if self.site.location.lat < 0. * u.deg : 
                self.icrs2topoxyz = np.array([[-1],[-1],[1]]) * self.icrs2topoxyz

#        print("leaving computesksy: self.celest = ", self.celest)


    def computebary(self) :  # this is a bit expensive so split off.
        """computes and set bary corrections using astropy routines.  Slow; nearly all users will want 'computequickbary'.
        """
        self.baryvcorr = self.celest.radial_velocity_correction(obstime = self.t,
            location = self.site.location).to(u.km / u.second) 
#        print "baryvcorr: ", self.baryvcorr.to(u.km/u.s)
#        print type(self.baryvcorr.to(u.km/u.s))
        self.barytcorr = self.t.light_travel_time(self.celnow, kind='barycentric',
           location = self.site.location, ephemeris = 'builtin')
        self.tbary = self.t + self.barytcorr
#        print "barytcorr: ", self.barytcorr.to(u.s)
#        print type(self.barytcorr.to(u.s))

    def computequickbary(self) :  
        """compute and set barycentric corrections to 0.1s and 10m/s. Reasonably fast.
        """

        # The call to computebary is very expensive in
        # order to get ultimate accuracy.  When less precision is 
        # required it should  be possible to make it much faster by 
        # computing the corrections directly.

        # The goal is to get something like 0.1 sec in light time, and 
        # 10 meters per second in velocity.  

        earthpos = get_body_barycentric('earth',self.t)

        # astropy doesn't seem to offer a velocity function.  
        # For moderate precision, the simplest numerical derivative works:

        x1 = get_body_barycentric('earth', self.t - (500. * u.s))
        x2 = get_body_barycentric('earth', self.t + (500. * u.s))

        earthvel = (x2 - x1) / (1000. * u.s)

        # ignore geographic correction to light time delay -- it's
        # 21 msec for a source at the zenith.

        # nice to have a builtin dot product function for this!
        extra = earthpos.dot(self.cel_J2000.cartesian)  # extra distance to barycenter
        tcorr = extra / ((299792.458 * u.km) / (1.0 * u.s))
        self.barytcorr = TimeDelta(tcorr)   # needs to be a TimeDelta

        self.tbary = self.t + self.barytcorr  # this is the right sign!
        
        #  earth-rotation part of velocity correction 

        diurnalvel = CartesianRepresentation(  
            -1. * self.site.diurnalspeed * np.sin(self.lst),
                  self.site.diurnalspeed * np.cos(self.lst),
                  self.site.diurnalspeed * 0.)

        # totalvel should be earth orbit + earth rotation velocity wrt barycenter
        
        totalvel = earthvel + diurnalvel 

        self.baryvcorr = totalvel.dot(self.cel_J2000.cartesian).to(u.km/u.s)

        # print("self.barytcorr ",self.barytcorr,"self.baryvcorr",self.baryvcorr)
 
    def computesunmoon(self) :
        """Compute and set the many moon and sun quantities."""
        self.lst = lpsidereal(self.t, self.site.location)
        self.moonpos, self.moondist = accumoon(self.t,self.site.location)  
        self.moonha = self.lst - self.moonpos.ra
        (self.moonaltit, self.moonaz, parang) = altazparang(self.moonpos.dec, self.moonha, self.site.location.lat)
        self.sunpos = lpsun(self.t)  
        self.sunha = self.lst - self.sunpos.ra
        (self.sunaltit, self.sunaz, parang) = altazparang(self.sunpos.dec, self.sunha, self.site.location.lat)
        self.twi = ztwilight(self.sunaltit)

        self.sunmoonang = self.sunpos.separation(self.moonpos) 
        self.moonillumfrac = 0.5 * (1. - np.cos(self.sunmoonang))
        self.moonobjang = self.celnow.separation(self.moonpos)
        (self.moonphasedescr, self.lunage, self.lunation) = phase_descr(self.t.jd)
        # print "age %f lunation %d" % (self.lunage, self.lunation)

#        print "moon illum frac:",self.moonillumfrac
#        print "moon-obj ang:", self.moonobjang
#        print "moon altit",self.moonaltit,"obj altit",self.altit

        self.lunsky = lunskybright(self.sunmoonang,self.moonobjang,thorconsts.KZEN,self.moonaltit,
             self.altit,self.moondist, self.sunaltit)

#        print "lunsky: ",self.lunsky

#        print "lst",self.lst
#        print "moon",self.moonpos
#        print "moon ha, alt, az", self.moonha, self.moonaltit, self.moonaz
#        print "sun",self.sunpos
#        print "sun ha, alt, az", self.sunha, self.sunaltit, self.sunaz
#        print "twilight diff: " , self.twi

    def computeplanets(self) : 
        """compute and store planetary positions and magnitudes."""

#        print "starting planets"
        planetlist = ['mercury','venus','mars','jupiter','saturn','uranus','neptune']

        # to get magnitudes need to get sun and earth positions too

        sunbary = get_body_barycentric('sun',self.t)
#        print 'sunbary:' 
#        print sunbary
        earthbary = get_body_barycentric('earth',self.t)
#        print 'earthbary:'
#        print earthbary

        for p in planetlist : 
            # get celestial position of planet
            self.planetdict[p] = get_body(p,self.t)
        #    print "get_body gives",p, self.planetdict[p]
        #  This returns a position in "GCRS: which is earth-centered.

            # now get sun-centered location (slightly different from bary) to get
            # magnitude.

            pbary = get_body_barycentric(p,self.t)
            psun =  sunbary - pbary # vector from sun to planet
            psundist = np.sqrt(psun.x ** 2 + psun.y ** 2 + psun.z ** 2) # modulus

            pearth = earthbary - pbary  # vector from planet to earth
            pearthdist = np.sqrt(pearth.x ** 2 + pearth.y ** 2 + pearth.z ** 2)
   
            # for inner planets, use polynomials for the phase angle dependence,
            # which is not at all trivial.

            if p == 'mercury' or p == 'venus' or p == 'mars' : 
                # angle between sun and earth as viewed from planet
                phaseang = angle_between(psun,pearth)
                # print "phaseang:",phaseang
                phasefac = np.polyval(thorconsts.PLANETPHASECOEFS[p],phaseang.value)
                # print "phasefac:",phasefac
                self.planetmags[p] = phasefac + 5. * np.log10(psundist.to(u.AU).value * pearthdist.to(u.AU).value)
                # print "mag:",self.planetmags[p]

            # outer planets are close enough to phase zero all the time to ignore the phase angle.
            else : 
                phasefac = thorconsts.PLANETPHASECOEFS[p]
                self.planetmags[p] = phasefac + 5. * np.log10(psundist.to(u.AU).value * pearthdist.to(u.AU).value)
                # print "mag:",self.planetmags[p]
            # saturn will not be good because there's no ring-tilt dependence factored in.
 
        fr = currentgeocentframe(self.t)
        #print "frame attributes:",fr.get_frame_attr_names()
#        print """
#after transformation:
#   """

        # we want these in equinox of date for plotting.  They'll be 
        # converted back to J2000 for the table when the coordinates are loaded
        # into the observation instance.

        for p in planetlist : 
            self.planetdict[p] = self.planetdict[p].transform_to(fr)
#            print p, self.planetdict[p].to_string('hmsdms')
#        print "ending planets"
        

    def setnightevents(self) : 
        """Compute the events (sunset etc) for a single night.  """

        # self.midnight also computed in computesky, but it's cheap.

        self.midnight = local_midnight_Time(self.t,self.site.localtz) 
      
        self.lstmid = lpsidereal(self.midnight,self.site.location)

        # if you're looking at the same night, from the same site,
        # lst mid will be the same.   Don't bother with the calculation

        if self.previouslstmid != None :
           if abs(self.previouslstmid - self.lstmid) < 0.001 * u.deg :
               # print "no night event calc'n -- same."
               return

        sunmid = lpsun(self.midnight) 
   
        # sunrise and sunset altitudes are complicated and initialized
        # with the site for efficiency's sake.  Twilight altitude is
        # fixed at -18 degrees, so not costly to set it here.
 
        twialt = Angle(-18.,unit=u.deg)     
        twialt12 = Angle(-12.,unit=u.deg)

        # Compute sunset, sunrise, and twilight times.
        # Start by computing the approximate hour angles at which these
        # occur, for the dec that the sun has at midnight.

        # for this purpose, wrap hour angles at noon so that all for a given night
        # are positive

        # Find hour angle at which the dec of the sun (evaluated at
        # local midnight) rises or sets
        sunsetha = ha_alt(sunmid.dec,self.site.location.lat,self.site.setalt)  
        # If the dec of the sun never rises or sets -- possible in the arctic -- set 
        # flags for later use; "0" is normal, it rises and sets, "1" is it's always
        # up (midnight sun) and "-1" is it never rises.
        if sunsetha > (500. * u.rad) : 
            self.sunsetflag = 1   # sun always up
        elif sunsetha < (-500. * u.rad) :
            self.sunsetflag = -1   # sun never rises or sets 
        else : self.sunsetflag = 0 
        sunriseha = Angle(2. * np.pi, unit = u.rad) - ha_alt(sunmid.dec,self.site.location.lat,self.site.risealt)  

        twilightha = ha_alt(sunmid.dec,self.site.location.lat,twialt)  # positive, correct for evening
#        print "sunsetha, sunriseha, twilightha",sunsetha,sunriseha,twilightha
        # Again, set flag in case twilight never ends (high-latitude summer) or the sun doesn't get
        # higher than -18 degrees 
        if twilightha > (500. * u.rad) : 
            self.twilightflag = 1   # never gets dark
        elif twilightha < (-500. * u.rad) :
            self.twilightflag = -1   # fully dark all night (only happens near pole and solstice)
        else : self.twilightflag = 0 

        # do the same for 12-degree twilight; with the sun between 12 and 18 degrees below the horizon the
        # sky is fairly dark and one can work on brighter object, standard stars and so on.

        twilight12ha = ha_alt(sunmid.dec,self.site.location.lat,twialt12)  # positive, correct for evening
        if twilight12ha > (500. * u.rad) : 
            self.twilight12flag = 1   # never gets dark
        elif twilight12ha < (-500. * u.rad) :
            self.twilight12flag = -1   # fully dark all night (only happens near pole and solstice)
        else : self.twilight12flag = 0 
 
        hasunmid = (self.lstmid - sunmid.ra).wrap_at(24. * u.hour)
        #print "hasunmid:",hasunmid
        #print "midnight",self.midnight

        self.tnightcenter = self.midnight - TimeDelta(hasunmid.hour / 24. - 0.5, format = 'jd')
        #self.lstnightcenter = lpsidereal(self.tnightcenter,self.site.location)

        #print 'night center',self.nightcenter

        if self.sunsetflag == 0 :  # if dec of sun is such that sunset and sunrise occur
            sunsetguess = hasunmid - sunsetha    # hour angle difference from sun's posn at midnight
            sunriseguess =  sunriseha - hasunmid
        if self.twilightflag == 0 : 
            evetwiguess = hasunmid - twilightha
            morntwiguess = Angle(2.*np.pi, unit=u.rad) - twilightha - hasunmid
        if self.twilight12flag == 0 : 
            evetwi12guess = hasunmid - twilight12ha
            morntwi12guess = Angle(2.*np.pi, unit=u.rad) - twilight12ha - hasunmid

        #print "sunsetguess, sunriseguess",sunsetguess,sunriseguess.hour
        #print "eve, morn",evetwiguess,morntwiguess.hour

        # convert to time differences

        if self.sunsetflag == 0 :
            TDsunset = TimeDelta(sunsetguess.hour / 24., format = 'jd')
            TDsunrise = TimeDelta(sunriseguess.hour / 24., format = 'jd')

        #print "tdsunset, tdsunrise",TDsunset,TDsunrise

        if self.twilightflag == 0 :
            TDevetwi = TimeDelta(evetwiguess.hour / 24., format = 'jd')
            TDmorntwi = TimeDelta(morntwiguess.hour / 24., format = 'jd')
        #print "TDeve, TDmorn",TDevetwi,TDmorntwi
        if self.twilight12flag == 0 :
            TDevetwi12 = TimeDelta(evetwi12guess.hour / 24., format = 'jd')
            TDmorntwi12 = TimeDelta(morntwi12guess.hour / 24., format = 'jd')

        # form into times and iterate to accurate answers.
     
        if self.sunsetflag == 0 :
            self.tsunset = self.midnight - TDsunset  # first approx
            #print "first approx",self.tsunset
            self.tsunset = jd_sun_alt(self.site.setalt, self.tsunset, self.site.location)
     
            self.tsunrise = self.midnight + TDsunrise  # first approx
            #print "first approx",self.tsunrise
            self.tsunrise = jd_sun_alt(self.site.risealt, self.tsunrise, self.site.location)

        if self.twilightflag == 0 :
     
            self.tevetwi = self.midnight - TDevetwi
            self.tevetwi = jd_sun_alt(twialt, self.tevetwi, self.site.location)
     
            self.tmorntwi = self.midnight + TDmorntwi
            self.tmorntwi = jd_sun_alt(twialt, self.tmorntwi, self.site.location)

        if self.twilight12flag == 0 :
     
            self.tevetwi12 = self.midnight - TDevetwi12
            self.tevetwi12 = jd_sun_alt(twialt12, self.tevetwi12, self.site.location)
            #self.lsteve12 = lpsidereal(self.tevetwi12,self.site.location)
     
            self.tmorntwi12 = self.midnight + TDmorntwi12
            self.tmorntwi12 = jd_sun_alt(twialt12, self.tmorntwi12, self.site.location)
            #self.lstmorn12 = lpsidereal(self.tmorntwi12,self.site.location)

        # and, moonrise and set times for that night.

        moonmid = lpmoon(self.midnight, self.site.location)
        hamoonmid = self.lstmid - moonmid.ra
        hamoonmid.wrap_at(12. * u.hour, inplace = True)

        #print "moon at midnight",moonmid
        #print "hamoonmid: ",hamoonmid.hour, 'hr'

        roughlunarday = TimeDelta(1.0366, format = 'jd')

        moonsetha = ha_alt(moonmid.dec,self.site.location.lat,self.site.setalt)

        # Using the midnight position of the moon to assess whether it actually
        # rises or sets in a 12-hour window around that time is problematic, 
        # since the moon's dec can move pretty quickly.  This is a rare 'corner
        # case' that only matters at very high latitudes so I'm not going to worry 
        # about it too much.  
        if moonsetha > (500. * u.rad) : 
            self.moonsetflag = 1   # moon always up 
        elif moonsetha < (-500. * u.rad) :
            self.moonsetflag = -1   # moon always below horizon
        else : self.moonsetflag = 0  
 
        moonsetdiff = moonsetha - hamoonmid   # how far from setting at midnight
        # find nearest setting point 
        if moonsetdiff.hour >= 12. : moonsetdiff = moonsetdiff - Angle(24. * u.hour)
        if moonsetdiff.hour < -12. : moonsetdiff = moonsetdiff + Angle(24. * u.hour)
        TDmoonset = TimeDelta(moonsetdiff.hour / 24., format = 'jd')
        self.tmoonset = self.midnight + TDmoonset

        #print "moonset first approx:",self.tmoonset
        #print "aiming for set alt = ",self.site.setalt
        self.tmoonset = jd_moon_alt(self.site.setalt, self.tmoonset, self.site.location)
        #print "moonset: ",self.tmoonset # .to_datetime(timezone = localtzone)

        moonriseha = -1. * ha_alt(moonmid.dec,self.site.location.lat,self.site.risealt) # signed 
        moonrisediff = moonriseha - hamoonmid  # how far from riseting point at midn.
        # find nearest riseing point 
        if moonrisediff.hour >= 12. : moonrisediff = moonrisediff - Angle(24. * u.hour)
        if moonrisediff.hour < -12. : moonrisediff = moonrisediff + Angle(24. * u.hour)
        TDmoonrise = TimeDelta(moonrisediff.hour / 24., format = 'jd')
        self.tmoonrise = self.midnight + TDmoonrise
        #print "moonrise first approx:",self.tmoonrise
        #print "aiming for rise alt = ",self.site.risealt
        self.tmoonrise = jd_moon_alt(self.site.risealt, self.tmoonrise, self.site.location)
        #print "moonrise: ",self.tmoonrise # .to_datetime(timezone = localtzone)

        # Save this to avoid re-doing unnecessarily.  If lstmid is exactly the same,
        # then the night and location are almost certainly unchanged.

        

        self.previouslstmid = self.lstmid

    def compute_hours_up(self)  : # involves night events but is specific to this object.
        """compute how long an object is up during the night, defined by
        12-degree (or nautical) twilight.  This requires setnightevents 
        to have been run previously.  Computes how long the object is less
        that 3.0 airmass, 2.0 airmasses, and 1.5 airmasses."""

        # this requires setnightevents to have been run for the same night.

        minalt, maxalt = min_max_alt(self.site.location.lat, self.celnow.dec) 

        if self.twilight12flag == 0 : 

            self.ha_mid = (self.lstmid - self.celnow.ra).wrap_at(12. * u.hour)
            #print("self.ha_mid.value",self.ha_mid.value)
            #print("self.ha_mid.hourangle",self.ha_mid.hourangle)
            deltattran = TimeDelta(self.ha_mid.hourangle / 24., format = 'jd') /  1.0027379093
            self.ttransit = self.midnight - deltattran
            #print("self.midnight = ",self.midnight)
            #print("self.ttransit = ",self.ttransit)
            
            if minalt < thorconsts.ALT30 and maxalt > thorconsts.ALT30 :  
                # if this dec passes through 3 airmasses
                ha30 = ha_alt(self.celnow.dec,self.site.location.lat,thorconsts.ALT30)
                dt30 = TimeDelta(ha30.hourangle / 24., format = 'jd') / 1.0027379093
                jd30_1 = self.ttransit - dt30   # Time of rise through 3 airmasses
                jd30_2 = self.ttransit + dt30   # Time of set past 3 airmasses 
            #    print("jd30_1 = ",jd30_1)
            #    print("jd30_2 = ",jd30_2)
                self.uptime30 = hrs_up(jd30_1,jd30_2,self.tevetwi12,self.tmorntwi12)

            elif minalt > thorconsts.ALT30 : self.uptime30 = (self.tmorntwi12 - self.tevetwi12)
            elif maxalt < thorconsts.ALT30 : self.uptime30 = thorconsts.ZERO_TIMEDELTA
            #print("time above 3 airm", self.uptime30)

            if minalt < thorconsts.ALT20 and maxalt > thorconsts.ALT20 :  
                # if it passes through 2 airmass 
                ha20 = ha_alt(self.celnow.dec,self.site.location.lat,thorconsts.ALT20)
                dt20 = TimeDelta(ha20.hourangle / 24., format = 'jd') / 1.0027379093
                jd20_1 = self.ttransit - dt20
                jd20_2 = self.ttransit + dt20
                self.uptime20 = hrs_up(jd20_1,jd20_2,self.tevetwi12,self.tmorntwi12)
            elif minalt > thorconsts.ALT20 : self.uptime20 = (self.tmorntwi12 - self.tevetwi12)
            elif maxalt < thorconsts.ALT20 : self.uptime20 = thorconsts.ZERO_TIMEDELTA
            #print("time above 2 airm", self.uptime20)
 

            if minalt < thorconsts.ALT15 and maxalt > thorconsts.ALT15 :  
                # if it passes through 1.5 airmasses
                ha15 = ha_alt(self.celnow.dec,self.site.location.lat,thorconsts.ALT15)
                dt15 = TimeDelta(ha15.hourangle / 24., format = 'jd') / 1.0027379093
                jd15_1 = self.ttransit - dt15
                jd15_2 = self.ttransit + dt15
                self.uptime15 = hrs_up(jd15_1,jd15_2,self.tevetwi12,self.tmorntwi12)
            elif minalt > thorconsts.ALT15 : self.uptime15 = (self.tmorntwi12 - self.tevetwi12)
            elif maxalt < thorconsts.ALT15 : self.uptime15 = thorconsts.ZERO_TIMEDELTA
            #print("time above 1.5 airm", self.uptime15)

    def printnow(self, use_slow_bary = False) :
        """Compute the instantaneous circumstances to ensure they're current and 
        print out a nicely formatted display."""

        # first ensure that they're up to date ...

        self.computesky()
        self.computesunmoon()
        self.computeplanets()
        if use_slow_bary :
            self.computebary()
        else :
            self.computequickbary()

        print(" ")
        print("Site : %s;  E longit = %s, lat = %s" % (self.site.name,
              self.site.location.lon.to_string(unit = u.deg, sep=' '), 
              self.site.location.lat.to_string(unit = u.deg, sep=' ')))
        print(" ")
        print("     J2000:  %s  %s       (in %s)" % 
                 (self.cel_J2000.ra.to_string(unit = u.hourangle, sep = ' ',
                 precision=2, pad = True),
                 self.cel_J2000.dec.to_string(sep = ' ',precision = 1, pad = True,
                 alwayssign = True),
                 self.constel))
        eqout = self.t
        eqout.format = "jyear_str"
        print("%s :  %s  %s" % (eqout.value,
                 self.celnow.ra.to_string(unit = u.hourangle, sep = ' ',
                 precision=2, pad = True), self.celnow.dec.to_string(sep = ' ',precision = 1,
                              pad = True, alwayssign = True)))

        print(" ")
        ut = self.t.to_datetime()
        local = self.t.to_datetime(timezone = self.site.localtz)
#        localdow = dows[datetime.weekday()]
        print("UT date and time    : %s    JD %s " % (ut.strftime("%a %Y-%m-%d %H:%M:%S"), self.t.jd))
        print("local date and time : %s" % (local.strftime("%a %Y-%m-%d %H:%M:%S")))
        print(" ")
        print("Local mean sidereal time: %s " % (self.lst.to_string(unit = u.hourangle, sep = ' ',
            precision = 0)))
        print(" ")
        parang_opposite = self.parang + Angle(180 * u.deg) 
        parang_opposite.wrap_at(180. * u.deg)
        print("Hour angle: %s  AltAz: %5.1f, %6.1f  Parallactic: %4.1f [%4.1f]" % \
                        (self.hanow.to_string(unit = u.hourangle, sep = ' ', precision = 0, pad = True,
                         alwayssign = True),
                        self.altit.deg, self.az.deg, self.parang.deg, parang_opposite.deg))
        if self.altit < 0. : print("Below horizon.")
        elif self.airmass > 10. : print("Airmass > 10.")
        else : print("Airmass: %6.3f" % (self.airmass))
        print(" ")
        print("Moon: %s   Alt,Az %4.1f, %4.1f" % (self.moonphasedescr,
                      self.moonaltit.deg, self.moonaz.deg))
        if self.moonaltit > 0. :  # give more detail on the moon if it's up.
            print("Moon ra and dec:  %s  %s  (%s)" % (
                 self.moonpos.ra.to_string(unit = u.hourangle, sep = ' ',
                 precision=0,pad=True), self.moonpos.dec.to_string(sep = ' ',fields = 2,
                 pad = True, alwayssign = True), eqout.value))
            print("Illum. fract : %5.3f   Moon-obj ang: %5.1f deg" % (self.moonillumfrac, 
                           self.moonobjang.deg))
            if self.lunsky != 99. :
                print("Lunar sky brightness %4.1f mag/sq arcsec" % (self.lunsky))

        else : print("The moon is down. ") 
        print(" ")
 
        if self.sunaltit.deg < -18. : 
            print("The sun is down; there is no twilight.")
        elif self.sunaltit.deg < 0. :
            print("In twilight; sky %4.1f mag brighter than dark sky." % self.twi)
        else : 
            print("THE SUN IS UP.")
        print("Sun RA and dec:  %s  %s (%s);  AltAz %4.1f, %5.1f" % (
            self.sunpos.ra.to_string(unit = u.hourangle, sep = ' ',
               precision=1, pad = True), self.sunpos.dec.to_string(sep = ' ',precision = 0,
               pad = True, alwayssign = True),
               eqout.value, self.sunaltit.deg, self.sunaz.deg))
        print(" ") 
        print("Barycentric corrns: add %7.2f sec and %6.2f km/s to observed." % \
            (self.barytcorr.to(u.s).value,self.baryvcorr.to(u.km/u.s).value))
        print("Barycentric JD (UTC system): %14.5f." % (self.tbary.jd))

    def printnight(self, use_local_time = True) :
        """Print a nicely formatted display of the rise/set times etc.."""

        self.setnightevents() 
    
        if use_local_time : 
           tz = self.site.localtz
           print("Night events; times listed are local.\n")
        else :
           tz = None
           print("Night events; times listed are UT.\n")
        sunset = self.tsunset.to_datetime(timezone = tz)
        print("         Sunset:  %s" % (time_rounded_to_minute(sunset, incl_date = True, incl_day = True)))
        endtwi = self.tevetwi.to_datetime(timezone = tz)
        print("  Twilight Ends:  %s" % (time_rounded_to_minute(endtwi, incl_date = True, incl_day = True)))
        nghtctr = self.tnightcenter.to_datetime(timezone = tz)
        print("Center of Night:  %s" % (time_rounded_to_minute(nghtctr, incl_date = True, incl_day = True)))
        begtwi = self.tmorntwi.to_datetime(timezone = tz)
        print("Twilight Begins:  %s" % (time_rounded_to_minute(begtwi, incl_date = True, incl_day = True)))
        sunrise = self.tsunrise.to_datetime(timezone = tz)
        print("        Sunrise:  %s" % (time_rounded_to_minute(sunrise, incl_date = True, incl_day = True)))

        print(" ")
        moonrise = self.tmoonrise.to_datetime(timezone = tz)
        moonset = self.tmoonset.to_datetime(timezone = tz)

        if self.tmoonrise < self.tmoonset :
           print("       Moonrise:  %s" % (time_rounded_to_minute(moonrise, incl_date = True, incl_day = True)))
           print("        Moonset:  %s" % (time_rounded_to_minute(moonset, incl_date = True, incl_day = True)))
        else : 
           print("        Moonset:  %s" % (time_rounded_to_minute(moonset, incl_date = True, incl_day = True)))
           print("       Moonrise:  %s" % (time_rounded_to_minute(moonrise, incl_date = True, incl_day = True)))

    
        
if __name__ == "__main__" :

#    sitedict = get_sites()
#    obsite = sitedict['mdm']
   
#    obgeo = obsite.location.to_geodetic()
#    print obgeo

    # cel = SkyCoord("21:29:36.2 -47:04:08",unit = (u.hourangle, u.deg), frame = 'icrs')

#    print "year, month, day, hr, min"
#
#    x = raw_input().split()
#    year = int(x[0])
#    month  = int(x[1])
#    day = int(x[2])
#    hr = int(x[3])
#    minute = int(x[4])
#    
#    dt = obsite.localtz.localize(datetime(year,month,day,hr,minute,0))
#    t = Time(dt)
#    print t.jd
#
#    dt2 = obsite.localtz.localize(datetime(2000,1,1,0,0,0)) 
#    t2 = Time(dt2)
#    print t2.jd
#    
    print("ra dec: ")
    x = raw_input()
    cel = SkyCoord(x,unit = (u.hourangle, u.deg), frame = 'icrs')

    o = Observation(celest = cel, t = "2018-07-22T23:00:00", use_local_time = True, 
       site = None, default_site = 'keck') 

    obgeo = o.site.location.to_geodetic()
    print(obgeo)
    print( o.celest)
    print( o.t)
    print( "lst: ",o.lst)
    print( "celnow: ",o.celnow)
    print( "hanow : ",o.hanow)
    print( "alt, az, parang ",o.altit, o.az, o.parang)
          
    o.setnightevents()

    print("sunset:  ", o.tsunset.to_datetime(timezone = o.site.localtz))
    print("eve twi: ",o.tevetwi.to_datetime(timezone = o.site.localtz))
    print("night ctr:",o.tnightcenter.to_datetime(timezone = o.site.localtz))
    print("morn twi:",o.tmorntwi.to_datetime(timezone = o.site.localtz))
    print("sunrise: ",o.tsunrise.to_datetime(timezone = o.site.localtz))
    print("moonset: ",o.tmoonset.to_datetime(timezone = o.site.localtz))
    print("moonrise:",o.tmoonrise.to_datetime(timezone = o.site.localtz))

    o.computesunmoon()
    
    o.computeplanets()

    #print "c21:"
    #print c2l

    celnowxyz = cel2xyz(o.celnow)
    #print "celnowxyz:"
    #print celnowxyz

    topoxyz = o.current2topoxyz.dot(celnowxyz)
 
    #print "topoxyz:"
    #print topoxyz
    # print "topoxyz[0]",topoxyz[0]
    #az = np.arctan2(topoxyz[0],topoxyz[1]) * thorconsts.DEG_IN_RADIAN
    #alt = np.arcsin(topoxyz[2]) * thorconsts.DEG_IN_RADIAN
    #print "alt  az",alt,az

    #fullmat = c2l.dot(prec)   # both!
    #celxyz = cel2xyz(o.celest)   # icrs!
    #topo2 = fullmat.dot(celxyz)
    #print "topo2:"
    #print topo2
    #az2 = np.arctan2(topo2[0],topo2[1]) * thorconsts.DEG_IN_RADIAN
    #alt2 = np.arcsin(topo2[2]) * thorconsts.DEG_IN_RADIAN
    #print "alt2  az2",alt2,az2

    (bright2000, brightmags, brightcolors, brightnames) = getbrightest("cartesian_bright.dat")
 
    (projectedx,projectedy) = skyproject(o.icrs2topoxyz,bright2000)

    (objx,objy) = skyproject(o.current2topoxyz,celnowxyz)

    for i in range(0,len(brightmags)) : 
        size = (5. - 0.9 * brightmags[i]) 
        if size > 0.: 
            plt.plot(projectedx[i],projectedy[i],'bo',markersize = size)

    plt.plot(objx,objy,'ro')

    plt.xlim(-1.,1.)
    plt.ylim(-1.,1.)
    plt.show()
    
    
        
        
