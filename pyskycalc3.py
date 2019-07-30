#!/usr/bin/env python3

# Divorced from PGPLOT -- uses native Tkinter drawing commands
# to make picture of sky.  Allows mouse and keyboard input,
# as in the Java version.  

# Can open a socket to allow an external program to feed coords 
# etc. into skycalc.

# copyright John Thorstensen, Dartmouth College, 2005 January.
# 
# Anyone may use this program for non-commercial purposes, 
# and copy it and alter it as they like,  but this banner must 
# be preserved, and my name must be preserved in window titles,
# i.e. proper credit must be given. 

# import sys
# import _skysub
from datetime import datetime, timedelta
from thorsky.thorskyutil import flmoon, lunskybright, thorconsts
from thorsky.thorskyclasses3 import *
import numpy as np
from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter.font import nametofont
from copy import deepcopy 
from time import sleep   # for waiting when needed
#import threading
import urllib
#import panstarrs as ps

import os
import platform

# MacOS doesn't lift windows properly.  Found this bizarre fix suggested on stackexchange.
# Also inserted it after starting Tk and it still didn't work correctly.
#if platform.system() == "Darwin" :
#    os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')
 
#from PIL import Image as PyIm
#from PIL import ImageTk as PyImTk   # for auto-popping finding charts.  Needs to be 
# recent for interlaced png, and renamed to avoid name collisions.

#BRIGHTSTARFILE = "/home/thorsten/src/java/JSkyCalc/brightest.dat"

# The effect of the TINYPIX flag on MacOS is to make the sky display
# larger, without much affecting the font size.
# TINYPIX = TRUE  # flag to turn on scaling for Latitude E7450 screen
TINYPIX = False  # flag to turn on scaling for Latitude E7450 screen

# buttons optimized for Linux come off all scrunched on a Mac.  
# Let's try a fix here ... keep the name of the function short.

def butwid(linuxwid) :
    if platform.system() == "Darwin" :
        return int(1.5 * linuxwid)
    else : return linuxwid

# we'll need this later.
if platform.system() == "Darwin" :
    is_MacOS = True
else : is_MacOS = False

# On MacOS, Tkinter does not correctly pass back the event location of 
# a keystroke on a canvas.  This is documented here:

# https://stackoverflow.com/questions/21517114/python-tkinter-key-event-locations-lost-in-macos
# This implements the fix suggested in the answer on that page, which 
# appears to work well. 

def keylocation(event):
    # This is an attempt to work around a bug in MacOS Tkinter event key location.
    px,py = event.widget.winfo_pointerxy()
    rx,ry = (event.widget.winfo_rootx(), event.widget.winfo_rooty())
    cx,cy = (px-rx, py-ry)
    # msg = "event xy %d %d" % (cx, cy)
    # label1.config(text=msg)  
    return(cx,cy)
 
# import pdb
# pdb.set_trace()

# import matplotlib.pyplot as plt

def circlecenter(x1,y1,x2,y2,x3,y3) :
   # Doug McIlroy says that all great circles will project to circles
   # in a gnomic projection.
   # So, I'll compute three points on any weird circle and just 
   # plot the circle
   # formula for the center of a circle given three points.
   # http://mathforum.org/library/drmath/view/55239.html
 # pdb.set_trace()

#    ht = array([x1**2 + y1**2, y1, 1, x2**2 + y2**2 , y2, 1,
#               x3**2 + y3**2, y3, 1], shape = (3,3))
   ht = [[x1**2 + y1**2, y1, 1], [x2**2 + y2**2 , y2, 1],
               [x3**2 + y3**2, y3, 1]]
#    hb = array([x1, y1, 1, x2, y2, 1, x3, y3, 1],shape = (3,3))
   hb = [[x1, y1, 1],[x2, y2, 1],[x3, y3, 1]]
#    dhb = determinant(hb)
   dhb = np.linalg.det(hb)
#   h = determinant(ht) / (2. * dhb)
   h = np.linalg.det(ht) / (2. * dhb)
# 
#    kt = array([x1, x1**2 + y1**2, 1, x2, x2**2 + y2**2, 1, 
#                x3, x3**2 + y3**2, 1], shape=(3,3))
   kt = [[x1, x1**2 + y1**2, 1],[ x2, x2**2 + y2**2, 1], [x3, x3**2 + y3**2, 1]]
#    k = determinant(kt) / (2. * dhb)
   k = np.linalg.det(kt) / (2. * dhb)
# 
   radius = np.sqrt((x1 - h) ** 2 + (y1 - k) ** 2)

   return (h, k, radius)
 
# this was defined as a method of SkyDisplay but is re-used in the 
# airmass graph.

def mooncirc(obs) :

   # astropy version -- runs, not checked.
   
   # computes topocentric ra, dec, illumination angle, and PA of 
   # the cusps of the moon's illumination.

   sidt = obs.lst

   codecsun = Angle(90. * u.deg) - obs.sunpos.dec 
   codecmoon = Angle(90. * u.deg) - obs.moonpos.dec
   dra = obs.sunpos.ra - obs.moonpos.ra
   dra.wrap_at(180. * u.deg) 

   # spherical law of sines + law of cosines needed to get
   # sine and cosine of pa separately; this gets quadrant etc.!

   pasin = np.sin(dra) * np.sin(codecsun) / np.sin(obs.sunmoonang)
   pacos = (np.cos(codecsun) - np.cos(codecmoon) * np.cos(obs.sunmoonang)) / \
      (np.sin(codecmoon) * np.sin(obs.sunmoonang))
   pa = np.arctan2(pasin,pacos)  # aran2 is defined in numpy

#   print "pa of arc from moon to sun is %5.2f deg." % (pa * _skysub.DEG_IN_RADIAN)

   cusppa = pa - Angle(90. * u.deg)  # line of cusps ... 

   jdfull = flmoon(obs.lunation,2)  # jd of full on this lunation
#   print "lunation %d, jdfull %f" % (obs.lunation,jdfull)

   if obs.t.jd > jdfull :  # if past full, invert
       moonsun = obs.sunmoonang * -1.
   else : 
       moonsun = obs.sunmoonang

#   print "obs.sunmoonang ",obs.sunmoonang," moonsun ",moonsun

   moonsun.wrap_at(360. * u.deg, inplace = True)
   
   # print "moonsun after , cusp pa", moonsun, cusppa

   # silly to return moon coords already in obs, but other things expect it ...

   return [obs.moonpos,moonsun,cusppa,pa]
 
def twilightcolor(obs, moonzero = 23., moonfac = 8.)  :

   # sets background color used for
   # indicating twilight.   Routine hands back a tuple of 
   # (R,G,B) normalized 0 to 1.
   # The formulae below give a little step at 18 degree twilight,
   # followed by a steep ramp to 5 mag of twilight, followed by a
   # gentler ramp up to full daylight.

   if obs.sunaltit < Angle(-18. * u.deg) :  # total darkness 
     # UNLESS the moon is up
      val = 0   # R, G, and B will all be "val" 
      moonfac = moonfac # scales brightness of grey sky
      if(obs.moonaltit > Angle(0. * u.deg)) : 
#         [georam,geodm,geodism,toporam,topodecm,topodistm] = \
#              _skysub.accumoon(obs.jd,obs.lat,obs.sidereal.val,
#                        obs.elevsea)  # need topodistm 
         zenithlunsky = lunskybright(obs.sunmoonang, 
            Angle(90. * u.deg) - obs.moonaltit, thorconsts.KZEN, obs.moonaltit, 
            Angle(90. * u.deg), obs.moondist, obs.sunaltit)
         # print "zenithlunsky = ",zenithlunsky
         if zenithlunsky < moonzero :
            val = int(moonfac * (moonzero - zenithlunsky))
      return "#%02x%02x%02x" % (val,val,val)  # black
   if obs.sunaltit > Angle(-0.8 * u.deg) :  # sun is UP
      return "#%02x%02x%02x" % (61, 123, 141)   # bluish.
   elif obs.twi > 5. :
      fac = (obs.twi - 5) / 10.
      return "#%02x%02x%02x" % (int(51 + fac * 10) ,int(69 + fac * 54), int(113 + fac * 28))
   else :
      fac = (obs.twi+4) / 9. 
      return "#%02x%02x%02x" % (int(51 * fac), int(69 * fac), int(113 * fac))


class SkyDisplay(Toplevel) :

   sidt = 0.
#   xpixint = 800
#   ypixint = 700
   if TINYPIX :
      xpixint = 1000
      ypixint = 875
      bigfont = "helvetica, 16"
      midfont = "helvetica, 14"
      planetfont = "helvetica, 11"
      objfont = "helvetica, 10"
   else : 
      xpixint = 800
      ypixint = 700
      bigfont = "helvetica, 14"
      midfont = "helvetica, 12"
      planetfont = "helvetica, 9"
      objfont = "helvetica, 8"
   xpix = float(xpixint)
   ypix = float(ypixint)
   aspect = xpix / ypix
   xmid = 0.
   ymid = 0.
   halfwidthy = 0.88
   halfwidthx = halfwidthy * aspect
   pixperunit = ypix / (2. * halfwidthy)
   brlist = []
   magconst1 = 0.0015
   magslope = 0.0013
   magzpt = 4.8
   gridcolor = "#990000"
   xobjs = []
   yobjs = []
   isvisible = True

   def __init__(self, obs, objs) :

       # save default scale factors ... 
       self.halfwidthxfull = self.halfwidthx
       self.halfwidthyfull = self.halfwidthy
       self.pixperunitfull = self.pixperunit 

       # and away we go.
       Toplevel.__init__(self)
       self.w = Canvas(self, width=self.xpixint, height=self.ypixint)
       self.w.pack()
       self.readbright()
       self.redraw(obs, objs)
#       self.gethelp()
#       self.helpwin.withdraw()

   def gethelp(self) :
       self.helpwin = Toplevel(self)
       self.helpwin.title("Display commands:")
       f = Frame(self.helpwin)
       t = Text(f,width=40,height=15,foreground="green",background="black") 
     #  t,foregr="green"
     #  t.bckgr="black"
       coms = """Sky Display commands:
  left mouse: sel. nearest listed obj
  mid mouse : sel. these coords
        h,s : nearest plotted HR star 

          f : step forward in time
          b : step back in time
          n : set to now

        z,i : zoom in
          o : zoom out
        for i in range(1:5) : 
          p : pan to cursor
        r,1 : restore default scaling

          q : quit (hide this win).
"""
       t.insert('1.0',coms)
       t.pack()
       f.pack()
 
   def showhelp(self) :
       # print "lifting help window ... "
       self.helpwin.lift()

   def killhelp(self) :
       # print "hiding help window ... "
       self.helpwin.withdraw()

   def redraw(self, obs, objs) :

       if not self.isvisible : return

       sinlat = np.sin(obs.site.location.lat) 
       coslat = np.cos(obs.site.location.lat)

       skcolor = twilightcolor(obs,moonzero = 23., moonfac=8.)
       self.w.create_rectangle(0,0,self.xpixint,self.ypixint, fill=skcolor)  
                # fill with background sky color
       self.plotgrid(coslat,sinlat)
       self.drawclock(obs, self.xmid + 0.86 * self.halfwidthx, 
                           self.ymid + 0.80 * self.halfwidthy, 
                           0.1 * (self.pixperunitfull / self.pixperunit))
       self.putbanner(obs, self.xmid + 0.45 * self.halfwidthx, 
                           self.ymid + -0.9 * self.halfwidthy) 
       self.plotbright(obs)
       self.drawmoon(obs)

       self.drawsun(obs,coslat,sinlat)
       self.drawplanets(obs,coslat,sinlat) 
       self.drawobjs(obs,objs,coslat,sinlat)
       self.drawcurrent(obs,coslat,sinlat)
 
   def xytopix(self, x, y) :

      # This is bizarre.  These things seem to somehow promote to 
      # astropy quantities by the simple act of doing arithmetic.  Solution
      # appears to be to convert the dimensionless Quantities back to float.

      #print("self.xpix = ",self.xpix)
      #expr = 0.5 * self.xpix * (1. + (x - self.xmid) / self.halfwidthx)
      #print("xpix, type(xpix)",self.xpix,type(self.xpix))
      #print("xmid, type(xmid)",self.xmid,type(self.xmid))
      #print("halfwidthx, type(halfwidthx)",self.halfwidthx,type(self.halfwidthx))
      #print("expr: ", expr, "type(expr) ",type(expr))
      #print("float(expr): ", float(expr))
      #print("round(float(expr)):",round(float((expr))))

      xout = int(round(float(0.5 * self.xpix * (1. + (x - self.xmid) / self.halfwidthx))))
      yout = int(round(float(0.5 * self.ypix * (1. - (y - self.ymid) / self.halfwidthy))))
      xout = int(round(float(0.5 * self.xpix * (1. + (x - self.xmid) / self.halfwidthx))))
      yout = int(round(float(0.5 * self.ypix * (1. - (y - self.ymid) / self.halfwidthy))))
      return [xout, yout]

   def pixtoxy(self, xpixint, ypixint) :
      x = self.xmid + self.halfwidthx * (2 * xpixint / self.xpix - 1.)
      y = self.ymid + self.halfwidthy * (1. - 2. * ypixint / self.ypix) 
      return (x,y)

   def zoom(self, ev, zoomfac, obs, objs) :  # zoom in on cursor
      # needs to be fixed for MacOS Tkinter bug
      if is_MacOS :
          (cx, cy) = keylocation(ev) 
          xycent = self.pixtoxy(cx,cy)
      else :
          xycent = self.pixtoxy(ev.x,ev.y)
      self.xmid = xycent[0]
      self.ymid = xycent[1]
      self.halfwidthy = self.halfwidthy / zoomfac
      self.halfwidthx = self.halfwidthx / zoomfac
      self.pixperunit = self.pixperunit * zoomfac
      self.redraw(obs, objs)

   def zoomout(self, ev, zoomfac, obs, objs) :  # zoom out 

      if is_MacOS :
          (cx, cy) = keylocation(ev) 
          xycent = self.pixtoxy(cx,cy)
      else :
          xycent = self.pixtoxy(ev.x,ev.y)

      xycent = self.pixtoxy(ev.x,ev.y)
      self.xmid = xycent[0]
      self.ymid = xycent[1]
      self.halfwidthy = self.halfwidthy * zoomfac
      self.halfwidthx = self.halfwidthx * zoomfac
      self.pixperunit = self.pixperunit / zoomfac
      self.redraw(obs, objs)

   def zoomdef(self, ev, obs, objs) :  # zoom to original  
      self.xmid = 0.
      self.ymid = 0.
      self.halfwidthy = self.halfwidthyfull
      self.halfwidthx = self.halfwidthxfull
      self.pixperunit = self.pixperunitfull
      self.redraw(obs, objs)

   def pan(self, ev, obs, objs) :  # move center

      if is_MacOS :
          (cx, cy) = keylocation(ev) 
          xycent = self.pixtoxy(cx,cy)
      else :
          xycent = self.pixtoxy(ev.x,ev.y)

      self.xmid = xycent[0]
      self.ymid = xycent[1]
      self.redraw(obs, objs)
 
   def drawline(self, x1,y1,x2,y2,fillcolor="",lwidth=1.0) : # user coords
      xy1 = self.xytopix(x1,y1)
      xy2 = self.xytopix(x2,y2)
      # print "xy1, xy2",xy1,xy2
      self.w.create_line(xy1[0],xy1[1],xy2[0],xy2[1],fill=fillcolor,width=lwidth)

   def drawcircle(self,x,y,radius,fillcolor="",outlinecolor="black", dash = "") :

      xy = self.xytopix(x - radius, y + radius) # one corner
     # print "circle: xy = ",xy
      diam = 2. *  self.pixperunit * radius
      self.w.create_oval(xy[0],xy[1],xy[0]+diam,xy[1]+diam,
          fill=fillcolor,outline=outlinecolor, dash=dash)

   def puttext(self,x,y,text="blank",fillcolor="white",font="helvetica, 12") :
      xy = self.xytopix(x,y)
      self.w.create_text(xy[0],xy[1],text=text,fill=fillcolor,font=font)

   def putbanner(self, obs, x, y) :
      ut = obs.t.to_datetime()
      local = obs.t.to_datetime(timezone = obs.site.localtz)
      outstr =            "      Site : " + obs.site.name
      outstr = outstr + "\nLocal time : " + local.strftime("%a  %Y-%m-%d  %H:%M:%S")
      self.puttext(x,y,text=outstr,fillcolor="cyan",font=self.midfont)
 
   def drawclock(self, obs, x, y, radius) :
      
     # cald = jd2cal(obs.jd, obs.stdz, obs.use_dst) 
      cald = obs.t.to_datetime(timezone = obs.site.localtz)
    #  time = float(tp[3]) + float(tp[4]) / 60. + float(tp[5]) / 3600.
      time = cald.hour + cald.minute / 60. + cald.second / 3600.
    #  dst_inuse = False
    #  if obs.use_dst != 0 : 
    #     if abs(cald[7] / 24. - cald[9]) > 0.01 : dst_inuse = True

      if time >= 12. : pm = True 
      else : pm = False 
      
      while time > 12. :  time = time - 12.

      # Is it DST?  Python's datetime stuff makes it very difficult to find out.

      #  https://medium.com/@nqbao/python-timezone-and-daylight-savings-e511a0093d0

      dstoffset = obs.site.localtz.localize(datetime(cald.year,cald.month,cald.day,cald.hour,cald.minute,cald.second)).dst()
      if dstoffset == timedelta(hours=0) : dst = False
      else : dst = True

      # print "dstoffset, dst:",dstoffset,dst

      xy1 = self.xytopix(x - radius, y + radius)

      for i in range(0,12) : 
          ang = Angle((30. * i) * u.deg)
          # print "i, angle ",i,angle
          cosang = np.cos(ang)
          sinang = np.sin(ang)
          tickx1 = x + 0.9 * radius * cosang
          tickx2 = x + radius * cosang
          ticky1 = y + 0.9 * radius * sinang
          ticky2 = y + radius * sinang
          # print "tickx1 tickx2 ticky1 ticky2",tickx1,tickx2,ticky1,ticky2
          self.drawline(tickx1,ticky1,tickx2,ticky2,fillcolor="cyan")
      self.drawcircle(x,y,radius,outlinecolor="cyan")

      mins = 60. * (time - int(time)) 
      hrangle = 0.523599 * time
      minangle = 0.104719755 * mins

      self.drawline(x,y,
               x + 0.6 * radius * np.sin(hrangle),
               y + 0.6 * radius * np.cos(hrangle),fillcolor="cyan",lwidth=2.5)
      self.drawline(x,y,
               x + 0.83 * radius * np.sin(minangle),
               y + 0.83 * radius * np.cos(minangle),fillcolor="cyan")

      if pm == 0 : decstr = "AM " 
      else      : decstr = "PM " 

      if dst : decstr = decstr + obs.site.localtzdayabbrev
      else : decstr = decstr + obs.site.localtzabbrev

      self.puttext(x,y - 1.25 * radius,text=decstr,fillcolor="cyan")

   def skyproject(self, ha, dec, coslat, sinlat) :

       # supports legacy methods that require it.  Not expensive
       # for a few points, but the orthographic methods developed
       # in 2018 exploit numpy arrays for speed and elegance.

       # ha and dec are Angles.
       
       #print "dec = ",dec
       cosdec = np.cos(dec)
       x = cosdec * np.sin(ha)
       y = cosdec * np.cos(ha)
       z = np.sin(dec)

       #print "xyz in skyproject",x,y,z
 
       ypr = sinlat * y - coslat * z
       zpr = coslat * y + sinlat * z

       #print "zpr ypr",zpr,ypr
 
       zdist = np.arccos(zpr)
       r = np.tan(zdist / 2.)
       inground = np.sqrt(x**2 + ypr**2)

       #print "zdist, r, inground",zdist, r, inground
       if r == 0. : 
          return (0.,0.)
       elif sinlat > 0 :
          return (r * (x/inground), -1. * r * (ypr / inground))
       else :
          return (-1. * r * (x / inground), r * (ypr / inground))
 
#    REPLACED this pixtocelest with a matrix-based version that should
#    be faster and is certainly cleaner.
#    def pixtocelest(self, obs, xpixin, ypixin) :
#       x = self.xmid + self.halfwidthx * (2. * xpixin / self.xpix - 1.)
#       y = self.ymid + self.halfwidthy * (1 - 2. * ypixin / self.ypix)
#       mod = np.sqrt(x**2 + y**2) 
#       if obs.site.location.lat < 0. * u.deg : 
#          x = -1 * x
#          y = -1 * y
#       
#       alt = (np.pi / 2) - 2 * np.arctan(mod)
#       az = np.arctan2(-1. * x, y)
#       zt = np.sin(alt)
#       xt = -1 * np.cos(alt) * np.sin(az)
#       yt = -1 * np.cos(alt) * np.cos(az)
#       coslat = np.cos(obs.site.location.lat) # / cooconsts.DEG_IN_RADIAN)
#       sinlat = np.sin(obs.site.location.lat) # / cooconsts.DEG_IN_RADIAN)
#       yc = yt * sinlat + zt * coslat
#       zc = -1 * yt * coslat + zt * sinlat
# 
#       mod = np.sqrt(xt ** 2 + yc ** 2 + zc ** 2) 
#       xt = xt / mod
#       yc = yc / mod
#       zc = zc / mod 
# 
#       xy = np.sqrt(xt**2 + yc**2) 
#       pseudoha = np.arctan2(yc, xt) # * cooconsts.HRS_IN_RADIAN
# 
#       ra = pseudoha - Angle(6. * u.hour) + obs.lst
#       dec = np.arcsin(zc) # * cooconsts.DEG_IN_RADIAN
# 
#       # cel = celest([ra,dec,obs.equinox])
#       # cel.quickpr()
# 
#       celraw = SkyCoord(ra,dec,frame = obs.nowfr) 
#       if obs.nowfr.name == "precessedgeocentric" :
#           cel = celraw.transform_to('gcrs')
#       else :
#           cel = celraw.transform_to('icrs')
#       return cel
  
   def pixtocelest(self, obs, xpixin, ypixin, return_xyz = False) :

      # From the pixels, get the orthographic coordinates.

      # print("xpixin, ypixin ... ",xpixin,ypixin)

      # this established that tkinter on macos does no return the 
      # event coordinates correctly if you type a letter.

      x = self.xmid + self.halfwidthx * (2. * xpixin / self.xpix - 1.)
      y = self.ymid + self.halfwidthy * (1 - 2. * ypixin / self.ypix)

      # find the topocentric XYZ coordinates corresponding.
      rad = np.sqrt(x ** 2 + y ** 2)   # this is also tan(zdist/2)
      zdist = 2 * np.arctan(rad)       # zenith distance.
      #print "xpixin, ypxin, zdist",xpixin,ypixin,zdist
      vert = np.cos(zdist)  
      east = -1. * x * (1 + vert)   # names are correct only in N hemi.
      south = -1. * y * (1 + vert) 
      #print "vert south east ",vert,south,east
      #print "modulus:", np.sqrt(vert ** 2 + south ** 2 + east ** 2)
      topovec = np.array([[south],[east],[vert]])  

      # now rotate from topocentric back to icrs by multiplying by
      # inverse transformation matrix!  The southern hemisphere signs
      # are included in the original icrs2topoxyz, and the resulting
      # celvec is correct.

      celvec = np.linalg.inv(obs.icrs2topoxyz).dot(topovec)
      # print "celvec:",celvec

      if return_xyz : return celvec

      ratmp = np.arctan2(celvec[1,0],celvec[0,0])
      dectmp = np.arcsin(celvec[2,0])

     # testcel = SkyCoord(ratmp,dectmp,unit=(u.radian,u.radian),frame='icrs')
     # print("celest: ",testcel.ra.to_string(unit=u.hourangle),testcel.dec.to_string(unit=u.deg))

      return SkyCoord(ratmp,dectmp,unit=(u.radian,u.radian),frame='icrs')

   def readbright(self) : 
       (self.bright2000, self.brightmags, self.brightcolors, 
             self.brightnames) = getbrightest()
       # should automatically read from the file included from the package,
       # so file name is obsolete.
       #      self.brightnames) = getbrightest("cartesian_bright.dat")
       # print "read %d bright stars." % (len(self.brightcolors)) 
   
   def plotbright(self, obs) :
       
       (projectedx, projectedy) = skyproject(obs.icrs2topoxyz, self.bright2000)
       # isnorth = obs.site.location.lat > 0. * u.deg
       for i in range(0,len(self.brightmags)) : 
           dotrad = self.magconst1 + self.magslope * (self.magzpt - self.brightmags[i])
           # rgb = "#%02x%02x%02x" % (self.brightcolors[i][0],self.brightcolors[i][1],self.brightcolors[i][2])
        #   if isnorth :
           # flipping of the display is now included in the icrs2topoxyz matrix.
           self.drawcircle(projectedx[i],projectedy[i],dotrad,fillcolor=self.brightcolors[i],
                 outlinecolor = "") 
        #   else :   # flip map south of equator
        #      self.drawcircle(-1. * projectedx[i],-1. * projectedy[i],dotrad,fillcolor=self.brightcolors[i],
        #         outlinecolor = "") 
 
   def nearestbright(self,xyzin) :  # here b/c brlist is local.
       # xyzin is the cartesian unit vector of the icrs coordinates the user marked.

       # insanely clever use of numpy to find index of the bright star that
       # is closest in 3-d vector space; numpy just 'does the right thing'.
       # runtime is imperceptible.

       diffvecs = self.bright2000 - xyzin     # subtracts columnwise
       squares = diffvecs * diffvecs          # squares
       sumsq = np.sum(squares, axis=0)        # sums along columns
       # nearestargument = sumsq.argmin()
       # print(self.bright2000[0][nearestargument],self.bright2000[1][nearestargument],self.bright2000[2][nearestargument])
       return sumsq.argmin()                  # finds index of minimum!
 
   def plotgrid(self, coslat, sinlat) :
       # plots the ha grid and equator.
       # first the equator

       ang0 = Angle(0. * u.hour)
       hr6 = Angle(6. * u.hour)
       min6h = -1. * hr6 
       min20deg = Angle(-20. * u.deg)
       deg10 = Angle(10. * u.deg)
       deg40 = Angle(40. * u.deg)
       
       pt0 = self.skyproject(min6h,ang0,coslat,sinlat)
       pt1 = self.skyproject(ang0,ang0,coslat,sinlat) # ha = 0, dec = 0
       pt2 = self.skyproject(hr6,ang0,coslat,sinlat)
       ctr = circlecenter(pt0[0],pt0[1],pt1[0],pt1[1],pt2[0],pt2[1])
       self.drawcircle(ctr[0],ctr[1],ctr[2],outlinecolor=self.gridcolor)
       self.drawcircle(0.,0.,1.,outlinecolor=self.gridcolor) # horizon
       for i in (-4.,-2.,2.,4.,6.) : # every two hrs of HA
           hacirc = Angle(i * u.hour)
           pt0 = self.skyproject(hacirc,min20deg,coslat,sinlat)
      #     self.drawcircle(pt0[0],pt0[1],0.004,fillcolor="red")
           pt1 = self.skyproject(hacirc,deg10,coslat,sinlat)
      #     self.drawcircle(pt1[0],pt1[1],0.004,fillcolor="red")
           pt2 = self.skyproject(hacirc,deg40,coslat,sinlat)
      #     self.drawcircle(pt2[0],pt2[1],0.004,fillcolor="red")
           ctr = circlecenter(pt0[0],pt0[1],pt1[0],pt1[1],pt2[0],pt2[1])
           self.drawcircle(ctr[0],ctr[1],ctr[2],outlinecolor=self.gridcolor)
       self.drawline(0.,1.,0.,-1.,fillcolor=self.gridcolor) # meridian           
       if sinlat > 0. : 
          self.puttext(0.,0.84,text="N",fillcolor=self.gridcolor,font=self.bigfont)
          self.puttext(0.,-0.84,text="S",fillcolor=self.gridcolor,font=self.bigfont)
          self.puttext(-0.95,0.,text="E",fillcolor=self.gridcolor,font=self.bigfont)
          self.puttext(0.95,0.,text="W",fillcolor=self.gridcolor,font=self.bigfont)
       else :
          self.puttext(0.,0.84,text="S",fillcolor=self.gridcolor,font=self.bigfont)
          self.puttext(0.,-0.84,text="N",fillcolor=self.gridcolor,font=self.bigfont)
          self.puttext(-0.95,0.,text="W",fillcolor=self.gridcolor,font=self.bigfont)
          self.puttext(0.95,0.,text="E",fillcolor=self.gridcolor,font=self.bigfont)
       # dashed circles at 2, 3, and 4 airmasses - true aimasses.
       self.drawcircle(0.,0.,0.57875,outlinecolor=self.gridcolor,dash=(2,4))
       self.drawcircle(0.,0.,0.70963,outlinecolor=self.gridcolor,dash=(2,4))
       self.drawcircle(0.,0.,0.77876,outlinecolor=self.gridcolor,dash=(2,4))
  
# moon and planets.
 
   def roughradectoxy(self,center,x,y,scale) :
     # takes a central SkyCoord, and a set of x y offsets along
     # ra and dec respectively, together with a scale factor (how many 
     # degrees to a unit), and delivers a list of ras and decs corresponding
     # to x and y.  Uses crude "graph paper" approximation, not tangent plane.
    
        cosdec = np.cos(center.dec) #  / cooconsts.DEG_IN_RADIAN) 
        raout = []
        decout = []
        for i in range(0,len(x)) :
            raout.append(center.ra.hour + scale * x[i] / (15. * cosdec))
            decout.append(center.dec.deg + scale * y[i])
      
        return [raout,decout]
 
       
   def moonedges(self,obs) :
      
     # delivers a list of two lists -- one is ra,decs of a set of points on the
     # illuminated lunar limb, the other is the ra, decs of a set of points on 
     # the terminator. 
     # revision uses a different approach -- compute
     # the limb and terminator in a system which is aligned with the
     # line of cusps, and then rotate the resulting arrays into 
     # position using matrix multiplication.  Much cleaner.  
    
    
         mooninfo = mooncirc(obs) 

         #print "mooninfo back from mooncirc:" 
         # for m in mooninfo : print m
       
         moonsun = mooninfo[1]
#         while moonsun < 0. :
#            moonsun = moonsun + cooconsts.TWOPI
         #print "moonsun = %f" % moonsun.deg
         #print "back from printing moonsun."
    
         limbx = []
         limby = []
         termx = []
         termy = []
       
         drange = np.pi / 10.  
         for p in np.arange(0.,np.pi+0.001,drange) :
#          # set up limb and terminator with x-axis along cusp line,
#          # limb of moon in top half-plane
             limbx.append(np.cos(p))
             limby.append(np.sin(p))
             termx.append(np.cos(p))  # need a 2nd copy later
             termy.append(np.sin(p) * np.cos(moonsun))
#          # cos(moonsun) takes care of phase angle 
       
         pa = mooninfo[3]  # pa from moon to sun 
#    #   print "moonsun %f = %f deg" % (moonsun, moonsun * _skysub.DEG_IN_RADIAN)
#    #   print "pa = %f = %f deg" % (pa, pa * _skysub.DEG_IN_RADIAN)
         turnmoon = [[np.cos(pa),np.sin(pa)],[-1. * np.sin(pa),np.cos(pa)]]
#    #   print "turnmoon = ", turnmoon
# 
#       # rotation matrix to turn moon to appropriate pa
#    
         limb = [limbx,limby]
         #print "limb before = ",limb
# 
#       # this is easy!  Just lay the x and y on top of each other in a 
#       # matrix, and ... 
#    
#       # multiply them and 
         limb = np.matrix(turnmoon).dot(np.matrix(limb))  # translated orig to matrix mult
#       #print "limb after = ",limb
  #       # strip out the x and y as separate arrays and 
          # limbx = limb[0]
         limbx = limb.A[0:]   # A is the "array" method
         limby = limb.A[1:]

#       # do the same for the terminator, and finally 
         term = [termx,termy]
         term = np.matrix(turnmoon).dot(np.matrix(term))
#       # print "term = ",term
         termx = term.A[0:]
         termy = term.A[1:]
# 
        # print 'termx ',termx
#    
#       # Now, how large to draw the moon symbol?  Unfortunately,
#       # scaling this with zenith distance requires us to know the 
#       # zenith dist of the moon right here ... 
#    
         coszover2 = np.cos((Angle(90. * u.deg) - obs.moonaltit) / 2.) 
         moonscale = 3. * coszover2
#    #   print "scale ... %f" % moonscale
#        
#       #print "limbx = ",limbx
#       #print "limby = ",limby
#       #print "termx = ",termx
#       #print "termy = ",termy
#       
         limbradec = self.roughradectoxy(mooninfo[0],limbx[0],limby[0],moonscale)
         termradec = self.roughradectoxy(mooninfo[0],termx[0],termy[0],moonscale)
# 
        # print "limbradec = ",limbradec
        # print "termradec = ",termradec
# 
         xlimb = []
         ylimb = []
         xterm = []  # x of terminator, not the computer terminal type!
         yterm = [] 
# 
         coslat = np.cos(obs.site.location.lat) 
         sinlat = np.sin(obs.site.location.lat) 
#      
         for i in range(0,len(limbradec[0])) :
              ha = obs.lst - Angle(limbradec[0][i] * u.hour)
#          #print 'ha = ',ha
#          #print 'i ',i,' limbradec[1][i] ',limbradec[1][i]
              xy = self.skyproject(ha,(limbradec[1][i] * u.deg),coslat,sinlat) 
              xlimb.append(xy[0])
              ylimb.append(xy[1])
         for i in range(0,len(termradec[0])) :
              ha = obs.lst - Angle(termradec[0][i] * u.hour)
              xy = self.skyproject(ha,Angle(termradec[1][i] * u.deg),coslat,sinlat) 
              xterm.append(xy[0])
              yterm.append(xy[1])
#    
#       #print "returning."
         return [xlimb,ylimb,xterm,yterm]
   
   def drawmoon(self, obs) :
       stuff = self.moonedges(obs)
       for i in range(0,len(stuff[0])-1) :
           self.drawline(stuff[0][i],stuff[1][i],
                             stuff[0][i+1],stuff[1][i+1],
                              fillcolor="yellow")
          # print i,stuff[0][i],stuff[1][i],"->",stuff[0][i+1],stuff[1][i+1]
       for i in range(0,len(stuff[2])-1) :     
          self.drawline(stuff[2][i],stuff[3][i],
                             stuff[2][i+1],stuff[3][i+1],
                            fillcolor="yellow")
 
   def drawsun(self, obs, coslat, sinlat) :  # sun symbol.
       suncent = self.skyproject(obs.sunha,obs.sunpos.dec,coslat,sinlat)
       self.drawcircle(suncent[0],suncent[1],radius=0.004,fillcolor="yellow",outlinecolor="yellow")
       self.drawcircle(suncent[0],suncent[1],radius=0.025,outlinecolor="yellow")

   def drawcurrent(self,obs,coslat,sinlat) :
#       print "drawcurrent - obs.celnow = ",obs.celnow
#       print "obs.hanow, obs.celnow.dec", obs.hanow, obs.celnow.dec
       xy = self.skyproject(obs.hanow, obs.celnow.dec,coslat,sinlat)
#       print "xy",xy
       xyp = self.xytopix(xy[0],xy[1])
       self.w.create_rectangle(xyp[0]+8,xyp[1]+8,xyp[0]-8,xyp[1]-8,outline="#55AAFF")
# #
   def drawplanets(self,obs,coslat,sinlat) :
   
       planx = []
       plany = []    # to hand back for a cursor find ...
       for p in obs.planetdict.keys() :  
           o = obs.planetdict[p]
           ha = obs.lst - o.ra
           xy = self.skyproject(ha,o.dec,coslat,sinlat)
#           print p
#           print "ha, o.dec",ha,o.dec
#           print "xy",xy
           planx.append(xy[0].value) 
           plany.append(xy[1].value) 
           dotrad = self.magconst1 + self.magslope * (self.magzpt - obs.planetmags[p])
           self.drawcircle(xy[0].value,xy[1].value,dotrad,fillcolor="yellow",outlinecolor="yellow")
           self.puttext(xy[0].value,xy[1].value - 0.028, text = p.capitalize(), 
                 fillcolor="yellow", font = self.planetfont)
       return [planx,plany]
 
   def drawobjs(self, obs, objs, coslat, sinlat) :
       for o in objs.keys() :
          # ugly ugly.  Should store Cartesians with the object list.
          celnow = objs[o].transform_to(obs.nowfr)
          ha = obs.lst - celnow.ra
          xy = self.skyproject(ha,objs[o].dec,coslat,sinlat)
          self.puttext(xy[0],xy[1],text=o,fillcolor="#00FF00",font=self.objfont)

 
# Airmass window part 

class AirmDisplay(Toplevel) :

   # adapted from the Java "AirmassDisplay".

   # This draws a map of the sky without using any graphics package --
   # instead it uses primitive python drawing tools.  It works nicely,
   # though.

   def __init__(self,obs) : #  , objs) :
      if TINYPIX :
         self.xpix = 960
         self.ypix = 600
      else : 
         self.xpix = 800
         self.ypix = 500
      self.gridcolor = "#990000"
      self.isvisible = True

      self.jdstart = 2450000.  # useful to hand out.

      self.xlobord = 0.
      self.xhibord = 1.
      self.ylobord = 0.
      self.yhibord = 1.   # user coordinates of full frame edge

      self.xlo = 0.       # decimal hours at start of night
      self.xhi = 1.
      self.ylo = 3.2
      self.yhi = 0.9
      self.loxvppix = 100.  # pixel value of left viewport edge
 
      # and away we go.
      Toplevel.__init__(self)
      self.w = Canvas(self, width=self.xpix, height=self.ypix)
      self.w.pack()

      self.previousobs = None

   def redraw(self, obs, objname = "", is_visible = True) : # , objs) :
       # date = _skysub.new_date_time()

       # skip the redraw if things haven't changed.  It's expensive.

       if self.previousobs != None : 
           # don't re-draw the airmass window unless it's up and has changed.
           # The winddow is not redrawn if the object is in the same position, and 
           # the local # sidereal time at midnight to very high accuracy.  That
           # will basically only occur if it's the same night and the same timee.
           if obs.cel_J2000.separation(self.previousobs.cel_J2000) < 0.01 * u.deg  and \
              abs((obs.lstmid - self.previousobs.lstmid).deg) < 0.0001 :
               # print "No airmass redraw - same."
               return 

       # start by painting the the main underlying rectangle -- black background.

       self.w.create_rectangle(0,0,self.xpix,self.ypix, fill="black")  

       # and then quit if the window isn't visible.  Not redrawing saves a lot
       # of execution time, and by blanking the window first, avoids any impression 
       # that the info in the window is up-to-date when you pop it. 

       if not self.isvisible : return

       ao = deepcopy(obs)  # ao = "airmass observation"

       xwidth = float(self.xpix)
       yheight = float(self.ypix)
       xvplo, xvphi, yvplo, yvphi = (0.08,0.88,0.15,0.90) 
       self.loxvppix = xwidth * xvplo  # pixel value of left viewport edge
 
       #print "ping 0"
       
       # get jdstart and jdend (sunset and sunrise)
       self.jdstart = ao.tsunset.jd
       jdend = ao.tsunrise.jd  # these are floats.
       self.timespan = (jdend - self.jdstart)  * 24.
       if ao.twilightflag == 0 : 
           darkhrs = (ao.tmorntwi.jd - ao.tevetwi.jd) * 24.
       else : 
           darkhrs = 0.

       # print 'jdstart, jdend, timespan',self.jdstart,jdend,self.timespan
       # print "darkhrs = ",darkhrs

       ao.settime(ao.tsunset) 
       ao.computesky()  # this does the "redo_coords" precession etc 

       # compute time 
       # calao = ao.caldat(stdz=ao.stdz,use_dst=ao.use_dst)  # should be local time
       # print calao

       # xlo is decimal hours local time of day at sunset
       dt = ao.tsunset.to_datetime(timezone = ao.site.localtz)
       self.xlo = dt.hour + dt.minute / 60. + dt.second / 3600.
       #print "decimal local sunset time is ",self.xlo
       self.xhi = self.xlo + self.timespan
       #print "self.xhi = ",self.xhi

       # find and save ut offset 
#       utoffset = int(math.floor(ao.decimalUT() - self.xlo + 0.0001))
       deltat = dt.utcoffset()   # a deltatime
       utoffset = 24. * deltat.days + deltat.seconds / 3600.  # in hours

       #print "ping1"

       xendevetwi = self.xlo + 24. * (ao.tevetwi.jd - ao.tsunset.jd) 
       xbegmorntwi = self.xlo + 24. * (ao.tmorntwi.jd - ao.tsunset.jd)
       fracx = xvphi - xvplo
       self.xlobord = self.xlo - xvplo * self.timespan / fracx
       self.xhibord = self.xhi + (1 - xvphi) * self.timespan / fracx

       span = self.yhi - self.ylo
       fracy = yvphi - yvplo
       self.ylobord = self.ylo - yvplo * span / fracy
       self.yhibord = self.yhi + (1. - yvphi) * span / fracy

       #print "xlobord, xhibord ",self.xlobord,self.xhibord," yloboard,yhibord ",self.ylobord,self.yhibord

       # fill in the background rectangles -- have to approximate a gradient
       # for twilight since canvas doesn't have it ...

       twilight_dt = 10. / 1440.   # 10-minute steps.

       jdlo = self.jdstart
       jdhi = jdlo + twilight_dt
       x_low = self.xlo + (jdlo - self.jdstart) * 24.  # low edge in hrs
       while jdhi < ao.tevetwi.jd :
         jdmid = (jdlo + jdhi) / 2.
         ao.settime(jdmid) 
         ao.computesky(redo_coords = False)
         ao.computesunmoon() 
         twicolor = twilightcolor(ao)
         x_high = self.xlo + (jdhi - self.jdstart) * 24.  # high edge in hrs
         xyll = self.xytopix(x_low,self.ylo) # lower left
         xyur = self.xytopix(x_high,self.yhi) # upper right
         self.w.create_rectangle(xyll[0],xyll[1],xyur[0],xyur[1], fill = twicolor, outline=twicolor)
         
         x_low = x_high   # advance the loop
         jdlo = jdhi
         jdhi = jdhi + twilight_dt

       #print "eve done"

       # now let's make grey rectangles for lunar sky brightness! 
       # Commented out for now -- leaving the dark part black.

       # does the moon rise or set during the night?  Find which occurs closest
       # to midnight
# 
#        moon_dt = 30. / 1440.   # half-hour
# 
#        rise_gap = abs(ao.tmoonrise - ao.tcent)
#        set_gap = abs(ao.jdmoonset - ao.jdcent)
#        if rise_gap < set_gap : 
#          #  print "start coloring at moonrise."
#           if ao.tmoonrise > ao.tevetwi :
#              jdlo = ao.tmoonrise.jd
#           else :
#              jdlo = ao.tevetwi.jd
#        else :
#          # print "start coloring at evening twilight."
#           jdlo = ao.jdevetwi
# 
#        jdhi = jdlo + moon_dt
#        x_low = self.xlo + (jdlo - self.jdstart) * 24.
#
#       while jdhi < ao.tmoonset.jd and jdlo < ao.tmorntwi.jd :
#         jdmid = (jdlo + jdhi) / 2.
#         ao.setut(jdmid) 
#         ao.computesky()
#         ao.computesunmoon() 
#         twicolor = twilightcolor(ao,moonzero=23, moonfac = 10.)
#         x_high = self.xlo + (jdhi - self.jdstart) * 24.  # high edge in hrs
#         xyll = self.xytopix(x_low,self.ylo) # lower left
#         xyur = self.xytopix(x_high,self.yhi) # upper right
#         self.w.create_rectangle(xyll[0],xyll[1],xyur[0],xyur[1], fill = twicolor, outline=twicolor)
#         
#         x_low = x_high   # advance the loop
#         jdlo = jdhi
#         jdhi = jdhi + moon_dt 



       # draw the morning twilight

       jdlo = ao.tmorntwi.jd
       jdhi = jdlo + twilight_dt
       x_low = self.xlo + (jdlo - self.jdstart) * 24.

       while jdhi < ao.tsunrise.jd :
         jdmid = (jdlo + jdhi) / 2.
         ao.settime(jdmid) 
         ao.computesky(redo_coords = False)
         ao.computesunmoon() 
         twicolor = twilightcolor(ao)
         x_high = self.xlo + (jdhi - self.jdstart) * 24.  # high edge in hrs
         xyll = self.xytopix(x_low,self.ylo) # lower left
         xyur = self.xytopix(x_high,self.yhi) # upper right
         self.w.create_rectangle(xyll[0],xyll[1],xyur[0],xyur[1], fill = twicolor, outline=twicolor)
         
         x_low = x_high
         jdlo = jdhi
         jdhi = jdhi + twilight_dt

       #print "ping2"

       # draw the frame

       self.drawline(self.xlo,self.ylo,self.xhi,self.ylo,fillcolor="white")
       self.drawline(self.xhi,self.ylo,self.xhi,self.yhi,fillcolor="white")
       self.drawline(self.xhi,self.yhi,self.xlo,self.yhi,fillcolor="white")
       self.drawline(self.xlo,self.yhi,self.xlo,self.ylo,fillcolor="white")

       # draw tick marks and grid along X, for both local and UT, and label.

       hrint = int(self.xlo)   # integer hour where tick mark falls.

       majorxticksize = 0.04
       
       localtimelabely = 3.3
       utlabely = 3.43

       timelabelx = self.xlo - 0.05 * self.timespan

       self.puttext(timelabelx, localtimelabely,text = "Local:",fillcolor="white") 
       self.puttext(timelabelx, utlabely,text = "UT:",fillcolor="white") 

       while hrint < self.xhi :
          if hrint > self.xlo :
             # grid line 
             self.drawline(hrint,self.ylo,hrint,self.yhi,fillcolor=self.gridcolor,dash=(2,3))
             # ticks
             self.drawline(hrint,self.ylo,hrint,self.ylo - majorxticksize,fillcolor="white")
             self.drawline(hrint,self.yhi,hrint,self.yhi + majorxticksize,fillcolor="white")
             utint = hrint - utoffset
             while utint >= 24 : utint = utint - 24
             while utint < 0 : utint = utint + 24
             if hrint >= 24 : localout = hrint - 24
             else : localout = hrint 
             self.puttext(hrint, localtimelabely, text = "%02d" % localout, fillcolor="white")
             self.puttext(hrint, utlabely, text = "%02d" % utint, fillcolor="white")

          hrint = hrint + 1

       # draw tick marks on Y, and grid, and label.

       airmlabelx = self.xlo - 0.03 * self.timespan
       xmajorticklength = 0.01 * self.timespan

       airm = 1.0
       delta_airm = 0.5
       while airm < self.ylo :
          self.drawline(self.xlo,airm,self.xhi,airm,fillcolor=self.gridcolor,dash=(2,3))
          self.drawline(self.xlo,airm,self.xlo+xmajorticklength,airm,fillcolor="white")
          self.drawline(self.xhi,airm,self.xhi-xmajorticklength,airm,fillcolor="white")
          self.puttext(airmlabelx,airm,text="%3.1f" % airm,fillcolor="white")
          airm = airm + delta_airm 

       # And now that axes are all set up, plot the object's airmass 
       # through the night ... 

       #print "ping3"

       self.plotairmass(ao)

       # add a title banner

       bannery = self.yhi - 0.08
       bannerx = (self.xhi + self.xlo) / 2.
       if objname != '' :
          self.puttext(bannerx, bannery, "%s :  %s  %s  (J2000)" % (objname.strip(), 
               obs.cel_J2000.ra.to_string(unit = u.hourangle, sep = ':',precision = 2, pad = True),
               obs.cel_J2000.dec.to_string(sep = ':',precision = 1, pad = True, alwayssign = True)))
       else :
          self.puttext(bannerx, bannery, "%s  %s  (J2000)" % (
               obs.cel_J2000.ra.to_string(unit = u.hourangle, sep = ':',precision = 2, pad = True),
               obs.cel_J2000.dec.to_string(sep = ':',precision = 1, pad = True, alwayssign = True)))
 
       # And add a vertical line at the program time.  "obs" has been left alone through 
       # this whole process.
       # On second thought, DON'T do this because we're only redrawing the window if the
       # night and object have changed.  Redrawing the window when the time has shifted by
       # (say) an hour is expensive.
       #  xobs = self.xlo + (obs.t.jd - self.jdstart) * 24.
       #  self.drawline(xobs,self.ylo,xobs,self.yhi,fillcolor="grey",dash=(2,2))

       #print "ping4"

       self.previousobs = deepcopy(obs)   # last time it was redrawn.

       #print "ping5"

   def xytopix(self, x, y) :   # user coordinates to pixel coordinates
      retx = self.xpix * (x - self.xlobord) / (self.xhibord - self.xlobord)
      rety = self.ypix * (1. - (y - self.ylobord) / (self.yhibord - self.ylobord))
   #   print "ypix, y, ylobord, yhibord",self.ypix, y, self.ylobord, self.yhibord
   #   print "frac: ",(y - self.ylobord) / (self.ylobord - self.yhibord)
   #   print x, y, " -> ",retx,rety
      return [retx, rety]

   def xpixtojd(self,xevent) :  
      hr = float(xevent - self.loxvppix) * (self.xhibord - self.xlobord)/ float(self.xpix) 
      jdout = self.jdstart + hr / 24.
      return jdout
 
   def drawline(self, x1, y1, x2, y2,fillcolor="",lwidth=1.0, dash = "") : # user coords 
      xy1 = self.xytopix(x1,y1)
      xy2 = self.xytopix(x2,y2)
      # print xy1,xy2
      self.w.create_line(xy1[0],xy1[1],xy2[0],xy2[1],fill=fillcolor,width=lwidth,dash = dash)

   # def puttext(self,x,y,text="blank",fillcolor="white",font="helvetica, 14") :
   def puttext(self,x,y,text="blank",fillcolor="white",font="helvetica, 12") :
      # print "text ",text,"x,y ",x,y," -> ",
      xy = self.xytopix(x,y)
      # print xy[0],xy[1]
      self.w.create_text(xy[0],xy[1],text=text,fill=fillcolor,font=font)

   def plotairmass(self,o,fillcolor="white", largehacolor="red") : 

      dt = 0.0139 # timestep in days; about 20 minutes.
      jd = o.tsunset.jd

      o.settime(jd)
      o.computesky()

      x1 = self.xlo # + (o.t.jd - o.tsunset.jd) * 24. 
      y1 = o.airmass
      xy1 = self.xytopix(x1,y1)

#      while o.airmass > 0. and o.airmass < self.ylo and abs(o.hanow.hour) < 6. and jd < o.tsunrise.jd : 
#      print("in plotairmass, o.tsunset.jd, o.tsunrise.jd:",o.tsunset.jd,o.tsunrise.jd)
      while jd < o.tsunrise.jd : 
         jd = jd + dt 
         o.settime(jd)
         o.computesky(redo_coords = False)  # don't precess every time, for speed.
#         print("in plotairmass loop, jd  airmass",jd,o.airmass)
         x2 = self.xlo + (jd - o.tsunset.jd) * 24.
         y2 = o.airmass
         xy2 = self.xytopix(x2,y2)
         if o.airmass > 0. and o.airmass < 3. : # and abs(o.hanow.hour) < 6. :
            if abs(o.hanow.hour) < 6. :
                self.w.create_line(xy1[0],xy1[1],xy2[0],xy2[1],fill=fillcolor)
            else :
                self.w.create_line(xy1[0],xy1[1],xy2[0],xy2[1],fill=largehacolor)
         xy1 = xy2[:]

       #  print jd, o.airmass,xy1
      
# I copied all this stuff into this method to plot the moon, but I haven't implemented it yet.
# This needs to be put in a standalone function to be called by both skydisplay and airmass.
 
#    def moonedges(self,obs) :
#      
#     # delivers a list of two lists -- one is ra,decs of a set of points on the
#     # illuminated lunar limb, the other is the ra, decs of a set of points on 
#     # the terminator. 
#    # revision uses a different approach -- compute
#    # the limb and terminator in a system which is aligned with the
#    # line of cusps, and then rotate the resulting arrays into 
#    # position using matrix multiplication.  Much cleaner.  
#    
#    
#       mooninfo = mooncirc(obs) 
# 
#    #   print "ok, mooninfo = ", mooninfo
#       
#    #   pgbeg("/xterm")
#    #   pgsvp(0.1,0.9,0.1,0.9)
#    #   pgwnad(1.,-1.,-1.,1.)
#       
#       moonsun = mooninfo[2]
#       while moonsun < 0. :
#          moonsun = moonsun + _skysub.TWOPI
#    #   print "moonsun = %f" % moonsun
#    
#       limbx = []
#       limby = []
#       termx = []
#       termy = []
#       
#       drange = _skysub.PI / 10.  
#       for p in numpy.arange(0.,_skysub.PI+0.001,drange) :
#          # set up limb and terminator with x-axis along cusp line,
#          # limb of moon in top half-plane
#          limbx.append(cos(p))
#          limby.append(sin(p))
#          termx.append(cos(p))  # need a 2nd copy later
#          termy.append(sin(p) * cos(moonsun))
#          # cos(moonsun) takes care of phase angle 
#       
#       pa = mooninfo[4]  # pa from moon to sun 
#    #   print "moonsun %f = %f deg" % (moonsun, moonsun * _skysub.DEG_IN_RADIAN)
#    #   print "pa = %f = %f deg" % (pa, pa * _skysub.DEG_IN_RADIAN)
#       turnmoon = [[cos(pa),sin(pa)],[-1. * sin(pa),cos(pa)]]
#    #   print "turnmoon = ", turnmoon
# 
#       # rotation matrix to turn moon to appropriate pa
#    
#       limb = [limbx,limby]
#       #print "limb before = ",limb
# 
#       # this is easy!  Just lay the x and y on top of each other in a 
#       # matrix, and ... 
#    
#       # multiply them and 
#       limb = numpy.matrix(turnmoon) * numpy.matrix(limb)
#       #print "limb after = ",limb
#       # strip out the x and y as separate arrays and 
#       # limbx = limb[0]
#       limbx = limb.A[0:]   # A is the "array" method
#       limby = limb.A[1:]
#    
#       # do the same for the terminator, and finally 
#       term = [termx,termy]
#       term = numpy.matrix(turnmoon) * numpy.matrix(term)
#       # print "term = ",term
#       termx = term.A[0:]
#       termy = term.A[1:]
# 
#       # print 'termx ',termx
#    
#       # Now, how large to draw the moon symbol?  Unfortunately,
#       # scaling this with zenith distance requires us to know the 
#       # zenith dist of the moon right here ... 
#    
#       coszover2 = cos((90. - obs.altmoon) / (2. * _skysub.DEG_IN_RADIAN))
#       moonscale = 3. * coszover2
#    #   print "scale ... %f" % moonscale
#        
#       #print "limbx = ",limbx
#       #print "limby = ",limby
#       #print "termx = ",termx
#       #print "termy = ",termy
#       
#       limbradec = self.roughradectoxy(mooninfo[0],mooninfo[1],limbx[0],limby[0],moonscale)
#       termradec = self.roughradectoxy(mooninfo[0],mooninfo[1],termx[0],termy[0],moonscale)
# 
#       #print "limbradec = ",limbradec
#       #print "termradec = ",termradec
# 
#       xlimb = []
#       ylimb = []
#       xterm = []  # x of terminator, not the computer terminal type!
#       yterm = [] 
# 
#       coslat = cos(obs.lat / _skysub.DEG_IN_RADIAN)
#       sinlat = sin(obs.lat / _skysub.DEG_IN_RADIAN) # wasteful but trivial
#      
#       for i in range(0,len(limbradec[0])) :
#          ha = obs.sidereal.val - limbradec[0][i]
#          #print 'ha = ',ha
#          #print 'i ',i,' limbradec[1][i] ',limbradec[1][i]
#          xy = self.skyproject(ha,limbradec[1][i],coslat,sinlat) 
#          xlimb.append(xy[0])
#          ylimb.append(xy[1])
#       for i in range(0,len(termradec[0])) :
#          ha = obs.sidereal.val - termradec[0][i]
#          xy = self.skyproject(ha,termradec[1][i],coslat,sinlat) 
#          xterm.append(xy[0])
#          yterm.append(xy[1])
#    
#       #print "returning."
#       return [xlimb,ylimb,xterm,yterm]
#   
#    def drawmoon(self, obs) :
#       stuff = self.moonedges(obs)
#       for i in range(0,len(stuff[0])-1) :
#          self.drawline(stuff[0][i],stuff[1][i],
#                             stuff[0][i+1],stuff[1][i+1],
#                              fillcolor="yellow")
#          # print i,stuff[0][i],stuff[1][i],"->",stuff[0][i+1],stuff[1][i+1]
#       for i in range(0,len(stuff[2])-1) :     
#          self.drawline(stuff[2][i],stuff[3][i],
#                             stuff[2][i+1],stuff[3][i+1],
#                            fillcolor="yellow")
# 

# A couple of classes needed for the auto-finding chart popper.

class UI(Label):

    def __init__(self, master, im):

        if im.mode == "1":
            # bitmap image
            self.image = PyImTk.BitmapImage(im, foreground="white")
            Label.__init__(self, master, image=self.image, bg="black", bd=0)

        else:
            # photo image
            self.image = PyImTk.PhotoImage(im)
            Label.__init__(self, master, image=self.image, bd=0)

# my own class to do the display

class picturewindow(Toplevel) :
   def __init__(self,imagepath,rotangle = 0.) :  # imagepath is the image file.
       # rotangle is clockwise in degrees
       rotrad = rotangle / cooconsts.DEG_IN_RADIAN
       Toplevel.__init__(self)
       # self.geometry("600x600+100-100")
       self.title(imagepath)
       im = PyIm.open(imagepath)

       (width,height) = im.size

       if TINYPIX and height < 800. : 
          factor = 800. / height
          width = width * factor
          height = height * factor 
          im = im.resize((int(width), int(height)))
#          print "image resized to ",width, height 

       if height > 900. :
          factor = 900. / height
          width = width * factor
          height = height * factor 
          im = im.resize((int(width), int(height)))
          print("image resized to ",width, height)

       if rotangle != 0. : 
#       if True : 
      #    print "width,height ",width,height
          rotated_width = int(width * abs(cos(rotrad)) + height * abs(sin(rotrad)))
          rotated_height = int(height * abs(cos(rotrad)) + width * abs(sin(rotrad)))
      #    print "rotated width and height are ",rotated_width, rotated_height
          imrot = im.rotate(rotangle,expand=1).resize((rotated_width, rotated_height))
          UI(self,imrot).pack()
       else : 
          UI(self,im).pack()

       self.bind('q',(lambda event: self.destroy()))

### START start Start of CHART chart Chart section ####

# Take charts out from under DScat section ...

#def feedchart(vars) :
##   for v in vars:
##     print v.get()
#   win_raise(chartwin)
#   tempcoo = celest([0.,0.,2000.])
##   print 'vars[0].get() = ',vars[0].get()
#   tempcoo.ra.val = getradec(vars[3].get(), "h")
#   tempcoo.dec.val = getradec(vars[4].get(), "d")
##   print "going for ",tempcoo.ra.val,tempcoo.dec.val
#   chartwin.getcharts(tempcoo.ra.val,tempcoo.dec.val)
# 
class chartselwin(Toplevel) :

    def __init__(self) :
        Toplevel.__init__(self)
        self.title('Finder charts')
        self.topbox = Frame(self)
        self.topbox.pack()
        self.chartvars = self.makeform(self.topbox)
        self.textfr = Frame(self)
        self.textwin = ScrolledText(self.textfr,width=70,height=18,bckgr="#16520d",
          foregr="white")
        self.textfr.pack(expand = YES, fill = BOTH)
        self.textwin.pack(expand = YES, fill = BOTH)
        Button(self.topbox, text = "Erase", command = 
          (lambda : self.textwin.erasetext()),width=3).pack(side=LEFT,
            anchor=E, fill=Y)
        Button(self.topbox,text="Hide",command = 
          (lambda : self.withdraw()),width=3).pack(side=LEFT,anchor=E, fill=Y)
 
        self.directions = """
 This window is a tool for displaying jpg or png finding charts 
 that are in a directory.  Charts are assumed N-top E-left in 
 original form, and can be flipped and rotated as desired (options
 for this are tuned to the MDM Obsevatory MIS) .
 
 You need a master list in which charts are listed one per line
 with J2000 coordinates, e.g.:
 
 myobj.jpg   22:33:44  55:00:11  
 
 "Read chart list" - reads the master list
 "List charts for this location" - searches for charts within a few
     arcmin of the present ra and dec (read from main window)
 "Display highlighted chart" - You highlight a chart using the left
   mouse button, then this command converts chart into "temp.fits" and 
   displays it using DS9, flipped and rotated if desired.
 """
        Button(self.topbox,text="Help",command = 
                (lambda d = self.directions : self.textwin.appendline(d)),width=3).pack(side=LEFT,
                anchor=E,fill=Y)
        Button(self.topbox,text="Read chart\nlist",command =
          (lambda : self.readchartlist()),width=7).pack(side=LEFT,
                anchor=E,fill=Y)
        Button(self.topbox,text="List charts for\nthis location",command=
          (lambda i = invars: self.getcharts(i)),width=10).pack(side=LEFT,
                anchor=E,fill=Y)
        Button(self.topbox,text="Display high-\nlighted chart",command =
          (lambda : self.displaychart()),width=9).pack(side=LEFT,
                anchor=E,fill=Y)

#       self.directions = """
#This window is a tool for displaying jpg or png finding charts 
#that are in a directory.  Charts are assumed N-top E-left in 
#original form, and can be flipped and rotated as desired (options
#for this are tuned to the MDM Obsevatory MIS) .
#
#You need a master list in which charts are listed one per line
#with J2000 coordinates, e.g.:
#
#myobj.jpg   22:33:44  55:00:11  
#
#"Read chart list" - reads the master list
#"List charts for this location" - searches for charts within a few
#    arcmin of the present ra and dec (read from main window)
#"Display highlighted chart" - makes a "temp.fits" file and 
#    displays it using ds9, flipped and rotated if desired.
#"""
#
    def makeform(self,topbox) :
        chartfields = ['Slit?','Rotator','Path to charts','Chart key file']
        chartvars = []
        entframe = Frame(topbox)
        for f in chartfields :
           row = Frame(topbox)
           lab = Label(row,width=11,text=f)
           var = StringVar()
           if f == 'Slit?' :
              var.set('y')
              lab.pack(side=LEFT)
              row.pack(side=TOP,fill=X,expand=YES)
              but = Radiobutton(row,text="yes",variable=var,value='y')
              but.pack(side=RIGHT,expand=YES,fill=X)
              but = Radiobutton(row,text="no",variable=var,value='n')
              but.pack(side=RIGHT,expand=YES,fill=X)
           else :
              if f == 'Rotator' :
                 var.set('0.')
              if f == 'Path to charts' :
                 var.set('/home/thorsten/iraf/charts/')
              if f == 'Chart key file' :
                 var.set('keyfile')
              lab.pack(side=LEFT)
              row.pack(side=TOP,fill=X,expand=YES)
              ent = Entry(row,width=10)
              ent.pack(side=RIGHT,expand=YES,fill=X)
              ent.config(textvariable = var)

           chartvars.append(var)

        return(chartvars)

    def readchartlist(self) :
 
        # reads the chart listing.

                path = self.chartvars[2].get()
                path = path.strip()
                if path[-1] != '/' :
                        path = path + '/'
                file = self.chartvars[3].get()
                file = file.strip()     
                try :
                        chartkey = open(path+file,"r")
                except :
                        self.textwin.appendline("Can't get the chart key list, '%s'.  No dice." % (path+file))
                        return -1
        
                for l in chartkey.readlines() :
                        try :
                                x = l.split()
                                ra = todeci(x[1])
                                dec = todeci(x[2])
                                name = path + x[0]
                                if ra != None and dec != None : 
                                        charts.append([name,ra,dec] + x[3:])
                                else :
                                        self.textwin.appendline("Line didn't parse: %s" % (l.strip()))
                        except :
                                self.textwin.appendline("Line didn't parse: %s" % (l.strip()))
                        
                self.textwin.appendline("Loaded information on %d charts." % len(charts))
                return(len(charts))
        
    def getcharts(self,vars) :  # feeding with vars from main window.

        tempcoo = celest([0.,0.,2000.])
        tempcoo.ra.val = getradec(vars[1].get(),"h")
        tempcoo.dec.val = getradec(vars[2].get(),"d")
        eqstr = string.upper(vars[3].get())
        if eqstr == "DATE" :
                eqtmp = julian_ep(obs.jd)       
        else : eqtmp = float(vars[3].get())

        tempcoo.selfprecess(2000.)
        ra =  tempcoo.ra.val
        dec = tempcoo.dec.val
        
        # found = []
        rastr = vars[1].get
        nfound = 0
        # print "ra, dec = ",ra,dec
        if len(charts) == 0 :
           self.textwin.appendline("Oops! You must load the chart list first.")
           
        for c in charts :
           dist = _skysub.subtend(ra,dec,c[1],c[2]) * \
                cooconsts.ARCSEC_IN_RADIAN
           if dist < 240. :   # 4 arcmin -- generous ... 
              self.textwin.appendline(c[0])
              nfound = nfound + 1
        if nfound == 0 :
           self.textwin.appendline("Sorry, no charts catalogued at that position.")

    # for some reason ImageMagick dropped the fits converter.  So display with "display" 
    # -- i.e. imagemagick -- which is more direct in any case.

    def displaychart(self) :
        sel = self.textwin.selection_get()
        sel.replace(" ","")  # patch at ends in case of sloppy sel

        slitview = (self.chartvars[0].get() == 'y')
        rotangle = float(self.chartvars[1].get())
#        if ok == 0 :
        for duh in range(0,1) :  # dummy to preserve indentation in csae of reneg later.
           if slitview :
# With Andor and Yorke Brown slitviewer, doing a flip in XY yields N at top, 
# E to left, with 18 degree offset
        #       imageangle = 18. + rotangle # for ImageMagick
        # replaced ImageMagick with PyIm 2015 March -- ImageMagick was bailing on
        # usable images which PyIm managed fine with.  HOWEVER I was "on" the 1.3m 
        # so COULD NOT TEST ROTATION.
               imageangle = -18. - rotangle # for PyIm  
#              os.system("display -rotate %s %s &" % (imageangle, sel))  # ImageMagick
               pw = picturewindow(sel, rotangle = imageangle)            # my PyIm class
           else :
#              os.system("display %s &" % (sel))
               pw = picturewindow(sel)

### END end End of Chart CHART chart section ##########

### START of DSS image server section ### 

class ImageWindow(Toplevel) :

        # Little window for driving ds9.

        def __init__(self) :
                Toplevel.__init__(self)
                self.title('Image server config')
                self.topbox = Frame(self)
                self.vars = self.makeform(self.topbox)
                Button(self.topbox,text = "Hide",command = 
                        (lambda : self.withdraw())).pack(side=LEFT,anchor=E,fill=Y)
                Button(self.topbox,text="Grab\nDSS",command = 
                        (lambda i = invars: self.grabDSS(i))).pack(side=LEFT,anchor=E,fill=Y)
                Button(self.topbox,text="Mark chart",command = 
                        (lambda i = invars: self.chartmarkup(i))).pack(side=LEFT,anchor=E,fill=Y)
                Button(self.topbox,text="Quick\nblink",command = 
                        (lambda i = invars: self.quickblink(i))).pack(side=LEFT,anchor=E,fill=Y)
                Button(self.topbox,text="PAN\nSTARRS",command = 
                        (lambda i = invars: self.grabpanstarrs(i))).pack(side=LEFT,anchor=E,fill=Y)
        def makeform(self,topbox) :
                variables = {}
                imagefields = ['Server','Survey','size (arcmin)','parity','rotation angle']
                for f in imagefields :
                        row = Frame(topbox)
                        row.pack(side=TOP,anchor = NW)
                        lab = Label(row,width=15,text=f)
                        ent = Entry(row, width=10)
                        var = StringVar()
                        if f == "Server" :
                                var.set('stsci')
                                lab.pack(side=LEFT)
                                but=Radiobutton(row,text="STScI",variable=var,value='stsci').pack(side=RIGHT,expand=YES,fill=X) 
                                but=Radiobutton(row,text="CfA",variable=var,value='sao').pack(side=RIGHT,expand=YES,fill=X)     
                        elif f == "Survey" :
                                var.set('dss2red')
                                lab.pack(side=LEFT)
                                but=Radiobutton(row,text="DSS2-red",variable=var,value='dss2red').pack(side=RIGHT,expand=YES,fill=X)    
                                but=Radiobutton(row,text="best available",variable=var,value='ALL').pack(side=RIGHT,expand=YES,fill=X)  
                                but=Radiobutton(row,text="DSS1",variable=var,value='dss').pack(side=RIGHT,expand=YES,fill=X)    
                        elif f == "parity" :
                                var.set('standard')
                                lab.pack(side=LEFT)
                        else  :
                                ent.config(textvariable = var)
                                lab.pack(side=LEFT,anchor = NW)
                                ent.pack(side=LEFT,anchor = NW)
                        if f == 'size (arcmin)' :
                                var.set('8.') 
                        if f == 'rotation angle' :
                                var.set('0.')
                        variables[f] = var 
                topbox.pack()

                return variables

        def ds9_is_open(self) :
                is_open = True
                tmptmp = os.popen("xpaaccess ds9")  # returns a file-like object of output 
                for l in tmptmp.readlines() :
                        if l.find("no") > -1 :
                                is_open = False
                tmptmp.close()
                # print "ds9_is open returning ",is_open
                return is_open
                                
        def grabDSS(self,inv) :  # inv is the invars 
                
                rotation = float(self.vars['rotation angle'].get())
                if self.vars['parity'].get() == 'standard' :
                        parity = 1.
                else :
                        parity = -1.
                server = self.vars['Server'].get()   
                survey = self.vars['Survey'].get() 
                print("server = ",server, "survey = ",survey)
                if server == 'sao' :
                        survey = 'dss1'
                        self.vars['Survey'].set(survey)   # force, CfA only has DSS1.
                size = float(self.vars['size (arcmin)'].get())

                rastr = inv['RA'].get()
                decstr = inv['dec'].get()
                eqstr = inv['equinox'].get()
   
                eq = float(eqstr)
                if eq == 2000. :
                    tempcoo = SkyCoord(rastr,decstr,frame='icrs',unit=(u.hour,u.deg))
                else : 
                    eqstr = 'J' + eqstr
                    fr = FK5(eqstr)
                    tempcoo = SkyCoord(rastr,decstr,frame = fr,unit=(u.hour,u.deg))
        
                rastr = tempcoo.ra.to_string(unit=u.hour,sep=":")
                decstr = tempcoo.dec.to_string(sep = ":")

        # 2008 Jan -- let's rewrite this using xpa access points and a more sophisticated everything.

                ntries = 0 
 
                if not self.ds9_is_open() :   # there is no ds9 up ...
                        # os.spawnlp(os.P_DETACH,"ds9")
                        if parity == 1 :
                                os.system("ds9 -wcs align yes -rotate %8.0f &" % (rotation))
                        else :
                                os.system("ds9 -wcs align yes -rotate %8.0f -orient x &" % (rotation))
                        sleep(1.0)
                        while (not self.ds9_is_open()) and ntries < 10:
                                # print("waiting for ds9 ...")
                                sleep(1.)  # have to wait for process to start before commanding it ... this is a guess
                                ntries = ntries + 1     
                # print("We are still here.")
                if self.ds9_is_open() :
                        os.system('xpaset -p ds9 cmap invert yes')
                        os.system('xpaset -p ds9 scale mode 97')
                        os.system('xpaset -p ds9 height 680')
                        os.system("xpaset -p ds9 dssstsci size %4.1f %4.1f" % (size,size))
                        if server == "stsci" :
                                if survey == "dss2red" :
                                        os.system("xpaset -p ds9 dssstsci survey poss2ukstu_red")
                                        os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
                                elif survey == "ALL" :  # literally all, not best of all.
                                        print("Getting all DSS images from STScI.")
                                        os.system("xpaset -p ds9 dssstsci survey poss2ukstu_red")
                                        os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
                                        os.system("xpaset -p ds9 dssstsci survey poss2ukstu_blue")
                                        os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
                                        os.system("xpaset -p ds9 dssstsci survey poss2ukstu_ir")
                                        os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
                                        os.system("xpaset -p ds9 dssstsci survey poss1_red")
                                        os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
                                        os.system("xpaset -p ds9 dssstsci survey poss1_blue")
                                        os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
                                        os.system("xpaset -p ds9 dssstsci survey quickv")
                                        os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
                                        os.system("xpaset -p ds9 dssstsci survey phase2_gsc2")
                                        os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
                                        os.system("xpaset -p ds9 dssstsci survey phase2_gsc1")
                                        os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
                                else :
                                        os.system("xpaset -p ds9 dssstsci survey all")
                                        os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
                        else :
                                os.system("xpaset -p ds9 dsssao survey all")   
                                os.system("xpaset -p ds9 dsssao coord %s %s" % (rastr,decstr))
                        os.system("xpaset -p ds9 regions show yes")  # seems to turn off for some reason.
                        os.system('xpaset -p ds9 zoom to fit')
                                        
                else :
                        print("DS9 failed to start!")

#    Old grab_dss code worked by constructing an elaborate string and doing it all on the command line.  
#               if parity == 1. :
#                       comstring = "ds9 -dss coord %s %s -dss server %s -dss survey %s -dss size %4.1f %4.1f -zscale -rotate %8.0f -wcs align yes &" % (rastr,decstr,server,survey,size,size,rotation)
#               else :
#                       comstring = "ds9 -dss coord %s %s -dss server %s -dss survey %s -dss size %4.1f %4.1f -zscale -rotate %8.0f -orient x -wcs align yes &" % (rastr,decstr,server,survey,size,size,rotation)
        #       print "woulda said:"
#               print comstring
#               os.system(comstring)
#
#        def grabpanstarrs(self,inv) :  # inv is the invars 
#            tempcoo = celest([0.,0.,2000.])
#            tempcoo.ra.val = getradec(inv[1].get(),"h")
#            tempcoo.dec.val = getradec(inv[2].get(),"d")
#            size = float(self.vars[2].get())
#            ps.query_cutout(tempcoo,arcminsize = size)
# 
#        def chartmarkup(self, inv) :  # Marks up the (already displayed) image on ds9.   
#
#                tempcoo = celest([0.,0.,2000.])
#                tempcoo.ra.val = getradec(inv[1].get(),"h")
#                tempcoo.dec.val = getradec(inv[2].get(),"d")
#                eqstr = string.upper(inv[3].get())
#                if eqstr == "DATE" :
#                        tempcoo.equinox = julian_ep(obs.jd)     
#                else : tempcoo.equinox = float(invars[3].get())
#                tempcoo.selfprecess(2000.)
#        
#                rastr = tempcoo.ra.tripletstring(delin = ":",places = 2, showsign = 0)
#                decstr = tempcoo.dec.tripletstring(delin = ":",places = 1, showsign = 1)
#
#                cosdec = cos(tempcoo.dec.val / cooconsts.DEG_IN_RADIAN)
#                size = float(self.vars[2].get())
#                rahalfwidth = (size * 2.5 / (3600. * cosdec)) # hrs
#                dechalfwidth = size / 120.  # deg
#                namelabel = celest([tempcoo.ra.val + rahalfwidth * 0.6,tempcoo.dec.val + dechalfwidth * 0.9,2000.])
#                ralabel = celest([tempcoo.ra.val - rahalfwidth * 0.4,tempcoo.dec.val + dechalfwidth * 0.9, 2000.])
#                declabel = celest([tempcoo.ra.val - rahalfwidth * 0.4,tempcoo.dec.val + dechalfwidth * 0.85, 2000.])
#                compass = celest([tempcoo.ra.val - rahalfwidth * 0.7, tempcoo.dec.val + dechalfwidth * 0.8, 2000.])
#                comment1 = celest([tempcoo.ra.val + rahalfwidth * 0.6,tempcoo.dec.val - dechalfwidth * 0.92,2000.])
#                linestart = celest([tempcoo.ra.val - rahalfwidth * 0.7, tempcoo.dec.val - dechalfwidth * 0.5, 2000.])
#                linestop = celest([tempcoo.ra.val - rahalfwidth * 0.7, linestart.dec.val + (1./30.), 2000.])
#                linelabel = celest([tempcoo.ra.val - rahalfwidth * 0.67, linestart.dec.val + (1./60.), 2000.])
#                os.system('xpaset -p ds9 cmap invert yes')
#                os.system('xpaset -p ds9 scale mode 99')
#                # with vertical mode, can fit a sizeable window on laptop screen!
#                os.system('xpaset -p ds9 width 600')
#                os.system('xpaset -p ds9 height 600')
#                os.system('xpaset -p ds9 zoom to fit')
#                os.system('xpaset -p ds9 regions format ds9')
#                os.system('xpaset -p ds9 regions system wcs')
#                os.system('xpaset -p ds9 regions sky icrs')
#                os.system('xpaset -p ds9 regions skyformat sexagesimal')
#        #       os.system('xpaset -p ds9 regions wcs yes')
#                os.system('xpaset -p ds9 zoom to fit')
##               os.system('echo "global font=\"helvetica 14 normal\"" | xpaset ds9 regions')
#                #os.system('echo "icrs; circle %s %s 0.15\' # color=red" | xpaset ds9 regions' % 
#                #               (tempcoo.ra.tripletstring(delin = ":",showsign="no"), 
#                #               (tempcoo.dec.tripletstring(delin = ":",showsign="yes"))))
##               print "1"
#                os.system('echo "icrs; text %s %s # textangle=0 text={%s} font=\\\"helvetica 12 normal\\\"" | xpaset ds9 regions' % 
#                        (namelabel.ra.tripletstring(delin=":",showsign="no"),
#                                        namelabel.dec.tripletstring(delin=":",showsign="yes"),
#                                        inv[0].get()))
##               print "2"
#                os.system('echo "icrs; text %s %s # textangle=0 text={%s} font=\\\"helvetica 12 normal\\\"" | xpaset ds9 regions' 
#                                        %  (declabel.ra.tripletstring(delin=":",showsign="no"),
#                                        declabel.dec.tripletstring(delin=":",showsign="yes"),
#                                        decstr))
##               print "3"
#                os.system('echo "icrs; text %s %s # textangle=0 text={%s} font=\\\"helvetica 12 normal\\\"" | xpaset ds9 regions' 
#                                        % (ralabel.ra.tripletstring(delin=":",showsign="no"),
#                                        ralabel.dec.tripletstring(delin=":",showsign="yes"),
#                                        rastr))
##               print "4"
#                os.system('echo "icrs; compass %s %s 0.6\' # compass=icrs \\\"N\\\" \\\"E\\\" 1 1" | xpaset ds9 regions' 
#                                        %  (compass.ra.tripletstring(delin=":",showsign="no"),
#                                        compass.dec.tripletstring(delin=":",showsign="yes")))
#                os.system('echo "icrs; line %s %s %s %s" | xpaset ds9 regions' % 
#                                       (linestart.ra.tripletstring(delin=":",showsign="no"),
#                                        linestart.dec.tripletstring(delin=":",showsign="yes"),
#                                        linestop.ra.tripletstring(delin=":",showsign="no"),
#                                        linestop.dec.tripletstring(delin=":",showsign="yes")))
#
#                os.system('echo "icrs; text %s %s # textangle=90 text={2 arcmin}" | xpaset ds9 regions' 
#                                        % (linelabel.ra.tripletstring(delin=":",showsign="no"),
#                                        linelabel.dec.tripletstring(delin=":",showsign="yes")))
##               print "5"
#                os.system('xpaset -p ds9 regions select all')
#                os.system('xpaset -p ds9 regions color black')
#                os.system('xpaset -p ds9 regions width 2')
#                # Finally figured out how to get the font by adding an explicit specification of the 
#                # font to each command -- it involves nested escapes giving \\\"
#                os.system('xpaset -p ds9 regions select none')
#                # draw position marker here to keep it thin and red.
#                os.system('echo "icrs; circle %s %s 0.15\' # color=red" | xpaset ds9 regions' % 
#                                (tempcoo.ra.tripletstring(delin = ":",showsign="no"), 
#                                (tempcoo.dec.tripletstring(delin = ":",showsign="yes"))))
#
#        def quickblink(self, inv) :  # grabs all STScI DSS images and blinks 'em.
#
#            tempcoo = celest([0.,0.,2000.])
#            tempcoo.ra.val = getradec(inv[1].get(),"h")
#            tempcoo.dec.val = getradec(inv[2].get(),"d")
#            eqstr = string.upper(inv[3].get())
#            if eqstr == "DATE" :
#                tempcoo.equinox = julian_ep(obs.jd)     
#            else : tempcoo.equinox = float(invars[3].get())
#            tempcoo.selfprecess(2000.)
#            
#            rastr = tempcoo.ra.tripletstring(delin = ":",places = 3, showsign = 0)
#            decstr = tempcoo.dec.tripletstring(delin = ":",places = 2, showsign = 1)
#
#            if not self.ds9_is_open() :
#                os.system("ds9 -wcs align yes &")
#            sleep(1.0)
#            ntries = 0
#            while(not self.ds9_is_open()) and ntries < 10 :
#               sleep(0.5)
#               ntries = ntries + 1
#
#            if ntries == 10 :
#               print("ds9 didn't open.")
#               return
#
#            os.system('xpaset -p ds9 cmap invert yes')
#            os.system('xpaset -p ds9 scale mode 97')
#            os.system('xpaset -p ds9 width 600')
#            os.system('xpaset -p ds9 height 600')
#            os.system('xpaset -p ds9 dssstsci size 5. 5.')
#            os.system("xpaset -p ds9 dssstsci survey poss2ukstu_red")
#            os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
#            os.system("xpaset -p ds9 dssstsci survey poss2ukstu_blue")
#            os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
#            os.system("xpaset -p ds9 dssstsci survey poss2ukstu_ir")
#            os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
#            os.system("xpaset -p ds9 dssstsci survey poss1_red")
#            os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
#            os.system("xpaset -p ds9 dssstsci survey poss1_blue")
#            os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
#            os.system("xpaset -p ds9 dssstsci survey quickv")
#            os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
#            os.system("xpaset -p ds9 dssstsci survey phase2_gsc2")
#            os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
#            os.system("xpaset -p ds9 dssstsci survey phase2_gsc1")
#            os.system("xpaset -p ds9 dssstsci coord %s %s" % (rastr,decstr))
#
#            os.system("xpaset -p ds9 regions show yes")
#            os.system("xpaset -p ds9 tile no")
#            os.system("xpaset -p ds9 regions format ds9")
#            os.system("xpaset -p ds9 regions system wcs")
#            os.system("xpaset -p ds9 regions sky icrs")
#            os.system("xpaset -p ds9 regions skyformat sexagesimal")
#            os.system("xpaset -p ds9 regions wcs yes")
#            os.system('echo "icrs; circle %s %s 0.15\'" | xpaset ds9 regions' %
#                         (tempcoo.ra.tripletstring(delin = ":",showsign="no"),
#                         (tempcoo.dec.tripletstring(delin = ":",showsign="yes"))))
#
#            os.system('xpaset -p ds9 regions select all')
#            os.system('xpaset -p ds9 regions color black')
#            os.system('xpaset -p ds9 regions width 2')
#            os.system('xpaset -p ds9 regions system wcs')
#            os.system('xpaset -p ds9 zoom to fit')
#            os.system("xpaset -p ds9 match frames wcs")
#            os.system('xpaset -p ds9 regions select none')
#            sleep(2.0)
#            os.system("xpaset -p ds9 blink yes")
#
### END of DSS image server section   ###

        
# Here comes the main MAIN Main program. 

# Some utility routines adapted from Mark Lutz' *Programming Python*,
# published by O'Reilly.

# Adapted from Mark Lutz' "Programming in Python", O'Reilly

class ScrolledList(Frame) : 
        def __init__(self, options, parent=None, kill=1) :
                Frame.__init__(self,parent)
                self.pack(expand=YES, fill=BOTH)
                self.makeWidgets(options)
                self.val = 0
                self.onepass = kill
        def handleList(self,event) :
                index = self.listbox.curselection()
                self.val = self.listbox.get(index)
                # if self.onepass == 1: 
                #       self.quit()  # kill box after any selection.
        def makeWidgets(self,options) :
                sbar = Scrollbar(self)
                list = Listbox(self,relief=SUNKEN)
                sbar.config(command=list.yview)
                list.config(yscrollcommand=sbar.set)
                sbar.pack(side=RIGHT, fill=Y)
                list.pack(side=LEFT,expand=YES, fill=BOTH)
                pos = 0
#               print "options %s" % options
                for label in options:
#                       print "label %s" % label
                        list.insert(pos,label)
                        pos = pos + 1
                list.bind('<Double-1>', self.handleList)
                self.listbox=list
                self.sb = sbar

# Largely from Lutz' Programming Python book, with some
# enhancements I added.



class ScrolledText(Frame) :
        def __init__(self, parent = None, text = '', file = None, width=80, height=24, bckgr = "#330066",foregr="yellow") :
                Frame.__init__(self, parent)   # make me expendable
                self.pack(expand = YES, fill=BOTH)
                self.bckgr = bckgr
                self.foregr = foregr

                self.makewidgets(width=width, height=height,bck=self.bckgr,
                                fgr=self.foregr)
                self.settext(text,file)
                 
        def makewidgets(self,width=80,height=24,bck="#330066",fgr="yellow") :
                sbar = Scrollbar(self)
                text = Text(self, relief = SUNKEN,width=width,height=height,
                        bg=bck,fg=fgr,font=("Courier",14,"bold"))
                sbar.config(command=text.yview)      # xlink sbar and text
                text.config(yscrollcommand=sbar.set) # each moves the other
                sbar.pack(side=LEFT, fill=Y)
                text.pack(side=LEFT, expand=YES, fill=BOTH)
                self.text = text
        def settext(self, text='', file = None) :
                if file:
                        text = open(file, 'r').read()
                self.text.delete('1.0',END)  # delete current text
                self.text.insert('1.0',text) # add at line 1, col 0
                self.text.mark_set(INSERT, '1.0')
                self.text.focus()            # save user a click
        def gettext(self) :
                return self.text.get('1.0', END+'-1c')  
        def getsel(self) :
                return self.text.selection_get()
        def appendline(self,line='\n',maxlines = 10000) :
                self.text.insert(END,line + "\n")
                if float(self.text.index(END)) > maxlines :
                        #print maxlines, self.text.index(END), \
                        #        self.text.index(END) > float(maxlines)
                        self.text.delete('1.0','2.0')
                self.text.see(END)
        def erasetext(self) :
                self.text.delete('1.0',END) 
        def dumptext(self,filename) :
                # outstuff = self.text.selection_get() 
                outf = open(filename,"a")
                outf.write(self.gettext())
                outf.close()

# quitter widget is a direct copy of Lutz' routine, except without the 
# __name__ == '__main__' trick.

#############################################
# a quit button that verifies exit requests;
# to reuse, attach an instance to other guis
#############################################

from tkinter import *                          # get widget classes
from tkinter import messagebox
#from tkMessageBox import asktokcancel           # get canned std dialog
#from tkMessageBox import asktokcancel           # get canned std dialog

class Quitter(Frame):                          # subclass our GUI
    def __init__(self, parent=None):           # constructor method
        Frame.__init__(self, parent)
        self.pack()
      #  widget = Button(self, text='Quit', command=self.quit, bg="#ffaa75")
      #  widget = Button(self, text='Quit', command=self.quit, width=butwid(3), bg="#ffaa75")
        if is_MacOS : 
            widget = Button(self, text='Quit', command=self.quit, width=butwid(3), highlightbackground="#ffaa75")
        else :
            widget = Button(self, text='Quit', command=self.quit, width=butwid(3), bg="#ffaa75")
        widget.pack(expand=YES, fill=BOTH, side=LEFT)
    def quit(self):
        ans = messagebox.askokcancel('Verify exit', "Really quit?")
        if ans: Frame.quit(self)

# Routines to create windows and so on ... 

#def makehelpframe(parentwin, intext = helptext) :
#
#        textfr = Frame(parentwin)
#        textfr.pack()
#        scrtxt = ScrolledText(textfr,text = intext)
#        
#        Button(textfr, text = 'Hide',command = (lambda: parentwin.withdraw())).pack(side = TOP)
        

#obs.setsite('mdm')   # THIS SETS THE DEFAULT SITE - to MDM on Kitt Peak.
                   # to change, refer to cooclasses.py for codes or to 
                   # customize for another site not already coded.

def moonwarningcolor(obs, objalt = 90.) :
        if obs.moonaltit < (0. * u.deg) or objalt < 0. : return normalentry  
        if obs.sunaltit > (-12. * u.deg) : return normalentry  # twilight dominates
        if obs.moonobjang > (25. * u.deg) :
                if obs.lunsky > 21.5 : return normalentry
                elif obs.lunsky > 19.5 : 
                        return "lightblue"
                else : return "#DD88FF"  # a light purple
        if obs.moonobjang < (10. * u.deg) : return "red"  # always flag this ...
        if obs.lunsky > 21.5 : return normalentry
        if obs.lunsky > 19.5 : return "yellow"
        if obs.lunsky > 18. : return "orange"
        else : return "red"

def twilightwarningcolor(obs) :
        if obs.sunaltit > (-0.8 * u.deg) : return "lightblue"
        if obs.sunaltit < (-18. * u.deg)  : return normalentry
        if obs.twi < 4. : return "yellow"
        if obs.twi < 8. : return "orange"
        return "red"

def airmasswarningcolor(obs,down_is_red = True) :
        
        if obs.altit < (0. * u.deg) : 
                if not down_is_red :
                        return "lightgreen"
                else :
                        return "red"
        if obs.airmass <= 2. : return normalentry
        if obs.airmass > 4. : return "red"
        if obs.airmass > 3. : return "orange"
        else : return "yellow"

def hawarningcolor(obs, limit) :
        
        if (obs.hanow > -1. * limit) and (obs.hanow < limit) : return normalentry
        else : return "orange"

def planetwarning(planets, ra, dec, tolerance) :

    return
#       
#       warning = ""
#       for i in range(1,10) :
#               if i != 3 :
#                       ang = subtendang(planets[i][1],planets[i][2],ra,dec) \
#                               * 57.2957795130823
#                       if ang < tolerance :
#                               warning = warning + " " + planets[i][0]
#       if warning == "" :
#               warning = " --- "
#               return [warning,normalentry]
#       else :
#               return [warning,"orange"]


def unpack_day_time(var1,var2,var3) :

# Takes the three text entry fields that give the date, time, and 
# UT/local flag, and updates the "class instant".

#       time1 = var1.get()    
        x = var1.split()
        #print "x:",x
        year = int(x[0])
        month = int(x[1])  # leave it numeric for now
        day = int(x[2]) 
#       time2 = var2.get()
        x = var2.split() 
        #print "x:",x
        hour = int(x[0])
        minute = int(x[1])
        try :   # allow for there to be no seconds field.
            floatsecond = float(x[2])
            second = int(floatsecond)
            microsecond = int((floatsecond - second) * 1000000.)
        except :
            second = 0
            microsecond = 0
        localfl = var3
        #print "year month day hour min sec",year,month,day,hour,minute,second
        #print "localfl ",localfl
        if localfl == 'y' or localfl == 'Y' :
                #obsin.settime(datetime(year,month,day,hour,minute,second,microsecond), 
                obs.settime((year,month,day,hour,minute,second,microsecond), 
                     use_local_time = True)
        else :
                #obsin.settime(datetime(year,month,day,hour,minute,second,microsecond)) 
                obs.settime((year,month,day,hour,minute,second,microsecond),
                     use_local_time = False) 

def equinox_to_float(framein) :
    # The Skycoord.frame.equinox is a weird time object with a represenation like
    # 'J2018.61'.  Unpacking the equinox as a floating-point number takes a bit of
    # doing.
    if framein.name == 'icrs' or framein.name == 'gcrs' : 
        return 2000.  
    # weirdly, 'fk5' needs to be lowercase for this test to work.
    elif framein.name == 'fk5' or framein.name == 'precessedgeocentric' :
        eqstr = framein.equinox.value    # which is a string
        return float(eqstr.replace('J',''))
    else : 
        print("unhandled frame type, ",framein.name)
        return None
        
# Circumstances computes and refreshes almost everything.

def circumstances(variables,outvars,havars,almvars,coovars,planetvars,
                boxes,haboxes,planetboxes, skdisp, airmdisp) :

        # Need to do things in the order site, time, then celestial position,
        # sine the equinox can be "NOW".

        sitecode = variables['sitecode'].get()   # site
        if sitecode == 'x' :   # specify whole site
                longitin = variables['longit'].get()
                latin = variables['lat'].get()
                stdz = variables['stdz'].get()
                stdtabbr = variables['std_abbrev'].get()
                dstabbr = variables['dst_abbrev'].get()
                elevsea = variables['elevsea'].get()
                eleveast = variables['eleveast'].get()
                elevwest = variables['elevwest'].get()
                obsname = variables['obs_name'].get()
                mysite = obs_site(name = obsname, 
                   loc = EarthLocation.from_geodetic(longitin,latin,elevsea),
                   tzstr = stdz, tzstdabbr = stdtabbr, tzdayabbr = dstabbr,
                   westheight = float(elevwest), eastheight = float(eleveast))
                obs.setsite(mysite)
        else : 
                obs.setsite(sitecode)
                variables['longit'].set(obs.site.location.lon.to_string(unit = u.deg, decimal = True, precision = 3))
                variables['lat'].set(obs.site.location.lat.to_string(unit = u.deg, decimal = True, precision = 3,pad = True, alwayssign = True))
                variables['stdz'].set(obs.site.tzstr)
                variables['elevsea'].set("%4.0f" % obs.site.location.height.value)
                variables['eleveast'].set("%4.0f" % obs.site.height_above_east.value)
                variables['elevwest'].set("%4.0f" % obs.site.height_above_west.value)
                variables['obs_name'].set(obs.site.name)
                variables['std_abbrev'].set(obs.site.localtzabbrev)
                variables['dst_abbrev'].set(obs.site.localtzdayabbrev)

        localflag = variables['Time is:'].get()

#        unpack_day_time(obs,variables['date'].get(),variables['time'].get(),variables['Time is:'].get()) 
        unpack_day_time(variables['date'].get(),variables['time'].get(),localflag) 

        outvars['jd'].set("%15.5f" % obs.t.jd)
        boxes['jd'].config(bg = incolor1)

        if localflag == 'y' :  
            tdout = obs.t.to_datetime(timezone = obs.site.localtz)
        else :
            tdout = obs.t.to_datetime()
        variables['date'].set(tdout.strftime("%Y %m %d  %a"))
        variables['time'].set(tdout.strftime("%H %M %S"))
        
        # print "day-time set up."

        rain = variables['RA'].get()
        decin = variables['dec'].get()
        eqin = variables['equinox'].get()

        # if setcelest is handed a tuple, it uses astropy to interpret the input
        # strings pretty flexibly.

        obs.setcelest((rain,decin,eqin))
        variables['RA'].set(obs.celest.ra.to_string(unit = u.hourangle,sep=' ',precision=2,
               pad = True)) 
        variables['dec'].set(obs.celest.dec.to_string(sep=' ',precision=1, pad = True,
               alwayssign = True)) 

        variables['equinox'].set("%7.2f" % equinox_to_float(obs.celest.frame))
        
        # print "*** in circumstances ... sky and moon coming ... "

        obs.computesky()

        #print "sky computed"
        obs.computesunmoon()
        #print "sunmoon too"

        if planetwin.state() == 'normal' :
           obs.computeplanets()
        else :
           # print "planetwin.state() = ",planetwin.state(), "no computation."
           pass
        # print "planets done"

        # obs.computebary()
        obs.computequickbary()
        # print "bary done"

        obs.setnightevents()

        obs.compute_hours_up()

#        print "*** in circumstances ... sky and moon finished. "

        # check for planet warnings ... 

        warning_planet = fillplform(planetvars,planetboxes)
      #  print "warning_planet = ",warning_planet
           
        outstr = ''
        if len(warning_planet) != 0 :
            for p in warning_planet : 
               outstr = outstr + ' ' + p.capitalize()
            boxes['planet_proxim'].configure(bg = 'red')
            boxes['planet_proxim'].configure(fg = 'white')
            outvars['planet_proxim'].set(outstr)
        else :
            outvars['planet_proxim'].set(' ')
            boxes['planet_proxim'].configure(bg = normalentry)
            boxes['planet_proxim'].configure(fg = 'black')

        # Need to explicitly reset the time variable to catch day-of-week changes
#       if variables['Time is:'].get() == 'y' :
#               calstr = obs.calstring(style = 1, stdz = obs.stdz, use_dst = obs.
#                       use_dst, print_day = 1, daycomma = 0,secdigits = 1) 
#       else :
#               calstr = obs.calstring(style = 1, print_day = 1, daycomma = 0,
#                       secdigits = 1)

#       [ymd,hms] = ymd_hms_string(calstr)
#       variables[4].set(ymd)
#       variables[5].set(hms)

        outvars['sidereal'].set(obs.lst.to_string(unit = u.hourangle, sep = ' ', 
                                precision = 0))
        hawarn = hawarningcolor(obs,Angle(6. * u.hour))
        boxes['ha'].configure(bg = hawarn)
        outvars['ha'].set(obs.hanow.to_string(unit = u.hourangle, sep = ' ', 
            precision = 0, alwayssign = True))

        airmasswarn = airmasswarningcolor(obs, down_is_red=True)

        boxes['airmass'].configure(bg = airmasswarn)
        if obs.altit > 0. :
                if obs.airmass < 10. :
                        outvars['airmass'].set("%6.3f" % obs.airmass)
                else :
                        outvars['airmass'].set("More than 10")
        else :
                outvars['airmass'].set("(Down.)")

        if obs.altit < 0. or obs.airmass > 4. :
                boxes['airmass'].configure(bg = "red")
        elif obs.airmass > 3. :
                boxes['airmass'].configure(bg = "orange")
        elif obs.airmass > 2. :
                boxes['airmass'].configure(bg = "yellow")
        else :  boxes['airmass'].configure(bg = normalentry)
#       if obs.altit > 0. :
#               outvars[3].set("%6.3f" % obs.secz)
#       else :
#               outvars[3].set("---")
        outvars['alt_az'].set("%5.2f    az %5.2f" % (obs.altit.deg,obs.az.deg))
        parang_opposite = obs.parang + Angle(180. * u.deg)
        parang_opposite.wrap_at(180. * u.deg)
        outvars['parallactic'].set("%5.1f [%5.1f]  deg" % (obs.parang.deg,parang_opposite.deg))

        outvars['sunradec'].set(obs.sunpos.ra.to_string(unit=u.hourangle, sep = ' ',precision=1,
                   pad = True) + '  ' + 
                obs.sunpos.dec.to_string(sep = ' ',precision=0,pad = True, alwayssign = True))
        outvars['sunaltaz'].set("%4.1f      az %5.1f" % (obs.sunaltit.deg,obs.sunaz.deg))

        twilightwarn = twilightwarningcolor(obs)
        boxes['ztwilight'].configure(bg = twilightwarn)
        if obs.sunaltit < (0. * u.deg) and obs.sunaltit >= (-18. * u.deg) :
                outvars['ztwilight'].set("%5.1f  mag (blue)" % obs.twi)
        elif obs.sunaltit >= (0. * u.deg) :
                outvars['ztwilight'].set("Daytime.")
        else : 
                outvars['ztwilight'].set("No twilight.")

        moonwarn = moonwarningcolor(obs) 
        outvars['moonphase'].set(obs.moonphasedescr)
        outvars['moonradec'].set(obs.moonpos.ra.to_string(unit=u.hourangle, sep = ' ',precision=0,
                   pad = True) + '  ' +
                   obs.moonpos.dec.to_string(sep = ' ',precision=0,pad = True, alwayssign = True,fields=2))
        outvars['moonaltaz'].set("%4.1f      az %5.1f" % (obs.moonaltit.deg,obs.moonaz.deg))
        if obs.moonaltit > (-1 * u.deg)  :
                outvars['illumfrac'].set("%5.3f" % obs.moonillumfrac)
        else :
                outvars['illumfrac'].set(" Moon is down. ")
        boxes['lunsky'].configure(bg = moonwarn)
        if obs.moonaltit > (-1 *u.deg) and obs.sunaltit < (-12. * u.deg) and obs.altit > (0. * u.deg) :
                outvars['lunsky'].set("%5.1f mag." % obs.lunsky)
        elif obs.moonaltit <= (-1 * u.deg) or obs.altit < (0. * u.deg) :        
                outvars['lunsky'].set("---")
        elif obs.sunaltit >= (-12 * u.deg) and obs.sunaltit < (0.* u.deg):
                outvars['lunsky'].set("(Bright twilight.)")
        else :
                outvars['lunsky'].set("(Daylight.)")
        boxes['moon-obj ang.'].configure(bg = moonwarn)
        outvars['moon-obj ang.'].set("%5.1f  deg." % obs.moonobjang.deg)
        outvars['baryjd'].set("%14.5f   [%+6.1f s]" % (obs.tbary.jd,obs.barytcorr.sec))
        outvars['baryvcor'].set("%6.2f  km/s" % obs.baryvcorr.value)
        outvars['constel'].set(" %s" % obs.constel)
#       fillhaform(havars,haboxes)
        fillalmform(almvars)
        fillcooform(coovars)
        
        # print "forms filled"
        skdisp.redraw(obs,objs2000)

        # print "skdisp redraw done"

        objn = variables['objname'].get()

        if objn.find('initializes at') > -1 or objn.find('null') > -1 :
            airmdisp.redraw(obs)
        else :
            airmdisp.redraw(obs, variables['objname'].get())

        # print "exiting."



def steptime(variables, outvars, havars, almvars, coovars, planetvars, 
        boxes, haboxes, planetboxes, skdisp, airmdisp, forward) :
        
#       print "stepping the time ... "
        obs.advancetime(variables['timestep'].get(),forward)
        if variables['Time is:'].get() == 'y' :
            localt = obs.t.to_datetime(timezone = obs.site.localtz)
            variables['date'].set(localt.strftime('%Y %m %d   %a'))  
            variables['time'].set(localt.strftime('%H %M %S'))
        else :
            ut = obs.t.to_datetime()
            variables['date'].set(ut.strftime('%Y %m %d   %a'))  
            variables['time'].set(ut.strftime('%H %M %S'))
        
        circumstances(variables,outvars,havars,almvars,coovars, planetvars,
                boxes, haboxes, planetboxes, skdisp, airmdisp)

def convert_time(variables, outvars, havars, almvars, coovars, planetvars,
        boxes, haboxes, planetboxes, skdisp, airmdisp) :
        
        if variables['Time is:'].get() == 'y' : # convert local to UT
                dt = obs.t.to_datetime()
                variables['Time is:'].set('n')
        else :
                dt = obs.t.to_datetime(timezone = obs.site.localtz)
                variables['Time is:'].set('y')
        variables['date'].set(dt.strftime("%Y %m %d   %a"))

        circumstances(variables,outvars,havars,almvars, coovars, planetvars,
                boxes, haboxes, planetboxes, skdisp, airmdisp)
        

def set_to_now(variables, outvars, havars, almvars,coovars, planetvars, 
        boxes, haboxes, planetboxes, skdisp, airmdisp) :
        # last variable; if < 5., just do it; otherwise, repeat periodically

        # print "variables[6].get() = ",variables[6].get()
#       if(variables[18].get() >= 5.) : repeat = True
#       else : repeat = False
#       proceed = True
        
#       while(proceed) :

#       obs.setlocal('NOW')

        obs.t = Time(ttime.time(),format='unix')
    
        if variables["Time is:"].get() == 'y' :
            dt = obs.t.to_datetime(timezone = obs.site.localtz)
        else :
            dt = obs.t.to_datetime()
        
        variables['date'].set(dt.strftime('%Y %m %d   %a'))
        variables['time'].set(dt.strftime('%H %M %S'))
        
        circumstances(variables,outvars,havars,almvars,coovars, planetvars,
                boxes, haboxes, planetboxes, skdisp, airmdisp)


# This has proven to be too prone to accidental mouse hits.

#def set_to_moused_time(event, variables, outvars, havars, almvars, coovars, planetvars,
#    boxes, haboxes, planetboxes, skdisp, airmdis) :
#    # for catching a mouse event in the airmass window and setting time to that.
#    
#    jdmarked = airmdisp.xpixtojd(event.x) 
#
#    obs.setut(instantinput = jdmarked, stdz = obs.stdz, use_dst = obs.use_dst)
#    if variables[6].get() == 'y' :
#        calstr = obs.calstring(style = 1, stdz = obs.stdz, \
#             use_dst = obs.use_dst, print_day = 1, daycomma = 0) 
#    else :
#        calstr = obs.calstring(style = 1, print_day = 1, daycomma = 0) 
#   
#    [ymd,hms] = ymd_hms_string(calstr)
#    variables[4].set(ymd)
#    variables[5].set(hms)
#   
#    circumstances(variables,outvars,havars,almvars,coovars, planetvars,
#                boxes, haboxes, planetboxes, skdisp, airmdisp)
#

def set_to_moused_obj(event, variables, outvars, havars, almvars,coovars, planetvars, 
        boxes, haboxes, planetboxes, skdisp, airmdisp) :

        # This is immune to the location bug in MacOS because it is a mouse button, not
        # a keystroke.

        markcel = skdisp.pixtocelest(obs,event.x,event.y)

        # this generates 2000 now.
        # print "markcel = ",markcel.quickpr()

        minsep = 100000000000000. * u.deg
        minkey = None 

        for k in objs2000.keys() :
        #       print "key = ",k,"obj",objs2000[k]
                sep = objs2000[k].separation(markcel)
                if sep < minsep :
                        minsep = sep
                        minkey = k

        if minkey == None : return
        
        # print "minkey = ",minkey,"objs2000 = ",objs2000[minkey]

        mincel = objs2000[minkey]

        variables['objname'].set(minkey)
        variables['RA'].set(mincel.ra.to_string(unit=u.hourangle,pad=True,sep=' '))
        variables['dec'].set(mincel.dec.to_string(pad=True,sep=' ',alwayssign = True))
        variables['equinox'].set("2000.")
        
        circumstances(variables,outvars,havars,almvars,coovars, planetvars,
                        boxes, haboxes, planetboxes, skdisp, airmdisp)

def set_to_moused_coord(event, variables, outvars, havars, almvars,coovars, planetvars, 
        boxes, haboxes, planetboxes, skdisp, airmdisp) :

        markcel = skdisp.pixtocelest(obs,event.x,event.y)

        variables['objname'].set("marked:")
        variables['RA'].set(markcel.ra.to_string(unit=u.hourangle,sep=' ',pad=True))
        variables['dec'].set(markcel.dec.to_string(sep = ' ',pad=True,alwayssign = True))
        variables['equinox'].set("2000.") # always hands back 2000 now.
        
        circumstances(variables,outvars,havars,almvars,coovars, planetvars,
                        boxes, haboxes, planetboxes, skdisp, airmdisp)

def set_to_moused_bright(event, variables, outvars, havars, almvars,coovars, planetvars, 
        boxes, haboxes, planetboxes, skdisp, airmdisp) :

        # https://stackoverflow.com/questions/21517114/python-tkinter-key-event-locations-lost-in-macos
        # I'm implementing the workaround in that page

        # print("responding to 'h'")
        # markcel = skdisp.pixtocelest(obs,event.x,event.y)

        if is_MacOS :
            (cx, cy) = keylocation(event) 
            markxyz = skdisp.pixtocelest(obs,cx,cy,return_xyz = True)
        else :
            markxyz = skdisp.pixtocelest(obs,event.x,event.y,return_xyz = True)
        # print("markxyz = ",markxyz)

        i = skdisp.nearestbright(markxyz)

        variables['objname'].set(skdisp.brightnames[i])

        ratmp = np.arctan2(skdisp.bright2000[1,i],skdisp.bright2000[0,i])
        dectmp = np.arcsin(skdisp.bright2000[2,i])
        
        starc = SkyCoord(ratmp,dectmp,unit=(u.radian,u.radian),frame='icrs')
        variables['RA'].set(starc.ra.to_string(unit=u.hourangle,sep = ' ',precision=2,pad=True))
        variables['dec'].set(starc.dec.to_string(sep=' ',precision=1,pad=True,alwayssign = True))
        variables['equinox'].set("2000.")

        circumstances(variables,outvars,havars,almvars,coovars, planetvars,
                        boxes, haboxes, planetboxes, skdisp, airmdisp)

#         # print "nearbr ",nearbr
#         if len(nearbr) > 3 :
#       
#               mincel = celest([nearbr[1],nearbr[2],2000.])
#       
#               variables[0].set(nearbr[0])
#               variables[1].set(mincel.ra.tripletstring(places=2,showsign=0))
#               variables[2].set(mincel.dec.tripletstring(places=1,showsign=0))
#               variables[3].set("%6.2f" % mincel.equinox)
#               
#               circumstances(variables,outvars,havars,almvars,coovars, planetvars,
#                               boxes, haboxes, planetboxes, skdisp, airmdisp)
#       
   
# Globally define two colors to flag inputs - one for site stuff.

incolor2 = "#efdee4" # very pale pinkish - accepts input only conditionally
incolor1 = "white"   # white fields accept input
normalentry = "#ebebeb"  # slightly grey - not used or input

# This define the main window and the input and output fields.

# These are the names of the fields one the left side of the gui;
# they are also used as dictionary keys.

# These are all in text-entry widgets EXCEPT: 
#  - 'Time is:' is controlled by a radio button
#  - site_abbrev is controlled by a separate radio-button panel

infields = ('objname','RA','dec','equinox','date','time','Time is:',\
    'timestep','obs_name','longit','lat', \
        'stdz', 'std_abbrev', 'dst_abbrev',  \
          'elevsea','eleveast', 'elevwest', 'sitecode') # , 'autoupdate')

def makeinform(root, fields) :

#       variables = []

        variables = {}    # change 2018 August - make 'variables' a directory since
                          # it'll have to be reorganized.

#       calstr = obs.calstring(style = 1, stdz = obs.stdz, \
#               use_dst = obs.use_dst, print_day = 1, daycomma = 0) 
#
#       [ymd, hms] = ymd_hms_string(calstr)

        localt = obs.t.to_datetime(timezone = obs.site.localtz)

        row = Frame(root)
        lab = Label(row, text = "Input Variables")
        row.pack(side = TOP)
        lab.pack(side = TOP)

        incolor = incolor1

        # paint in all the usual input variables in one color, then the
        # site variables in another.  
        
        for field in fields :
                row = Frame(root)
                lab = Label(row, width=11, text = field)
                var = StringVar()
                if field == 'sitecode' :
                        var.set('mdm')
                elif field == 'Time is:' :
                        row.pack(side = TOP, fill = X, expand = YES)
                        lab.pack(side = LEFT)
                        but = Radiobutton(row,text = "Local", variable = var, value = 'y')
                        but.pack(side = RIGHT, expand = YES, fill = X)
                        but = Radiobutton(row,text = "UT",variable = var, value = 'n')
                        but.pack(side = RIGHT, expand = YES, fill = X)
                        var.set('y')
                else :
                        if field == "obs_name" : # put in an extra label here ...
                                extra = Frame(root)
                                extralab = Label(extra, width = 20, height=2, text = "Site Parameters")
                                extra.pack(side = TOP, fill = X, expand = YES)
                                extralab.pack()
                                incolor = incolor2  # AND change the color ... 
                        row.pack(side = TOP, expand = YES)
                        lab.pack(side = LEFT, anchor = W)
                        ent = Entry(row, width=21)
                        ent.pack(side=RIGHT, expand = YES, fill = X)
                        ent.config(textvariable = var,bg = incolor)
 
                # Let's 'wake up' in J2000 coordinates.

                if field == 'RA' :
                        var.set(obs.cel_J2000.ra.to_string(unit = u.hourangle, sep = ' ',precision = 2, 
                             pad = True))
                elif field == 'dec' :
                        var.set(obs.cel_J2000.dec.to_string(sep = ' ',precision = 1, pad = True, 
                             alwayssign = True))
                elif field == 'equinox' :
                        var.set('%8.3f' % (equinox_to_float(obs.cel_J2000.frame)))

                # Fill in site parameters (likely the default at this point)

                elif field == 'longit' :
                        var.set(obs.site.location.lon.to_string(unit = u.deg, 
                           decimal = True, pad = True, precision = 3))
                elif field == 'lat' :
                        var.set(obs.site.location.lat.to_string(unit = u.deg, 
                           decimal = True, pad = True, alwayssign = True, 
                           precision = 3))
                elif field == 'stdz' :
                        var.set("%s" % obs.site.tzstr)
                elif field == 'elevsea' :
                        var.set("%4.0f" % obs.site.location.height.value)
                elif field == 'elevhoriz' :
                        var.set("%4.0f" % obs.site.elevhoriz.value)
                elif field == 'eleveast' :
                        var.set("%4.0f" % obs.site.height_above_east.value)
                elif field == 'elevwest' :
                        var.set("%4.0f" % obs.site.height_above_west.value)
                elif field == 'obs_name' :
                        var.set(obs.site.name)
                elif field == 'std_abbrev' :
                        var.set(obs.site.localtzabbrev)
                elif field == 'dst__abbrev' :
                        var.set(obs.site.localtzdayabbrev)
                elif field == 'date' :
                        var.set(localt.strftime("%Y %m %d   %a"))
                elif field == 'time' :
                        var.set(localt.strftime("%H %M %S"))
                elif field == 'timestep' :
                        var.set("1   h")
                        StepBox = ent
                elif field == 'objname' :
                        var.set('(initializes at zenith)')
                        ObjBox = ent
        
                # in 2018 August shifted from a list of fields to a dictionary
                # of fields, which is vastly more convenient.

                variables[field] = var

        return (variables,StepBox,ObjBox) 

def fillhaform(havars, haboxes)  :

        scratch = deepcopy(obs)   

        localtimestr = scratch.calstring(stdz = scratch.stdz, \
                use_dst = scratch.use_dst)
        
        x = string.split(localtimestr)
        ymd = x[0] + " " + x[1] + " " + x[2]
        if float(x[3]) >= 12. :
                midnstring = ymd + " 23 59 59.99"
        else :
                midnstring = ymd + " 0 0 0"
        # print "midnstring = ",midnstring,
        jdmid = time_to_jd(midnstring, stdz = scratch.stdz, \
                 use_dst = scratch.use_dst)
        # print "jdmid = ",jdmid

        one_hr = 1 / 24.
        one_sec = 1/86400.

        j = 0
        i = -10   # number of hours away from midnight  
        done = 0
        started = 0
        while i < 18 and j < 18 and done == 0 :
                scratch.jd = jdmid + i * one_hr + one_sec

                scratch.computesky()
                scratch.computesunmoon()

                if scratch.altsun < 8. :
                        # print "j = ",j,"ut = ",scratch.calstring()
                        for k in range(0,3) : # needs to be re-set if it had been blanked
                                haboxes[k][j].configure(bg = normalentry)  
                        if j == 0 :
                                started = 1
                        localdatestring = scratch.calstring(stdz = scratch.stdz,
                                use_dst = scratch.use_dst, style = 1, \
                                print_day = 1, daycomma = 0, dayabbrev = 1)
                        x = string.split(localdatestring)
                        localtimestr = x[6] + " " + x[1] \
                                + " " + x[2] + "   " +  x[3] + ":" + x[4]
                        havars[0][j].set(localtimestr)
        
                        utdatestring = scratch.calstring()
                        x = string.split(utdatestring)
                        uttimestr = x[3] + ":" + x[4]
                        havars[1][j].set(uttimestr)
                                
                        # print havars[0][j].get(),havars[1][j].get()
        
                        scratch.computesky()
                        scratch.computesunmoon()

                        hawarn = hawarningcolor(scratch, Angle(6.*u.hour))
                        moonwarn = moonwarningcolor(scratch.lunsky, scratch.altsun, scratch.obj_moon, \
                                scratch.altmoon, scratch.altit) 
                        twilightwarn = twilightwarningcolor(scratch.altsun, scratch.ztwilight)
                        airmasswarn = airmasswarningcolor(scratch, Angle(6.* u.hr))
                
                        
                        havars[2][j].set(scratch.sidereal.tripletstring(showsign = 0,places = -2, delin = ':'))
                        haboxes[3][j].configure(bg = hawarn)
                        havars[3][j].set(scratch.hanow.tripletstring(showsign = 1,places = -2, delin = ':'))
                        haboxes[4][j].configure(bg = airmasswarn)
                        if scratch.altit > 0. :
                                if scratch.airmass < 10. :
                                        havars[4][j].set("%6.3f" % scratch.airmass)
                                else :
                                        havars[4][j].set("> 10.")
                        else :
                                havars[4][j].set("(down)")
                        haboxes[5][j].configure(bg = moonwarn)
                        if(scratch.altmoon > -2.) :
                                havars[5][j].set("%5.1f" % scratch.altmoon)
                        else :
                                havars[5][j].set("---")
                        haboxes[6][j].configure(bg = twilightwarn)
                        if(scratch.altsun >= -18.) :
                                havars[6][j].set("%5.1f" % scratch.altsun)
                        else :
                                havars[6][j].set("---")
                        j = j + 1       
                elif started == 1 :
                        done = 1
                i = i + 1  
        while j < 18 :
                for k in range(0,7) :
                        havars[k][j].set("---")
                        haboxes[k][j].configure(bg = "darkgrey")
                j = j + 1

def printhatable(invars,havars)  :

        # Based closely on fillhaform, this version dumps the 
        # output to a file called "skycalc.out"

        outf = open("skycalc.out","a")

        scratch = deepcopy(obs)   # careful!!
        localtimestr = scratch.calstring(stdz = scratch.stdz, \
                use_dst = scratch.use_dst)
        
        x = string.split(localtimestr)
        ymd = x[0] + " " + x[1] + " " + x[2]
        if float(x[3]) >= 12. :
                midnstring = ymd + " 23 59 59.99"
        else :
                midnstring = ymd + " 0 0 0"
        # print "midnstring = ",midnstring,
        jdmid = time_to_jd(midnstring, stdz = scratch.stdz, \
                 use_dst = scratch.use_dst)
        # print "jdmid = ",jdmid

        scratch.jd = jdmid

        one_hr = 1 / 24.
        one_sec = 1/86400.

        outf.write("--- Hourly airmass for %s --- %s ---\n\n" % (invars[0].get(),ymd))
        outf.write("  Input coords: %s\n" % scratch.summarystring())
        outf.write("Current coords: %s\n\n" % scratch.precess(scratch.julian_epoch()).summarystring())

        outf.write("Local date & time    UT     LST      HA     airm   moonalt  sunalt\n\n")
        j = 0
        i = -10   # number of hours away from midnight  
        done = 0
        started = 0
        while i < 18 and done == 0 :
                scratch.jd = jdmid + i * one_hr + one_sec

                scratch.computesky()
                scratch.computesunmoon()

                if scratch.altsun < 8. :
                        for k in range(0,3) : # needs to be re-set if it had been blanked
                                haboxes[k][j].configure(bg = normalentry)  
                        if j == 0 :
                                started = 1
                        localtimestr = scratch.calstring(stdz = scratch.stdz, use_dst = scratch.use_dst, style = 2)
                        outf.write("%s  " % localtimestr)
        
                        utdatestring = scratch.calstring()
                        x = string.split(utdatestring)
                        uttimestr = x[3] + ":" + x[4]
                        outf.write(" %05s " % uttimestr)
                                
                        scratch.computesky()
                        scratch.computesunmoon()

                        outf.write(" %06s " % scratch.sidereal.tripletstring(showsign = 0,places = -2, delin = ':'))
                        outf.write(" %06s " % scratch.hanow.tripletstring(showsign = 1,places = -2, delin = ':'))
                        if scratch.altit > 0. :
                                if scratch.airmass < 10. :
                                        outf.write(" %6.3f " % scratch.airmass)
                                else :
                                        outf.write( "  > 10. ")

                        else :
                                outf.write(" (down) ")
                        if(scratch.altmoon > -2.) :
                                outf.write("   %5.1f " % scratch.altmoon)
                        else :
                                outf.write("     --- ")
                        if(scratch.altsun >= -18.) :
                                outf.write("  %5.1f " % scratch.altsun)
                        else :
                                outf.write("    --- ")
                        outf.write("\n")
                        j = j + 1       
                elif started == 1 :
                        done = 1
                i = i + 1  

        outf.write("\n\n")
        outf.close()

def makeplanetform(planetwin) :

        pl_name = {} # make the variables
        pl_ra   = {} 
        pl_dec  = {}
        pl_ha   = {}
        pl_airm = {}
        pl_prox = {}

        pl_name_box = {} # make the box variables
        pl_ra_box   = {}
        pl_dec_box  = {}
        pl_ha_box   = {}
        pl_airm_box = {}
        pl_prox_box = {}

        plpan = Frame(planetwin)
        plpan.pack(side = LEFT)
        row = Frame(plpan)
        row.pack(side = TOP)
        lab = Label(row,width=11,text = "Name   ")
        lab.pack(side = LEFT,fill=Y)
        lab = Label(row,width=8,text = "RA   ")
        lab.pack(side = LEFT,fill=Y)
        lab = Label(row,width=8,text = "Dec   ")
        lab.pack(side = LEFT,fill=Y)
        lab = Label(row,width=8,text = "HA   ")
        lab.pack(side = LEFT,fill=Y)
        lab = Label(row,width=8,text = "airmass ")
        lab.pack(side = LEFT,fill=Y)
        lab = Label(row,width=8,text = "proximity ")
        lab.pack(side = LEFT,fill=Y)

        for p in thorconsts.PLANETNAMES : 

            row = Frame(plpan)
            row.pack(side = TOP,fill=X,expand=YES)
            name = StringVar()
            ent = Entry(row,width = 9)
            ent.pack(side = LEFT, expand=YES, fill=X)
            ent.config(textvariable = name)
            pl_name[p] = name
            pl_name_box[p] = ent
                        
            ra = StringVar()
            ent = Entry(row,width = 9)
            ent.pack(side = LEFT, expand=YES, fill=X)
            ent.config(textvariable = ra)
            pl_ra[p] = ra 
            pl_ra_box[p] = ent 
        
            dec = StringVar()
            ent = Entry(row,width = 9)
            ent.pack(side = LEFT, expand=YES, fill=X)
            ent.config(textvariable = dec)
            pl_dec[p] = dec 
            pl_dec_box[p] = ent 
        
            ha = StringVar()
            ent = Entry(row,width = 9)
            ent.pack(side = LEFT, expand=YES, fill=X)
            ent.config(textvariable = ha)
            pl_ha[p] = ha
            pl_ha_box[p] = ent 
        
            airm = StringVar()
            ent = Entry(row,width = 9)
            ent.pack(side = LEFT, expand=YES, fill=X)
            ent.config(textvariable = airm)
            pl_airm[p] = airm 
            pl_airm_box[p] = ent 

            prox = StringVar()
            ent = Entry(row,width = 9)
            ent.pack(side = LEFT, expand=YES, fill=X)
            ent.config(textvariable = prox)
            pl_prox[p] = prox 
            pl_prox_box[p] = ent 
        
        row = Frame(plpan)
        row.pack(side=TOP)
        Button(row, text = 'Hide',command = (lambda: planetwin.withdraw())).pack(side = LEFT)

        plvars = [pl_name, pl_ra, pl_dec, pl_ha, pl_airm, pl_prox]
        plboxvars = [pl_name_box, pl_ra_box, pl_dec_box, pl_ha_box, 
                        pl_airm_box, pl_prox_box]
        return [plvars,plboxvars]

def fillplform(plvars, plboxvars) :
        
        # planets = computeplanets(obs.jd,obs.longit,obs.lat,0)
        
        plobs = deepcopy(obs)  # this gets site and date info ... 

        warning_planet = []

        # get a version of the current coordinates that is 
        # explicitly earth-centered, so it can be compared to the 
        # planetary positions without the God-awful origin shift
        # nonsense.

        geocentfrnow = currentgeocentframe(obs.t)
        testpos = obs.celest.transform_to(frame = geocentfrnow)
        
        for p in thorconsts.PLANETNAMES :
            #print "p, obs.planetdict[p]",p,obs.planetdict[p]
            plobs.setcelest(obs.planetdict[p])
            #print "plobs.celest", plobs.celest
            # these are both OK apparently.
            
            #print "p ra dec",p, plobs.cel_icrs.ra.to_string(unit = u.hourangle, 
            #   sep = ' ', precision = 1, pad = True), plobs.cel_icrs.dec.to_string(sep = ' ',
            #   precision = 0, pad = True, alwayssign = True)
            plobs.computesky()
            #print "now ra dec",p, plobs.celnow.ra.to_string(unit = u.hourangle, 
            #   sep = ' ', precision = 1, pad = True), plobs.celnow.dec.to_string(sep = ' ',
            #   precision = 0, pad = True, alwayssign = True)
            plvars[0][p].set(p.capitalize())
            plvars[1][p].set(plobs.cel_J2000.ra.to_string(unit = u.hourangle, 
               sep = ' ', precision = 1, pad = True)) 
            plvars[2][p].set(plobs.cel_J2000.dec.to_string(sep = ' ',
               precision = 0, pad = True, alwayssign = True))
            plvars[3][p].set(plobs.hanow.to_string(unit = u.hourangle,
               sep = ' ', precision = 0, alwayssign = True))
            if plobs.altit > (0. * u.deg)  :
                if plobs.airmass < 10. :
                    plvars[4][p].set("%7.2f" % plobs.airmass)
                else : plvars[4][p].set(" > 10.")
            else : plvars[4][p].set("(Down.)")

            prox = testpos.separation(plobs.celnow) 
#            print " ********* ", p.capitalize(), " ********* "
#            print "testpos",testpos
#            print "plobs.celnow",plobs.celnow
#            print "prox:",prox

            if prox < 1. * u.deg :
               for j in range(0,6) :
                   plboxvars[j][p].configure(bg = "orange")
               warning_planet.append(p)
            else :
               for j in range(0,6) :
                   plboxvars[j][p].configure(bg = normalentry)
            plvars[5][p].set("%6.1f" % prox.deg)

        return warning_planet

def makehaform(hawin) :

        ha_local = []  # make the variables
        ha_ut = []
        ha_lst = []
        ha_ha = []    # ha ha! 
        ha_airm = []
        ha_moonalt = []
        ha_sunalt = []

        ha_local_box = []  # make the variables
        ha_ut_box = []
        ha_lst_box = []
        ha_ha_box = []    # ha ha! 
        ha_airm_box = []
        ha_moonalt_box = []
        ha_sunalt_box = []

        
        hapan = Frame(hawin)
        hapan.pack(side = LEFT)
        row = Frame(hapan)
        row.pack(side = TOP)
        lab = Label(row,width=17,text = " Local")
        lab.pack(side = LEFT)
        lab = Label(row,width=6,text = " UT ")
        lab.pack(side = LEFT)
        lab = Label(row,width=6,text = " LST ")
        lab.pack(side = LEFT)
        lab = Label(row,width=6,text = " HA ")
        lab.pack(side = LEFT)
        lab = Label(row,width=7,text = " Airmass")
        lab.pack(side = LEFT)
        lab = Label(row,width=6,text = " moonalt")
        lab.pack(side = LEFT)
        lab = Label(row,width=6,text = " sunalt")
        lab.pack(side = LEFT)
        
        for i in range(0,18) :
                row = Frame(hapan)
                row.pack(side = TOP, fill = X, expand = YES)

                local = StringVar()
                ent = Entry(row, width = 17)
                ent.pack(side = LEFT, expand = YES, fill = X)
                ent.config(textvariable = local)
                ha_local.append(local) 
                ha_local_box.append(ent) 

                ut = StringVar()
                ent = Entry(row, width = 6)
                ent.pack(side = LEFT, expand = YES, fill = X)
                ent.config(textvariable = ut)
                ha_ut.append(ut) 
                ha_ut_box.append(ent)

                lst = StringVar()
                ent = Entry(row, width = 6)
                ent.pack(side = LEFT, expand = YES, fill = X)
                ent.config(textvariable = lst)
                ha_lst.append(lst) 
                ha_lst_box.append(ent)

                ha = StringVar()
                ent = Entry(row, width = 6)
                ent.pack(side = LEFT, expand = YES, fill = X)
                ent.config(textvariable = ha)
                ha_ha.append(ha)
                ha_ha_box.append(ent)

                airm = StringVar()
                ent = Entry(row, width = 7)
                ent.pack(side = LEFT, expand = YES, fill = X)
                ent.config(textvariable = airm)
                ha_airm.append(airm)
                ha_airm_box.append(ent)

                moonalt = StringVar()
                ent = Entry(row, width = 6)
                ent.pack(side = LEFT, expand = YES, fill = X)
                ent.config(textvariable = moonalt)
                ha_moonalt.append(moonalt)
                ha_moonalt_box.append(ent)

                sunalt = StringVar()
                ent = Entry(row, width = 6)
                ent.pack(side = LEFT, expand = YES, fill = X)
                ent.config(textvariable = sunalt)
                ha_sunalt.append(sunalt)
                ha_sunalt_box.append(ent)

        havars = [ha_local,ha_ut,ha_lst,ha_ha,ha_airm,ha_moonalt,ha_sunalt]
                
        row = Frame(hapan)
        row.pack(side = TOP)

        Button(row, text = 'Hide',command = (lambda: hawin.withdraw())).pack(side = LEFT)
        Button(row, text = 'Dump to file "skycalc.out"',command = (lambda i=invars,h=havars : \
                        printhatable(i,h))).pack(side = LEFT)

        return [havars, \
                [ha_local_box,ha_ut_box,ha_lst_box,ha_ha_box,ha_airm_box,ha_moonalt_box,ha_sunalt_box]] 


# Alt coords window for such things as galactic coordinates and 
# ecliptic coords, PM updated coordinates, J2000 coords in
# decimal degrees, etc.

coofields = ('Current RA:','Current dec:','equinox :',
  'Proper motion','PM is:','input epoch','RA (pm only)',
  'Dec (pm only)','equinox_:',
  'Ecliptic latitude','Ecliptic longitude','Galactic latitude',
  'Galactic longitude',
  'RA 2000 (deg)', 'dec 2000 (deg)')

# In astropy version omitted 'parallax factors' andd 'TAI (Barycentric)'

def makecooform(coowin,coofields) :
        
        coovars = {}
        cooboxes = {}
        
        coopan = Frame(coowin)
        coopan.pack(side = LEFT)

        inputfields = ["Proper motion","input epoch","Galactic longitude", \
                "Galactic latitude"]

        for field in coofields :
                row = Frame(coopan)
                lab = Label(row,width=18,text=field)
                var = StringVar()
                row.pack(side = TOP, fill = X, expand = YES)
                if field == "Proper motion" : var.set("0.  0.       [mas/yr]")
                if field == "input epoch" : var.set("2000.")
                if field == "PM is:" : 
                        lab.pack(side = LEFT)
                        but = Radiobutton(row,text="dX/dt",variable = var,
                                value = 'x')
                        but.pack(side = LEFT, expand = YES, fill = X)
                        but = Radiobutton(row,text= "d(alpha)/dt",variable=var,
                                value = 'a')
                        but.pack(side = LEFT, expand = YES, fill = X)
                        but = Radiobutton(row,text= "mu/theta",variable=var,
                                value = 'p')
                        but.pack(side = LEFT, expand = YES, fill = X)
                        var.set('x')
                        # print "var set to", var.get()
                else :
                        ent = Entry(row, width=21)
                        lab.pack(side = LEFT, anchor = W)
                        ent.config(textvariable = var)
                        try :
                                ind = inputfields.index(field)
                                ent.config(bg = incolor1)  # globally defined
                        except :
                                ent.config(bg = normalentry)
                        ent.pack(side = LEFT)
                        cooboxes[field] = ent

                coovars[field] = var 

                
        Button(coopan, text = 'Hide', command=(lambda: coowin.withdraw())).pack(side = TOP)

        # print("Alt coords: returning from makecooform.")

        return [coovars, cooboxes]

def fillcooform(coovars) :

        scratch = deepcopy(obs)
        currentep = scratch.julyear    

        pmoption = coovars['PM is:'].get()     
        # print "pmoption = ",pmoption
        inputep = float(coovars['input epoch'].get())
        
        delta_ep = currentep - inputep  # how many years of PM
        cosdelta = np.cos(scratch.celest.dec)
        # print "cosdelta = ",cosdelta
        mas_per_hr = 5.40e7
        mas_per_deg = 3.60e6

        x = coovars['Proper motion'].get().split()
        try :
                pm1 = float(x[0])
                pm2 = float(x[1])
        except :
                print("Cannot parse proper motion, two fields needed.")
                pm1 = 0.
                pm2 = 0.

        # print "pm1 pm2 ",pm1,pm2
        
        if (pm1 != 0. or pm2 != 0.) and cosdelta > 0.0001 :
                # parse out the various PM options
                if pmoption == 'x'  :  # dX/dt, dY/dt
                        dra = delta_ep * pm1 / (mas_per_hr * cosdelta)
                        ddec = delta_ep * pm2 / mas_per_deg
                elif pmoption == 'a' : # d(alpha)/dt, dY/dt
                        dra = delta_ep * pm1 / mas_per_hr
                        ddec = delta_ep * pm2 / mas_per_deg
                elif pmoption == 'p' :
                        theta = pm2 / 57.2957795130823
                        dra = delta_ep * pm1 * np.sin(theta) / (mas_per_hr * cosdelta)
                        ddec = delta_ep * pm1 * np.cos(theta) / mas_per_deg
                newra = scratch.celest.ra.hour + dra
                newdec = scratch.celest.dec.deg + ddec

        # Adjusting a SkyCoord is tricky.  See:
        # http://docs.astropy.org/en/stable/coordinates/inplace.html

                scratch.celest.cache.clear()     
                # and need to use actual attribute names rather than their
                # aliases 'ra' and 'dec' ...
                scratch.celest.data.lon[()] = newra * u.hourangle
                scratch.celest.data.lat[()] = newdec * u.deg

        #       print "adjusting dra ddec = ",dra,ddec

        if cosdelta == 0. :
                print("Can't adjust for proper motion, dec = +- 90!")
                
        # It's a safe assumption that obs.nowfr will have been updated
        # with a computesksy
        outcoo = scratch.celest.transform_to(obs.nowfr)
        coovars['Current RA:'].set(outcoo.ra.to_string(unit=u.hour,precision=2,sep=' ',pad=True))
        coovars['Current dec:'].set(outcoo.dec.to_string(precision=1,sep = ' ',pad = True,
                                   alwayssign=True)) 
        coovars['equinox :'].set("%s %s" % (obs.nowfr.name, obs.nowfr.equinox))

        # input equinox, but current epoch ... 
        coovars['RA (pm only)'].set(scratch.celest.ra.to_string(unit=u.hour,precision=2,sep=' ',
                                    pad = True))
        coovars['Dec (pm only)'].set(scratch.celest.dec.to_string(precision=1,sep=' ',pad=True,
                                    alwayssign = True))
        coovars['equinox_:'].set("%s" % scratch.celest.frame.name)


        eclipt = scratch.celest.transform_to('geocentrictrueecliptic')
        coovars['Ecliptic latitude'].set("%6.2f" % (eclipt.lat.deg))
        coovars['Ecliptic longitude'].set("%6.2f" % (eclipt.lon.deg))

        coovars['Galactic latitude'].set("%6.2f" % (scratch.celest.galactic.b.deg))
        coovars['Galactic longitude'].set("%6.2f" % (scratch.celest.galactic.l.deg))

        coovars['RA 2000 (deg)'].set("%10.5f" % (scratch.cel_J2000.ra.deg))
        coovars['dec 2000 (deg)'].set("%10.5f" % (scratch.cel_J2000.dec.deg))

        # print("Alt coords: returning from fillcooform.")

#        leapsecs = getleapsec(scratch.jd)
#        coovars[16].set("%12.5f  [+%5.2f]" % (scratch.baryjd + leapsecs / 86400., leapsecs))

########  Nightly almanac window.  ###########

# Same scheme.  Note we have both "Moon :", and "Moon_:" because we want to be
# able to swap the order of rise and set to keep them in time order.

almfields = ('Sunset','Twilight Ends','Night center','Twilight Begins','Sunrise',\
                'Moon :','Moon_:')

def makealmform(almwin, almfields) :
        
        almvars = {}
        
        almpan = Frame(almwin)
        almpan.pack(side = LEFT)
        
        for field in almfields :
                row = Frame(almpan)
                lab = Label(row,width=12,text=field)
                var = StringVar()
                ent = Entry(row, width=21)
                ent.configure(bg = normalentry)
                row.pack(side = TOP, expand = YES)
                lab.pack(side = LEFT, anchor = W)
                ent.config(textvariable = var)
                ent.pack(side = LEFT)
                almvars[field] = var   # 2018 revision - almvars is a dictionary
        
        Button(almpan, text = 'Hide', command=(lambda: almwin.withdraw())).pack(side = TOP)

        return almvars

        
# almfields = ('Sunset','Twilight Ends','Night center','Twilight Begins','Sunrise',\
#               'Moon :','Moon_:')
# 0 = sunset, 1 = eve twi, 3 = night center, 4 = morn. twilight, 5 = sunrise,
# 6, 7 = moon rise and/or set in time order.

def fillalmform(almvars) :

        if obs.sunsetflag == 0 :   # sunrise and set do occur
                almvars['Sunset'].set(
                    time_rounded_to_minute(obs.tsunset.to_datetime(timezone = obs.site.localtz),
                          incl_date = True, incl_day = True))
                almvars['Sunrise'].set(
                    time_rounded_to_minute(obs.tsunrise.to_datetime(timezone = obs.site.localtz),
                         incl_date = True, incl_day = True))
                almvars['Night center'].set(
                    time_rounded_to_minute(obs.tnightcenter.to_datetime(timezone = obs.site.localtz), incl_date = True,
                        incl_day = True))
        elif obs.sunsetflag > 0 :
                almvars['Sunset'].set("Sun up all night")
                almvars['Sunrise'].set("Sun up all night")
                almvars['Night center'].set("Sun up all night")
        else :
                almvars['Sunset'].set("Sun never rises.")
                almvars['Sunrise'].set("Sun never rises.")
                almvars['Night center'].set("Sun never rises.")
                
        if obs.twilightflag == 0 : # twilight does occur 
                almvars['Twilight Ends'].set(time_rounded_to_minute(
               obs.tevetwi.to_datetime(timezone = obs.site.localtz), incl_date = True, incl_day = True))
                almvars['Twilight Begins'].set(time_rounded_to_minute(
               obs.tmorntwi.to_datetime(timezone = obs.site.localtz), incl_date = True, incl_day = True))
        elif obs.twilightflag > 0 :
                almvars['Twilight Ends'].set("Twilight all night.")
                almvars['Twilight Begins'].set("Twilight all night.")
        else :
                almvars['Twilight Ends'].set("No twilight.")
                almvars['Twilight Begins'].set("No twilight.")

        moonrisestring = "Rise: " + time_rounded_to_minute(obs.tmoonrise.to_datetime(timezone = obs.site.localtz),
                incl_date = True, incl_day = True, incl_year = False)
        moonsetstring = " Set: " + time_rounded_to_minute(obs.tmoonset.to_datetime(timezone = obs.site.localtz),
                incl_date = True, incl_day = True, incl_year = False)

        if obs.moonsetflag == 0 :
                # list rise and set variables in order.  "Moon" and "Moon_" are
                # distinguished because "Moon" comes first on the display layout.
                if obs.tmoonrise < obs.tmoonset :
                        risekey = 'Moon :'
                        setkey = 'Moon_:' 
                else :
                        risekey = 'Moon_:' 
                        setkey = 'Moon :' 
                if abs(obs.tmoonrise - obs.tnightcenter) < 0.5 * u.day:
                        almvars[risekey].set(moonrisestring)
                else :
                        almvars[risekey].set("---")
                if abs(obs.tmoonset - obs.tnightcenter) < 0.5 * u.day:
                        almvars[setkey].set(moonsetstring)
                else :
                        almvars[setkey].set("---")
        else :
                if obs.moonsetflag > 0 :
                        moonrisestring = "Moon may stay up."
                elif obs.moonsetflag < 0 : 
                        moonrisestring = "Moon may stay down."
#               elif abs(obs.jdmoonrise - -100.) < 0.01 :
#                       moonrisestring = "Moonrise calcn error"
#               if abs(obs.jdmoonset -1000.) < 0.01 :
#                       moonsetstring = "Moon may stay up?"
#               elif abs(obs.jdmoonset - -1000.) < 0.01 :
#                       moonsetstring = "Moon may stay down?"
#               elif abs(obs.jdmoonset - -100.) < 0.01 :
#                       moonsetstring = "Moonset calcn error"
                almvars['Moon :'].set(moonrisestring)
                almvars['Moon_:'].set(moonsetstring)

def makesiteform(sitewin, sitevar) :

        sitepan = Frame(sitewin)
        sitepan.pack(side = LEFT)
        # root = Frame(win)
        # root.pack()
        row = Frame(sitepan)
        # lab = Label(row,text = "Site menu")
        row.pack(side = TOP)
#       lab = Label(row,text = "(Takes effect on refresh.)")
#       lab.pack(side = TOP)
        # lab.pack(side = TOP)

        for k in sitedict.keys() :
                but = Radiobutton(sitepan,text = sitedict[k].name, variable = sitevar,
                                        value = k)
                but.pack(side = TOP, anchor = W)
        but = Radiobutton(sitepan,text = "Other (allow user entries)", variable = sitevar,
                                        value = 'x')    
        but.pack(side = TOP, expand = YES, fill = X)

        Button(sitepan, text = 'Hide',command = (lambda: sitewin.withdraw())).pack(side = TOP)
        sitevar.set('mdm')

# This makes the right side of the main window.  Again, nearly everything
# is a text-entry field, which are saved as a dictionary keyed to the name
# (label) of the field.  Nearly everything here is an output variable, though
# JD can be force with a carriage return.

outfields = ('sidereal','ha','airmass','alt_az','parallactic','jd', \
        'sunradec','sunaltaz','ztwilight','moonphase','moonradec','moonaltaz', \
        'illumfrac','lunsky','moon-obj ang.','baryjd','baryvcor','constel',
        'planet_proxim')


def makeoutform(root, fields) :
        
#       variables = []

        variables = {}  # 2018 Aug -- change variables and entryboxes to dictionaries

        entryboxes = {}

#       calstr = obs.calstring(style = 1, stdz = obs.stdz, use_dst = obs.use_dst) 
#       x = string.split(calstr)
#       ymd = x[0] + " " + x[1] + " " + x[2]
#       hms = x[3] + " " + x[4] + " " + x[5]

        row = Frame(root)
        lab = Label(row, text = "Output Variables")
        row.pack(side = TOP)
        lab.pack(side = TOP)
        
        for field in fields :
                row = Frame(root)
                lab = Label(row, width=12, text = field)
                ent = Entry(row, width=23)
                row.pack(side = TOP, fill = X, expand = YES)
                lab.pack(side = LEFT)
                ent.pack(side=RIGHT, expand = YES, fill = X)
                if field == "jd" :  # signal it's an input ...
                        ent.configure(bg = incolor1)
                else : ent.configure(bg = normalentry)
                var = StringVar()
                ent.config(textvariable = var)
                # variables.append(var)
                variables[field] = var
                entryboxes[field] = ent 
        return (variables, entryboxes)

# This function is bound to a carriage-return in the JD box, which forces
# the JD to the value you enter.
        
def force_jd(variables,outvars,havars,almvars,coovars,planetvars,
                boxes,haboxes,planetboxes, skdisp, airmdisp) :

        testin = float(outvars['jd'].get())
#       if testin > 2488069. :
#               print "Bad JD input, past 2099."
#               return()
#       elif testin < 2415388. :
#               print "Bad JD input, before 1901."
#               return()

        obs.settime(float(testin)) 

        if variables['Time is:'].get() == 'y' :
            tdout = obs.t.to_datetime(timezone = obs.site.localtz)
        else :
            tdout = obs.t.to_datetime()

        variables['date'].set(tdout.strftime("%Y %m %d  %a"))
        variables['time'].set(tdout.strftime("%H %M %S"))

        circumstances(variables,outvars,havars,almvars,coovars,planetvars,
                boxes,haboxes,planetboxes, skdisp, airmdisp) 

# This function sets the coords in the main window to galactic coords
# entered in the coordinate window.  It will be bound to carriage 
# returns in either of those boxes.

def force_galact(variables,outvars,havars,almvars,coovars,planetvars,
                boxes,haboxes,planetboxes, skdisp, airmdisp) :

        glongin = float(coovars['Galactic longitude'].get())
        glatin = float(coovars['Galactic latitude'].get())

        inputcel = SkyCoord(glongin * u.deg, glatin * u.deg, frame = 'galactic')
        tempcel = inputcel.transform_to('icrs')
        obs.setcelest(tempcel)
        
        variables['RA'].set(obs.celest.ra.to_string(unit=u.hour,precision=2,pad=True,sep=' '))
        variables['dec'].set(obs.celest.dec.to_string(precision=1,pad=True,alwayssign=True,sep=' '))
        variables['equinox'].set("2000.")

        circumstances(variables,outvars,havars,almvars,coovars,planetvars,
                boxes,haboxes,planetboxes, skdisp, airmdisp) 


### OBJECT LIST SECTION ... looks primitive in 2018, high on the list to rewrite.


objs2000 = {}  # globally defined, for plotting

class ObjList(ScrolledList) :
        def __init__(self,names,objs,variables,parent) :
                self.lab = Label(parent,text = "Double click to select:")
                self.lab.pack(side = TOP)
                Frame.__init__(self,parent)
                self.parent = parent
                self.pack()
                self.makeWidgets(names)
                self.objs = objs
                self.objnames = names       # keep a list for sorting etc ... 
                self.objnames_orig = names[:]  # keep a copy in order
                self.variables = variables  # need to keep a copy with object ...
                self.pack()
                
                buttonpan = Frame(parent)
                buttonpan.pack(side=LEFT)
                row1 = Frame(buttonpan)
                row1.pack()
                self.but1 = Button(row1, text='Alphabetize', 
                        command=(lambda : self.alphabetical()))
                self.but1.pack(side = LEFT, fill=Y)
                self.but2 = Button(row1, text='List\nby RA', 
                        command=(lambda : self.by_ra()))
                self.but2.pack(side = LEFT)
                row2 = Frame(buttonpan)
                row2.pack()
                # proximity sorting almost never useful in practice
                #self.but3 = Button(row2, text='List by\nProximity', 
                #        command=(lambda : self.by_proximity()))
                # self.but3.pack(side = LEFT)
                self.but4 = Button(row2, text='Original\norder', 
                        command=(lambda : self.restore_order()))
                self.but4.pack(side = LEFT)
                row3 = Frame(buttonpan)
                row3.pack()
                self.but5 = Button(row3, text='Dismiss\nList', 
                        command=(lambda : self.deleteList()))
                self.but5.pack(side = LEFT)
        def handleList(self,event) :
                index = self.listbox.curselection()[0]
                label = self.listbox.get(index)
                self.variables['objname'].set(label)
                # print("index, label",index,label)
                # print("objs:",self.objs)
                # print("self.objs[label]",self.objs[label])
                self.variables['RA'].set(self.objs[label].ra.to_string(unit = u.hour, precision=2,
                  pad=True,sep = ' ')) 
                self.variables['dec'].set(self.objs[label].dec.to_string(precision=1,pad=True,
                  sep = ' ',alwayssign=True))
                self.variables['equinox'].set("%8.2f" % equinox_to_float(self.objs[label].frame))
                circumstances(self.variables, outvars, havars, almvars, 
                    coovars, planetvars, outboxes, haboxes, planetboxes, skdisp, airmdisp)
        def deleteList(self) :
                # may be more than one list loaded, so deleted objs2000 
                # selectively ... 
                # print objs2000
                for k in self.objs.keys() :
                        # print "about to del ... ",k,
                        try :
                                del objs2000[k]
                                # print "deleted."
                        except KeyError :   # if on overlapping list, it's already gone
                                 # print "skipped."
                                pass
                self.objs.clear()

                self.listbox.destroy()
                self.sb.destroy()
                self.but1.destroy()
                self.but2.destroy()
                #self.but3.destroy()
                self.but4.destroy()
                self.but5.destroy()
                self.lab.destroy()
                self.parent.destroy()
        def alphabetical(self) :
                self.objnames.sort()
                self.listbox.delete(0,len(self.objnames))
                pos = 0
                for label in self.objnames :
                        self.listbox.insert(pos,label)
                        pos = pos + 1
        def restore_order(self) :
                self.objnames = self.objnames_orig[:]
                self.listbox.delete(0,len(self.objnames))
                pos = 0
                for label in self.objnames :
                        self.listbox.insert(pos,label)
                        pos = pos + 1
        def by_ra(self) :  # doesn't precess, but shouldn't matter much here.
                self.listbox.delete(0,len(self.objnames))
                radict = {} 
                for n in self.objnames :
                        ra = self.objs.ra.hour
                        # print "ra = ",ra
                        if not ra in radict : radict[ra] = n
                        else :
                                while ra in radict :
                                        ra = ra + 0.000001
                                radict[ra] = n
                
                pos = 0
                keysort = radict.keys()[:]
                keysort.sort()
                for r in keysort :
                        self.listbox.insert(pos,radict[r])
                        pos = pos + 1


class Objlistwin(Toplevel) :

        # 2018 Aug - for now keep MDM style input, i.e. free-format, 
        # space-separated format with for example

        #  my_star  22 33 44.5   -0 14 55   2000      


        def __init__(self, invars) :
                objfile = askopenfilename()
                Toplevel.__init__(self)
                entries = Frame(self)
                objs = {}
                keylist = []
                for l in open(objfile,'r').readlines() :
                        x = l.split()
                        rastr = x[1] + ' ' + x[2] + ' ' + x[3]
                        decstr = x[4] + ' ' + x[5] + ' ' + x[6]
                        eqstr = x[7]

                        # keep a globally-defined copy in equinox 2000
                        if float(eqstr) == 2000. :
                            tempcel = SkyCoord(rastr,decstr,unit = (u.hour, u.deg))
                        else :
                            eq = "J"+eqstr
                            tempcel = SkyCoord(rastr,decstr,unit=(u.hour,u.deg),
                                frame = FK5(equinox = eq))
                            tempcel.tranform_to('icrs')
                        if not x[0] in objs2000 :  #  avoid overlapping keys ... 
                                objs2000[x[0]] = tempcel
                                objs[x[0]] = tempcel
                                keylist.append(x[0]) 
                        else :
                                trialkey = x[0]
                                while trialkey in objs2000 :
                                        trialkey = '|' + trialkey + '|'
                                objs2000[trialkey] = tempcel
                                objs[trialkey] = tempcel
                                keylist.append(trialkey) 
                                # print "| ",x[0],"| loaded as |",trialkey,"|"
                objframe = ObjList(keylist, objs, invars , self)
        
def objbyname(variables,outvars,havars,almvars,coovars,planetvars,
                boxes,haboxes,planetboxes, skdisp, airmdisp) :

        # given an object name typed into the "objname" field, looks to see if it's
        # on the list and then proceeds to load the object and update everything.

        name_in = variables['objname'].get()
        name_in = string.replace(name_in," ","")
        
        if len(objs2000) == 0 : pass
                # if there's no list open, a CR in the name field does 
                # nothing.  No need to scold the user in that case ...
        elif name_in in objs2000 :
                obs.setcelest(objs2000[name_in])
                variables['RA'].set(obs.celest.ra.to_string(unit=u.hour,sep=' ',
                   pad=True))
                variables['dec'].set(obs.celest.dec.to_string(unit=u.deg,sep=' ',
                   pad=True,alwayssign=True)) 
                variables['equinox'].set("%8.3f" % equinox_to_float(obs.celest.frame))
                circumstances(variables,outvars,havars,almvars,coovars,planetvars, 
                        boxes,haboxes,planetboxes, skdisp, airmdisp) 
        else :  # a list is open, but the name is wrong -- scold.
                variables['objname'].set("%s NOT FOUND!" % name_in)

def getsimbad(variables,outvars,havars,almvars,coovars,planetvars,
                boxes,haboxes,planetboxes, skdisp, airmdisp) :

        # given an object name typed into the "objname" field, looks it up in 
        # simbad and proceeds to load it and update everything.

        name_in = variables['objname'].get()

        
# have to quote this string to get the spaces etc encoded, but it doesn't work if you 
# encode /,:,?,and =, so the second "safe" argument mus tbe used.
        myquery = urllib.quote('http://simbad.u-strasbg.fr/simbad/sim-script?script=format object "%IDLIST(1) | %COO(A D)"\n' + name_in,'/:?=')
        urlf = urllib.urlopen(myquery)
        simdata = urlf.read()
        print(simdata)
        x = simdata.split("\n")
        data_found = False
        name_coords_found = False
        for xx in x : 
                # print "line : ",xx
                if xx.find("::data") > -1 : 
                        data_found = True
                        # print "data found!"
                if data_found :
                        if xx.find("|") > -1 : 
                                # print "name/coords found!"
                                yy = xx.split("|")
                                name = yy[0].replace("NAME ","")
                                coopieces = yy[1].split()
                                rastr = coopieces[0] + " " + coopieces[1] + " " + coopieces[2]
                                decstr = coopieces[3] + " " + coopieces[4] + " " + coopieces[5]
                                # print "name, rastr, decstr"
                                name_coords_found = True
        if name_coords_found :
                variables['objname'].set(name)
                variables['RA'].set(rastr + "     hms")
                variables['dec'].set(decstr + "     dms")
                variables['equinox'].set("2000.")
                circumstances(variables,outvars,havars,almvars,coovars,planetvars, 
                        boxes,haboxes,planetboxes, skdisp, airmdisp) 
        else :  # a list is open, but the name is wrong -- scold.
                variables['objname'].set("%s - simbad failed." % name_in)

def airmass_str(airmass, places = 2) :  # stupid little utility to format an airmass
        
        if airmass < 0. :
                return "(down)"
        elif airmass >= 10. :
                return " >10. "
        else :
                if places == 2 : return "%5.2f " % airmass
                if places > 2 : return "%6.3f" % airmass
                if places < 2 : return " %4.1f " % airmass
        
class text_table_win(Toplevel) :

        # transient window for computing predicted eclipse timings.

        def __init__(self, invars) :
                Toplevel.__init__(self)
                self.title('Skycalc: Text Tasks')
                entries = Frame(self)

                fields = ['start date (UT)','end date (UT)','HJD eclipse','sigma(ecl)',
                        'P [d]','sigma(P)', 'max sunalt','max airmass',
                        'time step']

                self.vars = {}
                kount = 2
                for f in fields :
                        # arrange in two columns.
                        if kount == 2 or f == fields[:-1] :
                                row = Frame(entries) 
                                row.pack(side = TOP, expand = YES)
                                kount = 0
                        lab = Label(row,width=14,text=f) 
                        lab.pack(side = LEFT)
                        var = StringVar()
                        ent = Entry(row,width=21)
                        ent.pack(side=LEFT,expand=YES,fill=X)
                        ent.config(textvariable = var, bg = incolor1)
                        self.vars[f] = var
                        kount = kount + 1
                        if f == 'HJD eclipse' :
                                var.set('2453000.')
                        if f == 'sigma(ecl)' :
                                var.set('0.0001')
                        if f == 'P [d]' :
                                var.set('0.12345')
                        if f == 'sigma(P)' :
                                var.set('0.00001')
                        if f == 'start date (UT)' :
                                var.set('2019-01-19') 
                        if f == 'end date (UT)' :
                                var.set('2019-01-22')
                        if f == 'max sunalt' :
                                var.set('-12.')
                        if f == 'max airmass' :
                                var.set('3.')
                        if f == 'time step' :
                                var.set('1 h')
                entries.pack()
                buttonpan = Frame(self)
                buttonpan.pack()

                Button(buttonpan,text = "What's\nThis?",
                        command = (lambda: self.showhelptext())).pack(side=LEFT,fill=Y)
                Button(buttonpan,text = "Seasonal\n observability",command = 
                        (lambda: self.compute_observability())).pack(side = LEFT)
                Button(buttonpan,text = "Compute\n times",command = 
                        (lambda: self.compute_eclipses())).pack(side = LEFT)
                Button(buttonpan,text = "Compute\n phases",command = 
                        (lambda: self.compute_phases())).pack(side = LEFT)
                Button(buttonpan,text = "Erase\n output",command = 
                        (lambda: self.outwin.erasetext())).pack(side=LEFT)
                Button(buttonpan,text = "Dump output\n to 'skycalc.out'",
                        command = 
                        (lambda: self.outwin.dumptext('skycalc.out'))).pack(side=LEFT)
                Button(buttonpan,text = "Dismiss\n window",command = 
                        (lambda: self.destroy())).pack(side = LEFT)
                entries.pack()

                self.outwin = ScrolledText(self)
                
        def compute_phases(self) :

                start = Time(self.vars['start date (UT)'].get())
                end = Time(self.vars['end date (UT)'].get())

                if start > end :
                    self.outwin.appendline("end date is earlier than start date!  No output.")
                    return

                T0 = Time(float(self.vars['HJD eclipse'].get()),format = 'jd')
                #sigT0 = TimeDelta(float(self.vars['sigma(ecl)'].get()),format = 'jd')
                sigT0 = float(self.vars['sigma(ecl)'].get())
                P = TimeDelta(float(self.vars['P [d]'].get()), format = 'jd')
                #sigP = TimeDelta(float(self.vars['sigma(P)'].get()), format = 'jd')
                sigP = float(self.vars['sigma(P)'].get())
                maxsunalt = Angle(float(self.vars['max sunalt'].get()) * u.deg)
                maxairmass = float(self.vars['max airmass'].get())

                x = self.vars['time step'].get().split()
                number = float(x[0])
                mult = 3600.   # default to hours.
                if len(x) > 1 :
                        unit = x[1].upper()
                        if unit.find('S') == 0 : mult = 1.
                        if unit.find('M') == 0 : mult = 60.
                        elif unit.find('H') == 0 : mult = 3600.
                        elif unit.find('D') == 0 : mult = 86400.
                        elif unit.find('W') == 0 : mult = 604800.
                        elif unit.find('Y') == 0 : mult = 31536000.

                timeincr = TimeDelta(number * mult, format = 'sec')
        
                localfl = invars['Time is:'].get()   # 'y' or 'Y' for time-is-local

                scratch = deepcopy(obs)

                self.outwin.appendline(" ") 
                self.outwin.appendline("**** Phases for %s ****" % (invars['objname'].get()))
                self.outwin.appendline(" ") 
                self.outwin.appendline("T0 = HJD %16.7f +- %11.7f" % (T0.jd,sigT0))
                self.outwin.appendline(" P = %16.9f +- %11.9f" % (P.value,sigP))
                self.outwin.appendline(" ") 
                self.outwin.appendline("Coordinates: %s  %s  %8.2f" % (scratch.celest.ra.to_string(unit=u.hour,sep=' ',
                        precision=2,pad=True),scratch.celest.dec.to_string(sep = ' ',pad=True, precision=1, alwayssign = True),
                        equinox_to_float(scratch.celest.frame)))
                self.outwin.appendline(" ")
                self.outwin.appendline("%s, longit. %s  lat. %s." % (invars['obs_name'].get(),
                       invars['longit'].get(), invars['lat'].get()))
                self.outwin.appendline(" ")

                if localfl == 'y' or localfl == 'Y' :
                        timetype = "Local "
                else :
                        timetype = "   UT "
                headerline = " " + timetype +  \
                        "Date and time -  phase del-ph airm.    HA   moon-[alt]-sun "
                self.outwin.appendline(headerline)
                self.outwin.appendline(" ")
                
                n = 0
                t = start
                last_bary = None 
                oneday = TimeDelta(1.,format='jd') 
                
                nprint = 0
                maxlines = 200
                last_was_printed = True # or too many spaces. 

                while t < end and nprint < maxlines:
                        t = start + n * timeincr
                        # print("start:",start,"t:",t)
                        scratch.settime(t)
                        scratch.computesky(redo_coords = False)
                        scratch.computesunmoon()
                        # computebary takes nearly a second on a fast machine.
                        # refresh only after a day.
                        # ... actually, use computequickbary which is plenty good enough
                        if last_bary != None :
                            if(t - last_bary > oneday) :
                                scratch.computequickbary()
                                last_bary = t
                        else : 
                            scratch.computequickbary()
                            last_bary = t
                        tb = scratch.t + scratch.barytcorr
                        cycles = (tb - T0) / P
                     #   print("scratch.t",scratch.t,"tbary:",scratch.tbary)
                     #   print("cycles",cycles)
                        phase = cycles - np.floor(cycles)
                        sig = np.sqrt(sigT0 ** 2 + (cycles * sigP) ** 2)
                        # print("xxx = ",xxx)
                        # sig = TimeDelta(xxx, format='jd')
                        
                        dphi = sig / P.value
                     #   print('dphi = ',dphi)

#                         print "%5.0f %12.3f airm %5.2f alts %8.2f hanow %8.2f" % \
#                        # (n, tgeo, scratch.airmass, scratch.altsun, scratch.hanow.val)

                        if maxairmass > 0. and scratch.altit > 0. * u.deg  \
                             and scratch.airmass < maxairmass :  
                                airm_print = True
                        elif maxairmass < 0. :  # print all
                                airm_print = True 
                        else : airm_print = False
                        if scratch.sunaltit < maxsunalt and airm_print :
                                if not last_was_printed :
                                        self.outwin.appendline(" ")     
                                last_was_printed = True 
                                hastr = scratch.hanow.to_string(pad = True, sep=':', precision=0, alwayssign=True,
                                       fields = 2)
                                if localfl == 'y' or localfl == 'Y' :
                                    tdout = scratch.t.to_datetime(timezone = obs.site.localtz)
                                else :
                                    tdout = scratch.t.to_datetime()
                                calout = tdout.strftime("%a  %Y-%m-%d  %H:%M")
                          #      print("calout",calout)
                          #      print("phase dphi",phase, dphi)
                                outline = "%s %6.3f %6.3f " % (calout, phase, dphi)
                                if scratch.airmass < 10. and scratch.airmass > 0. :
                                    outline = outline + "%5.2f " % scratch.airmass
                                elif airmass > 0. :
                                    outline = outline + " >10. "
                                else : 
                                    outline = outline + "down. "
                                outline = outline + " %6s  %6.1f %6.1f" % \
                                        (hastr,scratch.moonaltit.deg,scratch.sunaltit.deg)
                                self.outwin.appendline(outline)
                                nprint = nprint + 1

                        else : last_was_printed = 0
                        n = n + 1

                        if nprint >= maxlines :
                                self.outwin.appendline("Halted after %d lines ...." %
                                        (maxlines))
#
        def compute_eclipses(self) :

                start = Time(self.vars['start date (UT)'].get())
                end = Time(self.vars['end date (UT)'].get())

                if start > end :
                    self.outwin.appendline("end date is earlier than start date!  No output.")
                    return

                T0 = Time(float(self.vars['HJD eclipse'].get()),format = 'jd')
                sigT0 = float(self.vars['sigma(ecl)'].get())
                P = TimeDelta(float(self.vars['P [d]'].get()), format = 'jd')
                sigP = float(self.vars['sigma(P)'].get())
                maxsunalt = Angle(float(self.vars['max sunalt'].get()) * u.deg)
                maxairmass = float(self.vars['max airmass'].get())

                localfl = invars['Time is:'].get()   # 'y' or 'Y' for time-is-local

                scratch = deepcopy(obs)

                n = np.floor((start - T0) / P) + 1.
                t = T0 + n * P
                # print("n, t",n,t)

                # variables for seeing how stale the bary corr is
                last_bary = None
                oneday = TimeDelta(1.,format='jd') 

                self.outwin.appendline(" ") 
                self.outwin.appendline("**** Ephemeris for %s ****" % (invars['objname'].get()))
                self.outwin.appendline(" ") 
                self.outwin.appendline("T0 = HJD %16.7f +- %11.7f" % (T0.jd,sigT0))
                self.outwin.appendline(" P = %16.9f +- %11.9f" % (P.value,sigP))
                self.outwin.appendline(" ") 
                self.outwin.appendline("Coordinates: %s  %s  %8.2f" % (scratch.celest.ra.to_string(unit=u.hour,sep=' ',
                        precision=2,pad=True),scratch.celest.dec.to_string(sep = ' ',pad=True, precision=1, alwayssign = True),
                        equinox_to_float(scratch.celest.frame)))
                self.outwin.appendline(" ")
                self.outwin.appendline("%s, longit. %s  lat. %s." % (invars['obs_name'].get(),
                       invars['longit'].get(), invars['lat'].get()))
                self.outwin.appendline(" ")

                if localfl == 'y' or localfl == 'Y' :
                        timetype = "Local "
                else :
                        timetype = "   UT "
                headerline = "    N   --" + timetype +  \
                        " Date and time  del-t del-ph  airm.   HA   moon-[alt]-sun "
                self.outwin.appendline(headerline)
                self.outwin.appendline(" ")
                
                nprint = 0
                maxlines = 200 
                last_was_printed = 1 # or too many spaces. 
                while t < end and nprint < maxlines:
                        t = T0 + n * P  # barycentric jd of eclipse or whatever
                        # print("n, P, t",n,P,t) 

                        scratch.settime(t)  # this is time at solar system barycenter.

                        # computebary takes nearly a second on a fast machine.
                        # refresh only after a day, to make response faster.
                        if last_bary != None :
                            if(t - last_bary > oneday) :
                           #     self.outwin.appendline("Refreshing bary, tcorr = %s " % 
                           #            (scratch.barytcorr.sec))
                                scratch.computebary()
                                last_bary = t
                        else : 
                            scratch.computebary()
                            last_bary = t
                    
                        tgeo = t - scratch.barytcorr  # this is time on earth (inverted bary corr).  

                        # print("tcorr ",scratch.barytcorr)
                        # print("geo aftere tcorr subbed",t)
                        sig = np.sqrt(sigT0 ** 2 + (n * sigP) ** 2)
                        dphi = sig / P.value
                        dtmin = sig * 1440.
                        # print("dphi, dtmin",dphi,dtmin)

                        scratch.settime(tgeo)
                        scratch.computesky(redo_coords = False)
                        scratch.computesunmoon()
                        #print("%5.0f %12.3f airm %5.2f alts %8.2f hanow %8.2f" % \
                        #  (n, tgeo.jd, scratch.airmass, scratch.sunaltit.deg, scratch.hanow.hour))

                        if maxairmass > 0. and scratch.altit > 0. * u.deg  \
                             and scratch.airmass < maxairmass :  
                                airm_print = True
                        elif maxairmass < 0. :  # print all
                                airm_print = True 
                        else : airm_print = False

                        if scratch.sunaltit < maxsunalt and airm_print :
                                if not last_was_printed :
                                        self.outwin.appendline(" ")     
                                last_was_printed = True 
                                hastr = scratch.hanow.to_string(pad = True, sep=':', precision=0,alwayssign=True,
                                   fields=2)
                                if localfl == 'y' or localfl == 'Y' :
                                    tdout = scratch.t.to_datetime(timezone = obs.site.localtz)
                                else :
                                    tdout = scratch.t.to_datetime()
                                calout = tdout.strftime("%a  %Y-%m-%d  %H:%M")
                                # print("calout",calout)
                                outline = "%7.0f  %s %5.0fm %6.3f " % (n, calout, dtmin, dphi)
                                if scratch.airmass < 10. and scratch.airmass > 0. :
                                    outline = outline + "%5.2f " % scratch.airmass
                                elif airmass > 0. :
                                    outline = outline + " >10. "
                                else : 
                                    outline = outline + "down. "
                                outline = outline + " %6s  %6.1f %6.1f" % \
                                        (hastr,scratch.moonaltit.deg,scratch.sunaltit.deg)
                                self.outwin.appendline(outline)
                                nprint = nprint + 1

                        else : last_was_printed = 0
                        n = n + 1

                        if nprint >= maxlines :
                                self.outwin.appendline("Halted after %d lines ...." %
                                        (maxlines))


        def compute_observability(self) :

                start = Time(self.vars['start date (UT)'].get())
                end = Time(self.vars['end date (UT)'].get())

                if start > end :
                    self.outwin.appendline("end date is earlier than start date!  No output.")
                    return

                scratch = deepcopy(obs)

                self.outwin.appendline(" ") 
                self.outwin.appendline("**** Seasonal Observability for %s ****" % (invars['objname'].get()))
                self.outwin.appendline(" ") 
                self.outwin.appendline("Coordinates: %s  %s  %8.2f" % (scratch.celest.ra.to_string(unit=u.hour,sep=' ',
                        precision=2,pad=True),scratch.celest.dec.to_string(sep = ' ',pad=True, precision=1, alwayssign = True),
                        equinox_to_float(scratch.celest.frame)))
                self.outwin.appendline(" ")
                self.outwin.appendline("Note: Twilight used is nautical (sun altitude -12 deg) not astronomical (-18)")
                self.outwin.appendline(" ")
                self.outwin.appendline("%s, longit. %s  lat. %s." % (invars['obs_name'].get(),
                       invars['longit'].get(), invars['lat'].get()))
                self.outwin.appendline(" ")
                self.outwin.appendline("--Moon, date--     @eve twi:      @nght ctr:    @morn twi:    nght hrs @airm:")
                self.outwin.appendline("                   HA   airm      HA   airm      HA   airm     <3    <2  <1.5")
                self.outwin.appendline(" ") 

#                # get min and max alt for coords of middle of interval
#                meanjd = (start.jd + end.jd) / 2.
#                scratch.jd = meanjd
#                scratch.computesky()
#                [min_alt, max_alt] = _skysub.min_max_alt(scratch.lat, 
#                        scratch.CoordsOfDate.dec.val)

                scratch.settime(start) 

                # find first new or full after start.jd

                (phdescr, lunage, lunation) = phase_descr(start.jd)
                # print(start.jd, phdescr, lunage, lunation)

                if lunage > 14. : nph = 2  # start with full 
                else : nph = 0
 
                jdlist = flmoon(lunation, nph)
                # print( "jd = ",jdlist)
                scratch.settime(jdlist) 
                # print( scratch.t)
         
                # now loop through until we're off the end ... 
                while jdlist < end.jd :
                         nph = nph + 2
                         if nph == 4 :
                                 lunation = lunation + 1
                                 nph = 0
                                 outline = "N: "
                         else :
                                 outline = "F: "
                         jdlist = flmoon(lunation,nph)
                         if jdlist > end.jd : break
                 
                         scratch.settime(jdlist)
#                         self.outwin.appendline("lunation %d nph %d - %s" % \
#                                (lunation,nph,scratch.t.to_datetime(timezone=scratch.site.localtz).strftime("%Y %m %d")))
 
               # self.outwin.appendline("--Moon, date--     @eve twi:       @nght ctr:     @morn twi:   nght hrs @airm:")
               # self.outwin.appendline("                    HA   airm       HA   airm      HA   airm     <3    <2  <1.5")

                         outline = outline + scratch.t.to_datetime(timezone=scratch.site.localtz).strftime("%Y %m %d  ")
 
                         scratch.computesky(redo_coords = False)
                         scratch.computesunmoon()
                         scratch.setnightevents()
                         scratch.compute_hours_up()

                         # evening 12-degree twilight
                         if scratch.twilight12flag == 0 :
                             scratch.settime(scratch.tevetwi12)
                             scratch.computesky(redo_coords = False) 
                             # print("eve12:", scratch.lst,scratch.airmass)
                             hastr = scratch.hanow.to_string(sep=':',fields=2,pad=True,alwayssign=True)
                             if scratch.airmass < 0. or scratch.airmass > 9.99 : 
                                 outline = outline + "  %s ----  " % (hastr)
                             else :
                                 outline = outline + "  %s %4.1f  " % (hastr,scratch.airmass)
                         else :  # twilight doesn't begin or end
                             outline = outline + "  ------  ----  "

                         # night center
                         scratch.settime(scratch.tnightcenter) 
                         scratch.computesky(redo_coords = False) 
                         hastr = scratch.hanow.to_string(sep=':',fields=2,pad=True,alwayssign=True)
                         if scratch.airmass < 0. or scratch.airmass > 9.99 : 
                             outline = outline + "  %s ----  " % (hastr)
                         else :
                             outline = outline + "  %s %4.1f  "  % (hastr,scratch.airmass)
                        
                         # and morning 12-degree twilight.
                         if scratch.twilight12flag == 0 :
                             scratch.settime(scratch.tmorntwi12)
                             scratch.computesky(redo_coords = False) 
                             # print("morn12:", scratch.lst,scratch.airmass)
                             hastr = scratch.hanow.to_string(sep=':',fields=2,pad=True,alwayssign=True)
                             if scratch.airmass < 0. or scratch.airmass > 9.99 : 
                                 outline = outline + "  %s ----  " % (hastr)
                             else :
                                 outline = outline + "  %s %4.1f  " % (hastr,scratch.airmass)
                         else :  # twilight doesn't begin or end
                             outline = outline + "  -----  ----  "

                         outline = outline + "  %3.1f   %3.1f   %3.1f" % \
                           (scratch.uptime30.value * 24., scratch.uptime20.value * 24., scratch.uptime15.value * 24)

                         self.outwin.appendline(outline)
                         
                         

        def showhelptext(self) :
                helptext = """
This window is for computations which create tables of text output.
Once you're happy with the output, you can dump the contents 
to a text file "skycalc.out".  The screen is a primitive text
editor, so you can annotate your output directly if you like.

The calculations take the object coordinates and site info
from the main window, and (obviously) use other information from 
the fields on this window.

"Seasonal observability" tabulates circumstances for your object
at each new and full moon between the start and end dates.  Hour angle
and airmass are printed at three times of night: evening (18 degree) 
twilight, night center, and morning twilight.  After that comes 
the number of hours per night in which (a) there is no twilight, 
and (b) the object is less than a given airmass; tabulated
values are for 3, 2, and 1.5 airmasses.  This information is 
useful in planning observing time requests.  Of the information 
entered here, only start and end dates are used.

The "Compute times" option gives observed times of a strictly
periodic phenomeon, e.g. predicted binary star eclipse times.  
All the information in the form on this page is used, except for 
the time step.  The "Compute phases" option gives the phase at 
the interval specified by the time step.

To use these, enter the heliocentric (or more correctly the
BARYcentric) JD (e.g. of an observed eclipse) and its sigma, and 
the period and its sigma.  You can optionally filter output to 
pertain only to that visible from your site, by giving a maximum 
airmass (-1 is code for don't filter) and a maximum sun altitude 
(e.g. +90 doesn't filter, -18 filters to only events after 
twilight).  Times are local or UT depending on how the flag is set 
in the main window (but input ephemerides are always barycentric
julian dates).  The output is in the GEOCENTRIC time system, i.e. 
the time the signal arrives at earth, The period and JD errors 
are propagated to give phase and/or time uncertainties as appropriate. 
Don't neglect the uncertainties -- when they amount to an 
appreciable fraction of a cycle, the ephemeris is basically worthless.

[Note - barycentric and heliocentric times are the same within
a few seconds, because the solar system barycenter stays 
mostly within the body of the sun, which is about 2.3 light 
seconds in radius].
"""

                self.outwin.settext(helptext)

def showhelptext(parentframe = None) :
        helpwin = ScrolledText(parentframe)
        
root = Tk()

#if platform.system() == "Darwin" :
#    os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python3" to true' ''')

# sitedict is read in thorskyclasses, and is global there.
# This apparently works.

default_site = 'mdm'

obs = Observation(default_site = 'mdm')  # default site, right now, in zenith.

obs.computesky()
obs.computesunmoon()
obs.computeplanets()
obs.setnightevents()
obs.computebary()

if TINYPIX :
# scale fonts for Latitude E7450 screen ...
   default_font = nametofont("TkDefaultFont")
   default_font.configure(size=11)
   # put this in and font change propagates
   root.option_add("*Font", default_font)

# main window GUI

root.title("Sky Calculator by John Thorstensen, Dartmouth College")
toppan = Frame(root)
toppan.pack(side = TOP)
bottompan = Frame(root)
bottompan.pack(side = BOTTOM)
leftpan = Frame(toppan)
leftpan.pack(side = LEFT, expand = YES, fill=X, anchor=N)
(invars,StepBox,ObjBox) = makeinform(leftpan,infields)

rightpan = Frame(toppan)
rightpan.pack(side = LEFT, expand = YES, fill = X)
(outvars,outboxes) = makeoutform(rightpan,outfields)

# sky and airmass display windows

skdisp = SkyDisplay(obs,objs2000)
 
airmdisp = AirmDisplay(obs)
 
#hawin = Toplevel()
#hawin.title("Hourly Airmass")
#(havars,haboxes) = makehaform(hawin) 
# hawin.withdraw()

# other windows withdrawn after making

almwin = Toplevel()
almwin.title("Nightly Almanac")
almvars = makealmform(almwin,almfields)
if not is_MacOS :
    almwin.withdraw()

# print "almwin is made ... "

coowin = Toplevel()
coowin.title("Coordinates ...")
[coovars, cooboxes] = makecooform(coowin,coofields)
if not is_MacOS :
    coowin.withdraw()

planetwin = Toplevel()
planetwin.title("Planets (J2000)")
(planetvars, planetboxes) = makeplanetform(planetwin)
warning_planet = fillplform(planetvars, planetboxes)

# For some reason if you immediately withdraw a window
# on MacOS you can't get it back.

if not is_MacOS :
    planetwin.withdraw()

# UNIMPLEMENTED PIECES MOMENTARILY TURNED OFF FOR DEVEL.

havars = None
#coovars = None
haboxes = None

sitewin = Toplevel()
sitewin.title("Site menu")

makesiteform(sitewin,invars['sitecode'])

if not is_MacOS :  # this doesn't work (yet) on MacOS
    imwin = ImageWindow()   # interface to DS9
    if not is_MacOS :  # if and when this does get working, we will not 
                       # want to withdraw on MacOS since then it won't come back.
        imwin.withdraw()   

# print "calling circumstances ... "
#circumstances(invars,outvars,havars,almvars,coovars, planetvars, 
#        outboxes, haboxes, planetboxes, skdisp, airmdisp)

### BUTTONS, BINDINGS, AND A FEW WINDOW UTILITIES.

topbuttonpan = Frame(bottompan)
# buttonpan.pack(side=BOTTOM)
topbuttonpan.pack(side=TOP)
Button(topbuttonpan, text = 'Refresh\n output',
   command = (lambda v=invars,o=outvars,h=havars,a=almvars,c=coovars,
      p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp, ad=airmdisp: 
        circumstances(v,o,h,a,c,p,b,hb,pb,sd,ad)),width=butwid(5)).pack(side=LEFT)
Button(topbuttonpan, text = 'Step\n time',
   command = (lambda v=invars,o=outvars,h=havars,a=almvars,c=coovars,
      p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp: 
        steptime(v,o,h,a,c,p,b,hb,pb,sd,ad,True)),width=butwid(3)).pack(side=LEFT)
Button(topbuttonpan, text = 'Set to\n now',
   command = (lambda v=invars,o=outvars,h=havars,a=almvars,c=coovars,
      p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp: 
        set_to_now(v,o,h,a,c,p,b,hb,pb,sd,ad)),width=butwid(4)).pack(side=LEFT)

Button(topbuttonpan, text = 'UT <->\n local',
   command = (lambda v=invars,o=outvars,h=havars,a=almvars,c=coovars,
     p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp: 
        convert_time(v,o,h,a,c,p,b,hb,pb,sd,ad)),width=butwid(4)).pack(side=LEFT,fill=Y)

Button(topbuttonpan, text = 'Planets',
   command = (lambda w = planetwin : win_raise(w)),width=butwid(4)).pack(side=LEFT,fill=Y)

Button(topbuttonpan, text = 'Site\n menu',
   command = (lambda w = sitewin : win_raise(w)),width=butwid(3)).pack(side=LEFT)

Button(topbuttonpan, text = 'Sky\nDisplay',
   command = (lambda sd = skdisp, o=obs, obj=objs2000 : skydisplay_toggle(sd,o,obj)),width=butwid(4)).pack(side=LEFT,fill=Y)
Button(topbuttonpan, text = 'Airmass\nDisplay',
   command = (lambda ad = airmdisp, o=obs : airmdisplay_toggle(ad,o)),width=butwid(4)).pack(side=LEFT,fill=Y)
Button(topbuttonpan, text = 'Nightly\n Almanac',
   command = (lambda w = almwin: win_raise(w)),width=butwid(5)).pack(side=LEFT)

bottombuttonpan = Frame(bottompan)
bottombuttonpan.pack(side = TOP)
Button(bottombuttonpan, text = 'Alt.\nCoords',
   command = (lambda w = coowin : win_raise(w)),width=butwid(4)).pack(side=LEFT)

Button(bottombuttonpan, text='Get object\nlist', 
   command=(lambda v =invars : 
    Objlistwin(v)),width=butwid(7)).pack(side=LEFT)
6
Button(bottombuttonpan, text = 'Text\nTables',
   command = (lambda v = invars : text_table_win(v)),width=butwid(4)).pack(side=LEFT)

if not is_MacOS :  # Don't put up this button on Macs 'cause it doesn't work yet.
    Button(bottombuttonpan, text = 'DSS\n window',
       command = (lambda i = imwin : win_raise(i)),width=butwid(4)).pack(side=LEFT,anchor=W,fill=Y)

quitbutton = Quitter(bottombuttonpan).pack(side=LEFT,fill=Y)

def win_raise(window) : 
        # takes a pre-existing, possibly iconified window and either
        # de-iconifies it or raises it (the .lift() method) if it's 
        # hidden behind other windows.  This is a nice behavior.
        if window.state() == 'normal' : 
                window.lift()
                # Having trouble with this on MacOS; digging in stackexchange ... 
                # window.attributes('-topmost',True)
                # window.after_idle(root.attributes,'-topmost',False)
        else : 
                geotmp = window.geometry()  # it's saved when window is withdrawn
                window.deiconify()          # takes no arguments ... 
                if geotmp != '1x1+0+0' :    # never-popped windows are zero size
                        window.geometry(geotmp)     # put it back where it was!

       # print "window.state is now ",window.state()

def skydisplay_toggle(skd, obs, objs) :
    if skd.isvisible :
       skd.isvisible = False
       skd.withdraw()
    else :
       skd.isvisible = True
       win_raise(skd)
       skd.redraw(obs, objs)

def airmdisplay_toggle(ad, obs) :
    if ad.isvisible :
       ad.isvisible = False
       ad.withdraw()
    else :
       ad.isvisible = True
       win_raise(ad)
       ad.redraw(obs)

root.bind('<Return>',(lambda event,v=invars,o=outvars,h=havars,a=almvars,
  c=coovars,p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp: 
        circumstances(v,o,h,a,c,p,b,hb,pb,sd,ad)))

outboxes['jd'].bind('<Return>',(lambda event,v=invars,o=outvars,h=havars,a=almvars,
  c=coovars,p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp: 
        force_jd(v,o,h,a,c,p,b,hb,pb,sd,ad)))

# left mouse - select nearest object
skdisp.bind('<Button-1>',(lambda event,v=invars,o=outvars,h=havars,a=almvars,
  c=coovars,p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp:
        set_to_moused_obj(event,v,o,h,a,c,p,b,hb,pb,sd,ad)))

# center mouse - load coords of mouse itself
skdisp.bind('<Button-2>',(lambda event,v=invars,o=outvars,h=havars,a=almvars,
  c=coovars,p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp:
        set_to_moused_coord(event,v,o,h,a,c,p,b,hb,pb,sd,ad)))

# h (for "HR") or s (for star) - load nearest plotted bright star
skdisp.bind('h',(lambda event,v=invars,o=outvars,h=havars,a=almvars,
  c=coovars,p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp:
        set_to_moused_bright(event,v,o,h,a,c,p,b,hb,pb,sd,ad)))

# right mouse - help text - appears while right mouse button held down
skdisp.bind('<Button-3>',(lambda event : skdisp.gethelp()))

# right mouse - help text disappears when right mouse button released
skdisp.bind('<ButtonRelease-3>',(lambda event : skdisp.killhelp()))

# h (for "HR") or s (for star) - load nearest plotted bright star
skdisp.bind('h',(lambda event,v=invars,o=outvars,h=havars,a=almvars,
  c=coovars,p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp:
        set_to_moused_bright(event,v,o,h,a,c,p,b,hb,pb,sd,ad)))

skdisp.bind('s',(lambda event,v=invars,o=outvars,h=havars,a=almvars,
  c=coovars,p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp:
        set_to_moused_bright(event,v,o,h,a,c,p,b,hb,pb,sd,ad=airmdisp)))

# f - step forward in time
skdisp.bind('f',(lambda event,v=invars,o=outvars,h=havars,a=almvars,
  c=coovars,p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp: 
        steptime(v,o,h,a,c,p,b,hb,pb,sd,ad,True)))

# n - set to now
skdisp.bind('n',(lambda event,v=invars,o=outvars,h=havars,a=almvars,
  c=coovars,p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp: 
        set_to_now(v,o,h,a,c,p,b,hb,pb,sd,ad)))

# b - step backward in time
skdisp.bind('b',(lambda event,v=invars,o=outvars,h=havars,a=almvars,
  c=coovars,p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp: 
        steptime(v,o,h,a,c,p,b,hb,pb,sd,ad,False)))

# q - quit display (hide, it doesn't get destroyed)
skdisp.bind('q',(lambda event,sd = skdisp, os = obs, ob = objs2000 : 
        skydisplay_toggle(sd, os, ob)))

# Zoom -- zoom to mouse position (x 1.8); bound to 'z' and 'i'
skdisp.bind('z',(lambda event, os = obs, ob = objs2000 :
        skdisp.zoom(event, 1.8, os, ob)))
skdisp.bind('i',(lambda event, os = obs, ob = objs2000 : 
        skdisp.zoom(event, 1.8, os, ob)))

# o - zoom out
skdisp.bind('o',(lambda event, os = obs, ob = objs2000 :
        skdisp.zoomout(event, 1.8, os, ob)))

# r or 1 - restore 1-to-1
skdisp.bind('r',(lambda event, os = obs, ob = objs2000 : 
        skdisp.zoomdef(event, os, ob)))
skdisp.bind('1',(lambda event, os = obs, ob = objs2000 : 
        skdisp.zoomdef(event, os, ob)))

# p or c - pan center
skdisp.bind('p',(lambda event, os = obs, ob = objs2000 : 
        skdisp.pan(event, os, ob)))
skdisp.bind('c',(lambda event, os = obs, ob = objs2000 : 
        skdisp.pan(event, os, ob)))


airmdisp.bind('q',(lambda event,ad = airmdisp, os = obs : 
        airmdisplay_toggle(ad, os)))

#airmdisp.bind('<Button-1>',(lambda event,v=invars,o=outvars,h=havars,a=almvars,
#  c=coovars,p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp:
#        set_to_moused_time(event,v,o,h,a,c,p,b,hb,pb,sd,ad)))

# enable reverse transformation from galactic coords by binding carriage
# returns in those boxes.

cooboxes['Galactic longitude'].bind('<Return>',(lambda event,v=invars,o=outvars,h=havars,a=almvars,
  c=coovars,p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp: 
        force_galact(v,o,h,a,c,p,b,hb,pb,sd,ad)))

cooboxes['Galactic latitude'].bind('<Return>',(lambda event,v=invars,o=outvars,h=havars,a=almvars,
  c=coovars,p=planetvars,b=outboxes,hb=haboxes,pb=planetboxes,sd=skdisp,ad=airmdisp: 
        force_galact(v,o,h,a,c,p,b,hb,pb,sd,ad)))

# Finally, call circumstances once to fill in the GUI, and go to mainloop.

# print "calling circumstances ... "
circumstances(invars,outvars,havars,almvars,coovars, planetvars, 
        outboxes, haboxes, planetboxes, skdisp, airmdisp)

root.mainloop()
