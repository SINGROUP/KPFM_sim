#!/usr/bin/python3 -u
# ***** System information: ***** #
#
ppm_path = './PPM_master/'      # path (absolute or relative) to your PP-AFM code #
#
# ******************************* #

import os
import numpy as np
import matplotlib as mpl;  mpl.use('Agg'); print("plot WITHOUT Xserver"); # this makes it run without Xserver (e.g. on supercomputer) # see http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
import matplotlib.pyplot as plt
import sys
import glob
import shutil

sys.path.append(ppm_path)

import pyProbeParticle                as PPU     
import pyProbeParticle.GridUtils      as GU
import pyProbeParticle.PPPlot         as PPPlot
from   pyProbeParticle            import basUtils
from   pyProbeParticle            import elements 
#import pyProbeParticle.core           as PPC
import pyProbeParticle.HighLevel      as PPH
import pyProbeParticle.cpp_utils      as cpp_utils

# ***************** Description ******************************** #
#                                                                #
# this code/script lies on the edge of KPFM_sim and Probe        #
# Particle model (PP-AFM) - it takes the forces calculated as    #
# dE/dz from (scan of) CP2K calculations and save it as xsf/npy  #
# file as the PP-AFM is used to work with it (dF) or plot        #
#                                                                #
# ***************** PARAMETERS ********************************* #

debug = False

safe  = False # the safe procudere is checking and working only if all the points are calculated #
maxz  = None # the amount of z points/ None -- means all. If some points does not exits, they are Inf instead #

Q = 0.00 # for writing a new directory and params.ini #
K = 1.00 # for writing a new directory and params.ini #
A = 1.00 # just for F-df in params.ini #
data_format = "xsf" # xsf or npy #

#*****************    END     *********************************#

minf = 1000000 ;# minimum in all scans ???

#ls = os.listdir()
ls = [f for f in glob.glob("force_*_dE.txt")]
lsn = [f.replace("_"," ").replace("x","").replace("y","") for f in ls]
xys = np.array([[float(f.split()[1]),float(f.split()[2])] for f in lsn])



if debug:
    print ("debug:")
    print (xys)
    print (xys.shape)
    print

mixy = np.array([np.min(xys[:,0]),np.min(xys[:,1])])
maxy = np.array([np.max(xys[:,0]),np.max(xys[:,1])])
lxy  = maxy - mixy # I guess here we expect the same size of /x/and /y/
xu   = np.sort(np.unique(xys[:,0]))
yu   = np.sort(np.unique(xys[:,1]))
dx   = xu[1]-xu[0]
dy   = yu[1]-yu[0]

if debug:
    print (mixy,maxy,lxy)
    print (dx, dy)
print ("All the file-names; read overall ",xys.shape[0],"files with these unique X and Y values")
print (xu)
print (yu)
if safe: 
    if not (len(xu)*len(yu) == xys.shape[0]) : raise AssertionError ("the amount of unique X * Y values is different from the total amount of files", xys.shape[0])
else:
    print ("no assertion check in this file")
print ("continue to read the files")
    
tmpfilename="force_x"+str(xu[0])+"_y"+str(yu[0])+"_dE.txt"
if debug:
    print("tmpfilename",tmpfilename)
f0 = np.genfromtxt(tmpfilename)
if maxz is not None:
    f0 = f0[:maxz] # cutting to first maxz points #
maxc = len(f0)
if debug:
    print (f0) # - z descending
    print ("f0[0]",f0[0])
    print ("f0[-1]",f0[-1])
miz = np.min(f0[:,0])
maz = np.max(f0[:,0])
lz = maz - miz
dz = f0[0,0]-f0[1,0] # just be sure, that you have always the same interval between the points in your input files #
if debug:
    print (miz,maz,lz)
print (dz)

lvec = np.array([ [mixy[0],mixy[1],miz],
                  [lxy[0],0.,0.],
                  [0.,lxy[1],0.],
                  [0.,0.,lz]            ] )
ndim = np.array([len(f0),len(yu),len(xu)] ) # [z, y, x]

print ("lattice vector prepare")
print ("lvec:",lvec)
print ("ndim (z,y,x)):",ndim)
print ("going to read the files themselves")

print ("*** In any case the size of the whole field is based on the SIZE and AMOUNT of points from the FIRST calculated point !!! ***")

xsfout = np.zeros(ndim)

if debug:
    print (xsfout.shape)
    
for iy, y in enumerate(yu):
    for ix, x in enumerate(xu):
        if debug:
            print ("x,y:",x, y)
            print ("ix,iy:",ix, iy)
        tmpfilename="force_x"+str(x)+"_y"+str(y)+"_dE.txt"
        if safe : # there always are the files
            tmp = np.genfromtxt(tmpfilename)
        else:
            try:
                tmp = np.genfromtxt(tmpfilename)
                minf = len(tmp) if len(tmp) < minf else minf # maximal amount of points in each 
                print (ix,iy,"file exist, len(tmp):", len(tmp))
            except:
                tmp = np.zeros((maxc,2))
                tmp+= np.inf # putting there infinity
                print (ix,iy,"no file")	
        if len(tmp) == maxc:
            tmps = tmp[:maxz] # tmps - file to be saved
        else: # maxz is none -- using length of the first file 
            if debug:
                print("debug: maxc, len(tmp)",maxc,len(tmp) )
            tmps = np.zeros((maxc,2))
            tmps+= np.inf # putting there infinity
            ltmp = len(tmp) if len(tmp) <= maxc else maxc
            tmps[:ltmp] = tmp[:ltmp]
        xsfout[:,iy,ix]=tmps[::-1,1] # 

print ("all the forces read")

if debug:
    print (xsfout)
    print ("debug: minf",minf)
    
dirname = "Q%1.2fK%1.2f" %(Q,K)
print("saving the dE forces as relaxed_scan for ", dirname)
if not os.path.exists( dirname ):
    os.makedirs( dirname )
GU.save_scal_field( dirname+'/OutFz', xsfout, lvec, data_format=data_format )


if os.path.isfile("params.ini"):
    print("params.ini already exist in the directory -- moved to parans.ini-bak ")
    shutil.move("params.ini","params.ini-bak")
print ("forces saved, writing params.ini")
f=open("params.ini",'w')
print("charge         %1.2f                   # effective charge of probe particle [e]  #"                            %(Q),                   file=f)
print("stiffness      %1.2f  %1.2f    20.00   # [N/m] harmonic spring potential (x,y,R) components #"                 %(K,K),                 file=f)
print("scanMin        %2.6f  %2.6f    %2.6f   # start of scanning (x,y,z) {for tip, so PP is by r0Probe(z) lower}  #" %(mixy[0],mixy[1],miz), file=f)
print("scanMax        %2.6f  %2.6f    %2.6f   # end   of scanning (x,y,z) {for tip, so PP is by r0Probe(z) lower}  #" %(maxy[0],maxy[1],maz), file=f)
print("scanStep       %2.6f  %2.6f    %2.6f   # division of the scanning grid 0.1 is standart  #"                     %(dx,dy,dz),            file=f)
print("Amplitude      %1.2f                   # [Angstrom] peak-to-peak oscilation amplitude for conversion Fz->df #" %(A),                   file=f)
print("imageInterpolation None                # gives just full squares instead of interpolated values #" ,                                   file=f)
f.close()


# from Mathematica xsf file sorting definition - outer loop z - middle loop y - inner loop x (all ascending) #)
#For[i = 1, i <= n[[3]], i++, 
#  For[j = 1, j <= n[[2]], j++, 
#   For[k = 1, k <= n[[1]], 
#    k++, {xsf10[[i, j, k]] = xsfin[[maxat1 + a + 10, 1]], 
#     a = a + 1}]]];
# --- end --- #

print ("end")