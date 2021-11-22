#!/usr/bin/python3 -u

import os
import numpy as np
import sys
import glob
import shutil
from ase import atoms
from ase.io import *
import math

from   ctypes import c_int, c_double, c_char_p
import ctypes

lp = os.path.dirname(os.path.realpath(__file__))
print ("DEBUG: lib_path", lp)


# ***** System information: ***** #
# path (absolute or relative) to  #
# your PP-AFM code !!! Ommited    #
# ppm_path = './PPM_master/'      #
# - uncomment previous and        #
# - following if you want to save #
# - as xsf ....                   #
# ******************************* #
#sys.path.append(ppm_path)
#import pyProbeParticle                as PPU     
#import pyProbeParticle.GridUtils      as GU

# ***************** Description ******************************** #
#                                                                #
# This script will take a single geometry from a xxx.traj file   #
# and creates an (aproximate) potential between the metallic tip #
# and sample in the KPFM measurements.                           #
# Those here in-built functions can be part of outer scripts     #
#                                                                #
# ****************** constants ********************************* #
#                                                                #
R_BOHR = 0.529772 # Bohr-radius in angstroms                     #
eV_to_au = 0.03674932540261308                                   #
#                                                                #
lib_ext   ='_lib.so'                                             #
#                                                                #
# ***************** PARAMETERS ********************************* #
runable      = False ;# use false all over here, unless debugging#
debug        = False ;
save_pre_opt = False ; 
save_npy     = True  ;
# ----------------- Main parameters ---------------------------- #

g_file = 'geometry_x0.0_y0.0.traj'; # from where the geometry is taken - in the future this is supposed to be from the data base
idx = 0;

V_tip = 1.00 # at the moment this is tip-Voltage #

# !@@@ division of sample and tip based on the height - carefull here !!!!  #
#sample_atoms = 84
#sample_first = True
# the heigh is chosen since the atoms are cc on the edges of the cell, because of the boundary conditions #
cc = 11.0 # the maximum "height" of the bottom (sample) atoms ; because of xzy - height refers to y in xyz file

# defining the grid: --- see bellow --- #
# inx, ... - amount of voxels in each direction #
# ddd - differencing of the small cell ndx, ndy, ndz -- 3x3 matrix #
inx = 320 ; iny = 625 ; inz = 288 ;
ddd = np.array([[8.9118733969396402E-002,0.0,0.0],
                [0.0,9.0706854378510865E-002,0.0],
                [0.0,0.0,8.5754541751130356E-002]])
ddd *= R_BOHR # ddd is in BOHRs, while goemetry and this script anyway works in Angstroms @
# ideal is to run a first try run in CP2K and then get the necessary fields for this - see bellow:
'''
Restart from density | ERROR! | CUBE FILE NOT COINCIDENT WITH INTERNAL GRID            1
Restart from density |          151  DIFFERS FROM          320
Restart from density |   0.18864885321402539        0.0000000000000000        0.0000000000000000       DIFFERS FROM    8.9118733969396402E-002
Restart from density | ERROR! | CUBE FILE NOT COINCIDENT WITH INTERNAL GRID            2
Restart from density |          300  DIFFERS FROM          625
Restart from density |    0.0000000000000000       0.18876044789079074        0.0000000000000000       DIFFERS FROM    9.0706854378510865E-002
Restart from density | ERROR! | CUBE FILE NOT COINCIDENT WITH INTERNAL GRID            3
Restart from density |          131  DIFFERS FROM          288
Restart from density |    0.0000000000000000        0.0000000000000000       0.18831740147122483       DIFFERS FROM    8.5754541751130356E-002
'''
# furter technical parameters #
zer = 1E-7 # zero (precision) #
sh  = 0.0 # shift in the voxel centre ? #
r_a = 1.0 # atomic radius !!!  - size how big the same potential as on the atom centre #
n_add = 1; #5 # 20; # ammount of additional indices around the original cell #

#************* technical parameters ***************************#
precond = 0 ; # 1- true ; seems to be faster without it        #
inner_step = 1000; # amount of inner steps                     #
#*****************    END     *********************************#

cpp_name='ElField' # uncomment the following two if you want to always compile the C++ functions #
#cpp_utils.compile_lib( cpp_name  )
#cpp_utils.make( "all"  )
lib    = ctypes.CDLL( lp+"/" + cpp_name + lib_ext )    # load dynamic librady object using ctypes

# define used numpy array types for interfacing with C++ -----------------------------

array1i = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array1d = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')
array4d = np.ctypeslib.ndpointer(dtype=np.double, ndim=4, flags='CONTIGUOUS')

# python wrape to C++ functions ------------------------------------------------------

#void create_xyzGrid(double *xyzgrid, int * dims, double *vox_vec, double *g_null){
lib.create_xyzGrid.argtypes = [array4d,array1i,array2d,array1d]
lib.create_xyzGrid.restype = None
def create_xyz_Grid(xyzgrid,g_vec,g_0):
    # creates a grid (4D array - nx * ny * nz * 3 ) of the xyz positions for every voxel, of the cube that will be calculated #
    # xyzgrid - zeros in & xyz out ; g_vec - small voxel vectors ; g_0 - origin of the cube file #
    dims = np.array(xyzgrid.shape[0:3],dtype=np.int32)
    lib.create_xyzGrid(xyzgrid,dims,g_vec,g_0)

#    void fix_all(int * fixtop, int * fixbot, int *fixnum, double * xyzgrid, int * dims, int * nat, double * pos, int * top_at, int * top_fix_at, int * bot_at, int *  bot_fix_at)
lib.fix_all.argtypes = [array1i, array1i, array1i, array4d, array1i, array1i, array2d, array1i, array1i, array1i, array1i]
lib.fix_all.restype = None
def fix_all(xyzgrid,pos,top_atoms,fixed_top_atoms,bot_atoms,fixed_bot_atoms):
    # this function will create anything necessary for the preparation of the V calculations - it will create 3 arrays of indices: #
    # fixtop - indices which belongs to the tip (+0.5 Vtip) #
    # fixbot - indices which belongs to the sample (-0.5 Vtip) #
    # fixnum - indices for each /y/ (z) line for the better pre - opt which is done for the lines#
    # inputs are: xyzgrid - 4d grid of the positions of voxels; pos - n_at * 3 matrix of xyz positions of atoms; #
    # -||- : top_atoms - indices of atoms of the tip; fixed_top_atoms - indices of fixed atoms of the tip #
    # -||- : bot_atoms - indices of atoms of the sample; fixed_bot_atoms - indices of fixed atoms of the sample #
    dims = np.array(xyzgrid.shape[0:3],dtype=np.int32)
    fixtop = np.zeros((dims[0]*dims[1]*dims[2]),dtype=np.int32)
    fixbot = np.zeros((dims[0]*dims[1]*dims[2]),dtype=np.int32)
    fixnum = np.zeros((2+dims[0]*dims[2]*2),dtype=np.int32) # first 2 are lengths of fixtop and fixbot, others are z(y in reality) points for each x and y max(bot) and min(top) for the V ramp -- preconditioner #
    fixnum[3::2]=100000; # some insanely hight number for the initial comparison
    if debug:
         print("debug: fixtop.shape",fixtop.shape)
         print("debug: fixnum (in)",fixnum)
    nta  = len(top_atoms)
    ntfa = len(fixed_top_atoms)
    nba  = len(bot_atoms)
    nbfa = len(fixed_bot_atoms)
    atnum = np.array((nta,nba,ntfa,nbfa),dtype=np.int32)
    top_atoms=np.array(top_atoms,dtype=np.int32)
    fixed_top_atoms=np.array(fixed_top_atoms,dtype=np.int32)
    bot_atoms=np.array(bot_atoms,dtype=np.int32)
    fixed_bot_atoms=np.array(fixed_bot_atoms,dtype=np.int32)
    lib.fix_all(fixtop,fixbot,fixnum,xyzgrid,dims,atnum,pos,top_atoms,fixed_top_atoms,bot_atoms,fixed_bot_atoms)
    ifixt=fixnum[0];ifixb=fixnum[1];
    if debug:
        print("debug: ifixt, ifixb, yfixt, yfixb",fixnum)
        #print("debug: fixtop[ifixt-2..ifixt++]", fixtop[ifixt-2],fixtop[ifixt-1],fixtop[ifixt],fixtop[ifixt+1],fixtop[ifixt+2])
        #print("debug: fixbot[ifixt-2..ifixt++]", fixbot[ifixb-2],fixbot[ifixb-1],fixbot[ifixb],fixbot[ifixb+1],fixbot[ifixb+2])
        #np.savetxt("fix_top_indices.txt" ,fixindices, fmt='%i')
    fixtop = np.delete(fixtop,np.arange(ifixt,dims[0]*dims[1]*dims[2],dtype=np.int32)) # removing the not - necessary (null) parts of the array #
    fixbot = np.delete(fixbot,np.arange(ifixb,dims[0]*dims[1]*dims[2],dtype=np.int32)) # removing the not - necessary (null) parts of the array #
    if debug:
        print("debug: len(fixtop), len(fixbot)", len(fixtop), len(fixbot))
        print("debug: fixtop[-2..0]", fixtop[-2],fixtop[-1],fixtop[0])
        print("debug: fixbot[-2..0]", fixbot[-2],fixbot[-1],fixbot[0])
    return fixtop, fixbot, fixnum;

#    void prepare_V( int n_add,  int * dims_small, int * fixtop, int * fixbot, int *fixnum, int * optind,  double * V,  double V_tip){
lib.prepare_V.argtypes = [ c_int ,     array1i,      array1i,      array1i,       array1i,      array1i,     array3d ,     c_double ]
lib.prepare_V.restype = None
def prepare_V(fixtop,fixbot,fixnum, V_tip,n_add,dims):
    # this function creates a first guess for the potential (Voltage) - V - inside the metallic KPFM tip-sample system #
    # together with this it also creates array of indices - optind - which will be changed durign the optimization -- this speeds-up the calculations #
    # inputs: fixtop - indices of voxels bellonging to the tip; fixbot - indices of voxels bellonging to the sample; #
    # -||-: fixnum - top and bottom indices of last fixed voxel in each /y/ (z) line; n_add - number of added voxels in /x/ and /z/ direction #
    # -||-: dims - dimensions of the cube f voxels (nx, ny, nz) #
    dims = np.array(dims,dtype=np.int32)
    dims_large = dims.copy()  ; # print( "debug: dims_small" , dims       );
    dims_large[::2] += 2*n_add; # print( "debug: dims_large" , dims_large );
    V = np.zeros(dims_large)  ; # print( "debug: dims_large are not passed to C++" )
    optind=np.zeros(dims[0]*dims[1]*dims[2]+1,dtype=np.int32);
    if debug:
        print( "debug: dims_small" , dims       );
        print( "debug: dims_large" , dims_large );
        print( "debug: dims_large are not passed to C++" )
        print( "debug: n_add",n_add)
    lib.prepare_V( int(n_add), dims , fixtop,fixbot,fixnum, optind, V, V_tip)
    optind = np.delete(optind,np.arange(optind[0],dims[0]*dims[1]*dims[2]+1,dtype=np.int32));
    if debug:
        print("debug: len(optind), optind[0]", len(optind), optind[0])
    return V, optind;

#     void opt_V( int n_add, int * dims_small, int *optind, double * V, double * Vout, float prec, int precond, int inner_step){ // prec - precission
lib.opt_V.argtypes = [ c_int ,     array1i,      array1i,    array3d ,     array3d ,   c_double ,     c_int ,         c_int,  c_int]
lib.opt_V.restype = None
def opt_V(n_add, dims, optind, V, prec, precond,inner_step=1000, idx=0):
    # full C++ procedure for calculating the potential - via iterative solving of the Laplace equation in the optind voxels - indices of voxels, that are optimized #
    # n_add - amount of added voxels on /x/ and /z/ sides; dims (nx , ny, nz) - dimensions of the inner V field ; prec - maximal difference between (control) steps ... #
    # ... so the consistency is achieved (<1e-6 seems to be fine); precond - idea about faster preconditioner, no advantages in the control run ; #
    # inner_step ... amount of steps in between the control steps - 1000 seems to be fine, since the control step is much slower, than the normal step. #
    # idx - index -- rank of the cpu core, if mpi is used (otherwised 0) #
    dims = np.array(dims,dtype=np.int32);
    Vtmp = V.copy()
    if debug:
        print("debug: prec", prec)
    lib.opt_V( int(n_add), dims ,  optind, V, Vtmp, prec, precond, inner_step, idx)
    if debug:
        print("debug: Vout - V:", Vtmp-V)
    return V;

#    void print_cube( int n_at, int * mol_Z, double * mol_xyz, double * grid_origin, double * grid_vec , int * dims, double * data, char *fname){
#    // number of atoms ; Z - atoms n_at vector; x,y,z - atoms n_at x 3; 3x vector         ; 3x3 matrix       ;  3x vector ; nx*ny*nz - one dimesional data - already rescalled in python; every distances already in bohrs from python
lib.writeCube.argtypes = [ c_int ,     array1i,      array2d,    array1d ,     array2d ,   array1i ,     array3d ,         c_char_p]
lib.writeCube.restype = None
def save_cube_c(density, mol_xyz, grid_origin, grid_vec, file_path='pot.cube'):
    # faster function to save data through C++ - density - data in eV; mol_xyz: 0-x, 1-y, 2-z, 3-Z for all atoms n_at*4 matrix; grid_origin - origin of the cube file in Ang; #
    # grid_vec - the small differential vector in And; file_path - how to save#
    n_at  = len(mol_xyz)
    dims  = np.array(density.shape,dtype=np.int32)
    mol_Z = np.array(mol_xyz[:,3],dtype=np.int32).copy();
    g_v   = grid_vec.copy()
    g_o   = grid_origin.copy()
    xyz   = np.array(mol_xyz[:,:3]).copy();
    if debug:
        print ("mol_xyz:", mol_xyz)
        print ("g_o:", g_o);
        print ("g_v:", g_v);
        print ("xyz:", xyz);
    xyz /= R_BOHR # scalling everything back to BORHs - cube inner lenght scale #
    g_v /= R_BOHR
    g_o /= R_BOHR
    if debug:
        print ("g_o:", g_o);
        print ("g_v:", g_v);
        print ("xyz:", xyz);
        print ("mol_Z:", mol_Z)
    data = density.copy()
    data *= eV_to_au ; # scalling to atoic units #
    lib.writeCube(n_at, mol_Z, xyz, g_o, g_v , dims, data, file_path.encode());
    print("--- Cube file written ---")


# functions: ---------------------------------------------------------------------------

def write_E(x):
    out = ' '+"{:.5E}".format(x) if x >= 0 else "{:.5E}".format(x);
    return out;

def save_cube(density, mol_xyz, grid_origin, grid_vec, file_path='hartree.cube', cube_head=None):
    # original function from Niko Oinonen - density - data in eV; mol_xyz: 0-x, 1-y, 2-z, 3-Z for all atoms n_at*4 matrix; grid_origin - origin of the cube file in Ang; #
    # grid_vec - the small differential vector in And; file_path - how to save#
    N = len(mol_xyz)
    density = density.copy()
    density *= eV_to_au
    grid_shape = density.shape
    mol_xyz = mol_xyz.copy()
    with open(file_path, 'w') as f:
        f.write('Comment line\nComment line\n')
        if cube_head is None:
            f.write(f'{N:5d} {" ".join([str("{:11.6f}".format(o/R_BOHR)) for o in grid_origin])}\n')
            for i in range(3):
                f.write(f'{grid_shape[i]:5d} {" ".join([str("{:11.6f}".format(v)) for v in grid_vec[i]/R_BOHR])}\n')
        else:
            f.write(cube_head)
        for x in mol_xyz:
            x[:3] /= R_BOHR
            f.write(f'{int(x[-1]):5d} {0.0:11.6f} {x[0]:11.6f} {x[1]:11.6f} {x[2]:11.6f}\n')
        ai = 0;
        for i in range(grid_shape[0]):
            for j in range(grid_shape[1]):
                for k in range(grid_shape[2]):
                    f.write(f' {write_E(density[i, j, k])}')
                    ai += 1 ;
                    if (ai % 6) == 0:
                        f.write('\n')
                # f.write('\n') # !! This part is noy in CP2K !!#


def separate_top_bottom(z_pos,cc,fi): # z_pos - z(/y/ in xyz file) positions of all atoms ; cc -centre ; fi - fixed atoms indices #
    # function that creates 4 array - ti (top atom indices); bi (bottom atom indices); fit (indices of top fixed atoms); fib (indices of bottom fixed atoms). #
    # from the height of atoms (ff), position of the dividing centre (cc) and indices of the fixed atoms (fi). #
    nat=len(z_pos);nfat=len(fi)
    ti =[];
    bi =[];
    fit=[];
    fib=[];
    for i in range(nat):
        if z_pos[i] < cc:
            bi.append(i)
        else:
            ti.append(i)
    for i in range(nfat):
        if z_pos[fi[i]] < cc:
            fib.append(fi[i])
        else:
            fit.append(fi[i])
    return ti, bi, fit, fib;

def copy_arround_borders(pos,lvs,sd): # sd - safe distance
    # copy the atoms left, right, front and back at the borders, within the safe distance so the PBC are OK, for the fixing of voxels, pre-conditioner and everything #
    normlvs = np.array((lvs[0]/lvs[0, 0], [0,100.,0], lvs[2]/lvs[2, 2]),dtype=np.double)
    if debug:
        print("debug: normlvs",normlvs)
    #print("debug:", normlvs[0])
    pos2=[];
    for i in range(len(pos)):
        if np.linalg.norm(pos[i,0:3:2])<sd: # original point
            #print("debug: norm, pos", np.linalg.norm(pos[i,0:3:2]) ,pos[i,0:3:2])
            pos2.append(pos[i]+lvs[0]+lvs[2])
        if np.linalg.norm(pos[i,0:3:2]-lvs[0,0:3:2]-lvs[2,0:3:2])<sd: # the ending point of the cell
            #print("debug: norm, pos", np.linalg.norm(pos[i,0:3:2]-lvs[0,0:3:2]-lvs[2,0:3:2]),pos[i,0:3:2],-lvs[0,0:3:2]-lvs[2,0:3:2])
            pos2.append(pos[i]-lvs[0]-lvs[2])
        if np.linalg.norm(pos[i,0:3:2]-pos[i,0]*normlvs[0,0:3:2])<sd: # the front side of the cell
            #print("debug: norm, pos", np.linalg.norm(pos[i,0:3:2]-pos[i,0]*normlvs[0,0:3:2]),pos[i,0:3:2],-pos[i,0]*normlvs[0,0:3:2])
            pos2.append(pos[i]+lvs[2])
        if np.linalg.norm(pos[i,0:3:2]-pos[i,2]*normlvs[2,0:3:2])<sd: # the left side of the cell
            pos2.append(pos[i]+lvs[0])
        if np.linalg.norm(pos[i,0:3:2]-lvs[2,0:3:2]-(pos[i,0]-lvs[2,0])*normlvs[0,0:3:2])<sd: # the back side of the cell
            pos2.append(pos[i]-lvs[2])
        if np.linalg.norm(pos[i,0:3:2]-lvs[0,0:3:2]-(pos[i,2]-lvs[0,2])*normlvs[2,0:3:2])<sd: # the back side of the cell
            pos2.append(pos[i]-lvs[0])
    if debug:
        print("debug: pos",pos)
        print("debug: pos2",pos2)
    pos2=pos if pos2 == [] else np.concatenate((pos,np.array(pos2,dtype=np.double),pos))
    if debug:
        print("debug: len(pos)",len(pos))
        print("debug: len(pos2)",len(pos2))
    return pos2;


# ----  definition of functions here: -------
def create_biased_cube(geom,V_tip, final_pot_name='final_pot_opt_'+str(zer)+'.cube', cube_head=None, cc=cc, idx=0):
    '''
    the main function, that for given geometry in an ASE objec that is given in **geom** and prepare the electrostatic potential for the Tip-Sample system with tip voltage **Vtip** #
    the final cube file is written into the **final_pot_name** file (this way multiple functions can be runned in the same time); cube_head is important for the exact match with the cp2k code #
    **cc** is centre of the geometry - plane which decided what is __tip__ (above) and what is __sample__ (below);
    **idx*** -- rank of CPU, when MPI is used, otherwise 0 -> for finding out, from which core the output came from ;
    !!! beware CP2K code is super strict, only the python save_cube is working properly for writing the cube file as it seems. #
    '''
    print ("Going to create KPFM (metallic) tip-sample electrostatic field -- for given geometry and creates cube file:", final_pot_name);
    pos = geom.positions
    n_at = len(pos)
    mol_xyz = np.zeros((n_at,4))
    mol_xyz[:,:3]=pos
    mol_xyz[:,3] =geom.get_atomic_numbers()
    lvs = geom.get_cell()
    print("copying atoms around the cell, because of PBC")
    pos2= copy_arround_borders(pos,lvs,0.1)
    g_vec = ddd.copy()

    nx = inx; dx = lvs[0,:]/nx ;
    ny = iny; dy = lvs[1,:]/ny ;
    nz = inz; dz = lvs[2,:]/nz ;
    ndim=np.array([nx,ny,nz]) # now everything moved to x->y->z
    g_or = np.array([0.,0.,0.])
    for i in range(3):
        #g_vec[i] /= ndim[i] ; # already defined in the beginning; #
        print ("g_vec",i,":",g_vec[i])
    
    if debug:
        print ("debug: nx, ny, nz", nx, ny, nz)
        print ("debug: dx",dx)
        print ("debug: dy",dy)
        print ("debug: dz",dz)
    
    xyzarr = np.zeros((nx,ny,nz,3)); # now everything x->y->z # cube way
    print("creating xyz array in C++")
    create_xyz_Grid(xyzarr,g_vec,g_or)
    print("back in Python")
    if debug:
        print("xyzarr - only on z axis:")
        print(xyzarr[0,0,:,2])
        print("xyzarr - only on y axis:")
        print(xyzarr[0,:,0,1])
        print("xyzarr - only on x axis:")
        print(xyzarr[:,0,0,0])
    print("getting indices and separating tip and sample ; ")
    fi = geom.constraints[0].get_indices()
    if debug:
        print ("debug: fi",fi)
        #print ("debug: lvec",lvec)
        print ("debug: pos2[:,1]",pos2[:,1])
    ti, bi, fit, fib = separate_top_bottom(pos2[:,1],cc,fi)
    if debug:
        print ("debug: ti",ti)
        print ("debug: fit",fit)
        print ("debug: bi",bi)
        print ("debug: fib",fib)
    # ti - tip index ; bi - bottom (sample) index ; zmt - z min tip ; zmb - z max tip ;; Last 2 - from where the voltage is automatically fixed;
    print("fixing the indices in C++")
    fixtop, fixbot, fixnum=fix_all(xyzarr,pos2,ti,fit,bi,fib)
    if debug:
        print ("debug: fixtop",fixtop)
        print ("debug: fixbot",fixbot)
        tmp = fixtop - np.roll(fixtop,1) # just to check, if all the numbers in fixtop are raising higher #
        mask = tmp == 0
        print ("debug: mask",mask)
        print ("debug: fixtop[mask], fixtop[0]",fixtop[mask],fixtop[0])
        print ("debug: len(fixtop[mask])",len(fixtop[mask]))
        tmp = fixbot - np.roll(fixbot,1) # just to check, if all the numbers in fixbot are raising higher #
        mask = tmp == 0
        print ("debug: mask",mask)
        print ("debug: len(fixbot[mask])",len(fixbot[mask]))
        print("debug: pre - prepare_V (fixtop,fixbot,fixnum,V_tip,n_add,ndim)", fixtop, fixbot, fixnum, V_tip, n_add, ndim)
    # --- now going to C++ procedures ---- #
    print("preparing the field in C++")
    Varr, optind = prepare_V(fixtop,fixbot,fixnum,V_tip,n_add,ndim);
    if debug:
        Vtmp2 = Varr.copy()
    if save_pre_opt:
        print("saving pre-optimized potential small")
        save_cube(Varr[n_add:nx+n_add,:,n_add:nz+n_add],mol_xyz,g_or,g_vec,file_path='pot_test.cube', cube_head=cube_head) # y is no longer larger - not needed
    g_or2 = g_or - np.array([g_vec[0,0]*n_add+g_vec[2,0]*n_add,0.,g_vec[0,2]*n_add+g_vec[2,2]*n_add])
    if debug:
        print ("debug: g_or2", g_or2)
    if save_pre_opt:
        print("saving pre-optimized potential full")
        save_cube_c(Varr[:,:,:],mol_xyz,g_or2,g_vec,file_path='pot_full_test.cube');
    if debug:
        print ("debug: Vtmp2-Varr", Vtmp2-Varr)
        print ("debug: optind", optind)
    print("optimizing the field through iterations, all in C++")
    Varr = opt_V(n_add, ndim, optind, Varr, zer, precond,inner_step=inner_step, idx=idx)
    # !!!! still some weird behaviour on the edges .... !!!! #
    print ("fully optimized potential - going to saving;")
    if save_npy:
        np.save(final_pot_name+"npy",Varr[n_add:nx+n_add,:,n_add:nz+n_add])
    save_cube(Varr[n_add:nx+n_add,:,n_add:nz+n_add],mol_xyz,g_or,g_vec,file_path=final_pot_name, cube_head=cube_head) # y is no longer larger - not needed
    #GU.save_scal_field( 'V_out', Varr[n_add:nx+n_add,n_add:ny+n_add,n_add:nz+n_add], lvec, data_format="xsf" ,cube_head=cube_head)
    print ("Everything saved for given geometry and cube file:", final_pot_name);
    return ;

# --------- main run -------------------------

if runable:
    geom = read(g_file,index=idx );
    create_biased_cube(geom,V_tipfinal_pot_name=g_file+"_"+str(idx)+"_final_pot.cube",cc=cc);
    print ("--- done ---")

#
