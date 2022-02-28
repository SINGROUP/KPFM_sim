#!/bin/bash
# check modules.txt for updating following lines: # module purge
#module load gcc/9.3.0  openmpi/4.0.3 cp2k/7.1-omp 
export OMP_NUM_THREADS=1
ulimit -s unlimited
export ppath=./programs_KPFM # 
export PYTHONPATH=$PYTHONPATH:$ppath/CP2k_mtools3:$ppath/DFT_gridIO3:$ppath/KPFM_FEM3:$ppath/KPFM_sim3   #:$ppath/pycp2k
export PYTHONPATH=$PYTHONPATH:/u/88/krejcio1/unix/WORK/lib_cp2k_copy/
KPFM_GLOBAL_SCRIPTS=$ppath/KPFM_sim3/scripts
ORIG_DIR=`pwd`

## after adjust the modules above and PYTHONPATH, this script will allow to get energies and forces (and also atomic forces) for predefined voltages in Vf ##

if true ; then
  y=0.0  ## adjust /x/ and /y/ as well as voltages /Vf/ <- for chanrges those which are around 0 should be the most important (other left here for cross-checking) ## 
  x=0.0  ## look at wiki for more details ##
  Vf="1000.0"
  Vf+=" -0.31 -0.2 0.2 0.31 -100.31 -100.2 100.2 100.31"
  Vf+=" 0.1 0.05 0.01 -0.01 -0.05 -0.1 100.1 100.05 100.01 -100.01 -100.05 -100.1"
  Vf+=" 102.0 101.7 101.5 101.3 101.0 100.7 100.5 100.3 -100.3 -100.5 -100.7 -101.0 -101.3 -101.5 -101.7 -102.0"
  Vf+=" 2.0 1.7 1.5 1.3 1.0 0.7 0.5 0.3 -0.3 -0.5 -0.7 -1.0 -1.3 -1.5 -1.7 -2.0"
  for Vi in $Vf; do
    echo "x y V:" $x $y $Vi
    python3 $KPFM_GLOBAL_SCRIPTS/extract_descend_tip_data.py -i glob_res/afm.db -p $x $y -V $Vi -n -N # -n -- no atomic forces -N -- no geometry
  done
fi
