#!/bin/bash
module purge
module load cp2k/6.1-openmpi-scalapack-popt anaconda/2020-04-tf2
export OMP_NUM_THREADS=1
ulimit -s unlimited

#source /home/krejcio1/.bashrc
eval "$(conda shell.bash hook)"
conda activate KPFM1
progdir=/home/krejcio1/wrkdir/Programs/ana-20.04tf2_py3.7.7

export PYTHONPATH=$PYTHONPATH:$progdir/KPFM_sim3:$progdir/CP2k_mtools3:$progdir/DFT_gridIO3

echo "extractinc data from database"

# Usage: python extract_descend_tip_data.py result_db_file x_tip y_tip
python3 extract_descend_tip_data.py glob_res/afm_forces.db 0.0 0.0

# Usage: python extract_descend_tip_data.py result_db_file x_tip y_tip
python3 extract_descend_tip_data.py glob_res/afm_forces.db 0.0 2.86

echo "Done, done"

