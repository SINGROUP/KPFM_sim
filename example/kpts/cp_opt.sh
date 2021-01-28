#!/bin/bash
# after tunning run_opt.py

echo "glob res are not enough worker is copied for planning; also prepared for 2 results"

cp glob_res/opt.db glob_res/afm.db
cp glob_res/opt.db afm_workers.db
cp glob_res/opt.db afm_res.db
cp glob_res/opt.db afm_res2.db
## this is for 2nd worker 

echo "done!"

