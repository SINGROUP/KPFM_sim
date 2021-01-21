#!/bin/bash
# after tunning run_opt.py

echo "glob res are not enough, so this is for 2 workers"

cp glob_res/opt.db glob_res/afm.db
cp glob_res/opt.db afm_workers.db
cp glob_res/opt.db afm_res.db
# this is for 2nd worker
cp glob_res/opt.db afm_res2.db

echo "done!"

