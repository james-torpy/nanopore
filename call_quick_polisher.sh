#!/bin/bash

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/single_cell/nanopore"
script_dir="$project_dir/scripts"
in_dir="$project_dir/FM/subset"

# get number of jobs:
no_jobs=$(\ls -U $in_dir | wc -l)

# call canu job:
qsub -wd $in_dir -V -pe smp 4 \
-N FM_canu -S /bin/bash -t 1-$no_jobs -tc 45 \
-l mem_requested=32G,h_vmem=32G,tmp_requested=80G -o /dev/null \
-e /dev/null $script_dir/quick_polisher.sh