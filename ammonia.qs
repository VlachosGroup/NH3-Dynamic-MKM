#!/bin/bash -l
#
# Sections of this script that can/should be edited are delimited by a
# [EDIT] tag.  All Slurm job options are denoted by a line that starts
# with "#SBATCH " followed by flags that would otherwise be passed on
# the command line.  Slurm job options can easily be disabled in a
# script by inserting a space in the prefix, e.g. "# SLURM " and
# reenabled by deleting that space.
#
# This is a batch job template for a program using multiple processor
# cores/threads on a single node.  This includes programs with OpenMP
# parallelism or explicit threading via the pthreads library.
#
# Do not alter the --nodes/--ntasks options!
#SBATCH --nodes=1
#SBATCH --ntasks=1
#
# [EDIT] Indicate the number of processor cores/threads to be used
#        by the job:
#
#SBATCH --cpus-per-task=1
#
# [EDIT] All jobs have memory limits imposed.  The default is 1 GB per
#        CPU allocated to the job.  The default can be overridden either
#        with a per-node value (--mem) or a per-CPU value (--mem-per-cpu)
#        with unitless values in MB and the suffixes K|M|G|T denoting
#        kibi, mebi, gibi, and tebibyte units.  Delete the space between
#        the "#" and the word SBATCH to enable one of them:
#
#SBATCH --mem=64G
# SBATCH --mem-per-cpu=1024M
#
# [EDIT] Each node in the cluster has local scratch disk of some sort
#        that is always mounted as /tmp.  Per-job and per-step temporary
#        directories are automatically created and destroyed by the
#        auto_tmpdir plugin in the /tmp filesystem.  To ensure a minimum
#        amount of free space on /tmp when your job is scheduled, the
#        --tmp option can be used; it has the same behavior unit-wise as
#        --mem and --mem-per-cpu.  Delete the space between the "#" and the
#        word SBATCH to enable:
#
# SBATCH --tmp=24G
#
# [EDIT] It can be helpful to provide a descriptive (terse) name for
#        the job (be sure to use quotes if there's whitespace in the
#        name):
#
#SBATCH --job-name=amm_01
#
# [EDIT] The partition determines which nodes can be used and with what
#        maximum runtime limits, etc.  Partition limits can be displayed
#        with the "sinfo --summarize" command.
#
# SBATCH --partition=standard
#
#        To run with priority-access to resources owned by your workgroup,
#        use the "_workgroup_" partition:
#
#SBATCH --partition=ccei_biomass
#
# [EDIT] The maximum runtime for the job; a single integer is interpreted
#        as a number of minutes, otherwise use the format
#
#          d-hh:mm:ss
#
#        Jobs default to the default runtime limit of the chosen partition
#        if this option is omitted.
#
#SBATCH --time=2-00:00:00
#
#        You can also provide a minimum acceptable runtime so the scheduler
#        may be able to run your job sooner.  If you do not provide a
#        value, it will be set to match the maximum runtime limit (discussed
#        above).
#
# SBATCH --time-min=0-00:02:00
#
# [EDIT] By default SLURM sends the job's stdout to the file "slurm-<jobid>.out"
#        and the job's stderr to the file "slurm-<jobid>.err" in the working
#        directory.  Override by deleting the space between the "#" and the
#        word SBATCH on the following lines; see the man page for sbatch for
#        special tokens that can be used in the filenames:
#
# SBATCH --output=%x-%j.out
# SBATCH --error=%x-%j.out
#
# [EDIT] Slurm can send emails to you when a job transitions through various
#        states: NONE, BEGIN, END, FAIL, REQUEUE, ALL, TIME_LIMIT,
#        TIME_LIMIT_50, TIME_LIMIT_80, TIME_LIMIT_90, ARRAY_TASKS.  One or more
#        of these flags (separated by commas) are permissible for the
#        --mail-type flag.  You MUST set your mail address using --mail-user
#        for messages to get off the cluster.
#
# SBATCH --mail-user='wittregr@udel.edu'
# SBATCH --mail-type=END,FAIL,TIME_LIMIT_90
#
# [EDIT] By default we DO NOT want to send the job submission environment
#        to the compute node when the job runs.
#
#SBATCH --export=NONE
#

#
# [EDIT] If you're not interested in how the job environment gets setup,
#        uncomment the following.
#
#UD_QUIET_JOB_SETUP=YES

#
# Do standard OpenMP environment setup:
#
. /opt/shared/slurm/templates/libexec/openmp.sh

#
# [EDIT] Add chemkin to the environment:
#
vpkg_require matlab/r2020b

srun matlab -nodisplay -nosplash -nodesktop -singleCompThread -r 'amm_main4(593)'
