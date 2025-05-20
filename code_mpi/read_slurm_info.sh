#!/bin/bash



./cond_mass_funct_1024_256_mpi $((4 * $SLURM_ARRAY_TASK_ID + $SLURM_PROCID + 1))



