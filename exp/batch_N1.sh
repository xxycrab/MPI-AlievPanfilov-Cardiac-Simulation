#!/bin/bash

sbatch N1_comet-24.slurm
sbatch N1_comet-48.slurm
sbatch N1_comet-72.slurm
sbatch N1_comet-96.slurm

sbatch N1_comet-24-nocomm.slurm
sbatch N1_comet-48-nocomm.slurm
sbatch N1_comet-72-nocomm.slurm
sbatch N1_comet-96-nocomm.slurm

sbatch N1_comet-24-ref.slurm
sbatch N1_comet-48-ref.slurm
sbatch N1_comet-72-ref.slurm
sbatch N1_comet-96-ref.slurm
