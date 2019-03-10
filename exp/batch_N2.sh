#!/bin/bash

sbatch N2_comet-96.slurm
sbatch N2_comet-192.slurm
sbatch N2_comet-240.slurm
sbatch N2_comet-384.slurm
sbatch N2_comet-480.slurm

sbatch N2_comet-96-nocomm.slurm
sbatch N2_comet-192-nocomm.slurm
sbatch N2_comet-240-nocomm.slurm
sbatch N2_comet-384-nocomm.slurm
sbatch N2_comet-480-nocomm.slurm

sbatch N2_comet-96-ref.slurm
sbatch N2_comet-192-ref.slurm
sbatch N2_comet-240-ref.slurm
sbatch N2_comet-384-ref.slurm
sbatch N2_comet-480-ref.slurm
