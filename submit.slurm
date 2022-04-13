#!/bin/bash
#SBATCH -J ba-prob
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -p gpu
#SBATCH -C V100
#SBATCH --gres=gpu:1
#SBATCH -q wildfire
module load pgi/19x
./ba-probability ba-BA-sys3-5_EEE.txt --output-file badf-BA-sys3-5_EEE.txt --num_theta 250 --num_length 150 --d_range 1.5 3.5 --theta_range 0.0 180.0 --d_width_factor 3.0 --q_width_factor 5.0
