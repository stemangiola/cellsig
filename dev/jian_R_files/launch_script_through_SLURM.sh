#!/bin/bash
#SBATCH --job-name=singlecore_job 
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1G

$1