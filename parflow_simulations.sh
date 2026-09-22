#!/bin/bash
#
# author: Sami Miaari
# e-mail: s.miaari@fz-juelich.de

#SBATCH --partition=dc-cpu
#SBATCH --nodes=1
#SBATCH --account=jjsc39
#SBATCH --time=20:00:00
#SBATCH --job-name="residtol10-5"
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=s.miaari@fz-juelich.de
#SBATCH --mail-type=ALL
#

export PARFLOW_DIR=/p/project1/cslts/miaari1/v2026_parflow/ParFlow_scripts/stages2026/parflow_install
source /p/project1/cslts/miaari1/v2026_parflow/ParFlow_scripts/stages2026/env_2026_cpus.ini

# Loop over test_case200 .. test_case399
for i in {0..999}; do
	cd /p/scratch/cslts/miaari1/tolerance_sensitivity/residualtol10-5/test_case${i} || continue
	tclsh infiltration.tcl
	tclsh exfiltration.tcl
done
