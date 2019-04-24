#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=2g
module load R
cd model
Rscript remoteCalibration_loop.R $SSBcali_value $rankCali_value $Q_value $K_value $R_value $output_file_name


