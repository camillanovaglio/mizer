#Slurm
#
#Rscript remoteCalibration.R $SSBcali_value $rankCali_value $Q_value $K_value $R_value $output_file_name
#
base_file_name='optim_SEA_kappa'

for SSBcali_value in {0..1}; do  
   for rankCali_value in {0..1}; do
	    for Q_value in {0..1}; do
   		    for K_value in {0..1}; do
   		    	for R_value in {0..1}; do

                    # build a unique output file name based on the options that will be used.

                    base_file_name="optim_SEA_kappa"

                    if [[ SSBcali_value -eq 1 ]]; then
                        base_file_name=${base_file_name}_SSB
                    fi  
                    if [[ rankCali_value -eq 1 ]]; then
                        base_file_name=${base_file_name}_Cali
                    fi  
                    if [[ Q_value -eq 1 ]]; then
                        base_file_name=${base_file_name}_Q
                    fi 
                    if [[ K_value -eq 1 ]]; then
                        base_file_name=${base_file_name}_K
                    fi 
                    if [[ R_value -eq 1 ]]; then
                        base_file_name=${base_file_name}_R
                    fi 

                    echo $base_file_name
                    base_file_name=${base_file_name}.RData

      				sbatch --export SSBcali_value=$SSBcali_value,rankCali_value=$rankCali_value,Q_value=$Q_value,K_value=$K_value,R_value=$R_value,output_file_name=$base_file_name calibration.q
	
			    done
		    done
	    done	
   done
done
