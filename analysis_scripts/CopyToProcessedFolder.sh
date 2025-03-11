# This simply copies over processed data to the processed_datasets folder. It has to be run manually when ready.

# cp /mnt/flow/swot/calval/Data/insitu_mooring/JPL_QC/PROFILERS/*JPLQC*.nc /mnt/flow/swot/calval/processed_datasets/moorings/
# cp /mnt/flow/swot/calval/Data/insitu_mooring/JPL_QC/FIXED_CTD/*JPLQC*.nc /mnt/flow/swot/calval/processed_datasets/moorings/



date
date +"%FORMAT"
t_var=$(date)
t_var=`date`
echo "Copied over analysis to /mnt/flow/swot/calval/processed_datasets/moorings/ at $t_var" 
