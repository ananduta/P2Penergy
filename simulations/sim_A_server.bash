#/bin/sh
# Request one node with 8 free processor cores
#PBS -l nodes=1:ppn=8:typel
# Mail me when the job ends for any reason
#PBS -m ae
# This is my email address
#PBS -M w.ananduta@tudelft.nl


# Go to the directory where I entered the qsub command
cd $PBS_O_WORKDIR

# Activate the Matlab version I want
module load matlab/2020a

# Run my M file and don't even try to display graphics
matlab -nodisplay -r sim_A
