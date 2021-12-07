#/bin/sh
# Request one node with 4 free processor cores
#PBS -l nodes=1:ppn=4:typel
# Mail me when the job ends for any reason
#PBS -m ae
# This is my email address
#PBS -M w.ananduta@tudelft.nl

# Go to the directory where I entered the qsub command
cd $PBS_O_WORKDIR

module load matlab/2020a

matlab -nodisplay -r sim_A
