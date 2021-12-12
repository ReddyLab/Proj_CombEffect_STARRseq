### Check which duke server I am at and load the corresponding modules
if echo $(pwd -P) | grep -q "gpfs"; then
    #echo "You are on Duke Server: HARDAC"
    module load Anaconda3/2019.10-gcb02
    conda activate /data/common/shared_conda_envs/ucsc
fi

if echo $(pwd -P) | grep -q "hpc"; then
    #echo "You are on Duke Server: DCC"
    module load BigWig
fi
