### Check which duke server I am at and load the corresponding modules
if echo $(pwd) | grep -q "gpfs"; then
    #echo "You are on Duke Server: HARDAC"
    module load bedtools2
fi

if echo $(pwd) | grep -q "hpc"; then
    #echo "You are on Duke Server: DCC"
    module load Bedtools
fi
