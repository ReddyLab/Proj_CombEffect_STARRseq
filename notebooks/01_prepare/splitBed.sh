#############################################
### Usage ./splitBed.sh input.bed outdir
#############################################

### get directory from command line
FP_IN=$1
FD_OT=$2

### loop through the input bed file and 
### grep each chromosome to a output file
for chr in `cut -f 1 $FP_IN | sort | uniq`;
do
    echo $chr
    grep -w $chr $FP_IN > $FD_OT/$chr.bed
done
