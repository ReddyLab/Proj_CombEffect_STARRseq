### FD_WRK=/data/reddylab/Kuei
if echo $(pwd) | grep -q "gpfs"; then
    FD_WRK="/data/reddylab/Kuei/out/CombEffect_STARR"
    FD_BASE="/data/reddylab/Kuei"
fi
if echo $(pwd) | grep -q "hpc"; then
    FD_WRK="/work/kk319/out/CombEffect_STARR"
    FD_BASE="/work/kk319"
fi

#FD_SRC=${FD_WRK}/source
#FD_OUT=${FD_WRK}/out/CombEffect_STARR
#FD_DAT=/data/reddylab/gjohnson/whole_genome_STARRseq/wgss3/alignment_and_processing/alignments