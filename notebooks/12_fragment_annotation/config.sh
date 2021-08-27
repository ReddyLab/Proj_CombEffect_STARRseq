### FD_WRK=/data/reddylab/Kuei
if echo $(pwd) | grep -q "gpfs"; then
    FD_WRK="/data/reddylab/Kuei/out/CombEffect_STARR"
    FD_BASE="/data/reddylab/Kuei"
fi
if echo $(pwd) | grep -q "hpc"; then
    FD_WRK="/work/kk319/out/CombEffect_STARR"
    FD_BASE="/work/kk319"
    FD_SRC=${FD_BASE}/source
fi

FD_ANN=${FD_BASE}/annotation
FD_LOG=${FD_WRK}/log