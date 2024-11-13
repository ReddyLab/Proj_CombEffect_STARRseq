#!/bin/bash
source config_rsync.sh

### set path
SRC=${SSH_DCC}:${PREFIX_DCC}/${FDIRY_DCC}
DES=${PREFIX_HAD}/${FDIRY_HAD}

### rsync
./rsync_run.sh ${SRC}/ ${DES}/