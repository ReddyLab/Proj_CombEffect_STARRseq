#!/bin/bash
source config_rsync.sh

### set path
SRC=${PREFIX_HAD}/${FDIRY_HAD}
DES=${SSH_DCC}:${PREFIX_DCC}/${FDIRY_DCC}

### rsync
./rsync_run.sh ${SRC}/ ${DES}/