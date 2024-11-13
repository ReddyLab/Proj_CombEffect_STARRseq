#!/bin/bash

source config_rsync.sh
ssh ${SSH_DCC} ls -l ${PREFIX_DCC}/${FDIRY_DCC}/${1}
echo