#!/bin/bash

### set path
source config_duke.sh

SRC=${FD_RES}
DES=${FD_BAC}/proj_combeffect.tar.gz

### show msg
echo "+++++++++++++++++++++++++++"
echo "Backup Info"
echo "    Directory: " ${SRC}
echo "    Archive:   " ${DES}
echo
echo "+++++++++++++++++++++++++++"
echo "Backup (TAR + GZIP)"

### rync run
while true; do
    read -p "Are you sure to perform backup? " yn
    case $yn in
        [Yy]* ) echo -e "Backup starts:\n";  tar -cvzf ${DES} -C ${SRC} .; echo; break;;
        [Nn]* ) echo -e "Backup is not performed.\n"; exit;;
        * )     echo "Please answer yes or no.";;
    esac
done

### if you want to decompress the tar.gz file; run 
### mkdir ${SRC}
### tar -xvzf ${DES} -C ${SRC}