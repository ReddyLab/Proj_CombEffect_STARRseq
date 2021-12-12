### https://unix.stackexchange.com/questions/72661/show-sum-of-file-sizes-in-directory-listing
#dir () { ls -FaGlh "${@}" | awk '{ total += $4; print }; END { print total }'; }
dir() { 
    ls -lhaG --color=always "${@}" |\
    sed -re 's/^([^ ]* ){3}//' |\
    awk '{ total += $1; print }; END { print total }'
}

### https://unix.stackexchange.com/questions/27013/displaying-seconds-as-days-hours-mins-seconds
function displaytime {
  local T=$1
  local D=$((T/60/60/24))
  local H=$((T/60/60%24))
  local M=$((T/60%60))
  local S=$((T%60))
  (( $D > 0 )) && printf '%d days ' $D
  (( $H > 0 )) && printf '%d hours ' $H
  (( $M > 0 )) && printf '%d minutes ' $M
  (( $D > 0 || $H > 0 || $M > 0 )) && printf 'and '
  printf '%d seconds\n' $S
}