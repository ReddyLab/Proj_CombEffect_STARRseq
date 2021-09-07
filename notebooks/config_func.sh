### https://unix.stackexchange.com/questions/72661/show-sum-of-file-sizes-in-directory-listing
#dir () { ls -FaGlh "${@}" | awk '{ total += $4; print }; END { print total }'; }
dir() { 
    ls -lhaG --color=always "${@}" |\
    sed -re 's/^([^ ]* ){3}//' |\
    awk '{ total += $1; print }; END { print total }'
}