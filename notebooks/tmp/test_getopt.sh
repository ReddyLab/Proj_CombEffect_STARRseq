#!/bin/bash


### https://sookocheff.com/post/bash/parsing-bash-script-arguments-with-shopts/
while getopts ":t:" opt; do
  case ${opt} in
    t )
      target=$OPTARG
      echo "TARGET: ${target}"
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      ;;
  esac
done
shift $((OPTIND -1))