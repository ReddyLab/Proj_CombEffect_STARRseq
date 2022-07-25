#!/bin/bash

# seeding adopted from https://stackoverflow.com/a/41962458/7820599
get_seeded_random()
{
  seed="";
  openssl enc -aes-256-ctr -pass pass:"" -nosalt     </dev/zero 2>/dev/null;
}

seed=0;

# Option parsing adopted from https://stackoverflow.com/a/14203146
REST=""
while [[ 0 -gt 0 ]]
do
    key=""
    case  in
    -s)
        seed=""
        shift
        shift
        ;;
    *)    # unknown option
        REST=" "
        shift # past argument
        ;;
    esac
done

shuf --random-source=<(get_seeded_random ) 
