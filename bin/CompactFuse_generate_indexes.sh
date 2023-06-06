#!/bin/bash

##GENERATION OF INDEXES FOR THE CALLING TOOLS

#Tool installer.
# -tools a,f,s,j default all


# Default variable values
a=false
s=false
f=false
j=false

# Check if any flags are provided
if [[ $# -eq 0 ]]; then
  # No flags provided, set all variables to true
  a=true
  s=true
  f=true
  j=true
else
  # Flags are provided
  while [[ $# -gt 0 ]]; do
    flag="$1"
    case "$flag" in
      -tools)
        shift
        # Loop through each character in the tools flag argument
        for char in $(echo "$1" | fold -w 1); do
          case "$char" in
            a)
              a=true
              ;;
            s)
              s=true
              ;;
            f)
              f=true
              ;;
            j)
              j=true
              ;;
          esac
        done
        ;;
    esac
    shift
  done
fi

#ASK FOR INSTALLATIONS
echo -e "This script downloads directly the premade indexes for hg38 \n "
echo "If any other build is needed please refer to their respective documentations and build it manually and include them inside the index folder as is done in the automatic installation"
already=false

if $a ; then
   wget -O ../INDEXES/CTAT_genome/CTAT.tar.gz https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play.tar.gz
   already=true
fi

if $s; then
    if already = false; then
        wget -O ../INDEXES/CTAT_genome/CTAT.tar.gz https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play.tar.gz
    fi
fi

if $f; then
    cd ../FusionCatcher_INS/data
    ./download-human-db.sh
    cd ../../bin
fi

if $j; then
   wget -O ../INDEXES/JAFFA/JaffaIndex.tar.gz https://figshare.com/ndownloader/files/25410494
fi