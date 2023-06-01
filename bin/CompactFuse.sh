#!/bin/bash


#Basic usage:
# - Select tools (default all).
# - tools a, j, s, f
# - Select analisis (default complete)
# - analisis complete, prioritization, calling
# If only prioritization is selected: A .txt file has to be generated containing the names of all samples.
# - Select File
# - reads Abs/File (File Contains a subfile for each fasta File/S1/S1_1 File/S1/S1_2)
# - 
#

Main_path=$PWD
Base_name=$(basename $1)

mkdir $Main_path/../Call_Results/$Base_name
mkdir $Main_path/../Call_Results/$Base_name/Arriba
mkdir $Main_path/../Call_Results/$Base_name/StarFusion
mkdir $Main_path/../Call_Results/$Base_name/JAFFA
mkdir $Main_path/../Call_Results/$Base_name/FusionCatcher



for dir in "$1"/*
do
  if [ -d "$dir" ]; 
    then
    parent_dir=$(basename "$dir")
    echo "$parent_dir" >> "$1/SRA_list.txt"
    ls "$dir" | sed 's/_1.fastq.gz//g' | sed 's/_2.fastq.gz//g' | sort | uniq >> "$1/SRA_list.txt"
    fi
done

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

# #ARRIBA
# mv $1 ../arriba_v2.4.0n
# mv SCRIPTS/arriba_fullrun.py ../arriba_v2.4.0n
# cd ../arriba_v2.4.0n
# python3 arriba_fullrun.py $1
# mv arriba_fullrun.py  ../PIPELINE/SCRIPTS
# mv $1 ../PIPELINE/$1
# cd ../PIPELINE

# #STARFUSION
# mv $1 ../STAR-Fusion-v1.12.0
# mv SCRIPTS/STARfusion_fullrun.py ../STAR-Fusion-v1.12.0
# cd ../STAR-Fusion-v1.12.0
# python3 STARfusion_fullrun.py $1
# mv STARfusion_fullrun.py ../PIPELINE/SCRIPTS
# mv $1 ../PIPELINE/$1
# cd ../PIPELINE

# #JAFFA
# mv $1 ../JAFFA  
# mv SCRIPTS/jaffa_fullrun.py ../JAFFA  
# cd ../JAFFA  
# python3 jaffa_fullrun.py $1
# mv jaffa_fullrun.py ../PIPELINE/SCRIPTS
# mv $1 ../PIPELINE/$1
# cd ../PIPELINE

# #FUSION CATCHER
# source activate pit2
# mv $1 ../fusioncatcher
# cd ../fusioncatcher
# WK=$PWD
# cd $1
# lsa=$(echo */)
# cd ..
# for val in $lsa; do
# 	echo "Runing command for ${val}"
# 	python2 bin/fusioncatcher.py -d data/ -i $WK/$1/$val -o $val
# 	mv $val /media/sequentia/sdb1/visitor3/$1/FusionCatcher
# done
# mv $1 ../PIPELINE/$1
# cd ../PIPELINE



if $a ; then
    python3 arriba_fullrun.py $1
fi

if $s; then
    python3 STARfusion_fullrun.py $1
fi

if $f; then
    source activate pit2
    lsa=$(echo */)
    cd ..
    for val in $lsa; do
	    echo "Runing command for ${val}"
	    python2 bin/fusioncatcher.py -d data/ -i $WK/$1/$val -o $val
    done
fi

if $j; then
    python3 jaffa_fullrun.py $1
fi