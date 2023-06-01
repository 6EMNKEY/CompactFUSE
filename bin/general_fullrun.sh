#!/bin/bash

#RUN FROM PIPELINE
mkdir /media/sequentia/sdb1/visitor3/$1
mkdir /media/sequentia/sdb1/visitor3/$1/Arriba
mkdir /media/sequentia/sdb1/visitor3/$1/StarFusion
mkdir /media/sequentia/sdb1/visitor3/$1/JAFFA
mkdir /media/sequentia/sdb1/visitor3/$1/FusionCatcher


#ARRIBA
mv $1 ../arriba_v2.4.0n
mv SCRIPTS/arriba_fullrun.py ../arriba_v2.4.0n
cd ../arriba_v2.4.0n
python3 arriba_fullrun.py $1
mv arriba_fullrun.py  ../PIPELINE/SCRIPTS
mv $1 ../PIPELINE/$1
cd ../PIPELINE

#STARFUSION
mv $1 ../STAR-Fusion-v1.12.0
mv SCRIPTS/STARfusion_fullrun.py ../STAR-Fusion-v1.12.0
cd ../STAR-Fusion-v1.12.0
python3 STARfusion_fullrun.py $1
mv STARfusion_fullrun.py ../PIPELINE/SCRIPTS
mv $1 ../PIPELINE/$1
cd ../PIPELINE

#JAFFA
mv $1 ../JAFFA  
mv SCRIPTS/jaffa_fullrun.py ../JAFFA  
cd ../JAFFA  
python3 jaffa_fullrun.py $1
mv jaffa_fullrun.py ../PIPELINE/SCRIPTS
mv $1 ../PIPELINE/$1
cd ../PIPELINE

#FUSION CATCHER
source activate pit2
mv $1 ../fusioncatcher
cd ../fusioncatcher
WK=$PWD
cd $1
lsa=$(echo */)
cd ..
for val in $lsa; do
	echo "Runing command for ${val}"
	python2 bin/fusioncatcher.py -d data/ -i $WK/$1/$val -o $val
	mv $val /media/sequentia/sdb1/visitor3/$1/FusionCatcher
done
mv $1 ../PIPELINE/$1
cd ../PIPELINE

