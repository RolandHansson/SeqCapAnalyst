#!/bin/bash

log_directory="benchmark_test"
mkdir "$log_directory"
cd $log_directory
logfile="benchmark_log"
echo -n "" > $logfile

function pipeline_test ()
{
echo "starting pipeline_test! cores = $1, target = $2"
cores=$1
target=$2 
name="${target}_${cores}"
echo $name > $logfile
pipeline.sh -n $cores -f "$name" -d 1 ${target}_R1_001.fastq
extract_log $name "${target}_${cores}_cores"
rm -r "$name"
rm -r Consensus
}

function extract_log () 
{
folder=$1
name=$2
log=`cat $folder/log.txt | grep "took"`
touch $name 
[[ -f $name ]] || exit
echo -e "$log"
echo -e "$log" > ${name}.log
echo -e "$name\n" "$log" >> $logfile
echo -e "$log" | grep -oP "[0-9]+h [0-9]+m [0-9]+s" | awk '{ print (3600*$1+60*$2+$3) }' > ${name}.log.seconds
}

############################################3
####     MAIN      ####
#######################

pipeline_test 18 IMIN44_S12_L001
pipeline_test 18 BL37577_S10_L001
pipeline_test 18 BL37590_S4_L001
pipeline_test 24 BL37590_S4_L001
pipeline_test 24 IMIN44_S12_L001
pipeline_test 24 BL37577_S10_L001
pipeline_test 30 IMIN44_S12_L001
pipeline_test 36 IMIN44_S12_L001



pipeline_test 6 IMIN44_S12_L001
pipeline_test 12 IMIN44_S12_L001

pipeline_test 6 BL37577_S10_L001
pipeline_test 12 BL37577_S10_L001

pipeline_test 6 BL37590_S4_L001
pipeline_test 12 BL37590_S4_L001

pipeline_test 3 IMIN44_S12_L001 
pipeline_test 2 IMIN44_S12_L001 
pipeline_test 1 IMIN44_S12_L001 
pipeline_test 4 IMIN44_S12_L001
pipeline_test 5 IMIN44_S12_L001

pipeline_test 1 BL37577_S10_L001 
pipeline_test 3 BL37577_S10_L001

pipeline_test 1 BL37590_S4_L001 
pipeline_test 3 BL37590_S4_L001
