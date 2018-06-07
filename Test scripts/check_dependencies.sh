#!/bin/bash

clear=1

function check_command
{
command_word=$1
package=${2:-$1} 
echo -n "Checking $package... "
command -v "$1" >/dev/null 2>&1 && { echo "OK"; } || { echo >&2 "NOT FOUND"; clear=0; }
}

check_command samtools
check_command bedtools 
check_command vcftools
check_command vcfutils.pl vcfutils
check_command fastqc FastQC 
check_command bowtie2 
check_command picard.jar picard
check_command seqtk 
check_command trimmomatic-0.36.jar trimmomatic
check_command mafft
check_command something_invalid

[ $clear -eq 1 ] && { echo "Good to go!"; } || { echo "Missing programs!"; }



