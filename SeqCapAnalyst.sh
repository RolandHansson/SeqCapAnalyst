#!/bin/bash
MINARGS=1

# Use Trimmomatic, bowtie2 on a pair of files, 
# collect statistics on file quality with FastQC and FastQ_Screen,
# do variant calling, make consensus sequence (genome AND exons), 
# perform alignment of exons consensus from files run, and
# collect data about probe and exon coverage for statistical analysis.

### REQUIREMENTS: 
# Requires input fastq files to end in "_R[12]_001.fastq"
# Requires compatible files in /basefiles folder: 
#   Genome reference fasta file 
#   exon GFF file (with gene ID) 
#   probe file (csv or txt) 

# Picard (Calls on "java -jar bin/picard.jar")
# Bowtie2 (Calls on "bowtie2", "bowtie2-build")
# Trimmomatic (calls "java -jar ~/bin/trimmomatic-0.36.jar") <-- Will need fix
# SAMtools (calls on "samtools1.7")
# bcftools (calls on "bcftools1.6") 
# vcftools
# bedtools
# vcfutils.pl (usually included in above)
# seqtk 
# Fastqc (calls on "fastqc")
# cigar_parser.py (included) 

#########################################################


##################
### CONSTANTS: ###
##################

### Calculate default number of cores
CORES=`grep -c ^processor /proc/cpuinfo` #Count number of cores.
CORES=$(($CORES - 1)) #Save one core for other tasks.
[ $CORES -gt 10 ] && CORES=10 #Tests show neglible impact on performance when using more than 6-10 cores (depending on file size). Thus, default value is capped at 10. 
[ $CORES -lt 1 ] && CORES=1 #In case pipeline is run on PC with one processor. 

SKIP=1 #Skip steps if resulting files already exist
INTERM=1 #Use existing indices and intermediate file rather than create new. 
KEEP=0 #Keep temporary files
DELETE=0 #Delete everything but statistics afterwards
CALLING=0 #Determine method of variant calling (GATK/VCF). 
ARGUMENTS=$@
PIPELINE="$HOME/bin/SeqCapAnalyst" #folder of pipeline-related files
INTERMEDIATES="$PIPELINE/intermediates"
BASE_FILES="$PIPELINE/base_files"
FOLDER=""
FILENAME=""

#Shortcuts to files, in case you want to run with a specific version. 
#samtools="$PIPELINE/third_party_programs/samtools-1.7/samtools" 
samtools="samtools"
#bcftools="$PIPELINE/third_party_programs/bcftools-1.6/bcftools"
bcftools="bcftools"
trimmomatic="$PIPELINE/third_party_programs/Trimmomatic-0.36/trimmomatic-0.36.jar"
#picard="$PIPELINE/third_party_programs/picard/build/libs/picard.jar"
picard="picard"

# References to files with variable names
ADAPTER="$PIPELINE/adapters/TruSeq3-PE-2.fa"
if [[ ! `ls $BASE_FILES/*.[Ff][Aa]* | grep -vE "*.[Ff][Aa][Ii]" | wc -l` -eq 1 ]]
then
    echo -e >&2 "Error in base_files folder, cannot determine fasta file"
    exit
else
    #Only one fasta file in base_files folder (ignore index)
    GENOME_FILE=`ls $BASE_FILES/*.[Ff][Aa]* | grep -vE "*.[Ff][Aa][Ii]"` 
fi

if [[ ! `ls $BASE_FILES/*.[Gg][Ff][Ff] | wc -l` -eq 1 ]]
then
    echo -e >&2 "Error in base_files folder, cannot determine gff file"
    exit
else
    #Only one GFF file in base_files folder
    GFF_FILE=`ls $BASE_FILES/*.[Gg][Ff][Ff]`
fi

    #Only one txt/csv file in base_files folder
if [[ `ls $BASE_FILES/*.[Tt][Xx][Tt] | wc -l` -eq 1 ]]
then
    PROBE_FILE=`ls $BASE_FILES/*.[Tt][Xx][Tt]`
elif [[ `ls $BASE_FILES/*.[Cc][Ss][Vv] | wc -l` -eq 1 ]]
then
    PROBE_FILE=`ls $BASE_FILES/*.[Cc][Ss][Vv]`
else
  echo -e >&2 "Error in base_files folder, cannot determine probe file (must be txt/csv)"
  PROBE_FILE="undefined"
    exit
fi

FILTER_OUT="removed_reads_temp"
FASTQ_SCREEN_CONF="$PIPELINE/fastq_screen.conf"
FASTQ_SCREEN_CONF_BEGIN="$PIPELINE/fastq_screen_begin.conf"

SCRIPT_START_TIME=`date +%s`
LAST_ITERATION=0
STRINGENCY=4
TEST=0
HELP="
INPUT: 
FASTQ files you wish to decontaminate

OPTIONS:
-h help file 
-f [name]: specify folder name (Note: does not work in batch mode.)
-n [int]: specify number of cores to use. (Will auto-set to [cores available] -1, or '10' as upper cap.
-a [file: path] specify an adapter other than default. 
-b [int] (default 1): specify stringency of bowtie2 matching, from 1(permissive) to 5(draconic)
-s 0 (default 1) : set pipeline to overwrite any already-existing files in target folder, rather than assume they are correct. 
-i 0 (default 1) : set pipeline to overwrite any already-existing files in intermediates folder, rather than assume they are correct. 
"


#################
### FUNCTIONS ###
#################

function move_on_to_next_file () 
{

if [ $DELETE -eq 1 ] #If you only want the stats
then
echo -e "removing non-stat files"

#delete all folders except stats
rm -r $FOLDER/cap_consensus
rm -r $FOLDER/scaffolds
rm -r $FOLDER/consensus

#delete all files except stats
rm $FOLDER/${FILENAME}_*
rm $FOLDER/cap_consensus
rm $FOLDER/multifastas
rm $FOLDER/scaffolds
rm $FOLDER/consensus
fi

### Delete redundant files
if [ $KEEP -eq 0 ]
then

echo -e "removing redundant files"
if [ -f $(echo $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam) ] 
then
  [ -f $(echo $FOLDER/${FILENAME}_double_single_end_genome.sam) ] && rm $FOLDER/${FILENAME}_double_single_end_genome.sam 
  [ ! -f $(echo $FOLDER/${FILENAME}_double_single_end_genome.bam) ] && rm $FOLDER/${FILENAME}_double_single_end_genome.bam 
  [ ! -f $(echo $FOLDER/${FILENAME}_double_single_end_genome_sorted.bam) ] && rm $FOLDER/${FILENAME}_double_single_end_genome_sorted.bam
fi
rm $FOLDER/temp*
[ ! -f $(echo $FOLDER/${FILENAME}_*.sam) ] || rm $FOLDER/$FILTER_OUT 

#rm $FOLDER/${FILENAME}_trimmed_single_end_reverse.fastq
#rm $FOLDER/${FILENAME}_trimmed_single_end.fastq
else
echo Keeping temporary files is enabled, skipping deletion step.
fi

### If less than one parameter, quit
test=$(echo $others | wc -w)
if [ $test -lt 1 ]
then
    echo parameter count is now $test
    echo "exiting"
    exit
fi 

### If more than one parameter, call self with [other parameters except for first]  
echo parameters: "$0" " -n $CORES $others "
$0 -n $CORES $others 
}

function die_unless () 
{
  if [ ! -f $1 ]
  then  
    echo -e >&2 "Cannot find file '$1'"
    move_on_to_next_file 
    exit ${2:-1}
#  else 
#    echo -e >&2 "Checked that file '$1' exists"
  fi 
}

function time_since () 
{
  SEC=$[`date +%s` - $1]
  MIN=$[$SEC/60]
  HOUR=$[$MIN/60]
  SEC=$[$SEC-($MIN*60)]
  MIN=$[$MIN-($HOUR*60)]
  
  echo "${HOUR}h ${MIN}m ${SEC}s" 
}

function reset_file () 
{
  [[ -f $1 ]] && rm $1
  touch $1
}

function mkdir_if_needed ()
{
[[ ! -d $1 ]] && { mkdir $1; echo "Made directory '$1'"; }
}


function mark_depth () 
{
j=-1
#length1=`echo $1 | tr -d "\n" | wc -c`
#length2=`echo $2 | tr -d "\n" | wc -c`
READ_LIMIT=${3-5}
#if [ ! $length1 -eq $length2 ] 
#then 
#  echo -e >&2 "Code string and read_nr string not equally long!\n '$1'($length1)\n'$2'($length2)" 
#  echo -e "Code string and read_nr string not equally long!\n ($length1) vs ($length2)" >> $FOLDER/cap_consensus_log
#  exit 0 
#fi

#echo -e >&2 "'$1'($length1)\n'$2'($length2)\n"
  while read -n1 i 
    do 
    j=$[j+1]
    read_nr=${2:j:1}
    [[ -z "$read_nr" ]] && echo -e >&2 "i=$i, j=$j, read_nr=$read_nr, read_limit=$READ_LIMIT"
    if [ $read_nr -lt $READ_LIMIT ]
    then 
      echo -n $i | tr 'AGCTN' 'agctn' 
    else
      echo -n $i | tr 'agctn' 'AGCTN'
    fi

  done < <(echo -n "$1" | tr -d "\n")
}

FILTERED=0
function return_unless_DNA_match ()
{
  #sequence=$1 
  avg_len=$1
  while read line; do 
    read read name cigar query MD <<< `echo -e "$line" | cut -f1,3,6,10,12` 
    length=`echo "$query" | wc -c`
    if [ $length -lt $avg_len ] ; then 
	echo "$line"
	continue
    fi
    #echo -e >&2 "line: '$line'"
    #echo -e >&2 "name: '$name'"
    #echo -e >&2 "read: '$read'"
    #echo -e >&2 "cigar: '$cigar'"
    #echo -e >&2 "query: '$query'"
    #echo -e >&2 "length: '$length'"
    #echo -e >&2 "avg_len: '$avg_len'"
    
    #test=`echo $cigar | tr -d '[:digit:]'`
    if [[ "$cigar" =~ [0123456789]+M$ ]] && [[ "$MD" =~ ^MD:Z:[0-9\s]+$ ]]; then 
      echo -e ">$name\t$read ($length bases)" >> $FOLDER/$FILTER_OUT
    else 
    echo "$line"
    fi
  done
}

READ_MINIMUM=30 
MIN_LENGTH_TO_DISCARD=40
function subsample_bam () #Note: will add 'R' or 'F' to read-names in resulting SAM. 
{
  #Create variables
  samfile=$1
  exon=$2
  >&2 echo -e "\e[1A Subsampling ${exon}...                        "

#Outsourced to separate file 
  bash $PIPELINE/subsample_bam.sh $samfile $exon $READ_MINIMUM 
}


function SAM_to_fasta ()
{
while read line; do
  echo ">$line" | cut -f1
  echo "$line" | cut -f10
done
}

function rev_comp () 
{
  echo $1 | tr "ATUGCYRSWKMBDHVNatugcyrswkmbdhvn" "TAACGRYSWMKVHDBNtaacgryswmkvhdbn" | rev
}

function add_leading_zeroes () #Sorts the file after adding leading zeroes to the first-column genome positions.  
{
infile=$1
echo >&2 "Adding leading zeroes to file '$infile'."
max_number_size=`cat $1 | cut -d":" -f2 | cut -d "(" -f1 | tr "-" "\n" | awk '{if (length($1)>length(a)) a=$1;} END{print length(a)}'` #Gets the maximum number of digits, to keep track of how many leading zeroes to add. 

grep '[^[:blank:]]' $infile | gawk -v x=$max_number_size \
'{match($0, /(.*:)([0-9]+)-([0-9]+).([+-]).\s*(.*)/, a); \
b=sprintf("%0"x"d",a[2]); \
c=sprintf("%0"x"d",a[3]); \
print a[1] b"-"c"("a[4]")\t"a[5];}' | \
LC_COLLATE='C' sort -V
}

function get_gene_ID_from_exon ()
{
die_unless $INTERMEDIATES/exon_to_gene.list
exon=`echo $1 | tr -d '[:space:]'`
result=`grep "$exon" "$INTERMEDIATES/exon_to_gene.list" | cut -f2 | head -n1`
#echo >&2 "result is '$result'"
if [ ! -z $result ] 
then 
  echo "$result"; 
else
  echo >&2 "Could not find matching gene for '$exon' in '$INTERMEDIATES/exon_to_gene.list'."
  echo "undefined"
fi
}

function parse_fastq_screen ()
{
infile=$STATS_FOLDER/${FILENAME}_trimmed_single_end_min36_screen.txt
match_ref=0
match_both=0
match_contam=0
match_none=0
contam_total=0

while read line ; do

#match no-hits line
if [[ "$line" == %* ]]
then 
  match_none=`echo $line | cut -f2 -d ' '`
  continue
fi

#match reference line
if [[ $line == Malaria* ]] 
then
  read total a b c d <<< `echo $line | cut -d ' ' -f2,5,7,9,11`
  match_ref=$[ (a + b)*100/total ]
  match_both=$[ (c + d)*100/total ]
  continue
fi

#match other lines
if [[ `echo $line | wc -c` -gt 10 ]] 
then
  read total a b <<< `echo $line | cut -d ' ' -f2,5,7`
#  echo "line is '$line', contam_total is $contam_total, a is $a, b is $b"
  contam_total=$[ contam_total + a + b ]
  continue
fi
done < <( cat $infile | grep -v "^#" | grep -v "Genome" | grep ".." )

match_contam=$[contam_total*100/total]
echo "% of Reads matching to reference: ${match_ref}%"
echo "Reference AND contaminant match:  ${match_both}%"
echo "Reads matching one contamination: ${match_contam}%"
echo "Reads matching none of the above: ${match_none}%"
}


#####################
### GET ARGUMENTS ###
#####################

while [[ $* ]]
do
    OPTIND=1
#    echo $1
    if [[ $1 =~ ^- ]]
    then
        getopts :f:s:n:a:b:d:i:g:th parameter
        case $parameter in
            h)  echo "$HELP"
		exit
                ;;
            f)  FOLDER=$OPTARG #choose name of output folder
		#Currently only works on single sample. 
                echo "FOLDER=$FOLDER"
                shift
                ;;
            a)  ADAPTER=$OPTARG #set path to trimmomatic adapter file
                echo "ADAPTER=$ADAPTER"
                shift
                ;;
            n)  CORES=$OPTARG 
                echo "CORES=$CORES"
                shift
                ;;            
	    b)  STRINGENCY=$OPTARG #Choose Bowtie2 stringency setting (0-6)
                echo "STRINGENCY=$STRINGENCY"
                shift
                ;;
            s)  SKIP=$OPTARG #If SKIP=0, make folder files from scratch. 
                echo "SKIP=$SKIP"
                shift
                ;;
            i)  INTERM=$OPTARG #If INTERM=0, remake intermediate files from scratch. 
                echo "INTERM=$INTERM"
                shift
                ;;
            t)  TEST=1 #Test that command-line-option system works properly. 
                echo "TEST=$TEST"
                shift
                ;;
            g)  CALLING=$OPTARG #Choose variant calling method.
                echo "CALLING=$CALLING"
                shift
                ;;
            d)  DELETE=$OPTARG #Delete result files afterwards, except stats.
                echo "DELETE=$DELETE"
                shift
                ;;
            *) echo "This is an invalid argument: '$parameter'"
                ;;
        esac
    else
        other_arguments="$other_arguments $1"
    fi
    shift
done

#echo "arguments - options: $other_arguments"

########################################
### ELIMINATE DOUBLES FROM ARGUMENTS ###
########################################

#Removes any arguments not ending in .fastq 
#Replaces "_R2_001.fastq" with "_R1_001.fastq"
#Removes duplicate entries

other_arguments=`echo $other_arguments | tr " " "\n" | grep ".fastq" | sed s/_R2_001.fastq/_R1_001.fastq/ | sort | uniq | tr "\n" " "`
#echo "arguments after duplicate removal: $other_arguments"

########################################

IFS=' ' read -r first_word others <<< "$other_arguments"
#echo "first_word in 'others' is $first_word"
FILENAME=$( echo "$first_word" | sed "s/_R[12]_001.fastq.*//" )
FILENAME=$( echo "$FILENAME" | sed "s/.fastq.*//" ) #If infiles don't have standard naming or are zipped
if [ -z "$FILENAME" ]; then exit; fi
if [ -z "$FOLDER" ] ; then FOLDER=$FILENAME; fi
if [ -z "$SKIP"  ]; then SKIP=1; fi

#Tests whether this is the last iteration
if [ -z "${others// }" ]      #If other arguments doesn't contain anything (but spaces)
then 
    LAST_ITERATION=1                     #Then this is the last iteration of the program. 
fi



###########################################

#####################
### SCRIPT START: ###
#####################

intro="\
New pipeline started at `date` on sample $FILENAME\n
Arguments: '$ARGUMENTS'\n
Target files after duplicate removal: $other_arguments\n
Process ID: $BASHPID\n
FOLDER is $FOLDER\n
FILENAME = $FILENAME\n
PIPELINE = $PIPELINE\n
GFF_FILE = $GFF_FILE\n
GENOME_FILE = $GENOME_FILE\n
ADAPTER = $ADAPTER\n
INTERM = $INTERM\n
CORES = $CORES\n
SKIP = $SKIP\n
STRINGENCY = $STRINGENCY\n
LAST_ITERATION = $LAST_ITERATION\n"

echo -e $intro

### Make folder to put new files in. 
echo Creating folder: $FOLDER
mkdir_if_needed $FOLDER
STATS_FOLDER="$FOLDER/fastq_stats"

### Make intermediates folder
[[ ! -d $INTERMEDIATES ]] && { 
  mkdir $INTERMEDIATES; 
  INTERM=0; #Redoing intermediates
  echo "Made directory '$INTERMEDIATES'"; }
mkdir_if_needed $INTERMEDIATES

### Check for basic parameters Test. If so, EXIT
if [ $TEST -eq 1 ] 
then
echo "Halting script because parameter test [-t] option was set."
exit
fi


### Make logfile
LOG=$FOLDER/log.txt

if [ $SKIP -eq 0 ] || [ `ls $FOLDER | wc -l` -lt 3 ]
then
touch $LOG
echo -e "Settings for creating these files:\n" > $LOG
else
echo logfile already exists, skipping.
fi

### Note in log that pipeline has run
echo -e "Called pipeline at `date`.\n $intro" >> $LOG


### Check base files
if [ $SKIP -eq 0 ] || [ $INTERM -eq 0 ]
then
# Check GFF file 
dos2unix $GFF_FILE
cat $GFF_FILE | grep -e "CDS" -e "^##" | tr "\t" "," | sed -r "s/,$/,NA/g" | tr "," "\t" > $BASE_FILES/temp_gff
cat $BASE_FILES/temp_gff > $GFF_FILE
rm $BASE_FILES/temp_gff

# Check Genome file
dos2unix $GENOME_FILE
fi

#########################################################################################
#### MAKING INTERMEDIATE (REUSEABLE) FILES ####
#########################################################################################

echo -e "\nChecking prerequisite files"

infile=$GENOME_FILE
outfile=$INTERMEDIATES/HtExons.fasta
### Make HtExons.fasta if not present
if [ $INTERM -eq 0 ] || [ ! -f $(echo $outfile) ]
then
die_unless $GFF_FILE
die_unless $infile 
echo -e "\ngenerating HtExons.fasta"
bedtools getfasta -fi $infile -bed $GFF_FILE -fo $outfile -s

echo made HtExons.fasta from gff file
else
echo $outfile already exists, skipping.
fi

#################################################################
### Make HtExons.fasta.fai if not present
infile=$INTERMEDIATES/HtExons.fasta
outfile=$INTERMEDIATES/HtExons.fasta.fai
if [ $INTERM -eq 0 ] || [ ! -f $(echo $outfile) ]
then
die_unless $infile
echo -e "\ngenerating $outfile"

$samtools faidx $infile

echo made $outfile
else
echo $outfile already exists, skipping.
fi

#################################################################
### Make HtGenome.fasta.fai if not present
infile=$GENOME_FILE
outfile=$GENOME_FILE.fai #Note: decided by samtools
if [ $INTERM -eq 0 ] || [ ! -f $(echo $outfile) ]
then
die_unless $infile
echo -e "\ngenerating $outfile"

$samtools faidx $infile

echo made $outfile
else
echo $outfile already exists, skipping.
fi

#################################################################
### Make Gene-exon database list if not present
infile="$GFF_FILE"
outfile1=$INTERMEDIATES/exon_to_gene.list
outfile2=$INTERMEDIATES/exon_to_gene_sorted.list
if [ $INTERM -eq 0 ] || [ ! -f $(echo $outfile1) ] || [ ! -f $(echo $outfile2) ]
then
echo -e "\ngenerating $outfile" 
die_unless $infile

#Make list of exon name, gene id, and length
if grep -q "gene_id" "$GFF_FILE" #If gff format uses "gene_id"
then
  cat $GFF_FILE | gawk 'match($0, /gene_id "([^"]+)/, x) { print $1":"$4-1"-"$5"("$7")\t"x[1]"\t"$5-$4+1 }' | tr '\/' '||' | sort >  $outfile1

elif grep -q "ID=" "$GFF_FILE"  #If gff format uses "ID"
then
  cat $GFF_FILE | gawk 'BEGIN{gene=0} /CDS/ {gene_found = match($0, /ID=([^;]+)/, x)
 if(gene_found == 0) x[1]="undefined.gene."++gene; else test="3"; print $1":"$4-1"-"$5"("$7")\t"x[1]"\t"$5-$4+1 }' | tr "\/" "||" > $outfile1
else 
    echo -e >&2 "Unfamiliar GFF file format"
    exit ${2:-1}
fi

#Find the largest number in the file, measure its length
#max_number_size=`cat $outfile | cut -d":" -f2 | cut -d "(" -f1 | tr "-" "\n" | awk '{if (length($1)>length(a)) a=$1;} END{print length(a)}'`

echo "" > $outfile2 # Should start with a newline
add_leading_zeroes $outfile1 >> $outfile2

#Add leading zeroes to every length
#grep '[^[:blank:]]' $outfile | gawk -v x=$max_number_size \
#'{match($0, /(.*:)([0-9]+)-([0-9]+).([+-]).\s*(.*)/, a); \
#b=sprintf("%0"x"d",a[2]); \
#c=sprintf("%0"x"d",a[3]); \
#print a[1] b"-"c"("a[4]")\t"a[5];}' | \
#LC_COLLATE='C' sort -V >> $outfile2

echo "made $outfile1 and $outfile2"
else
echo $outfile already exists, skipping.
fi


#################################################################
### Make probe-to-exon-to-gene list for statistics purposes
infile_probes=$PROBE_FILE 
infile_exons=$INTERMEDIATES/exon_to_gene.list
outfile=$INTERMEDIATES/probe_to_exon_to_gene.list
if [ $INTERM -eq 0 ] || [ ! -f $(echo $outfile) ]
then

#Useful function; outsourced to separate script. 
bash $PIPELINE/probes_to_exons.sh "$infile_probes" "$infile_exons" > "$outfile" 

echo made $outfile
else
echo $outfile already exists, skipping.
fi

####################################################################
### Make exon fasta file with ONLY probe-targeted exons, for statistics purposes. 

infile=$INTERMEDIATES/probe_to_exon_to_gene.list
infile2=$INTERMEDIATES/HtExons.fasta
outfile=$INTERMEDIATES/relevant_exons.fasta
outfile2=$INTERMEDIATES/relevant_exons.bed
if [ $INTERM -eq 0 ] || [ ! -f $(echo $outfile) ] || [ ! -f $(echo $outfile2) ]
then
echo "making $outfile, $outfile2"
die_unless $infile
die_unless $infile2
reset_file $outfile
reset_file $outfile2

for line in `cat $infile | cut -f3 | uniq | cut -f1 -d "(" ` ; do 
    grep -A1 "$line" "$infile2" >> $outfile
    chr=`echo $line | cut -f1 -d ':'`
    pos=`echo $line | cut -f2 -d ':'`
    read start stop <<< `echo $pos | tr '-' ' '`
    echo -e "$chr\t$start\t$stop" >> $outfile2
done

echo made $outfile and $outfile2
else
echo $outfile and $outfile2 already exist, skipping.
fi


####################################################################
### Make list of exon gc count, if not present
infile=$INTERMEDIATES/HtExons.fasta
outfile=$INTERMEDIATES/exons_gc_count.list
outfile2=$INTERMEDIATES/exons_gc_count_sorted.list
if [ $INTERM -eq 0 ] || [ ! -f $(echo $outfile) ] || [ ! -f $(echo $outfile2) ]
then
echo -e "\ngenerating $outfile" 
die_unless $infile

#Make column with GC count, then sort it. 
echo "" > $outfile
cat $infile | paste - - | while read line 
  do read exon code <<< $line 
  #echo "exon='$exon', code='$code'"
  exon=`echo $exon | tr -d '>'`
  echo -e "$exon\t`echo $code | tr -d -c 'CGcg' | wc -c`" >> $outfile 
done

#Find the largest number in the file, and measure its length
#max_number_size=`cat $outfile | cut -d":" -f2 | cut -d "(" -f1 | tr "-" "\n" | awk '{if (length($1)>length(a)) a=$1;} END{print length(a)}'`

#echo "" > $outfile2 # Should start with a newline?
add_leading_zeroes $outfile > $outfile2
#Add leading zeroes to every number
#grep '[^[:blank:]]' $outfile | gawk -v x=$max_number_size \
#'{match($0, /(.*:)([0-9]+)-([0-9]+).([+-]).\s*(.*)/, a); \
#b=sprintf("%0"x"d",a[2]); c=sprintf("%0"x"d",a[3]); \
#print a[1] b"-"c"("a[4]")\t"a[5];}' | \
#LC_COLLATE='C' sort >> $outfile2
echo "made $outfile and $outfile2"
else
echo $outfile2 already exists, skipping.
fi

#################################################################
### Make detailed exons vs genes file with gc content. 
infile=$INTERMEDIATES/exon_to_gene_sorted.list
infile2=$INTERMEDIATES/exons_gc_count_sorted.list
outfile=$INTERMEDIATES/Exons_vs_Genes.tsv
outfile2=$INTERMEDIATES/gc_content.tsv

if [ $INTERM -eq 0 ] || [ ! -f $(echo $outfile) ] || [ ! -f $(echo $outfile2) ]
then
echo -e "\ngenerating $outfile" 
mkdir_if_needed Consensus
die_unless $infile
die_unless $infile2

#Add exon GC count to the exon-to-gene file
join -1 1 -2 1 -a1 -eERROR -o '0,1.2,1.3,2.2' $infile $infile2 | sort -k2,1 > $outfile
echo made $outfile

#Adding header to Exon-vs-Genes file
echo -e "\ngenerating $outfile2" 
echo -e "\t\t\t\t" > $outfile2
echo -e "#queryName\tqueryLength\tComment\thitDescription\thitName\t%GC_ref" >> $outfile2
cat $outfile >> $outfile2

echo made $outfile2
else
echo $outfile2 already exists, skipping.
fi

#################################################################
### Make GC percentage files for separate statistics calculations
infile=$INTERMEDIATES/Exons_vs_Genes.tsv
outfile=$INTERMEDIATES/exons_gc_content_sorted.tsv
outfile2=$INTERMEDIATES/gc_content_sorted.tsv
if [ $INTERM -eq 0 ] || [ ! -f $(echo $outfile) ] || [ ! -f $(echo $outfile2) ]
then
echo -e "\ngenerating $outfile" 
die_unless $infile
reset_file $outfile
reset_file $outfile2

#Compile GC percentages per exon
cat $infile | while read a b c d ; do echo -e "$a\t`echo "scale=1; 100 * $d / $c " | bc`" ; done > $outfile

#Compile exon-specific data into gene-specific data with GC percentages per gene
awk '{a[$2]+=$3; b[$2]+=$4}END{for (i in a){c=100*b[i]/a[i]; printf "%s\t%s\t%3.1f\n", i, a[i], c;}}' $infile | sort -k1 > $outfile2

echo made $outfile
else
echo $outfile already exists, skipping.
fi

      ################
############################  
###   Checking Indices   ###
############################
      ################

#################################################################
### Build genome index with bowtie2. 
infile=$GENOME_FILE
example=$INTERMEDIATES/HtIndex.1.bt2
echo -e "\nchecking indices"
if [ $INTERM -eq 0 ] || [ ! -f $(echo $example) ]
then
echo constructing bowtie2 index for genome
die_unless $infile
BOWTIE2_OPTIONS="--threads $CORES -f"
bowtie2-build $BOWTIE2_OPTIONS \
$infile \
$INTERMEDIATES/HtIndex
echo -e "Bowtie2-build settings, genome: $BOWTIE2_OPTIONS \n" >> $LOG
echo built index for exons with basename $INTERMEDIATES/HtIndex
else
echo $example already exists, skipping.
fi

#################################################################
### Build exon index with bowtie2. 
infile=$INTERMEDIATES/HtExons.fasta
example=$INTERMEDIATES/HtExonIndex.1.bt2
if [ $INTERM -eq 0 ] || [ ! -f $(echo $example) ]
then
echo constructing bowtie2 index for exons
die_unless $infile
BOWTIE2_OPTIONS="--threads $CORES -f"
bowtie2-build  $BOWTIE2_OPTIONS \
$infile \
$INTERMEDIATES/HtExonIndex
echo -e "Bowtie2-build settings, exons: $BOWTIE2_OPTIONS \n" >> $LOG
echo built index for exons with basename $INTERMEDIATES/HtExonIndex
else
echo $example already exists, skipping.
fi

#################################################################
### Make dictionary of genome fasta files
infile=$GENOME_FILE
outfile=$INTERMEDIATES/HtGenome.dict
if [ $INTERM -eq 0 ] || [ ! -f $(echo $outfile) ]
then
echo -e "\nMaking dictionary from $infile"
[[ -f $(echo $outfile) ]] && rm $outfile

java -jar $picard CreateSequenceDictionary R=$infile O=$outfile

else
echo $outfile already exists, skipping.
fi

#################################################################
### Make dictionary of exon fasta files
infile=$INTERMEDIATES/HtExons.fasta
outfile=$INTERMEDIATES/HtExons.dict
if [ $INTERM -eq 0 ] || [ ! -f $(echo $outfile) ]
then
echo -e "\nMaking dictionary from $infile"
[[ -f $(echo $outfile) ]] && rm $outfile

java -jar $picard CreateSequenceDictionary R=$infile O=$outfile

else
echo $outfile already exists, skipping.
fi


##################################################################
### Index fastas for fastq_screen contamination search
changes=0
while read line; do 
  name=`echo $line | cut -d'.' -f1`
  if [ ! -f ${name}_index.1.bt2 ] 
  then 
    echo indexing $name
    bowtie2-build --threads $CORES -f $line ${name}_index > /dev/null 2> /dev/null
    changes=1
  fi
done < <(find $PIPELINE/fastq_screen_fastas/ -iname "*.fasta" -o -iname "*.fa")

### If any new sequence was detected, make a new conf file.
if [ $INTERM -eq 0 ] || [ $changes -eq 1 ] || [ ! -f $(echo $FASTQ_SCREEN_CONF) ]
then 
echo "updating fastq_screen.conf" 
updated_list="##Malaria\nDATABASE\tMalaria\t$INTERMEDIATES/HtIndex"
while read line; do 
  path=`echo $line | cut -d'.' -f1`
  name=`echo ${path##*/}`
  updated_list="$updated_list\n##${name}\nDATABASE\t${name}\t${path}_index"
done < <(find $PIPELINE/fastq_screen_fastas/ -iname "*.fasta" -o -iname "*.fa")

#Concatenate database list to original conf file. 
cat $FASTQ_SCREEN_CONF_BEGIN > $FASTQ_SCREEN_CONF
echo -e $updated_list >> $FASTQ_SCREEN_CONF

else
echo fastq_screen fastas already indexed, skipping.
fi


###########################################################
### Make BED file 
infile=$INTERMEDIATES/HtExons.fasta
outfile=$INTERMEDIATES/HtExons.bed
if [ $INTERM -eq 0 ] || [ ! -f $(echo $outfile) ]
then 
echo "making bed file" 
die_unless $infile
reset_file $outfile
while read line
do
scaffold_name=`echo "$line" | cut -d ":" -f1`
interval=`echo "$line" | grep -oP "\d+-\d+"`
IFS='-'; read start end <<< "$interval"
orientation=`echo "$line" | grep -oP "\(.\)+" | tr -d "()"`
echo -e "$scaffold_name\t$start\t$end\t$line\t.\t$orientation" >> $outfile 
done < <( cat "$infile" | grep "^>" | tr -d ">" ) 
IFS=$' \t\n'
else
echo $outfile already exists, skipping.
fi
echo ""



################################################################
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
################################################################
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
################################################################


      ##############
##########################
###   Trimming Files   ###
##########################
      ##############


###################################################
################################
### Run Trimmomatic, single-end:
infile=${FILENAME}_R1_001.fastq

outfile1=$FOLDER/${FILENAME}_trimmed_single_end.fastq
outfile2=$FOLDER/${FILENAME}_trimmed_single_end_min36.fastq
if [ $SKIP -eq 0 ] || ([ ! -f $(echo $outfile1) ] && [ ! -f $(echo $outfile2) ])
then
echo -e "\nrunning trimmomatic, single-end"
START_TIME=`date +%s`

#Accept zipped files
[[ -f ${infile}.gz ]] && infile=${infile}.gz
[[ -f ${infile}.gz ]] && infile=${infile}.bz2
die_unless $infile

TRIM_OPTIONS="-phred33 -threads $CORES MINLEN:36 SLIDINGWINDOW:4:15 ILLUMINACLIP:${ADAPTER}:2:30:10"

echo -e "\nCommand: $infile java -jar $PIPELINE/third_party_programs/Trimmomatic-0.36/trimmomatic-0.36.jar SE $infile $outfile1 $TRIM_OPTIONS"
java -jar $trimmomatic SE $infile $outfile1 $TRIM_OPTIONS

echo -e "Trimmomatic settings: $TRIM_OPTIONS" >> $LOG
echo -e "Trimmomatic(forward) took" `time_since $START_TIME` "to run. \n" >> $LOG
else
echo $outfile1 already exists, skipping.
fi

###########################################################3
########################################
### Run Trimmomatic, single-end reverse:
infile=${FILENAME}_R2_001.fastq
outfile1=$FOLDER/${FILENAME}_trimmed_single_end_reverse.fastq
outfile2=$FOLDER/${FILENAME}_trimmed_single_end_reverse_min36.fastq
if [ $SKIP -eq 0 ] || ([ ! -f $(echo $outfile1) ] && [ ! -f $(echo $outfile2) ])
then
echo -e "\nrunning trimmomatic, single-end reverse"
START_TIME=`date +%s`

[[ -f ${infile}.gz ]] && infile=${infile}.gz
[[ -f ${infile}.gz ]] && infile=${infile}.bz2
die_unless $infile

TRIM_OPTIONS="-phred33 -threads $CORES MINLEN:36 SLIDINGWINDOW:4:15 ILLUMINACLIP:${ADAPTER}:2:30:10"
#java -jar ~/bin/trimmomatic-0.36.jar SE \
java -jar $trimmomatic SE $infile $outfile1 $TRIM_OPTIONS

echo -e "Trimmomatic settings: $TRIM_OPTIONS" >> $LOG
echo -e "Trimmomatic(reverse) took" `time_since $START_TIME` "to make. \n" >> $LOG
else
echo $outfile1 already exists, skipping.
fi

### Remove small reads from trimmed files (check whether necessary):
if [ $SKIP -eq 0 ] || [ ! -f $(echo $FOLDER/${FILENAME}_trimmed_single_end_min36.fastq) ]
then
echo removing small reads from single-end
START_TIME=`date +%s`
die_unless $FOLDER/${FILENAME}_trimmed_single_end.fastq
die_unless $FOLDER/${FILENAME}_trimmed_single_end_reverse.fastq

cat $FOLDER/${FILENAME}_trimmed_single_end.fastq \
| paste - - - - \
| awk  'length($3)>35' \
| tr '\t' '\n' \
> $FOLDER/${FILENAME}_trimmed_single_end_min36.fastq

# Also, mark reverse reads AS reverse.
cat $FOLDER/${FILENAME}_trimmed_single_end_reverse.fastq \
| sed '1~4s/ /R /' \
| paste - - - - \
| awk  'length($3)>35' \
| tr '\t' '\n' \
> $FOLDER/${FILENAME}_trimmed_single_end_reverse_min36.fastq

echo -e "Secondary trimmings took" `time_since $START_TIME` "to do." >> $LOG
else
echo $FOLDER/${FILENAME}_trimmed_single_end_min36.fastq already exists, skipping.
fi

      ###############
###########################
###   Quality Control   ###
###########################
      ###############


#######################################################################
###########################################
### Do fastqc on trimmed fastq files. 
infile1=$FOLDER/${FILENAME}_trimmed_single_end_min36.fastq
infile2=$FOLDER/${FILENAME}_trimmed_single_end_reverse_min36.fastq
outfile1=$STATS_FOLDER/${FILENAME}_trimmed_single_end_min36_fastqc.html
if [ $SKIP -eq 0 ] || [ ! -f $(echo $outfile1) ]
then
echo making fastqc file
START_TIME=`date +%s`
die_unless $infile1
die_unless $infile2
mkdir_if_needed $STATS_FOLDER

fastqc $infile1 --outdir=$STATS_FOLDER
fastqc $infile2 --outdir=$STATS_FOLDER
echo -e "FastQC took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $outfile1 already exists, skipping.
fi


#######################################################################
###########################################
### Use Fastq_screen on trimmed fastq files to check for mosquito, bird and human contamination. 
input_file=$FOLDER/${FILENAME}_trimmed_single_end_min36.fastq
output_file=$STATS_FOLDER/${FILENAME}_trimmed_single_end_min36_screen.txt
if [ $SKIP -eq 0 ] || [ ! -f $(echo $output_file) ]
then
echo -e "\nscreening for contamination"
START_TIME=`date +%s`
die_unless $input_file
mkdir_if_needed $STATS_FOLDER

fastq_screen $input_file --aligner bowtie2 --conf $FASTQ_SCREEN_CONF --outdir $STATS_FOLDER --threads $CORES --force

echo -e "fastq_screen took" `time_since $START_TIME` "to run.\n" >> $LOG

else
echo $output_file already exists, skipping.
fi


      ########
####################
###   Matching   ###
####################
      ########


#######################################################################
###########################################
##### Matching double single-end genome
infile1=$FOLDER/${FILENAME}_trimmed_single_end_min36.fastq
infile2=$FOLDER/${FILENAME}_trimmed_single_end_reverse_min36.fastq 
current_file=$FOLDER/${FILENAME}_double_single_end_genome.sam
descendant1=$FOLDER/${FILENAME}_double_single_end_genome_sorted.bam
descendant2=$FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam
if [ $SKIP -eq 0 ] || ( [ ! -f $(echo $current_file) ] && [ ! -f $(echo $descendant1) ] && [ ! -f $(echo $descendant2) ] ) 
then
echo -e "\nstarting bowtie2, single-end genome, $current_file"
START_TIME=`date +%s`
die_unless $infile1
die_unless $infile2

#Stringency
case ${STRINGENCY-4} in 
1) 
BOWTIE2_OPTIONS="--local --threads $CORES -q -D 22.9 -R 3 -N 1 -L 23 -i S,1,0.16"
;;
2) 
BOWTIE2_OPTIONS="--local --threads $CORES -q -D 20 -R 2 -N 1 -L 20 -i S,1,0.5"
;;
3) 
BOWTIE2_OPTIONS="--local -N 1 --threads $CORES -q "
;;
4) 
BOWTIE2_OPTIONS="--very-sensitive-local --threads $CORES -q "
;;
5) 
BOWTIE2_OPTIONS="--sensitive-local --threads $CORES -q "
;;
*) 
BOWTIE2_OPTIONS="--threads $CORES -q "
;;
esac

bowtie2 $BOWTIE2_OPTIONS \
-U $infile1 \
-U $infile1 \
-x $INTERMEDIATES/HtIndex \
-S $current_file

echo -e "Bowtie2 settings, exons: $BOWTIE2_OPTIONS" >> $LOG
echo -e "SAMfile(genome), $current_file, took" `time_since $START_TIME` "to make. \n" >> $LOG
else
echo $current_file already exists, skipping.
fi


    ##########
##################
### SAM TO BAM ###
##################
    ##########


#######################################################################
###########################################
### Convert genome SAMfile to BAM format, sort it, add read group and index it. 
input_file=$FOLDER/${FILENAME}_double_single_end_genome.sam
output_file=$FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam
fastq_file=$FOLDER/${FILENAME}_trimmed_single_end_min36.fastq
if [ $SKIP -eq 0 ] || ([ ! -f $(echo $output_file) ] && [ ! -f $(echo $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam) ])
then
echo -e "\nconverting genome match to BAM"
START_TIME=`date +%s`
die_unless $input_file
die_unless $fastq_file

### Code for this was outsourced to a separate file; it is convenient to be able to call the function on its own. 
# fastq_file is used to find information for read group assignment.
# (file name is backup for read_group data; a default name like '512022_S1_L001' contains data on sample name and lane nr.)

bash $PIPELINE/SAM_to_indexed_BAM.sh $input_file $output_file $fastq_file $FILENAME 


echo -e "SAM_to_indexed_BAM took" `time_since $START_TIME` "to run.\n" >> $LOG

## Delete redundant file
rm $FOLDER/${FILENAME}_double_single_end_genome.sam

else
echo $output_file already exists, skipping.
fi



      #############
#########################
#### VARIANT CALLING ####
#########################
      #############


#######################################################################
###########################################
### Variant Calling, samtools.
if [ $CALLING -eq 0 ]
then
if [ $SKIP -eq 0 ] || [ ! -f $(echo $FOLDER/${FILENAME}_variants.vcf) ]
then

echo -e "\nvariant calling, genome"
START_TIME=`date +%s`
die_unless $GENOME_FILE
die_unless $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam

PILEUP_OPTIONS="-g -f"
BCF_OPTIONS="-m -v"
$samtools mpileup $PILEUP_OPTIONS $GENOME_FILE $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam > $FOLDER/${FILENAME}_variants.bcf

die_unless $FOLDER/${FILENAME}_variants.bcf
$bcftools call $BCF_OPTIONS $FOLDER/${FILENAME}_variants.bcf > $FOLDER/${FILENAME}_variants.vcf

echo -e "$samtools settings, pile-up genome: $PILEUP_OPTIONS" >> $LOG
echo -e "$samtools settings, call genome: $BCF_OPTIONS" >> $LOG
echo -e "Variant calling(genome) took" `time_since $START_TIME` "to do. \n" >> $LOG
else
echo $FOLDER/${FILENAME}_variants.vcf already exists, skipping.
fi


#######################################################################
###########################################
### Sorting VCF file
if [ $SKIP -eq 0 ] || [ ! -f $(echo $FOLDER/${FILENAME}_variants.bcf.csi) ]
then
echo -e "\nsorting VCF file"
die_unless $FOLDER/${FILENAME}_variants.bcf

START_TIME=`date +%s`
$bcftools index $FOLDER/${FILENAME}_variants.bcf

echo -e "Sorting VCF took" `time_since $START_TIME` "to make." >> $LOG
else
echo $FOLDER/${FILENAME}_variants_bcf.csi already exists, skipping.
fi
fi

#######################################################################
###########################################
### GATK Variant Calling. 
if [ $CALLING -eq 1 ]
then
current_file=$FOLDER/${FILENAME}_genome_variants_gatk.vcf
if [ $SKIP -eq 0 ] || [ ! -f $(echo $current_file) ]
then
echo "variant calling, genome GATK"
die_unless $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam

START_TIME=`date +%s`
PILEUP_OPTIONS="--variant_index_type LINEAR --variant_index_parameter 128000" #<---Check that this number is reasonable 

java -jar $PIPELINE/third_party_programs/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R $GENOME_FILE \
-I $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam \
-o $FOLDER/${FILENAME}_genome_variants_gatk.vcf $PILEUP_OPTIONS

echo -e "Variant calling settings, GATK: $PILEUP_OPTIONS" >> $LOG
echo -e "Variant calling (GATK) took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $current_file already exists, skipping.
fi
fi

      #############
#########       ##########
#### CONSENSUS-MAKING ####
#########       ##########
      #############



#######################################################################
###########################################
### Make consensus sequence.
infile="$FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam"
req_file="$infile.bai"
ref_file="$GENOME_FILE"
outfile="$FOLDER/${FILENAME}_consensus.fasta"
if [ $SKIP -eq 0 ] || [ ! -f $(echo $outfile) ] || [ ! -d $(echo $FOLDER/scaffolds) ]
then
echo -e "\nmaking consensus sequence"
die_unless $infile
die_unless $ref_file
die_unless $req_file

START_TIME=`date +%s`
mkdir_if_needed $FOLDER/scaffolds
touch $outfile
echo -e "starting...\n\n"

#Get SAM files of exons from genome match

cat $ref_file | grep "^>" | \
while read line trash; do 
name=`echo $line | tr -d '>' | sed "s/\s//"`
if [ `$samtools view $infile $name -c` -eq 0 ]; then
#echo -en "\e[1A"; echo -e "\e[0K\rskipping $name"
continue 
fi
#echo -en "\e[1A"; echo -e "\e[0K\rFinding reads from <${name}>";\
scaffold_length=`cat $ref_file | grep $name -A1 | tail -n1 | tr -d "\n" | wc -c`
{
$samtools mpileup -uf $ref_file $infile -r $name -d 100000 | $bcftools call -c --ploidy 1 > $FOLDER/scaffolds/${name}_consensus.vcf;\
} > /dev/null 2> /dev/null #Silence this script, so it doesn't report every time it's run.
#echo did $FOLDER/scaffolds/${name}_consensus.fastq; \
scaffold_sequence=`cat $FOLDER/scaffolds/${name}_consensus.vcf | perl $PIPELINE/third_party_programs/vcfutils.pl vcf2fq | seqtk seq -A | tail -n1 | tr -d "\n"` 
length_difference=$[scaffold_length-`echo $scaffold_sequence | wc -c`] 
#echo "length=$scaffold_length, length_difference=$length_difference"
if [ $length_difference -gt 0 ] ; then # samtools mpileup mislays the end of the sequence unless it has reads. 
 # echo "padding sequence with extra 'n'"
  extra_n=`printf '%*s' "$length_difference" | tr ' ' "n"`
  scaffold_sequence=`echo "${scaffold_sequence}$extra_n"`
fi
echo -e ">$name\n$scaffold_sequence" > $FOLDER/scaffolds/${name}_consensus.fasta; \
echo -e ">$name\n$scaffold_sequence" >> $outfile; \
echo -en "\e[1A"; echo "created $FOLDER/scaffolds/${name}_consensus.fasta" 
done
echo -en "\e[1A"; echo -e "\e[0K\rDone";\

echo -e "Samtools pileup consensus file took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $outfile already exists, skipping.
fi


#######################################################################
###########################################
#  Create exon consensus file. 
infile=$FOLDER/${FILENAME}_consensus.fasta
infile2=$INTERMEDIATES/HtExons.bed
outfile=$FOLDER/${FILENAME}_exon_consensus2.fasta
if [ $SKIP -eq 0 ] || [ ! -f $outfile ] || [ ! -f $outfile2 ]
then
echo -e "Creating $outfile, exon consensus file"
die_unless $infile

START_TIME=`date +%s`

#Cut out exons from consensus
bedtools getfasta -fi $infile -bed $infile2 -name -fo $FOLDER/temp

#Make version where exons without coverage are excluded.
reset_file $outfile
while read name sequence
do
  [[ $sequence =~ [ATGCatgc] ]] || continue
  echo -e "$name\n$sequence" >> $outfile
done < <( cat "$FOLDER/temp" | paste - - ) 
rm $FOLDER/temp

echo -e "$outfile took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $outfile already exists, skipping.
fi

#######################################################################
###########################################
#  Create coverage stats file.
infile=$FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam
outfile=$FOLDER/exon_coverage.tsv
ref_file=$INTERMEDIATES/exon_to_gene.list
if [ $SKIP -eq 0 ] || [ ! -f $outfile ]
then
echo -e "Creating $outfile, exon coverage stats file"
die_unless $infile
die_unless $ref_file
touch $outfile
echo "starting..."
START_TIME=`date +%s`

cat $ref_file | while read exon gene length ; do 
  echo >&2 -en "\e[1A" "Processing $exon                           \n" 
  position=`echo $exon | tr -d ">" | cut -f1 -d '('`; 
  echo -en "$position\t"; 
  echo -en "$length\t"; 
  $samtools depth -r "$position" "$infile" | awk -v l=$length '{if($3>2) mapped+=1; else unmapped+=1}END{if (mapped>0) print mapped/(l+1); else print 0}' ; 
done > $outfile

echo -e "$outfile took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $outfile already exists, skipping.
fi


#######################################################################
###########################################
#  Change exon consensus file capitalization depending on depth of coverage.
infile=$FOLDER/${FILENAME}_exon_consensus2.fasta
outfile=$FOLDER/${FILENAME}_exon_consensus3.fasta
if [ $SKIP -eq 0 ] || [ ! -f $outfile ]
then
echo -e "Creating $outfile, depth-capitalized exon consensus file"
die_unless $infile
mkdir_if_needed $FOLDER/cap_consensus
mkdir_if_needed $FOLDER/cap_consensus/exons
touch $outfile
reset_file $FOLDER/cap_consensus_log
echo "starting..."
START_TIME=`date +%s`

export FOLDER
export FILENAME
export INTERMEDIATES
export READ_LIMIT
export -f mark_depth 
export infile
export samtools

cat $infile | grep "^>" | \
xargs -I{} --max-procs $CORES bash -c 'name="{}"; 
exon_name=`echo $name | tr -d ">" | cut -d "(" -f1`;
echo -en "\e[1A" "Processing $exon_name                    \n";
exon_code=`cat "$infile" | grep "$name" -A1 | tail -n1`; #In case of exons too long for xargs
#Ignore if no reads
read_nr=`$samtools view $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $exon_name | wc -l`;
if [[ $read_nr -eq 0 ]] ; then 
echo -e "\e[1A" "ERROR: No reads in $exon_name\n"
exit 0
fi

#Note: Maximum readcount capped at 9 
exon_read_counts=`$samtools depth -aa "$FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam" -r "$exon_name" | cut -f3 | sed -e "s/[0-9]\{2,\}/9/g" | tr -d "\n" | tail -c +2`;  
current_fasta="$FOLDER/cap_consensus/exons/$exon_name.fasta";
if [[ "$name" == *"(-)"* ]] 
then
exon_read_counts=`echo "$exon_read_counts" | rev`;
fi
if [[ $exon_code =~ [ATGCNatgcn-]+ ]] ; then
  length1=`echo $exon_code | tr -d "\n" | wc -c`
  length2=`echo $exon_read_counts | tr -d "\n" | wc -c`
  if [ $length1 -eq $length2 ] ; then
    cap_code=`mark_depth "$exon_code" "$exon_read_counts" | sed "s/n/-/g"`;
    echo -e "$name\n$cap_code" > $current_fasta;
  else 
    echo -e "\e[1A" "ERROR in $name: $length1 != $length2\n"
fi
else
echo -e "\e[1A" "ERROR: No sequence found for exon $exon_name!\n cap_code=<$cap_code>\n exon_code = <$exon_code>\n exon_read_counts=<$exon_read_counts>\n" >> $FOLDER/cap_consensus_log;
fi
'
cat $FOLDER/cap_consensus/exons/*.fasta > $outfile
echo -e "$outfile took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $outfile already exists, skipping.
fi


###########################################
###################
#  Make multifasta. 
infile=$FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam
ref_file=$INTERMEDIATES/HtExons.fasta
ref_file2=$INTERMEDIATES/probe_to_exon_to_gene.list
outfile=${FOLDER}/${FILENAME}_exon_consensus.fasta
if [ $SKIP -eq 0 ] || [ ! -d $(echo $FOLDER/consensus) ]
then
echo -e "making multifastas"
die_unless $infile 
die_unless $PIPELINE/cigar_parser.py
START_TIME=`date +%s`
mkdir_if_needed $FOLDER/multifastas
mkdir_if_needed $FOLDER/consensus
mkdir_if_needed Consensus
#echo -e "" > $FILTER_OUT
ADDRESS=`pwd`/$FOLDER 
AVG_LEN=`$samtools view -F 0x4 $infile | cut -f 10 | awk '{EXON_READ_NR+=1; TOTAL_LENGTH+=length($1)}END{print int(TOTAL_LENGTH/EXON_READ_NR+0.5) }'`
echo -n "" > $FOLDER/$FILTER_OUT
echo -n "" > ${FOLDER}/cigar_parser.log
reset_file $FOLDER/temp.txt
EXON_CONTIGS=`cat $ref_file2 | cut -f3 | sort | uniq | wc -l` #Exons actually targeted by probes

#Make variables and 'return_unless_DNA_match()' available for xargs to use
export -f return_unless_DNA_match
export -f time_since
export -f subsample_bam
export samtools
export ADDRESS
export FOLDER
export FILENAME
export PIPELINE
export INTERMEDIATES
export BASE_FILES
export GENOME_FILE
export AVG_LEN
export FILTERED
export FILTER_OUT
export EXON_CONTIGS
#Iterate over names
echo "Starting..."
cat $ref_file | grep -o "Ht[^\n]*" | \
xargs -I{} --max-procs $CORES bash -c 'line="{}"; 
name=`echo "$line" | sed "s/\s//"`;
name_stripped=`echo "$name" | sed "s/(.*//"`;
scaffold=`echo "$name" | sed "s/:.*//"`;
#echo starting $FOLDER/$name;
#echo line = "$line";
START=`date +%s.%3N`;
#echo "start = $START";
nr_of_reads=`$samtools view $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $name_stripped | wc -l`;
#echo "nr of reads: $nr_of_reads";
if [[ nr_of_reads -gt 1000 ]] ; then #If there are more than 1000 reads, subsample. 
#echo "subsampling";
seq_len=`cat $INTERMEDIATES/HtExons.fasta | grep $name_stripped -A1 | tail -n1 | tr -d "\n" | wc -c`;
out=`subsample_bam $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $name_stripped | sort -k4n | $PIPELINE/cigar_parser.py $GENOME_FILE $name --speedy`;
else
out=`$samtools view $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $name_stripped | return_unless_DNA_match "40" | $PIPELINE/cigar_parser.py $GENOME_FILE $name --speedy`;
fi;
#echo "out = $out";
COUNTER=`cat "$FOLDER/temp.txt" | tr -d "\n" | wc -c`;
if [ ${#out} -gt 10 ]; then
time_nr=`echo "$(date +%s.%3N) - $START" | bc`;
#MAKE MULTIFASTAS
echo "$out" > $FOLDER/multifastas/$name.fasta;
#MAKE CONSENSUS
echo ">$name" >> $FOLDER/consensus/$name.fasta;
echo "$out" | tail -n 6 | grep Cut_consensus -A1 | tail -n1 >> $FOLDER/consensus/$name.fasta;
#echo "name=<$name>, sequence=<$sequence> "
size_nr=`cat $FOLDER/multifastas/$name.fasta | wc -c`;
line_nr=`$samtools view $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $name_stripped | wc -l`;
echo -e "$COUNTER\t$NUMBER\t$name\t$line_nr\t$size_nr\t$time_nr" >> ${FOLDER}/cigar_parser.log;
fi;
echo -n "." >> "$FOLDER/temp.txt";
if [ $[ COUNTER % 10 ] == 0 ] ; then 
echo -e "\e[1A" "$[ 100 * COUNTER / EXON_CONTIGS ] %    ($COUNTER of $EXON_CONTIGS)      ";
fi
exit 1;'
echo -en "\e[1A"; echo "Done!                                "
cat ${FOLDER}/consensus/*.fasta > $outfile
echo -e "Total reads processed: TOTAL_READS=$[(`cat ${FILENAME}_R[12]_001.fastq | wc -l`)/4]" >>$LOG 
echo -e "Total reads filtered out: $FILTERED" >>$LOG
echo -e "Multifastas/consensus took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $FOLDER/consensus already exists, skipping.
fi


#######################################################################
###########################################
#  Put consensus into consensus folder.

current_file=Consensus/cigar_parser/added_samples.list
if [ $SKIP -eq 0 ] || [ ! -f $(echo $current_file) ] || [ -z $(grep "$FOLDER" "$current_file") ]

then
echo -e "compiling scaffold multifastas"
START_TIME=`date +%s`
die_unless $FOLDER/${FILENAME}_consensus.fasta
mkdir_if_needed Consensus
mkdir_if_needed Consensus/cigar_parser
cp $FOLDER/${FILENAME}_consensus.fasta Consensus/${FILENAME}_consensus.fasta
cp $FOLDER/${FILENAME}_exon_consensus2.fasta Consensus/${FILENAME}_exon_consensus.fasta

sample=$FOLDER
cat $FOLDER/${FILENAME}_exon_consensus.fasta | paste - - | while read line ; do \
read exon sequence <<< `echo -e $line | sed "s/>//" | sed "s/([-+])//g"`
#echo "exon: $exon, sample: $sample, sequence: $sequence"
gene_ID=`get_gene_ID_from_exon $exon`
seq_dashes=`echo $sequence | tr 'n' '-'`
[[ ! $sequence =~ [ATGCatgc]+ ]] || echo -e ">${sample}_$exon\n$seq_dashes" >> Consensus/cigar_parser/${gene_ID}_$exon.fasta
done 

echo $sample >> $current_file #Mark that the job is done.
echo -e "multispecies multifastas took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo Consensus/${FILENAME}_consensus.fasta already exists, skipping.
fi


#######################################################################
###########################################
#  Compile multi-species exon multifastas.
current_file=Consensus/exons/added_samples.list
if [ $SKIP -eq 0 ] || [ ! -f $(echo $current_file) ] || [ -z $(grep "$FOLDER" "$current_file") ]

then
echo -e "generating multi-species exon multifastas."
START_TIME=`date +%s`
die_unless $FOLDER/${FILENAME}_exon_consensus2.fasta
mkdir_if_needed Consensus/exons 

sample=$FOLDER
cat $FOLDER/${FILENAME}_exon_consensus2.fasta | paste - - | while read line; do \
read exon sequence <<< `echo -e $line | sed "s/>//" | sed "s/([-+])//g"` 
gene_ID=`get_gene_ID_from_exon $exon`
seq_dashes=`echo $sequence | tr 'n' '-'`
#echo "exon: $exon, sample: $sample, sequence: $sequence"
[[ ! $sequence =~ [ATGCatgc]+ ]] || echo -e ">${sample}_$exon\n$seq_dashes" >> Consensus/exons/${gene_ID}_$exon.fasta; 
done

echo "$sample" >> $current_file #Mark that the job is done.
echo -e "multi-species exon multifastas took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $current_file already exists, skipping.
fi


#######################################################################
###########################################
#  Compile multi-species exon multifastas with capitalization, for comparison.
current_file=Consensus/cap_consensus/added_samples.list
if [ $SKIP -eq 0 ] || [ ! -f $(echo $current_file) ] || [ -z $(grep "$FOLDER" "$current_file") ]
then
echo -e "generating capitalized multi-species exon multifastas.\n"
START_TIME=`date +%s`
die_unless $FOLDER/${FILENAME}_exon_consensus3.fasta
mkdir_if_needed Consensus/cap_consensus
sample=$FOLDER

cat $FOLDER/${FILENAME}_exon_consensus3.fasta | paste - - | while read line 
do 
read exon sequence <<< `echo -e $line | sed "s/>//" | sed "s/([-+])//g"`; #Remove '>', '(+)' and '(-)'
echo -en "\e[1A"; echo "Processing $exon                                              "
gene_ID=`get_gene_ID_from_exon $exon`
#echo "exon: $exon, sample: $sample, sequence: $sequence"
if [[ $sequence =~ [ATGCatgc]+ ]] 
then 
seq_dashes=`echo $sequence | tr 'n' '-'`
echo -e ">${sample}_$exon\n$seq_dashes" >> "Consensus/cap_consensus/${gene_ID}_$exon.fasta"; 
fi
done

echo $sample >> $current_file #Mark that the job is done.
echo -e "multi-species exon capitalized multifastas took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $current_file already exists, skipping.
fi


      ###############################
##########################################
##### ALIGN MULTI-SPECIES FASTA FILES ####
##########################################
      ###############################


#######################################################################
###########################################
#  Align cigar_parser-derived multi-species exon multifastas
current_file=Consensus/cigar_parser/aligned_samples.list
previous_file=Consensus/cigar_parser/added_samples.list
if [ $LAST_ITERATION -eq 1 ] && ! cmp --silent $current_file $previous_file
then
echo -e "Aligning cigar-parser-derived exon consensus files\n"
START_TIME=`date +%s`

ls Consensus/cigar_parser | grep HtScaffold | grep -v aligned | while read name
do 
  new_name=`echo $name | sed "s/.fasta//"`_aligned.fasta 
  echo >&2 -en "\e[1A"; 
  echo >&2 "Aligning $new_name             "
  mafft --preservecase --thread $CORES --quiet "Consensus/cigar_parser/$name" | tr "\n" "\t" | sed "s/\t>/\n>/g" | sed "s/\t/\n/" | tr -d "\t" > "Consensus/cigar_parser/$new_name"
done

cat $previous_file > $current_file #Mark that the job is done.
echo -e "cigar_parser-derived multi-species exon multifastas took" `time_since $START_TIME` "to align.\n" >> $LOG
else
echo Skipping sequence alignment.
fi 

#######################################################################
###########################################
#  Align multi-species exon multifastas
current_file=Consensus/exons/aligned_samples.list
previous_file=Consensus/exons/added_samples.list
if [ $LAST_ITERATION -eq 1 ] && ! cmp --silent $current_file $previous_file
then
echo -e "Aligning exon consensus files\n"
START_TIME=`date +%s`

ls Consensus/exons | grep HtScaffold | grep -v aligned | while read name
do 
  new_name=`echo $name | sed "s/.fasta//"`_aligned.fasta 
  echo >&2 -en "\e[1A"; 
  echo >&2 "Aligning $new_name             "
  mafft --preservecase --thread $CORES --quiet "Consensus/exons/$name" | tr "\n" "\t" | sed "s/\t>/\n>/g" | sed "s/\t/\n/" | tr -d "\t" > "Consensus/exons/$new_name"
done

cat $previous_file > $current_file #Mark that the job is done.
echo -e "multi-species exon multifastas took" `time_since $START_TIME` "to align.\n" >> $LOG
else
echo Skipping sequence alignment.
fi 

#######################################################################
###########################################
#  Align capitalized multi-species exon multifastas
current_file=Consensus/cap_consensus/aligned_samples.list
previous_file=Consensus/cap_consensus/added_samples.list
if [ $LAST_ITERATION -eq 1 ] && ! cmp --silent $current_file $previous_file
then
echo -e "Aligning capitalized exon consensus files\n"
START_TIME=`date +%s`

ls Consensus/cap_consensus | grep HtScaffold | grep -v aligned | while read name
do 
  new_name=`echo $name | sed "s/.fasta//"`_aligned.fasta 
  echo >&2 -en "\e[1A"; 
  echo >&2 "Aligning $new_name             "
#  {
  mafft --preservecase --thread $CORES --quiet --retree 2 --maxiterate 0 "Consensus/cap_consensus/$name" | tr "\n" "\t" | sed "s/\t>/\n>/g" | sed "s/\t/\n/" | tr -d "\t" > "Consensus/cap_consensus/$new_name"
done

cat $previous_file > $current_file #Mark that the job is done.
echo -e "capitalized multi-species exon multifastas took" `time_since $START_TIME` "to align.\n" >> $LOG
else
echo Skipping sequence alignment.
fi 

         ############
  ##########################
#####                    #####
###       STATISTICS       ##################################################
#####                    #####
  ##########################
         ############



             #############
#########################################
##### LENGTH AND IDENTITY STATISTICS ####
#########################################
             #############


#######################################################################
######################################################
### Make length/identity statistics for each gene. 

infile=$FOLDER/${FILENAME}_exon_consensus3.fasta
ref_file=$INTERMEDIATES/exon_to_gene.list
outfile_gene1=$FOLDER/${FILENAME}_gene_length_and_id.tsv
outfile_gene2=$FOLDER/${FILENAME}_gene_length_and_id_test.tsv
outfile_gene3=Consensus/${FOLDER}_${FILENAME}_gene_length_and_id.tsv
outfile_exon1=$FOLDER/${FILENAME}_exon_length_and_id.tsv
outfile_exon2=$FOLDER/${FILENAME}_exon_length_and_id_test.tsv
outfile_exon3=$FOLDER/${FILENAME}_exon_length_and_id_test_sorted.tsv
if [ $SKIP -eq 0 ] || [ ! -f $(echo $outfile_gene1) ] || [ ! -f $(echo $outfile_gene2) ] || [ ! -f $(echo $outfile_exon1) ] || [ ! -f $(echo $outfile_exon2) ]
then
echo -e "\ngenerating $outfile_gene1, $outfile_exon1" 
START_TIME=`date +%s`
die_unless $infile
reset_file $outfile_gene1
reset_file $outfile_gene2
reset_file $outfile_exon1
reset_file $outfile_exon2
echo "" > $FOLDER/temp4

echo >&2 "Making Hash" 

#Make Gene-vs-Exon hash 
declare -A MOIN_REVISIONS exon_gene #Allow strings with parentheses etc to be keys.
declare -A MOIN_REVISIONS exon_length  
while read exon gene length 
do
    exon_gene[$exon]=$gene
    exon_length[$exon]=$length
done < $ref_file

echo >&2 "Looping over exons, making alignments" 
echo >&2 "Starting..."
#For each exon with reads:
cat $infile | paste - - | tr -d ">" | while read name code ; 
do 
echo -en >&2 "\e[1A Processing <${name}>                  \n";\
cat "$INTERMEDIATES/HtExons.fasta" | grep ">${name}" -A1 > $FOLDER/temp2
coverage=`echo $code | tr -cd "AGCTagct" | wc -c`
certain_coverage=`echo $code | tr -cd "AGCT" | wc -c`
certain_code=`echo $code | tr -c "AGCT" "N"`
echo -e ">$name\n$certain_code" > $FOLDER/temp
{
stretcher $FOLDER/temp2 $FOLDER/temp $FOLDER/temp3 #Needleman-Wunsch rapid global alignment of two sequences
} > /dev/null 2> /dev/null #Silence this script, so it doesn't report every time it's run.
alignment=`cat $FOLDER/temp3 | grep "Identity" | grep -oP "\d.*"`
gene=${exon_gene[$name]-"undefined"}
length=${exon_length[$name]-0}
if [ "$gene" == "undefined" ] 
then 
  length=`echo $code | tr -d "\n" | wc -c`
  echo >&2 "Exon $name has no connection to a gene; it may lack a Gene ID in the GFF file."
fi
#echo -e "$name\t$gene\t$certain_coverage\t$matches\t$identity\t$percent%"
matches=`echo "$alignment" | cut -d '/' -f1` 
misses=`bc <<< "$certain_coverage - $matches"` 
if [ $certain_coverage -eq 0 ] 
then 
  certain_identity=0
else 
  certain_identity=`bc <<< "scale=1; 100*$matches/$certain_coverage"` 	
fi
if [ $length -eq 0 ] 
then 
  percent_exon=0
  coverage_percent=0
else 
  percent_exon=`bc <<< "scale=1; 100*$coverage/$length"` 	
  coverage_percent=`bc <<< "scale=1; 100*$certain_coverage/$length"` 	
fi
 
echo -e "$name\t$gene\t$coverage\t$matches\t$alignment\t$percent_exon%" >> $outfile_exon1
echo -e "$name\t$gene\t$certain_coverage\t$matches\t$certain_identity\t$coverage_percent\t$alignment\t$misses" >> $outfile_exon2
cat $FOLDER/temp3 >> $FOLDER/temp4 
done

add_leading_zeroes $outfile_exon2 > $outfile_exon3 #Making exon stat file with leading zeroes. 

echo >&2 "Compiling $outfile_gene1 from exon data" 

#Compile exon data into gene data 
gawk '{a[$2]+=$3; b[$2]+=$4}END{for (i in a){if(a[i]==0) c=0; else c=100*b[i]/a[i]; printf "%s\t%s\t%3.1f\n", i, a[i], c;}}' $outfile_exon1 | sort -k1 > $outfile_gene1 

echo >&2 "Compiling $outfile_gene2" 

gawk '{a[$2]+=$3; b[$2]+=$4}END{for (i in a){if(a[i]==0) c=0; else c=100*b[i]/a[i]; printf "%s\t%s\t%3.1f\n", i, a[i], c;}}' $outfile_exon2 | sort -k1 > $outfile_gene2

cp $outfile_gene1 $outfile_gene3

echo made $outfile 
echo -e "$outfile took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $outfile already exists, skipping. 
[[ ! -f $(echo $outfile_gene3) ]] && cp $outfile_gene1 $outfile_gene3

fi


########################################################################
################################################################
### Record Probe data

infile=${FOLDER}/${FILENAME}_double_single_end_genome_sorted_rg.bam
outfile=${FOLDER}/${FILENAME}_probe_stats.tsv
outfile2=Consensus/${FOLDER}_${FILENAME}_probe_stats.tsv
ref_file=$PROBE_FILE
if [ -f $(echo $ref_file) ]
then
if [ $SKIP -eq 0 ] || [ ! -f $(echo $outfile) ]
then
echo -e "generating $outfile."
START_TIME=`date +%s`
die_unless $infile
die_unless $ref_file
if grep -q "Start" "$ref_file" #If probe file format uses Start and Stop positions
then
echo -e "probe\tGC\t${FILENAME}_reads\t${FILENAME}_coverage" > $outfile
grep -v "^#" "$ref_file" | grep . | tail -n +2 | while read scaffold probe seq repl strand chr start stop  ;
do 
    #echo "scaffold: $scaffold, probe: $probe, seq: $seq, repl: $repl, strand: $strand, chr: $chr, start: $start, stop: $stop"; 
    stop=$[stop-1]; #fix fencepost error
    pos="${chr}:${start}-${stop}" ; 
    echo -en "\e[1A Processing <${pos}>\n                      "
    length=`echo $seq | tr -cd 'AGCT' | wc -c`;
    GC=`echo $seq | tr -cd 'GC' | wc -c` ; 
    gc_content=`bc <<< "scale=1; 100*$GC/$length"` ; 
    echo -ne "$pos\t$gc_content\t" >> $outfile ; 
    reads=`$samtools view "$infile" "$pos" -c` ;
    echo -ne "$reads\t" >> $outfile ;
    if [ $reads -eq 0 ] 
    then
	echo "0" >> $outfile
    else
	$samtools depth "$infile" -r "$pos" | awk -v l=$length '{if($3>2) mapped+=1; else unmapped+=1}END{if (mapped>0) print (100*mapped/l); else print 0}' >> $outfile;
	fi
done
elif grep -q "Coordinates" "$ref_file"
then  
echo -e "probe\tGC\t${FILENAME}_reads\t${FILENAME}_coverage" > $outfile
tail -n +2 $ref_file | while read scaffold probe seq repl strand position ;
do 
	pos=`echo $position | sed "s/^chr//"` ; 
	echo -en "\e[1A Processing <${pos}>\n                      "
	length=`echo $seq | tr -cd 'AGCT' | wc -c`;
	GC=`echo $seq | tr -cd 'GC' | wc -c` ; 
	gc_content=`bc <<< "scale=1; 100*$GC/$length"` ; 
	echo -ne "$pos\t$gc_content\t" >> $outfile ; 
	reads=`$samtools view "$infile" "$pos" -c` ;
	echo -ne "$reads\t" >> $outfile ;
	if [ $reads -eq 0 ] 
	then
		echo "0" >> $outfile
	else
	$samtools depth "$infile" -r "$pos" | awk -v l=$length '{if($3>2) mapped+=1; else unmapped+=1}END{if (mapped>0) print (100*mapped/l); else print 0}' >> $outfile;
	fi
done
else echo "Could not find compatible algorithm for $ref_file!"
fi

cp $outfile $outfile2

echo made $outfile 
echo -e "$outfile took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $outfile already exists, skipping. 
[[ ! -f $(echo $outfile2) ]] && cp $outfile $outfile2
fi
else
echo "Skipping $outfile - probe file '$ref_file' not present"
fi

#######################################################################
###############################################################
### Compile gene length/identity statistics in consensus folder. 
infile=$INTERMEDIATES/gc_content_sorted.tsv
infile2=$FOLDER/${FILENAME}_gene_length_and_id_test.tsv
outfile=Consensus/${FOLDER}_${FILENAME}_gene_length.tsv
outfile2=Consensus/${FOLDER}_${FILENAME}_gene_identity.tsv
outfile3=Consensus/${FOLDER}_${FILENAME}_gene_coverage.tsv
if [ $SKIP -eq 0 ] || [ ! -f $(echo $outfile) ] 
then

START_TIME=`date +%s`
die_unless $FOLDER/${FILENAME}_gene_length_and_id.tsv
die_unless $FOLDER/${FILENAME}_gene_length_and_id_test.tsv

lineage="" # Will be filled in by hand; researcher will know better than BLAST search, especially with novel strains not present in DB.

echo -e "generating $outfile."
#Add length column to length file in consensus folder
echo -e "$lineage\n${FILENAME}_length" > $outfile
join -1 1 -2 1 -a1 -e'0' -o '2.2' $infile $infile2 >> $outfile

echo -e "generating $outfile2."
#Add identity column as identity file in consensus folder
echo -e "$lineage\n${FILENAME}_identity" > $outfile2
join -1 1 -2 1 -a1 -e'0' -o '2.3' $infile $infile2  >> $outfile2

echo -e "generating $outfile3."
#Calculate coverage in percent and add file to consensus folder
echo -e "$lineage\n${FILENAME}_coverage" > $outfile3
join -1 1 -2 1 -a1 -e'0' -o '1.2,2.2' $infile $infile2 | gawk '{cov=100*$2/$1; printf "%3.1f\n", cov;}' >> $outfile3
echo -e "$outfile took" `time_since $START_TIME` "to make.\n" >> $LOG
else

echo $outfile, $outfile2 and $outfile3 already exist, skipping.
fi

#######################################################################
###############################################################
### Make Consensus summary of gene stats. 

outfile="Consensus/gene_length_and_id.tsv"
#Join GC/length file, coverage file and identity file into one new file.  
paste $INTERMEDIATES/gc_content.tsv Consensus/*_gene_length.tsv Consensus/*_gene_identity.tsv Consensus/*_gene_coverage.tsv > $outfile


#######################################################################
###############################################################
### Compile exon length/identity statistics in consensus folder. 

infile=$FOLDER/${FILENAME}_exon_length_and_id_test_sorted.tsv
infile2=$INTERMEDIATES/exons_gc_count_sorted.list
outfile=Consensus/${FOLDER}_${FILENAME}_exon_identity.tsv
outfile2=Consensus/${FOLDER}_${FILENAME}_exon_coverage.tsv
outfile3=Consensus/${FOLDER}_${FILENAME}_exon_alignment.tsv
outfile4=Consensus/${FOLDER}_${FILENAME}_exon_mismatches.tsv
if [ $SKIP -eq 0 ] || [ ! -f $(echo $outfile) ] || [ ! -f $(echo $outfile2) ] || [ ! -f $(echo $outfile3) ] || [ ! -f $(echo $outfile4) ] 
then

START_TIME=`date +%s`
die_unless $infile
die_unless $infile2

lineage="" # Will be filled in by hand; researcher will know better than BLAST search, especially with novel strains not present in DB.

echo -e "generating $outfile."
#Add identity column to identity file in consensus folder
echo -e "$lineage\n${FILENAME}_identity" > $outfile
join -1 1 -2 1 -a1 -e'0' -o '2.5' -t $'\t' $infile2 $infile >> $outfile

echo -e "generating $outfile2."
#Add coverage column to coverage file in consensus folder
echo -e "$lineage\n${FILENAME}_coverage" > $outfile2
join -1 1 -2 1 -a1 -e'0' -o '2.6' -t $'\t' $infile2 $infile >> $outfile2

echo -e "generating $outfile3."
#Add alignment column to alignment file in consensus folder
echo -e "$lineage\n${FILENAME}_alignment" > $outfile3
join -1 1 -2 1 -a1 -e'0' -o '2.7' -t $'\t' $infile2 $infile >> $outfile3

echo -e "generating $outfile4."
#Add mismatches column to mismatches file in consensus folder
echo -e "$lineage\n${FILENAME}_mismatches" > $outfile4
join -1 1 -2 1 -a1 -e'0' -o '2.8' -t $'\t' $infile2 $infile >> $outfile4

echo -e "$outfile took" `time_since $START_TIME` "to make.\n" >> $LOG
else

echo $outfile, $outfile2, $outfile3 and $outfile4 already exists, skipping.
fi

###############################
########################
## Make exon GC file

infile=$INTERMEDIATES/HtExons.fasta
outfile=Consensus/exon_gc_content.tsv
outfile2=Consensus/exon_gc_content_labeled.tsv
if [ $SKIP -eq 0 ] || [ ! -f $(echo $outfile) ] || [ ! -f $(echo $outfile2) ]
then
echo -e "generating $outfile."
START_TIME=`date +%s`
die_unless $infile
reset_file $outfile

cat $infile | paste - - | tr -d ">" | while read name seq ; do 
  length=`echo $seq | tr -cd 'AGCT' | wc -c`;
  GC=`echo $seq | tr -cd 'GC' | wc -c` ; 
  gc_content=`bc <<< "scale=1; 100*$GC/$length"` 
  echo -e "$name\t$length\t$gc_content" >> $outfile
  done
echo -e "\nexon\tlength\tgc%" > $outfile2
cat $outfile >> $outfile2
echo -e "$outfile took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $outfile already exists, skipping.
fi

#######################################################################
###############################################################
### Make Consensus summary of exon stats. 

#Join exon GC/length file, coverage file and identity file into one new file.  
paste Consensus/exon_gc_content_labeled.tsv Consensus/*_exon_coverage.tsv Consensus/*_exon_identity.tsv > Consensus/exon_length_and_id.tsv



            ##############
######################################
##### SEQUENCE QUALITY STATISTICS ####
######################################
            ##############


################################################################
#################################################
### Make statistics for double single end. 

fastqc_file1=$FOLDER/fastq_stats/${FILENAME}_trimmed_single_end_min36_fastqc.html
fastqc_file2=$FOLDER/fastq_stats/${FILENAME}_trimmed_single_end_reverse_min36_fastqc.html
infile1=`ls ${FILENAME}_R1_001.fastq* | head -n1`
infile2=`ls ${FILENAME}_R2_001.fastq* | head -n1`
infile3=$FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam
infile4=$FOLDER/${FILENAME}_probe_stats.tsv
ref_file1="$GENOME_FILE"
ref_file2="$INTERMEDIATES/HtExons.fasta"
ref_file3="$INTERMEDIATES/probe_to_exon_to_gene.list"
ref_file4="$INTERMEDIATES/relevant_exons.fasta"
ref_file5="$INTERMEDIATES/relevant_exons.bed"
if [ $SKIP -eq 0 ] || [ ! -f $(echo $FOLDER/statistics_double_single_end.txt) ]
then
echo -e "\ncalculating basic statistics for $FILENAME, double single end \n"

die_unless $ref_file1
die_unless $ref_file2
die_unless $ref_file3
die_unless $ref_file4
die_unless $ref_file5
die_unless $FOLDER/$FILTER_OUT
die_unless $infile3
die_unless $infile4

STAT=$FOLDER/statistics_double_single_end.txt
touch $STAT
START_TIME=`date +%s`
echo -e "\nStatistics for $FILENAME, double single end: \n" > $STAT

GENOME_SIZE=`cat $ref_file1 | grep -v "^>" | tr -dc "ATGCatgc" | wc -c`
echo -e "Time: " $[`date +%s` - START_TIME] ": Size of reference genome (mapped): $GENOME_SIZE bases"
GENOME_CONTIGS=`cat $ref_file1 | grep "^>" | wc -l`
echo -e "Time: " $[`date +%s` - START_TIME] ": Number of contigs in reference counted: $GENOME_CONTIGS"
RELEVANT_EXONS_SIZE=`cat $ref_file4 | grep -v "^>" | tr -dc "ATGCatgc" | wc -c`
echo -e "Time: " $[`date +%s` - START_TIME] ": Total size of targeted exons: $RELEVANT_EXONS_SIZE bases"
EXON_MAX=`cat $ref_file2 | grep "^>" | wc -l` 
PROBE_NUMBER=`cat $ref_file3 | wc -l` 
EXON_CONTIGS=`cat $ref_file3 | cut -f3 | sort | uniq | wc -l` 
echo -e "Time: " $[`date +%s` - START_TIME] ": Number of targeted exons in reference counted: $EXON_CONTIGS out of $EXON_MAX"

if [[ "$infile1" == *.fastq ]] ; then
TOTAL_READS_FORWARD=$[`cat $infile1 | wc -l`/4]
else 
TOTAL_READS_FORWARD=$[`zcat $infile1 | wc -l`/4]
fi
echo -e "Time: " $[`date +%s` - START_TIME] ": Reads in forward file counted: $TOTAL_READS_FORWARD"
if [[ "$infile2" == *.fastq ]] ; then
TOTAL_READS_REVERSE=$[`cat $infile2 | wc -l`/4]
else 
TOTAL_READS_REVERSE=$[`zcat $infile2 | wc -l`/4]
fi
echo -e "Time: " $[`date +%s` - START_TIME] ": Reads in reverse file counted: $TOTAL_READS_REVERSE"
TOTAL_READS=$[$TOTAL_READS_FORWARD+$TOTAL_READS_REVERSE]
echo -e "Time: " $[`date +%s` - START_TIME] ": Total reads: $TOTAL_READS"


USED_READS_FORWARD=`grep -oP "<td>Total Sequences</td><td>\d+</td>" "$fastqc_file1" | grep -oP "\d+"`
echo -e "Time: " $[`date +%s` - START_TIME] ": Used reads in forward file counted: $USED_READS_FORWARD"
USED_READS_REVERSE=`grep -oP "<td>Total Sequences</td><td>\d+</td>" "$fastqc_file2" | grep -oP "\d+"`
echo -e "Time: " $[`date +%s` - START_TIME] ": Used reads in reverse file counted: $USED_READS_REVERSE"
READS_SINGLE=$[$USED_READS_FORWARD+$USED_READS_REVERSE]
echo -e "Time: " $[`date +%s` - START_TIME] ": Total used reads: $READS_SINGLE"

READS_FILTERED_OUT=`cat ${FOLDER}/$FILTER_OUT | wc -l`
cat ${FOLDER}/$FILTER_OUT | sort | tr "\t" "\n" | awk '!a[$0]++' > ${FOLDER}/removed_reads #Sort the (asynchronously generated) file under scaffold headers
echo -e "Time: " $[`date +%s` - START_TIME] ": Sorted removed reads: $READS_FILTERED_OUT"
TIME_READS=`date +%s`

MAPPED_READS_GENOME=`$samtools view -F 0x4 $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam | cut -f 1 | sort | uniq | wc -l `
TIME_MAPPED=`date +%s`
echo -e "Time: " $[`date +%s` - START_TIME] ": Mappings counted: $MAPPED_READS_GENOME"

read EXONS_HIT EXONS_HIT_TEN EXONS_HIT_HUNDRED<<< $(for line in `cat $ref_file3 | cut -f3 | uniq | cut -f1 -d "(" `; do $samtools view $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $line | wc -l; done | awk '{if($0>0) nr1+=1; if($0>9) nr10+=1; if($0>99) nr100+=1}END{print nr1" "nr10" "nr100}')

echo -e "Time: " $[`date +%s` - START_TIME] ": Hits in targeted exons counted. $EXONS_HIT exons with at least one hit, $EXONS_HIT_TEN with at least ten, $EXONS_HIT_HUNDRED with at least 100."

GENOME_BASES_HIT=0
GENOME_BASES_MISSED=0

read GENOME_BASES_HIT GENOME_BASES_MISSED <<< $($samtools depth -aa $infile3 | awk '{if($3>2) total+=1; else unmapped+=1}END{print total" "unmapped}')
GENOME_BASES_MISSED=${GENOME_BASES_MISSED-0}
read EXON_BASES_HIT EXON_BASES_MISSED <<< $($samtools depth -aa $infile3 -b $ref_file5 | awk '{if($3>2) total+=1; else unmapped+=1}END{print total" "unmapped}')
read PROBES_HIT PROBES_HIT_WELL <<< $(cat $infile4 | awk  '{if($4>0) hit+=1; if($4>95) hit_well+=1; }END{print hit" "hit_well}')

END_TIME=`date +%s`
echo -e "Time: " $[`date +%s` - START_TIME] ": Base hits counted: $GENOME_BASES_HIT"
RUNTIME=$(($END_TIME-$START_TIME))

echo -e "Number of reads, including bad ones: $TOTAL_READS.  \n" >> $STAT
echo -e "Number of used reads: $READS_SINGLE; $[100*$READS_SINGLE/$TOTAL_READS]% \n" >> $STAT
echo -e "Number of reference-identical reads filtered out: $READS_FILTERED_OUT; `echo "scale=3; 100*$READS_FILTERED_OUT/$READS_SINGLE" | bc`% of used reads\n" >> $STAT
echo -e "Reads matching genome: $MAPPED_READS_GENOME; `echo "scale=1; 100*$MAPPED_READS_GENOME/$READS_SINGLE" | bc`% of used reads \n" >> $STAT
echo -e "Number of genome contigs: $GENOME_CONTIGS \n" >> $STAT
echo -e "Number of exons in reference files: $EXON_MAX \n" >> $STAT
echo -e "Number of probes in probe file: $PROBE_NUMBER \n" >> $STAT
echo -e "Number of probes with coverage above 0% : $PROBES_HIT; `echo "scale=1; 100*$PROBES_HIT/$PROBE_NUMBER" | bc`% of probes \n" >> $STAT
echo -e "Number of probes with coverage above 95%: $PROBES_HIT_WELL; `echo "scale=1; 100*$PROBES_HIT_WELL/$PROBE_NUMBER" | bc`% of probes \n" >> $STAT
echo -e "Number of exons targeted by probes: $EXON_CONTIGS; `echo "scale=1; 100*$EXON_CONTIGS/$EXON_MAX" | bc`%  \n" >> $STAT
echo -e "Total bases in exons targeted by probes: $RELEVANT_EXONS_SIZE; `echo "scale=1; 100*$RELEVANT_EXONS_SIZE/$GENOME_SIZE" | bc`% of genome \n" >> $STAT
echo -e "Targeted exons with at least one match: $EXONS_HIT; `echo "scale=1; 100*$EXONS_HIT/$EXON_CONTIGS" | bc`%  \n" >> $STAT 
echo -e "Targeted exons with at least ten matches: $EXONS_HIT_TEN; `echo "scale=1; 100*$EXONS_HIT_TEN/$EXON_CONTIGS" | bc`%  \n" >> $STAT
echo -e "Targeted exons with at least a hundred matches: $EXONS_HIT_HUNDRED; `echo "scale=1; 100*$EXONS_HIT_HUNDRED/$EXON_CONTIGS" | bc`%  \n" >> $STAT 

echo -e "Coverage breadth: \n
Genome bases covered by at least 3 reads: $GENOME_BASES_HIT; `echo "scale=1; 100*$GENOME_BASES_HIT/($GENOME_BASES_HIT+$GENOME_BASES_MISSED)" | bc`%%  \n 
Genome bases missed: $GENOME_BASES_MISSED bases; `echo "scale=1; 100*$GENOME_BASES_MISSED/($GENOME_BASES_HIT+$GENOME_BASES_MISSED)" | bc`%%  \n
Targeted exon bases covered by at least 3 reads: $EXON_BASES_HIT; `echo "scale=1; 100*$EXON_BASES_HIT/$RELEVANT_EXONS_SIZE" | bc`%%  \n " >> $STAT

printf "AGCT bases counted: $GENOME_SIZE\nAll bases counted: $[$GENOME_BASES_HIT+$GENOME_BASES_MISSED]" >> $STAT

parse_fastq_screen >> $STAT

echo -e "Total Script runtime:" `time_since $SCRIPT_START_TIME` >> $STAT

echo -e "Calculating statistics took" `time_since START_TIME` "to do." >>$LOG 
echo -e "Program as a whole took:" `time_since $SCRIPT_START_TIME` >>$LOG
echo -e "Statistics recorded in file."
cat $STAT
echo -e "These statistics took $RUNTIME seconds to calculate. Counting reads: $[$TIME_READS-$START_TIME], counting mapped reads: $[$TIME_MAPPED-$TIME_READS], counting exons hit: $[$END_TIME-$TIME_MAPPED]\n\n"
else
echo statistics_double_single_end.txt already exists, skipping.
fi


#######################################################################
###########################################
#  Compile match statistics

outfile=$FOLDER/match_stats.txt
if [ $LAST_ITERATION -eq 1 ]
then
echo -e "Compiling statistics for all folders"
START_TIME=`date +%s`

echo "Statistics compiled at `date`" > $outfile
ls */stat* | while read name 
do 
  filename=`echo "$name" | cut -d '/' -f1`
  total=`grep "Number of reads" $name | grep -oP "[0-9]*"`
  matching=`grep "Reads matching genome" $name | grep -oP "[0-9\.%]*" | tr "\n" " "`
  echo -e "$filename $total $matching" >> $outfile
done

echo -e "$outfile took" `time_since $START_TIME` "to calculate.\n" >> $LOG
else
echo Skipping sequence alignment.
fi 

#######################################################################
###########################################
#  Compile probe statistics
infile=$FOLDER/${FILENAME}_probe_stats.tsv
if [ $LAST_ITERATION -eq 1 ]
then
echo -e "Compiling probe statistics for all folders"
die_unless $infile
START_TIME=`date +%s`

cat $infile | cut -f1,2 > Consensus/probes_stats.tsv
ls Consensus/*probe_stats.tsv | while read line ; 
do 
  paste Consensus/probes_stats.tsv <( cat $line | cut -f4 ) > temp ; 
  cat temp > Consensus/probes_stats.tsv ; 
done


echo -e "Compiling probe statistics took" `time_since $START_TIME` "to align.\n" >> $LOG
else
echo Skipping sequence alignment.
fi 

##################################
### GO ON TO NEXT FILE IN LIST ###
##################################

move_on_to_next_file
