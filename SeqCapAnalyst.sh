#!/bin/bash
MINARGS=1

# Use Trimmomatic, bowtie2 on a pair of files: 

### REQUIREMENTS: 
# Requires input fastq files to end in "_R[12]_001.fastq"

# Picard (Calls on "java -jar bin/picard.jar")
# Bowtie2 (Calls on "bowtie2", "bowtie2-build")
# Trimmomatic (calls "java -jar ~/bin/trimmomatic-0.36.jar") <-- Will need fix
# SAMtools (calls on "samtools")
# bcftools
# vcftools
# bedtools
# vcfutils.pl (usually included in above)
# seqtk 
# Fastqc (calls on "fastqc")
# cigar_parser.py (included) 

#########################################################

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
  [ ! -f $(echo $FOLDER/${FILENAME}_double_single_end_genome.bam) ] || rm $FOLDER/${FILENAME}_double_single_end_genome.bam 
  [ ! -f $(echo $FOLDER/${FILENAME}_double_single_end_genome_sorted.bam) ] || rm $FOLDER/${FILENAME}_double_single_end_genome_sorted.bam
fi

[ ! -f $(echo $FOLDER/${FILENAME}_*.sam) ] || rm $FOLDER/$FILTER_OUT || rm $FOLDER/$FILTER_OUT

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

READ_LIMIT=5
function mark_depth () 
{
j=-1
length1=`echo $1 | tr -d "\n" | wc -c`
length2=`echo $2 | tr -d "\n" | wc -c`
if [ $length1 -gt $length2 ] 
then 
  echo -e >&2 "Code string and read_nr string not equally long!\n '$1'($length1)\n'$2'($length2)"
  exit 0
fi

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

FILTER_OUT=""
FILTERED=0
function return_unless_DNA_match ()
{
  #sequence=$1 
  avg_len=$1
  while read line; do 
    query=`echo "$line" | cut -f10`
    read=`echo "$line" | cut -f1`
    name=`echo "$line" | cut -f3`
    length=`echo "$query" | wc -c`
    MD=`echo "$line" | cut -f12 | cut -c 6-`
    #echo -e >&2 "Seq: '$sequence'"
    #echo -e >&2 "avg_len: '$avg_len'"
    #echo -e >&2 "line: '$line'"
    #echo -e >&2 "query: '$query'"
    #echo -e >&2 "length: '$length'"
    #echo -e >&2 "MD: '$MD'"
    if [ $length -lt $avg_len ] ; then echo "$line" #; echo -e >&2 "kept because too short"
    else
#    [[ "$sequence" =~ "$query" ]] || echo "$line"
    if [[ "$MD" =~ ^[0-9\s]+$ ]] ; then 
      echo -e ">$name\t$read ($length bases)" >> ${FOLDER}/$FILTER_OUT
      export FILTERED=$((FILTERED+1))  
    else 
    echo "$line"
    fi
    fi
  done
}

READ_MINIMUM=10 
MIN_LENGTH_TO_DISCARD=40
function subsample_bam () #Note: will add 'R' or 'F' to read-names in resulting SAM. 
{
  #Create variables
  samfile=$1
  exon=$2
  >&2 echo -e "\e[1A Subsampling ${exon}...     "
  [[ VERBOSE -eq 1 ]] && >&2 echo "samfile=$samfile, exon=$exon"
  declare -A MOIN_REVISIONS endpoints #Array for end points
  declare -A MOIN_REVISIONS lines #Array to store SAM lines
  scaffold=`echo "$exon" | sed "s/:.*//"`;
  interval=`echo "$exon" | grep -oP "\d+-\d+"`
  IFS='-'; read start end <<< "$interval"
  IFS=$' \t\n'
  #Run loop for positions in exon
  [[ VERBOSE -eq 1 ]] && >&2 echo "start=$start, end=$end, IFS='$IFS'"
  for i in $(seq $start $end) ; do 

    #Remove entries from array
    for x in ${!endpoints[@]} ; do
      if [ ${endpoints[$x]} -lt $i ] ; then 
        #print to SAM and remove from arrays 
        [[ VERBOSE -eq 1 ]] && >&2 echo "endpoint=${endpoints[$x]}, i=$i"
        echo "${lines[$x]}"
        [[ VERBOSE -eq 1 ]] && >&2 echo -n "i=$i, removing $x (ends at ${endpoints[$x]}), " 
        unset 'endpoints[$x]' 
        unset 'lines[$x]'
        [[ VERBOSE -eq 1 ]] && >&2 echo "coverage is now ${#endpoints[@]}"
      fi 
    done

    coverage=${#endpoints[@]} # nr of entries in array 
    #echo "coverage='$coverage'"
    [[ coverage -lt $READ_MINIMUM ]] || continue #If enough reads, continue. 

    difference=$[READ_MINIMUM - coverage] 
    #Otherwise, load [min_reads] more files, if available, and add until [min_reads] is reached. (Note: some may have been added previously) 
    #echo "start=$start, end=$end, i='$i', difference=$difference"
    [[ VERBOSE -eq 1 ]] && >&2 echo "searching for '$samfile' '${scaffold}:${i}-$end'"
    list=`samtools view "$samfile" "${scaffold}:${i}-$end" | return_unless_DNA_match "$MIN_LENGTH_TO_DISCARD" | head -n $[ 5 + READ_MINIMUM * 2 ]`
    [[ VERBOSE -eq 1 ]] && >&2 echo "list=$list"
    OLDIFS=$IFS
    IFS=$'\n'; 
    for line in $list ; do
      [[ VERBOSE -eq 1 ]] && >&2 echo "line='$line'"
      name=`echo -e "$line" | cut -f1` 
      flag=`echo -e "$line" | cut -f2` 
      direction=`((($flag&16)>0)) && echo 'R' || echo 'F'`
      name="${name}$direction"
      if test "${endpoints[${name}]+isset}" ; then continue 
      #else 
      #  echo "'$name' not in list. List:"
      #  for y in ${!endpoints[@]} ; do echo -n "$y, " ; done ; echo ""
      fi
      if test "${endpoints[${name}]+isset}" ; then continue ; fi
     # [[ -z "${endpoints[${name}]}" ]] && continue #If already in list, skip.
      start_pos=`echo -e "$line" | cut -f4`
      code=`echo -e "$line" | cut -f10`
      cigar=`echo -e "$line" | cut -f6`

      [[ VERBOSE -eq 1 ]] && >&2 echo "name=$name, start_pos=$start_pos, cigar=$cigar, code=$code"
      
      #Consider cigar string. Length is # of 'M' + # of 'D', ignoring # of 'I'.
      length=$[0 `echo $cigar | sed -r "s/([0-9]+)[MDS]/+\1 /g" | sed -r "s/([0-9]+)I//g"`] #too complex?
       
      #Add name and code to list
      endpoints[$name]=$[start_pos+length-1]
      #echo "name=$name, start_pos=$start_pos, length=$length, endpoint=${endpoints[$name]}, i=$i"
      lines[$name]=$line
      #echo "new line = '${lines[$name]}'"
      coverage=${#endpoints[@]}
      #>&2 echo "$i, added $name, coverage is now $coverage" 
      [[ $coverage -lt $READ_MINIMUM ]] || break
    done
    IFS=$OLDIFS
  done 
  
#Remove entries remaining in array
for x in ${!endpoints[@]} ; do
    #print to SAM and remove from arrays
    echo "${lines[$x]}" 
    #echo "unsetting $x"
    unset 'endpoints[$x]' 
    unset 'lines[$x]'
done
}

function SAM_to_fasta ()
{
while read line; do
  echo ">$line" | cut -f1
  echo "$line" | cut -f10
done
}

function get_gene_ID_from_exon ()
{
die_unless $PIPELINE/exon_to_gene.list
exon=`echo $1 | tr -d '[:space:]'`
result=`grep "$exon" "$PIPELINE/exon_to_gene.list" | cut -f2 | head -n1`
#echo >&2 "result is '$result'"
if [ ! -z $result ] 
then 
  echo "$result"; 
else
  echo >&2 "Could not find matching gene for '$exon' in '$PIPELINE/exon_to_gene.list'."
  echo "unknown"
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


##################
### CONSTANTS: ###
##################

### Calculate default number of cores
CORES=`grep -c ^processor /proc/cpuinfo` #Count number of cores.
CORES=$(($CORES - 1)) #Save one core for other tasks.
[ $CORES -gt 10 ] && CORES=10 #Tests show neglible impact on performance when using more than 6-10 cores (depending on file size). Thus, default value is capped at 10. 
[ $CORES -lt 1 ] && CORES=1 #In case pipeline is run on PC with one processor. 

SKIP=1 #Skip steps if resulting files already exist
KEEP=0 #Keep intermediary files
DELETE=0 #Delete everything but statistics afterwards
ARGUMENTS=$@
BASE=""
PIPELINE="$HOME/bin/SeqCapAnalyst" #folder of pipeline-related files
FOLDER=""
FILENAME=""

# Names of files
REFERENCE="$PIPELINE/probes.fasta"
ADAPTER="$PIPELINE/adapters/TruSeq3-PE-2.fa"
GFF_FILE="$PIPELINE/HtGenomeGFF_CDS1000.gff"
MALAVI_FILE="$PIPELINE/MalAvi_MS_name_20170326_224422.fasta"
FILTER_OUT="removed_reads_temp"
FASTQ_SCREEN_CONF="$PIPELINE/fastq_screen.conf"
FASTQ_SCREEN_CONF_BEGIN="$PIPELINE/fastq_screen_begin.conf"

# Numbers for statistics
EXON_CONTIGS=2569
GENOME_CONTIGS=2987
SCAFFOLD_CONTIGS=700
EXONS_SIZE=2293615
GENOME_SIZE=23218265
SCAFFOLDS_SIZE=12400015

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
-s 0 (default 1) : set pipeline to overwrite any already-existing files, rather than assume they are correct. 

INPUT: 
OUTPUT: 

"

while [[ $* ]]
do
    OPTIND=1
#    echo $1
    if [[ $1 =~ ^- ]]
    then
        getopts :f:s:n:a:b:d:th parameter
        case $parameter in
            h)  echo "$HELP"
		exit
                ;;
            f)  FOLDER=$OPTARG #Currently only works on single file. 
                echo "FOLDER=$FOLDER"
                shift
                ;;
            a)  ADAPTER=$OPTARG 
                echo "ADAPTER=$ADAPTER"
                shift
                ;;
            n)  CORES=$OPTARG 
                echo "CORES=$CORES"
                shift
                ;;            
	    b)  STRINGENCY=$OPTARG 
                echo "STRINGENCY=$STRINGENCY"
                shift
                ;;
            s)  SKIP=$OPTARG
                echo "SKIP=$SKIP"
                shift
                ;;
            t)  TEST=1
                echo "TEST=$TEST"
                shift
                ;;
            d)  DELETE=$OPTARG
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
FILENAME=$( echo "$first_word" | sed s/_R[12]_001.fastq$// )
FILENAME=$( echo "$FILENAME" | sed s/.fastq$// ) #If infiles don't have standard naming
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
ADAPTER = $ADAPTER\n
CORES = $CORES\n
SKIP = $SKIP\n
STRINGENCY = $STRINGENCY\n
LAST_ITERATION = $LAST_ITERATION\n"

echo -e $intro

### Make folder to put new files in. 
echo Creating folder: $FOLDER
mkdir_if_needed $FOLDER
STATS_FOLDER="$FOLDER/fastq_stats"

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



#########################################################################################
#### INDEX BUILDING IN PIPELINE FOLDER ####
#########################################################################################

echo -e "\nChecking prerequisite files"

### Make HtExons.fasta if not present
if [ $SKIP -eq 0 ] || [ ! -f $(echo $PIPELINE/HtExons.fasta) ]
then
die_unless $GFF_FILE
echo -e "\ngenerating HtExons.fasta"
bedtools getfasta -fi $PIPELINE/HtGenome.fasta -bed $GFF_FILE -fo $PIPELINE/HtExons.fasta -s

echo made HtExons.fasta from gff file
else
echo $PIPELINE/HtExons.fasta already exists, skipping.
fi


### Make HtExons.fasta.fai if not present
infile=$PIPELINE/HtExons.fasta
outfile=$PIPELINE/HtExons.fasta.fai
if [ $SKIP -eq 0 ] || [ ! -f $(echo $outfile) ]
then
die_unless $infile
echo -e "\ngenerating $outfile"

samtools faidx $infile

echo made $outfile
else
echo $outfile already exists, skipping.
fi


### Make HtGenome.fasta.fai if not present
infile=$PIPELINE/HtGenome.fasta
outfile=$PIPELINE/HtGenome.fasta.fai
if [ $SKIP -eq 0 ] || [ ! -f $(echo $outfile) ]
then
die_unless $infile
echo -e "\ngenerating $outfile"

samtools faidx $infile

echo made $outfile
else
echo $outfile already exists, skipping.
fi

### Make exon_scaffolds.fasta if not present
if [ $SKIP -eq 0 ] || [ ! -f $(echo $PIPELINE/exon_scaffolds.fasta) ]
then
die_unless $PIPELINE/HtExons.fasta
echo -e "\ngenerating exon_scaffolds.fasta"

touch $PIPELINE/exon_scaffolds.fasta 
cat $PIPELINE/HtExons.fasta | grep -o "Ht[^\:]*" | sort | uniq | while read line; do cat $PIPELINE/HtGenome.fasta | grep $line -A 1 >>  $PIPELINE/exon_scaffolds.fasta; done

echo made exon_scaffolds.fasta from gff file
else
echo $PIPELINE/exon_scaffolds.fasta already exists, skipping.
fi


### Make Gene-exon database if not present
current_file=$PIPELINE/gc_content.tsv
if [ $SKIP -eq 0 ] || [ ! -f $(echo $current_file) ]
then
echo -e "\ngenerating $current_file" 
die_unless $GFF_FILE
mkdir_if_needed Consensus

#Make list of exon name, gene id, and length
cat $GFF_FILE | gawk 'match($0, /gene_id "([^"]+)/, x) { print $1":"$4-1"-"$5"("$7")\t"x[1]"\t"$5-$4+1 }' | sort > $PIPELINE/exon_to_gene.list
#Add column with GC content 
cat $PIPELINE/HtExons.fasta | paste - - | while read line 
do read exon code <<< $line 
#echo "exon='$exon', code='$code'"
exon=`echo $exon | tr -d '>'`
echo -e "$exon\t`echo $code | tr -d -c 'CGcg' | wc -c`" >> $PIPELINE/exons_gc_count.list 
done
sort $PIPELINE/exons_gc_count.list > $PIPELINE/exons_gc_count_sorted.list

join -1 1 -2 1 -a1 -eERROR -o '0,1.2,1.3,2.2' $PIPELINE/exon_to_gene.list $PIPELINE/exons_gc_count_sorted.list | sort -k2,1 > $PIPELINE/Exons_vs_Genes.tsv

#Compile exon-specific data into gene-specific data
#echo -e "\t\t" > $PIPELINE/gc_content.tsv
#echo -e "#Gene_Name\tGene_Length\t%GC_ref" >> $PIPELINE/gc_content.tsv
awk '{a[$2]+=$3; b[$2]+=$4}END{for (i in a){c=100*b[i]/a[i]; printf "%s\t%s\t%3.1f\n", i, a[i], c;}}' $PIPELINE/Exons_vs_Genes.tsv | sort -k1 >> $PIPELINE/gc_content_sorted.tsv

join -1 1 -2 1 -a1 -t $'\t' -e'NULL' -o '0,2.2,2.3,2.4,2.5,1.3' $PIPELINE/gc_content.tsv $PIPELINE/sheet_values.tsv | less > $PIPELINE/Exons_vs_Genes2.tsv

#Make exon-specific GC content for statistics program
reset_file $PIPELINE/exons_gc_content.tsv
cat $PIPELINE/Exons_vs_Genes.tsv | while read a b c d ; do echo -e "$a\t`echo "scale=1; 100 * $d / $c " | bc`" ; done >> $PIPELINE/exons_gc_content.tsv

echo -e "\t\t\t\t" > $current_file
echo -e "#queryName\tqueryLength\tComment\thitDescription\thitName\t%GC_ref" >> $current_file
cat $PIPELINE/Exons_vs_Genes2.tsv >> $current_file

echo made $current_file from gff file 
else
echo $current_file already exists, skipping.
fi

      ################
############################  
###   Checking Indices   ###
############################
      ################


### Build genome index with bowtie2. 
echo -e "\nchecking indices"
if [ $SKIP -eq 0 ] || [ ! -f $(echo $PIPELINE/HtIndex.1.bt2) ]
then
echo constructing bowtie2 index for genome
BOWTIE2_OPTIONS="--threads $CORES -f"
bowtie2-build $BOWTIE2_OPTIONS \
$PIPELINE/HtGenome.fasta \
$PIPELINE/HtIndex
echo -e "Bowtie2-build settings, genome: $BOWTIE2_OPTIONS \n" >> $LOG
echo built index for exons with basename $PIPELINE/HtIndex
else
echo $PIPELINE/HtIndex.1.bt2 already exists, skipping.
fi


### Build exon index with bowtie2. 
if [ $SKIP -eq 0 ] || [ ! -f $(echo $PIPELINE/HtExonIndex.1.bt2) ]
then
echo constructing bowtie2 index for exons
BOWTIE2_OPTIONS="--threads $CORES -f"
bowtie2-build  $BOWTIE2_OPTIONS \
$PIPELINE/HtExons.fasta \
$PIPELINE/HtExonIndex
echo -e "Bowtie2-build settings, exons: $BOWTIE2_OPTIONS \n" >> $LOG
echo built index for exons with basename $PIPELINE/HtExonIndex
else
echo $PIPELINE/HtExonIndex.1.bt2 already exists, skipping.
fi

### Make separate folder of exons
if [ $SKIP -eq 0 ] || [ ! -d $(echo $PIPELINE/exons) ]
then
echo constructing separate folder of exons
mkdir_if_needed $PIPELINE/exons
cat ${PIPELINE}/HtExons.fasta | paste - - | while read line; do NAME=`echo $line | grep -oP "(?<=>)[^\s]+"`; CODE=`echo $line | grep -oP "(?<=\s).*"`; echo -e ">${NAME}\n${CODE}" > $PIPELINE/exons/${NAME}.fasta; done
#Alternate method (note: compare times)
#cat $PIPELINE/HtExons.fasta | paste - - | tr -d ">" | while read name code ; do echo -e ">${name}\n${code}" > $PIPELINE/exons/${name}.fasta; done 

else
echo $PIPELINE/exons already exists, skipping.
fi

infile=$PIPELINE/exon_scaffolds.fasta  #THIS ONE MAY NEED TO BE DELETED
outfile=$PIPELINE/exon_scaffolds.dict
### Make dictionary of fasta files
if [ $SKIP -eq 0 ] || [ ! -f $(echo $outfile) ]
then
echo Making dictionary from $infile
java -jar $PIPELINE/third_party_programs/picard/build/libs/picard.jar CreateSequenceDictionary R=$infile O=$outfile
else
echo $outfile already exists, skipping.
fi

infile=$PIPELINE/HtGenome.fasta
outfile=$PIPELINE/HtGenome.dict
### Make dictionary of fasta files
if [ $SKIP -eq 0 ] || [ ! -f $(echo $outfile) ]
then
echo Making dictionary from $infile
java -jar $PIPELINE/third_party_programs/picard/build/libs/picard.jar CreateSequenceDictionary R=$infile O=$outfile
else
echo $outfile already exists, skipping.
fi

infile=$PIPELINE/HtExons.fasta
outfile=$PIPELINE/HtExons.dict
### Make dictionary of fasta files
if [ $SKIP -eq 0 ] || [ ! -f $(echo $outfile) ]
then
echo Making dictionary from $infile
java -jar $PIPELINE/third_party_programs/picard/build/libs/picard.jar CreateSequenceDictionary R=$infile O=$outfile
else
echo $outfile already exists, skipping.
fi

##################################################################
######################################################
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

#If any new sequence was detected, make a new conf file.
if [ $SKIP -eq 0 ] || [ $changes -eq 1 ] || [ ! -f $(echo $FASTQ_SCREEN_CONF) ]
then 
echo "updating fastq_screen.conf" 
updated_list="##Malaria\nDATABASE\tMalaria\t$PIPELINE/HtIndex"
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
if [ $SKIP -eq 0 ] || ([ ! -f $(echo $FOLDER/${FILENAME}_trimmed_single_end.fastq) ] && [ ! -f $(echo $FOLDER/${FILENAME}_trimmed_single_end_reverse_min36.fastq) ])
then
echo -e "\nrunning trimmomatic, single-end"
START_TIME=`date +%s`
die_unless ${FILENAME}_R1_001.fastq

TRIM_OPTIONS="-phred33 -threads $CORES MINLEN:36 SLIDINGWINDOW:4:15 ILLUMINACLIP:${ADAPTER}:2:30:10"
#java -jar ~/bin/trimmomatic-0.36.jar SE \
java -jar $PIPELINE/third_party_programs/Trimmomatic-0.36/trimmomatic-0.36.jar SE \
${FILENAME}_R1_001.fastq $FOLDER/${FILENAME}_trimmed_single_end.fastq \
$TRIM_OPTIONS
echo -e "Trimmomatic settings: $TRIM_OPTIONS" >> $LOG
echo -e "Trimmomatic(forward) took" `time_since $START_TIME` "to run. \n" >> $LOG
else
echo $FOLDER/${FILENAME}_trimmed_single_end.fastq already exists, skipping.
fi

###########################################################3
########################################
### Run Trimmomatic, single-end reverse:
if [ $SKIP -eq 0 ] || ([ ! -f $(echo $FOLDER/${FILENAME}_trimmed_single_end_reverse.fastq) ] && [ ! -f $(echo $FOLDER/${FILENAME}_trimmed_single_end_reverse_min36.fastq) ])
then
echo running trimmomatic, single-end reverse
START_TIME=`date +%s`
die_unless ${FILENAME}_R2_001.fastq

TRIM_OPTIONS="-phred33 -threads $CORES MINLEN:36 SLIDINGWINDOW:4:15 ILLUMINACLIP:${ADAPTER}:2:30:10"
#java -jar ~/bin/trimmomatic-0.36.jar SE \
java -jar $PIPELINE/third_party_programs/Trimmomatic-0.36/trimmomatic-0.36.jar SE \
${FILENAME}_R2_001.fastq $FOLDER/${FILENAME}_trimmed_single_end_reverse.fastq \
$TRIM_OPTIONS
echo -e "Trimmomatic settings: $TRIM_OPTIONS" >> $LOG
echo -e "Trimmomatic(reverse) took" `time_since $START_TIME` "to make. \n" >> $LOG
else
echo $FOLDER/${FILENAME}_trimmed_single_end_reverse.fastq already exists, skipping.
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
current_file=$STATS_FOLDER/${FILENAME}_trimmed_single_end_min36_fastqc.html
if [ $SKIP -eq 0 ] || [ ! -f $(echo $current_file) ]
then
echo making fastqc file
START_TIME=`date +%s`
die_unless $FOLDER/${FILENAME}_trimmed_single_end_min36.fastq 
die_unless $FOLDER/${FILENAME}_trimmed_single_end_reverse_min36.fastq
mkdir_if_needed $STATS_FOLDER

fastqc $FOLDER/${FILENAME}_trimmed_single_end_min36.fastq --outdir=$STATS_FOLDER
fastqc $FOLDER/${FILENAME}_trimmed_single_end_reverse_min36.fastq --outdir=$STATS_FOLDER
echo -e "FastQC took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $current_file already exists, skipping.
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
current_file=$FOLDER/${FILENAME}_double_single_end_genome.sam
if [ $SKIP -eq 0 ] || ([ ! -f $(echo $current_file) ] && [ ! -f $(echo $FOLDER/${FILENAME}_double_single_end_genome_sorted.bam) ]&& [ ! -f $(echo $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam) ] ) 
then
echo -e "\nstarting bowtie2, single-end genome, $current_file"
START_TIME=`date +%s`
die_unless $FOLDER/${FILENAME}_trimmed_single_end_min36.fastq
die_unless $FOLDER/${FILENAME}_trimmed_single_end_reverse_min36.fastq 

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
-U $FOLDER/${FILENAME}_trimmed_single_end_min36.fastq \
-U $FOLDER/${FILENAME}_trimmed_single_end_reverse_min36.fastq \
-x $PIPELINE/HtIndex \
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
# (file name is backup for read_group data; a default name like '512022_S1_L001' contains data on sample name and lane nr .)

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
### Variant Calling, genome. 
if [ $SKIP -eq 0 ] || [ ! -f $(echo $FOLDER/${FILENAME}_variants.vcf) ]
then

echo -e "\nvariant calling, genome"
START_TIME=`date +%s`
die_unless $PIPELINE/HtGenome.fasta
die_unless $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam

PILEUP_OPTIONS="-g -f"
BCF_OPTIONS="-m -v"
samtools mpileup $PILEUP_OPTIONS $PIPELINE/HtGenome.fasta $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam > $FOLDER/${FILENAME}_variants.bcf

die_unless $FOLDER/${FILENAME}_variants.bcf
bcftools call $BCF_OPTIONS $FOLDER/${FILENAME}_variants.bcf > $FOLDER/${FILENAME}_variants.vcf

echo -e "SAMtools settings, pile-up genome: $PILEUP_OPTIONS" >> $LOG
echo -e "SAMtools settings, call genome: $BCF_OPTIONS" >> $LOG
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
bcftools index $FOLDER/${FILENAME}_variants.bcf

echo -e "Sorting VCF took" `time_since $START_TIME` "to make." >> $LOG
else
echo $FOLDER/${FILENAME}_variants_bcf.csi already exists, skipping.
fi


#######################################################################
###########################################
### GATK Variant Calling. 
current_file=$FOLDER/${FILENAME}_genome_variants_gatk.vcf
if [ $SKIP -eq 0 ] || [ ! -f $(echo $current_file) ]
then
echo "variant calling, genome GATK"
die_unless $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam

START_TIME=`date +%s`
PILEUP_OPTIONS="--variant_index_type LINEAR --variant_index_parameter 128000" #<---Check that this number is reasonable 

java -jar $PIPELINE/third_party_programs/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R $PIPELINE/HtGenome.fasta \
-I $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam \
-o $FOLDER/${FILENAME}_genome_variants_gatk.vcf $PILEUP_OPTIONS

echo -e "Variant calling settings, GATK: $PILEUP_OPTIONS" >> $LOG
echo -e "Variant calling (GATK) took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $current_file already exists, skipping.
fi


      #############
#########       ##########
#### CONSENSUS-MAKING ####
#########       ##########
      #############



#######################################################################
###########################################
### Make consensus sequence.

if [ $SKIP -eq 0 ] || [ ! -f $(echo $FOLDER/${FILENAME}_consensus.fasta) ] || [ ! -d $(echo $FOLDER/scaffolds) ]
then
echo -e "\nmaking consensus sequence"
die_unless $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam
die_unless $PIPELINE/HtGenome.fasta

START_TIME=`date +%s`
mkdir_if_needed $FOLDER/scaffolds
touch $FOLDER/${FILENAME}_consensus.fasta
echo -e "starting...\n\n"

#Get SAM files of exons from genome match

cat $PIPELINE/HtGenome.fasta | grep -o "Ht[^\n]*" | \
while read line; do 
name=`echo $line | sed "s/\s//"`
if [ `samtools view $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $name -c` -eq 0 ]; then
#echo -en "\e[1A"; echo -e "\e[0K\rskipping $name"
continue 
fi
#echo -en "\e[1A"; echo -e "\e[0K\rFinding reads from <${name}>";\
scaffold_length=`cat $PIPELINE/HtGenome.fasta | grep $name -A1 | tail -n1 | tr -d "\n" | wc -c`
{
samtools mpileup -uf $PIPELINE/HtGenome.fasta $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam -r $name -d 100000 | bcftools call -c --ploidy 1 > $FOLDER/scaffolds/${name}_consensus.vcf;\
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
echo -e ">$name\n$scaffold_sequence" >> $FOLDER/${FILENAME}_consensus.fasta ; \
echo -en "\e[1A"; echo "created $FOLDER/scaffolds/${name}_consensus.fasta" 
done
echo -en "\e[1A"; echo -e "\e[0K\rDone";\

echo -e "Samtools pileup consensus file took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $FOLDER/${FILENAME}_consensus.fasta already exists, skipping.
fi


#######################################################################
###########################################
#  Create exon consensus file. 
current_file=$FOLDER/${FILENAME}_exon_consensus2.fasta
if [ $SKIP -eq 0 ] || [ ! -f $current_file ]
then
echo -e "Creating $current_file, exon consensus file"
die_unless $PIPELINE/HtExons.fasta

START_TIME=`date +%s`
reset_file $current_file

#CUT OUT EXONS AND MAKE EXON CONSENSUS
cat $PIPELINE/HtExons.fasta | grep -o "Ht[^\n]*" | while read line; do 
scaffold_name=`echo "$line" | grep -oP "Ht[^:]*"`
[[ -f $( echo "$FOLDER/scaffolds/${scaffold_name}_consensus.fasta" ) ]] || continue #Skip if file doesn't exist
interval=`echo "$line" | grep -oP "\d+-\d+"`
#echo "$interval"
IFS='-'; read start end <<< "$interval"
length=$[end-start]
scaffold_code=`cat "$FOLDER/scaffolds/${scaffold_name}_consensus.fasta" | grep -A1 "$scaffold_name" | grep -v "^>"`
sequence=${scaffold_code:start:length} #cut out exon substring from scaffold
#echo $scaffold_name
#echo $sequence
[[ $sequence =~ [ATGCatgc]+ ]] || continue  #Skip if sequence is just 'N' or '-'
echo -e ">$line\n$sequence" >> "$current_file"
done

echo -e "$current_file took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $current_file already exists, skipping.
fi

#######################################################################
###########################################
#  Create coverage stats file.
infile=$FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam
outfile=$FOLDER/exon_coverage.tsv
ref_file=$PIPELINE/exon_to_gene.list
if [ $SKIP -eq 0 ] || [ ! -f $outfile ]
then
echo -e "Creating $outfile, exon coverage stats file"
die_unless $infile
die_unless $ref_file
touch $outfile
echo "starting..."
START_TIME=`date +%s`

cat $ref_file | while read line gene length ; do name=`echo $line | grep -oP Ht[^\(]+`; echo -en "$name\t"; echo -en "$length\t"; samtools depth -aa -r "$name" $infile | awk -v l=$length '{if($3>2) mapped+=1; else unmapped+=1}END{if (mapped>0) print mapped/(l+1); else print 0}' ; done > $outfile

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
echo "starting..."
START_TIME=`date +%s`

export FOLDER
export FILENAME
export READ_LIMIT
export -f mark_depth 

cat $infile | paste - - | \
xargs -I{} --max-procs $CORES bash -c 'line="{}"; 
read name exon_code <<< $line; 
exon_name=`echo $name | tr -d ">" | cut -d "(" -f1`;
echo -e "\e[1A" "Processing $exon_name                    ";

#Ignore if no reads
read_nr=`samtools view $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $exon_name | wc -l`
if [[ $read_nr -eq 0 ]] ; then 
echo -e "ERROR: No reads in $exon_name\n";
exit 0
fi

#Fast-process if all reads
#grep $exon_name 

#Note: Maximum readcount capped at 9 
exon_read_counts=`samtools depth -aa "$FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam" -r "$exon_name" | cut -f3 | sed -e "s/[0-9]\{2,\}/9/g" | tr -d "\n" | tail -c +2`;  
current_fasta="$FOLDER/cap_consensus/exons/$exon_name.fasta";
touch "$current_fasta";

cap_code=`mark_depth "$exon_code" "$exon_read_counts" | sed "s/n/-/g"` 
if [[ $cap_code =~ [ATGCNatgcn-]+ ]] ; then
echo -e "$name\n$cap_code" > $current_fasta;
else
echo -e "\e[1A" "ERROR: No sequence found for exon $name!\cap_code=<$cap_code>\n exon_code = <$exon_code>\n exon_read_counts=<$exon_read_counts>\n"
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

if [ $SKIP -eq 0 ] || [ ! -d $(echo $FOLDER/consensus) ]
then
echo -e "making multifastas"
die_unless $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam #THIS MAY FAIL
die_unless $PIPELINE/cigar_parser.py
START_TIME=`date +%s`
mkdir_if_needed $FOLDER/multifastas
mkdir_if_needed $FOLDER/consensus
mkdir_if_needed Consensus
echo -e "" > $FILTER_OUT
ADDRESS=`pwd`/$FOLDER 
AVG_LEN=`samtools view -F 0x4 $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam | cut -f 10 | awk '{EXON_READ_NR+=1; TOTAL_LENGTH+=length($1)}END{print int(TOTAL_LENGTH/EXON_READ_NR+0.5) }'`
echo -n "" > ${FOLDER}/$FILTER_OUT
echo -n "" > ${FOLDER}/cigar_parser.log
reset_file $FOLDER/temp.txt
#Make variables and 'return_unless_DNA_match()' available for xargs to use
export -f return_unless_DNA_match
export -f time_since
export -f subsample_bam
export ADDRESS
export FOLDER
export FILENAME
export PIPELINE
export AVG_LEN
export FILTERED
export FILTER_OUT
export EXON_CONTIGS
#Iterate over names
echo "Starting..."
cat $PIPELINE/HtExons.fasta | grep -o "Ht[^\n]*" | \
xargs -I{} --max-procs $CORES bash -c 'line="{}"; 
name=`echo "$line" | sed "s/\s//"`;
name_stripped=`echo "$name" | sed "s/(.*//"`;
scaffold=`echo "$name" | sed "s/:.*//"`;
#echo starting $FOLDER/$name
#echo line = "$line"
START=`date +%s.%3N`;
nr_of_reads=`samtools view $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $name_stripped | wc -l`;
if [[ nr_of_reads -gt 1000 ]] ; then #If there are more than 1000 reads, subsample. 
seq_len=`cat $PIPELINE/HtExons.fasta | grep $name_stripped -A1 | tail -n1 | tr -d "\n" | wc -c`;
out=`subsample_bam $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $name_stripped | sort -k4n | $PIPELINE/cigar_parser.py $PIPELINE/HtGenome.fasta $name --speedy`;
else
out=`samtools view $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $name_stripped | return_unless_DNA_match "40" | $PIPELINE/cigar_parser.py $PIPELINE/HtGenome.fasta $name --speedy`;
fi;
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
line_nr=`samtools view $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $name_stripped | wc -l`;
echo -e "$COUNTER\t$NUMBER\t$name\t$line_nr\t$size_nr\t$time_nr" >> ${FOLDER}/cigar_parser.log;
fi;
echo -n "." >> "$FOLDER/temp.txt";
if [ $[ COUNTER % 10 ] == 0 ] ; then 
echo -en "\e[1A"; echo -e "$[ 100 * COUNTER / EXON_CONTIGS ] %    ($COUNTER of $EXON_CONTIGS)      ";
fi
exit 1;'
echo -en "\e[1A"; echo "Done!                                "
cat ${FOLDER}/consensus/*.fasta > ${FOLDER}/${FILENAME}_exon_consensus.fasta
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
cp $FOLDER/${FILENAME}_exon_consensus.fasta Consensus/${FILENAME}_exon_consensus.fasta

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

current_file=$FOLDER/${FILENAME}_gene_length_and_id.tsv
if [ $SKIP -eq 0 ] || [ ! -f $(echo $current_file) ]
then
echo -e "\ngenerating $current_file" 
START_TIME=`date +%s`
die_unless $FOLDER/${FILENAME}_exon_consensus3.fasta
reset_file $FOLDER/${FILENAME}_exon_length_and_id.tsv
echo "" > $FOLDER/temp4

#Make Gene-vs-Exon hash 
declare -A MOIN_REVISIONS exon_gene #Allow strings with parentheses etc to be keys.  
while read exon gene length # $length not used, but having it gets rid of that column. 
do
    exon_gene[$exon]=$gene
done < $PIPELINE/exon_to_gene.list

#For each exon:
cat $FOLDER/${FILENAME}_exon_consensus3.fasta | paste - - | tr -d ">" | while read name code ; 
do 
echo -e ">$name\n$code" > $FOLDER/temp
cat $PIPELINE/exons/${name}.fasta > $FOLDER/temp2
coverage=`echo $code | tr -cd "AGCTagct" | wc -c`
[ $coverage -eq 0 ] && continue #'Keeps script from dividing by zero if something goes wrong. 
certain_coverage=`echo $code | tr -cd "AGCT" | wc -c`
{
stretcher $FOLDER/temp2 $FOLDER/temp $FOLDER/temp3 #Needleman-Wunsch rapid global alignment of two sequences
} > /dev/null 2> /dev/null #Silence this script, so it doesn't report every time it's run.
identity=`cat $FOLDER/temp3 | grep "Identity" | grep -oP "\d.*"`
gene=${exon_gene[$name]}
#echo -e "$name\t$gene\t$certain_coverage\t$matches\t$identity\t$percent%"
matches=`echo "$identity" | cut -d '/' -f1` 
percent=`bc <<< "scale=1; 100*$matches/$coverage"` 
echo -e "$name\t$gene\t$coverage\t$matches\t$identity\t$percent%" >> $FOLDER/${FILENAME}_exon_length_and_id.tsv
echo -e "$name\t$gene\t$certain_coverage\t$matches\t$identity\t$percent%" >> $FOLDER/${FILENAME}_exon_length_and_id_test.tsv
cat $FOLDER/temp3 >> $FOLDER/temp4 
done

#Compile exon data into gene data 
gawk '{a[$2]+=$3; b[$2]+=$4}END{for (i in a){c=100*b[i]/a[i]; printf "%s\t%s\t%3.1f\n", i, a[i], c;}}' $FOLDER/${FILENAME}_exon_length_and_id.tsv | sort -k1 > $current_file 

echo made $current_file 
echo -e "$current_file took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $current_file already exists, skipping. 
fi


########################################################################
################################################################
### Record Probe data

infile=${FOLDER}/${FILENAME}_double_single_end_genome_sorted_rg.bam
outfile=${FOLDER}/${FILENAME}_probe_stats.tsv
outfile2=Consensus/${FOLDER}_${FILENAME}_probe_stats.tsv
ref_file=$PIPELINE/Ht_probes.csv
if [ $SKIP -eq 0 ] || [ ! -f $(echo $outfile) ] 
then
echo -e "generating $outfile."
START_TIME=`date +%s`
die_unless $infile
die_unless $ref_file
  
echo -e "probe\tGC\t${FILENAME}_reads\t${FILENAME}_coverage" > $outfile
tail -n +2 $ref_file | while read scaffold probe seq repl strand position ;
do 
	pos=`echo $position | sed "s/^chr//"` ; 
	length=`echo $seq | tr -cd 'AGCT' | wc -c`;
	GC=`echo $seq | tr -cd 'GC' | wc -c` ; 
	gc_content=`bc <<< "scale=1; 100*$GC/$length"` ; 
	echo -ne "$pos\t$gc_content\t" >> $outfile ; 
	reads=`samtools view "$infile" "$pos" -c` ;
	echo -ne "$reads\t" >> $outfile ;
	if [ $reads -eq 0 ] 
	then
		echo "0" >> $outfile
	else
	samtools depth -aa "$infile" -r "$pos" | awk -v l=$length '{if($3>2) mapped+=1; else unmapped+=1}END{if (mapped>0) print (100*mapped/l); else print 0}' >> $outfile;
	fi
done
cp $outfile $outfile2

echo made $outfile 
echo -e "$outfile took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $outfile already exists, skipping. 
fi

#######################################################################
###############################################################
### Compile length/identity statistics in consensus folder. 

current_file=Consensus/${FOLDER}_${FILENAME}_gene_length_and_id.tsv
if [ $SKIP -eq 0 ] || [ ! -f $(echo $current_file) ] 
then
echo -e "generating $current_file."
START_TIME=`date +%s`
die_unless $FOLDER/${FILENAME}_gene_length_and_id.tsv

lineage="" # Will be filled in by hand; researcher will know better than BLAST search, especially with novel strains not present in DB.

#Add coverage column to coverage file in consensus folder
echo -e "$lineage\n${FILENAME}_length" > Consensus/${FOLDER}_${FILENAME}_gene_coverage.tsv
join -1 1 -2 1 -a1 -e'0' -o '2.2' $PIPELINE/gc_content_sorted.tsv $FOLDER/${FILENAME}_gene_length_and_id_test.tsv >> Consensus/${FOLDER}_${FILENAME}_gene_coverage.tsv

#Add identity column to identity file in consensus folder
echo -e "$lineage\n${FILENAME}_identity" > Consensus/${FOLDER}_${FILENAME}_gene_identity.tsv
join -1 1 -2 1 -a1 -e'0' -o '2.3' $PIPELINE/gc_content_sorted.tsv $FOLDER/${FILENAME}_gene_length_and_id.tsv  >> Consensus/${FOLDER}_${FILENAME}_gene_identity.tsv


cp $FOLDER/${FILENAME}_gene_length_and_id.tsv Consensus/${FOLDER}_${FILENAME}_gene_length_and_id.tsv

echo -e "$current_file took" `time_since $START_TIME` "to make.\n" >> $LOG
else
echo $current_file already exists, skipping.
fi

#######################################################################
###############################################################
### Make summary of gc_content. 

#Join GC/length file, coverage file and identity file into one new file.  
paste $PIPELINE/gc_content.tsv Consensus/*_gene_coverage.tsv Consensus/*_gene_identity.tsv > Consensus/gene_length_and_id.tsv

            ##############
######################################
##### SEQUENCE QUALITY STATISTICS ####
######################################
            ##############


################################################################
#################################################
### Make statistics for double single end. 

if [ $SKIP -eq 0 ] || [ ! -f $(echo $FOLDER/statistics_double_single_end.txt) ]
then
echo -e "\ncalculating basic statistics for $FILENAME, double single end \n"
die_unless $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam
die_unless $FOLDER/${FILENAME}_trimmed_single_end_min36.fastq
die_unless $FOLDER/${FILENAME}_trimmed_single_end_reverse_min36.fastq
die_unless $PIPELINE/HtExons.fasta
die_unless ${FOLDER}/$FILTER_OUT

STAT=$FOLDER/statistics_double_single_end.txt
touch $STAT
START_TIME=`date +%s`
echo -e "\nStatistics for $FILENAME, double single end: \n" > $STAT
TOTAL_READS=$[(`cat ${FILENAME}_R[12]_001.fastq | wc -l`)/4]
READS_FILTERED_OUT=`cat ${FOLDER}/$FILTER_OUT | wc -l`
cat ${FOLDER}/$FILTER_OUT | sort | tr "\t" "\n" | awk '!a[$0]++' > ${FOLDER}/removed_reads #Sort the (asynchronously generated) file under scaffold headers
READS_SINGLE=$[(`cat $FOLDER/${FILENAME}_trimmed_single_end_min36.fastq | wc -l`)/4 + (`cat $FOLDER/${FILENAME}_trimmed_single_end_reverse_min36.fastq | wc -l`)/4]
TIME_READS=`date +%s`
echo -e "Time: " $[`date +%s` - START_TIME] ": Reads counted."
MAPPED_READS_GENOME=`samtools view -F 0x4 $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam | cut -f 1 | sort | uniq | wc -l `
TIME_MAPPED=`date +%s`
echo -e "Time: " $[`date +%s` - START_TIME] ": Mappings counted."

read EXONS_HIT <<< $(for line in `cat $PIPELINE/HtExons.fasta | grep -o "Ht[^\(]*"`; do samtools view $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $line | wc -l; done | awk '{if($0>0) nr+=1}END{print nr}')
read EXONS_HIT_TEN <<< $(for line in `cat $PIPELINE/HtExons.fasta | grep -o "Ht[^\(]*"`; do samtools view $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $line | wc -l; done | awk '{if($0>9) nr+=1}END{print nr}')
read EXONS_HIT_HUNDRED <<< $(for line in `cat $PIPELINE/HtExons.fasta | grep -o "Ht[^\(]*"`; do samtools view $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam $line | wc -l; done | awk '{if($0>99) nr+=1}END{print nr}')

echo -e "Time: " $[`date +%s` - START_TIME] ": Exon hits counted."

GENOME_BASES_HIT=0
GENOME_BASES_MISSED=0

read GENOME_BASES_HIT GENOME_BASES_MISSED <<< $(samtools depth -aa $FOLDER/${FILENAME}_double_single_end_genome_sorted_rg.bam | awk '{if($3>2) total+=1; else unmapped+=1}END{print total" "unmapped}')

END_TIME=`date +%s`
echo -e "Time: " $[`date +%s` - START_TIME] ": Base hits counted."
RUNTIME=$(($END_TIME-$START_TIME))

echo -e "Number of reads, including bad ones: $TOTAL_READS.  \n" >> $STAT
echo -e "Number of used reads: $READS_SINGLE; $[100*$READS_SINGLE/$TOTAL_READS]% \n" >> $STAT
echo -e "Number of reference-identical reads filtered out: $READS_FILTERED_OUT; `echo "scale=3; 100*$READS_FILTERED_OUT/$READS_SINGLE" | bc`% of used reads\n" >> $STAT
echo -e "Reads matching genome: $MAPPED_READS_GENOME; `echo "scale=1; 100*$MAPPED_READS_GENOME/$READS_SINGLE" | bc`% of used reads \n" >> $STAT
echo -e "Number of genome contigs: $GENOME_CONTIGS \n" >> $STAT
echo -e "Exons with at least one match: $EXONS_HIT; `echo "scale=1; 100*$EXONS_HIT/$EXON_CONTIGS" | bc`%  \n" >> $STAT 
echo -e "Exons with at least ten matches: $EXONS_HIT_TEN; `echo "scale=1; 100*$EXONS_HIT_TEN/$EXON_CONTIGS" | bc`%  \n" >> $STAT
echo -e "Exons with at least a hundred matches: $EXONS_HIT_HUNDRED; `echo "scale=1; 100*$EXONS_HIT_HUNDRED/$EXON_CONTIGS" | bc`%  \n" >> $STAT 

printf "Coverage breadth: \n
Genome bases covered by at least 3 reads: $GENOME_BASES_HIT; `echo "scale=1; 100*$GENOME_BASES_HIT/$GENOME_SIZE" | bc`%%  \n 
Genome bases missed: $GENOME_BASES_MISSED bases; `echo "scale=1; 100*$GENOME_BASES_MISSED/$GENOME_SIZE" | bc`%%  \n"

printf "Bases accounted for: \nGenome: $[$GENOME_BASES_HIT+$GENOME_BASES_MISSED] of $GENOME_SIZE; $[$GENOME_SIZE - ($GENOME_BASES_HIT+$GENOME_BASES_MISSED)] missed\n" >> $STAT

parse_fastq_screen >> $STAT

echo -e "Total Script runtime:" `time_since $SCRIPT_START_TIME`

echo -e "Calculating statistics took" `time_since START_TIME` >>$LOG "to do."
echo -e "Program as a whole took:" `time_since $SCRIPT_START_TIME` >>$LOG
echo -e "Time: " $[`date +%s` - START_TIME] ": Stats recorded in file."
cat $STAT
echo -e "These statistics took $RUNTIME seconds to calculate. Counting reads: $[$TIME_READS-$START_TIME], counting mapped reads: $[$TIME_MAPPED-$TIME_READS], counting exons hit: $[$END_TIME-$TIME_MAPPED]\n\n"
else
echo statistics_double_single_end.txt already exists, skipping.
fi


#######################################################################
###########################################
#  Compile match statistics

if [ $LAST_ITERATION -eq 1 ]
then
echo -e "Compiling statistics for all folders"
START_TIME=`date +%s`

echo "Statistics compiled at `date`" > match_stats.txt
ls */stat* | while read name 
do 
  filename=`echo "$name" | cut -d '/' -f1`
  total=`grep "Number of reads" $name | grep -oP "[0-9]*"`
  matching=`grep "Reads matching genome" $name | grep -oP "[0-9\.%]*" | tr "\n" " "`
  echo -e "$filename $total $matching" >> match_stats.txt 
done

echo -e "Compiling match statistics took" `time_since $START_TIME` "to align.\n" >> $LOG
else
echo Skipping sequence alignment.
fi 

#######################################################################
###########################################
#  Compile probe statistics

if [ $LAST_ITERATION -eq 1 ]
then
echo -e "Compiling probe statistics for all folders"
START_TIME=`date +%s`

cat $FOLDER/${FILENAME}_probe_stats.tsv | cut -f1,2 > Consensus/probes_stats.tsv
ls Consensus/*probe_stats.tsv | while read line ; 
do 
  paste Consensus/probes_stats.tsv <( cat $line | cut -f3 ) > temp ; 
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
