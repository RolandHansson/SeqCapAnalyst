MIN_LENGTH_TO_DISCARD=40
FOLDER=.
FILTER_OUT=filter_out.sam


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


function subsample_bam () #Note: will add 'R' or 'F' to read-names in resulting SAM. 
{
  #Create variables
  samfile=$1
  exon=$2
  READ_MINIMUM=$3
  temp_file=$(mktemp /tmp/foo.XXXXXXXXX)
  VERBOSE=0

  >&2 echo -e "\e[1A Subsampling ${exon}...                        "
  [[ $VERBOSE -eq 1 ]] && >&2 echo "samfile=$samfile, exon=$exon"
  declare -A MOIN_REVISIONS startpoints #Array for end points
  declare -A MOIN_REVISIONS lines #Array to store SAM lines
  scaffold=`echo "$exon" | sed "s/:.*//"`;
  interval=`echo "$exon" | grep -oP "\d+-\d+"`
  IFS='-'; read start end <<< "$interval"
  IFS=$' \t\n'
  #Run loop for positions in exon
  [[ $VERBOSE -eq 1 ]] && >&2 echo "start=$start, end=$end, IFS='$IFS'"
  for i in $(seq $end -1 $start) ; do #In reverse order
    [[ $i%30 -eq 1 ]] || [[ $i -eq $start ]] || [[ $i -eq $end ]] || continue  #*To save time, do this only every N positions in code
    #Remove entries from array
    for x in ${!startpoints[@]} ; do
      if [ ${startpoints[$x]} -gt $i ] ; then 
        #print to SAM and remove from arrays 
        [[ $VERBOSE -eq 1 ]] && >&2 echo "startpoint=${startpoints[$x]}, i=$i"
        echo "${lines[$x]}"
        [[ $VERBOSE -eq 1 ]] && >&2 echo -n "i=$i, removing $x (starts at ${startpoints[$x]}), " 
        unset 'startpoints[$x]' 
        unset 'lines[$x]'
        [[ $VERBOSE -eq 1 ]] && >&2 echo "coverage is now ${#startpoints[@]}"
      fi 
    done

    coverage=${#startpoints[@]} # nr of entries in array 
    [[ coverage -lt $READ_MINIMUM ]] || continue #If enough reads, continue. 

    difference=$[READ_MINIMUM - coverage] 
    #Otherwise, load [min_reads] more files, if available, and add until [min_reads] is reached. 
    #(Note: some may have been added previously) 
    [[ $VERBOSE -eq 1 ]] && >&2 echo "searching for '$samfile' '${scaffold}:${i}-$end'"
    #>&2 echo -e "\e[1A Subsampling ${exon}; at pos ${i} of $end                     "
    end1=$[$i+200] 
    j=$[i+30]
samtools1.7 view "$samfile" "${scaffold}:${i}-${j}" | return_unless_DNA_match "$MIN_LENGTH_TO_DISCARD" | head -n $[ 5 + READ_MINIMUM * 2 ] > $temp_file

    while read line; do
      IFS=$' \t\n'
      [ -z "$line" ] && continue #Skip empty lines
      read name flag t3 start_pos t5 cigar t7 t8 t9 code trash10 <<< "$line"
      [[ $VERBOSE -eq 1 ]] && >&2 echo -e "line='$line' \nname:'$name' flag:'$flag' '$t3' start:'$start_pos' '$t5' cigar_'$cigar' '$t7' '$t8' '$t9' code:'$code' '$trash10' \n"
      #>&2 echo -e "\nflag: '$flag', name:'$name', line: '$line'\n"

      direction=`((($flag&16)>0)) && echo 'R' || echo 'F'`
      name="${name}$direction"
      if test "${startpoints[${name}]+isset}" ; then continue 
      #else 
      #  echo "'$name' not in list. List:"
      #  for y in ${!startpoints[@]} ; do echo -n "$y, " ; done ; echo ""
      fi
      if test "${startpoints[${name}]+isset}" ; then continue ; fi #If already in list, skip.
      
      #Consider cigar string. Length is # of 'M' + # of 'D', ignoring # of 'I'.
      length=$[0 `echo $cigar | sed -r "s/([0-9]+)[MDS]/+\1 /g" | sed -r "s/([0-9]+)I//g"`] #too complex? Takes too long? 
      [[ $VERBOSE -eq 1 ]] && >&2 echo "name='$name', start_pos='$start_pos', cigar='$cigar', code='$code'"
       
      #Add name and code to list
      startpoints[$name]=$[start_pos]
      lines[$name]=$line
      #echo "new line = '${lines[$name]}'"
      coverage=${#startpoints[@]}
      #>&2 echo "$i, added $name, coverage is now $coverage" 
      [[ $coverage -lt $READ_MINIMUM ]] || break
    done < <( cat $temp_file )
    rm "$temp_file"
    IFS=$OLDIFS
  done 
  
#Remove entries remaining in array
for x in ${!startpoints[@]} ; do
    #print to SAM and remove from arrays
    echo "${lines[$x]}" 
    #echo "unsetting $x"
    unset 'startpoints[$x]' 
    unset 'lines[$x]'
done
}

subsample_bam $1 $2 ${3-10} 
