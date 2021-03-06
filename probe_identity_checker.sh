#This file matches probes to their corresponding exons and genes. 
#Declare variables
probe_to_exon_file="$1"
reference_file="$2"
consensus_file="$3"
verbose=0
exon_sequence=""

#Read in all exons
#Sort?

#Store gene_name, exon_name, Scaffold, start and end
#sort?

#Read in all probes, Store Scaffold, start and end.
#For each probe,

while true
do
  read -r name probe_pos exon_pos gene_name <&3 || break
  if [ $exon_sequence -eq "" ] #first time
  then #load first exon sequence
    read -r pos <&4 || break
    read -r exon_sequence <&4 || break
  fi

#Scroll down to the line that contains the exon

#Find common part of sequence

#Get identity 

#Print identity

  read sca start end <<< `echo "$pos" | tr ":-" "  "` 
  pos=`echo $pos | sed "s/chrHt/Ht/"`  
#  echo -en "pos=$pos\n start=$start\n end=$end"
  #Go to correct scaffold
#  >&2 echo "comparing $scaffold_nr vs $exon_scaffold_nr"
  until [ "$scaffold_nr" -le "$exon_scaffold_nr" ] 
  do
    [[ $verbose -eq 1 ]] && >&2 echo "$scaffold_nr > $exon_scaffold_nr, loading new exon"
    read -r exon_pos gene_name gene_length <&4 || break
    read exon_scaffold exon_start exon_end trash <<< `echo "$exon_pos" | cut -d '(' -f1 | tr ":-" "  "`
    exon_scaffold_nr=`echo "$exon_scaffold" | sed "s/.*fold//"`
    #>&2 echo "new exon: $exon_pos, nr $exon_scaffold_nr"
#    echo "exon_scaffold=$exon_scaffold, exon_start=$exon_start, exon_end=$exon_end"
#    echo "checked $scaffold vs $exon_scaffold"
  done
   
  #Go to correct exon 
  #Unless probe overlaps exon start, proceed to next exon.
  until ( [ $start -ge $exon_start -a $start -le $exon_end ] || [ $exon_start -gt $end ] || [ $end -ge $exon_start -a $end -le $exon_end ] || [ $exon_start -ge $start -a $exon_start -le $end ] || [ $exon_end -ge $start -a $exon_end -le $end ])
  do
    [[ $verbose -eq 1 ]] && >&2 echo " ( $pos ) does not overlap ( $exon_pos ). Loading new exon"
    read -r exon_pos gene_name gene_length <&4 || break
    read exon_scaffold exon_start exon_end <<< `echo "$exon_pos" | cut -d '(' -f1 | tr ":-" "  "`
    exon_scaffold_nr=`echo "$exon_scaffold" | sed "s/.*fold//"` 
    
#    >&2 echo "$scaffold_nr <-> $exon_scaffold_nr"
    if [ $exon_scaffold_nr -gt $scaffold_nr ] ; then 
      [[ $verbose -eq 1 ]] && >&2 echo "SCAFFOLD DOESN'T MATCH!" 
      break
    elif [ $exon_start -gt $end ] ; then
      [[ $verbose -eq 1 ]] && >&2 echo "EXON TOO FAR AHEAD!" 
      break    
    fi
#    echo "checked $start vs $exon_start"
  done

  #>&2 echo "FOUND: $pos corresponds to $exon_pos"
  if [ "$scaffold" != "$exon_scaffold" ] ; then
    echo -e "$name\t$pos\t\t"
    >&2 echo -e "WARNING: probe ${scaffold}:$start-$end did not overlap with any exon. Loop is at exon $exon_pos"
  elif ( [ $start -ge $exon_start ] && [ $end -le $exon_end ] ) ; then #print line
    echo -e "$name\t$pos\t$exon_pos\t$gene_name"
  else #not an error
    echo -e "$name\t$pos\t$exon_pos\t$gene_name"
#    >&2 echo -e "Note: probe ${scaffold}:$start-$end is out of bounds of exon $exon_pos"
  fi       

    #Load next probe and continue (they are ordered) 
done 3< $exons_file 4< $probe_to_exon_file
