#This file matches probes to their corresponding exons and genes. 
#Declare variables
probes_file=""
exons_file=""
exon_scaffold="" #initial value

#Read in all exons
#Sort?

#Store gene_name, exon_name, Scaffold, start and end
#sort?

#Read in all probes, Store Scaffold, start and end.
#For each probe,

cat $probes_file | while read scaffold start end ; 
while true
do
  read -r scaffold name seq rep str pos <&3 || break
  read sca start end <<< `echo "$pos" | tr ":-" "  "`   
  
  #Go to correct scaffold
  until [ $scaffold -eq $exon_scaffold ] 
  do
    read -r exon_pos gene_name gene_length <&4 || break
    read exon_scaffold exon_start exon_end <<< `echo "$exon_pos" | cut -d '(' -f1 | tr ":-" "  "`
  done
  if [ $scaffold -gt $exon_scaffold ] then
    echo >2 "ERROR: probe ${scaffold}:$start-$end did not find a match!"
  fi
  
  #Go to correct exon 
  #Unless probe start greater or equal to exon start, continue
  [[ $start -lt $exon_start ]] && continue
  #Double check: probe end lesser or equal to exon end.
  #Otherwise, proceed to next exon.
  until [ $end -gt $exon_end ] 
  do 
    read -r exon_pos gene_name gene_length <&4 || break
    read exon_scaffold exon_start exon_end <<< `echo "$exon_pos" | cut -d '(' -f1 | tr ":-" "  "`
  done
  
  if ( [ ! $scaffold -eq $exon_scaffold ] && [ ! $start -lt $exon_start ] && [ ! $end -gt $exon_end ]) 
  then #print line
    echo -e "$probe_name $probe_gc $probe_length $probe_identity $exon_line"
  else #report error
    echo >2 "ERROR: probe ${scaffold}:$start-$end did not find a match!"
  fi       
    #Load next probe and continue (they are ordered) 
  done
 

done
