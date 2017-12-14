#!/usr/bin/python
from __future__ import print_function #allows colored text?
# This script takes a SAM file and makes it into an aligned multifasta.
import sys
import re
import warnings
import argparse 
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

#Skip traceback and just print the warning
def custom_formatwarning(msg, *a): 
  if debug:
    print("\nwarning: " + str(msg) + '\n')
  return "warning: " + str(msg) + '\n'
warnings.formatwarning = custom_formatwarning

""" 
        Compact Idiosyncratic Gapped Alignment Report -> CIGAR
        CIGAR regex: ([0-9]+[MIDNSHPX=])+
        M alignment match (can be a sequence match or mismatch) <- DONE
        I insertion to the reference                            <- DONE
        D deletion from the reference                           <- DONE
        N skipped region from the reference                     -------
        S soft clipping (clipped sequences present in SEQ)      -------
        H hard clipping (clipped sequences NOT present in SEQ)  -------
        P padding (silent deletion from padded reference)       -------
        = sequence match                                        <- DONE
        X sequence mismatch                                     <- DONE

"""

########################
#   COMMAND OPTIONS    #
########################

parser = argparse.ArgumentParser()
parser.add_argument("reference_file", type=argparse.FileType('r'), help="fasta file to extract reference sequence from.")
parser.add_argument("reference_name", help="name of reference sequence.")
parser.add_argument("-i", "--input",  type=argparse.FileType('r'), help="fasta file to extract sequences from.")
parser.add_argument("-n", "--name", help="find specific sequence name")
parser.add_argument("-s", "--consensus_file", type=argparse.FileType('r'), help="(optional) external consensus file")
parser.add_argument("-b", "--blast_file", type=argparse.FileType('r'), help="(optional) external blast file, to help find start/end")
parser.add_argument("-g", "--gaps", help="Add gaps every 10 bases", action="store_true")
parser.add_argument("-f", "--full", help="Include tangential hits on reference sequence", action="store_true")
parser.add_argument("-c", "--colour", "--color", action="store_true", help="Colour output in the prompt. (Note: may render output unusable by other programs.)")
parser.add_argument("-d", "--debug", action="store_true", help="Adds debugging messages.")
parser.add_argument("-p", "--positioning", action="store_true", help="Adds position info to every base in sequence. Also turns on colour.")
parser.add_argument("--cigar", action="store_true", help="Adds CIGAR after name.")
parser.add_argument("--speedy", action="store_true", help="Skips the 'difference-from-consensus' part and outputs stream of fastas, followed by consensus.")
parser.add_argument("--skip", action="store_true", help="Skips the raw output fasta stream.")

#TO DO:
parser.add_argument("-o", "--consensus", action="store_true", help="Removes consensus sequence part.")
parser.add_argument("-r", "--reference", action="store_true", help="Removes reference sequence part.")
parser.add_argument("-m", "--multifasta", action="store_true", help="Removes multifasta part.")
args=parser.parse_args()

#########################
#     DEBUG OPTIONS     #
#########################

debug=0 #Debug messages
gaps_without_inserts=0 #Gap every 10 bases. Note: INSERTED bases do not count towards base count.
gaps_with_inserts=0 #Gap every 10 bases. Note: INSERTED bases DO count towards base count.
full=0 #'Full' alignment also includes matches with x<10 bases in common with reference. 
colour=0 #To do: switch off/on colours for insertions etc
positioning=0
find_name=0

if args.debug: debug=1
if args.gaps: gaps_with_inserts=1
if args.full: full=1
if args.colour: colour=1
if args.positioning: positioning=1
if args.name: find_name=args.name

################################
##   GET REFERENCE SEQUENCE   ##
################################

reference_file=args.reference_file
reference_name=args.reference_name
matchObject=re.match('(.*):(\d+)-(\d+)', reference_name)
if debug:
  print(args)
  print("reference_name = '{}' ".format(reference_name))
  print("match = '{} {} {}'".format(matchObject.group(1), matchObject.group(2), matchObject.group(3)))
reference_scaffold=matchObject.group(1)
reference_start=int(matchObject.group(2))
reference_end=int(matchObject.group(3))
if reference_end<reference_start: 
  reference_start, reference_end = reference_end, reference_start
if debug:
  print("Reference name: '{}' \nreference file: {}".format(reference_name, reference_file))
#reference_start-=1 #Not sure whether this is the right option
file_object=args.reference_file

for line in file_object: 
  if line.find(reference_scaffold)>-1:
    sequence = file_object.next().rstrip()
    reference = sequence[reference_start:reference_end]
    reference_length=len(reference)
    if debug: 
      print("Found reference scaffold: '{}'".format(line).rstrip())
      print("Full sequence: {}".format(sequence))
      print("Partial sequence:{}-{}: '{}'".format(reference_start, reference_end, reference))
    break
file_object.close()

########################################
######   GET EXTERNAL CONSENSUS   ######
########################################

if args.consensus_file:
  file_object=args.consensus_file
  reference_name=args.reference_name
  matchObject=re.match('(.*):(\d+)-(\d+)', reference_name)
  reference_scaffold=matchObject.group(1)
  ext_reference_start=int(matchObject.group(2))
  ext_reference_end=int(matchObject.group(3))

        
  for line in file_object: 
    if line.find(reference_scaffold)>-1:
      external_consensus = file_object.next().rstrip()
      external_consensus_exon = external_consensus[ext_reference_start:ext_reference_end]
      if debug: 
        print("Found external consensus: '{}'".format(line).rstrip())
        print("Full sequence: {}".format(external_consensus))
        print("Partial Sequence:{}-{}: '{}'".format(ext_reference_start, ext_reference_end, external_consensus_exon))
      break
  file_object.close()

########################################
######   GET EXTERNAL START/END   ######
########################################

#if args.blast_file:
#  file_object=args.blast_file
#  reference_name=args.reference_name  
#  for line in file_object: 
#    if line.find(reference_name)>-1:
#      line_split = split('\t', line)
#      blast_start = line_split[0]
#      matchObject=re.match('(.*)\s([.\d]+)\s([e-.\d]+)\s([.\d]+)\s([.\d]+)\s([.\d]+)\s(.*)', line)

######################################
######   SET GLOBAL VARIABLES   ######
######################################

#Make empty list of strings.
output=[]
#Make list of sequences
sequences=[]
startpoints =[]
endpoints=[]
cigars=[]
names=[]

names_dict={} #I don't know why I can't get output order otherwise... 
sequence_nr=0 # This seems so inefficient

insertions={} #dictionary 
aligned_bases={} #dictionary of finished alignments
aligned_colours={}#dictionary of finished alignment colours
count_bases={}#dictionary counting letters printed in particular sequence
colours={'yellow':'\033[93m', 'red':'\033[91m', 'blue':'\033[94m', 'green':'\033[92m', 'cyan':'\033[96m', 'magenta':'\033[95m', 'end':'\033[0m' }
col_end='\033[0m'
if not colour: 
  colours={'yellow':'', 'red':'', 'blue':'', 'green':'', 'cyan':'', 'magenta':'', 'end':'' }
  col_end=''
base_count=0

max_pos=reference_end
min_pos=reference_start
pattern = re.compile('([MIDNSHPX=])')

#Checking values
#print("length: {}, start: {}, end: {}".format(reference_length, reference_start, reference_end))

####################################
####         FUNCTIONS      #################################################
####################################

##################
### FUNCTION PRINT_BASE

def print_base(txt, txt_colour=0, base=1, name=0, no_print=0, consensus=1):
  global sequence_nr
  txt=str(txt)
  if not colour or not txt_colour:
    col=''
  else: 
    col=colours[txt_colour]
  if name not in aligned_bases:
    aligned_bases[name]=""
  if args.skip: no_print=1

  if base==0:
    if no_print: return 
    print(col, end='') #start color printing
    print(txt, end='')
    print(col_end, end='') #end color printing
    return
	  			

  if not no_print: print(col, end='') #start color printing
# Record data for making of consensus sequence:
  for i in txt:
    aligned_bases[name]+=i
    if not no_print: print(i, end='')
    try:
      print_base.base_count += 1
    except AttributeError:
      print_base.base_count = 0
    pos = len(aligned_bases[name])

    if consensus:
      try:
        aligned_bases_count[pos][i]+=1
      except AttributeError:
        warnings.warn("Unexpected base: {}".format(i))
        aligned_bases_count[pos][i]=1
      except KeyError: 
        warnings.warn("Unexpected pos: {} (base:'{}')".format(pos, i))
        aligned_bases_count[pos]={'A':0,'T':0,'G':0,'C':0,'N':0,'-':0, '+':0}
        aligned_bases_count[pos][i]=1

  if not no_print: print(col_end, end='')
  return
print_base.base_count=0


##################
### FUNCTION PRINT_SEQUENCE

def print_sequence(name, bases, start, end, cigar=0, no_print=0, consensus=1):
  global sequence_nr
  if args.skip: no_print=1
  if debug:
    print("print_sequence(name={}, start={}, end={}, cigar={}, bases={})".format(name, start, end, cigar, bases))
  if cigar and cigar.count('M')==len(cigar): 
    if debug: 
      print("CIGAR {} is all Ms".format(cigar))
    cigar=0
  add_cigar=''
  if args.cigar: 
    add_cigar+=' ' + str(cigar)
  if name not in names_dict.values():			
    names_dict[sequence_nr]=name #Making sure input order = output order.
    sequence_nr+=1
    #print("sequence_nr={}".format(sequence_nr))  
# Start printing
  aligned_bases[name]=""
  aligned_colours[name]=""
  print_base(">{}{}\n".format(name, add_cigar), base=0, name=name, consensus=consensus)

  insertions_done=[]
  pos=-1
  cigar_pos=-1
  j=min_pos-1
  while j<max_pos-1:
    j+=1

#Make counting easier for debugger: 
    if gaps_with_inserts and j!=min_pos:
      if (j+min_pos) % 10 == 0: print(' ', end='')

#Show position in alignment, for error finding
    if positioning:
      print_base(j , 'green', base=0, name=name, consensus=consensus)

#Account for beginning and end:

#'-' if before current sequence.
    if (start>j):  
      if (j in insertions.keys()):
        print_base(insertions[j]*'-' , 'yellow', name=name, consensus=consensus)
      print_base('-', name=name, consensus=consensus)  
      continue

#'-' if after current sequence.
    if ((1+j)>end):  
      if (j in insertions.keys()):
        print_base(insertions[j]*'-' , 'yellow', name=name, consensus=consensus)
      print_base('-', name=name, consensus=consensus)  
      continue

#Print letter from current sequence

#Establish position
    pos+=1
    cigar_pos+=1

#Take care of index error, if any    
    if pos>(len(bases)-1):
      warnings.warn("position is out of index for current sequence! \nlen:{}, pos:{}, j:{}, end:{}, calculated end:{}".format(reference_length, pos, j, end, start+len(str(cigar))))
      print_base('*', 'yellow', name=name, consensus=consensus)
      continue

#Establish current letter

    if not cigar:
      c='M'
    else: 
      c=cigar[cigar_pos]
      if args.cigar and not c=='I':
        print_base(str(cigar_pos)+c, 'cyan', name=name, base=0, consensus=consensus)



      #'-' if current one has an 'I' in the cigar string. 
    if (c=='I'):
      if j in insertions.keys():
        ins_len=insertions[j]
        insertions_done.append(j)
      else: 
        w = "\n\033[93mCIGAR doesn't match insertions record in {}! \n(Ref: {}) {}/{}\n".format(name, reference_name, j, len(names))
        w+= "CIGAR STRING: "+cigar+"\n"
        w+= "DNA SEQUENCE: "+bases+"\n"
        w+= "j:{}, pos:{}, calculated j:{}, start={}, end={}\n{} Insertions:{}".format(j, pos, pos+start, start, end, col_end, insertions)
        warnings.warn(w)
        ins_len=0
      insert=""
      cigar_pos-=1
      pos-=1
      while cigar[cigar_pos+1]=='I': #-->at least once, or this if-construct wouldn't trigger 
        pos+=1
        cigar_pos+=1
        if args.cigar:
          print_base(str(cigar_pos) + cigar[cigar_pos], 'cyan', name=name, base=0, consensus=consensus)
        if positioning:
          print_base(pos, 'blue', name=name, base=0, consensus=consensus )
        print_base(bases[pos], 'yellow', name=name, consensus=consensus )
        ins_len-=1
      print_base('-'*ins_len, 'yellow', name=name, consensus=consensus )
      j-=1
      continue

#'-' if any sequence - but not current one - has an 'I' in the cigar string on this position.
    if (j in insertions.keys() and j not in insertions_done):
      print_base(insertions[j]*'+','yellow', name=name, consensus=consensus)
           
#'-' if 'D' in the cigar string.
    if not (cigar==0) and (c=='D'): 
      print_base('-' , 'red', name=name, consensus=consensus)
      pos-=1
      continue  

#[letter], if 'M' OR '=' OR 'X' in the cigar string. 
    if pos>len(bases)-1: 
      warnings.warn("Position is out of index!\nPos: {}, length: {}, cigar_length: {}, c: {}".format(pos, len(bases), len(cigar), c))
      continue
    if (c=='M' or c=='=' or c=='X'):
      if positioning: 
        print_base(pos, 'blue', base=0, consensus=consensus)
      print_base(bases[pos], name=name, consensus=consensus)
      continue

#'-' if 'S' in the cigar string. 
    if (c=='S'): 
      #if positioning: 
      #  print_base(pos, 'magenta', base=0, consensus=consensus)
      #print_base('-', 'magenta', name=name, consensus=consensus)
      continue

#Otherwise, report error. 
    warnings.warn("encountered unhandled CIGAR character: '{}' in cigar string: '{}'".format(c, cigar))
  else: #Add a rowbreak after sequence is finished
    if not no_print: print()
    print_base.base_count=0

######################
########    AMBIGUITY RESOLVER    

def choose_base(pos, dic):
  highest=0
  total=0
  con_base='N'
  #print("pos: {}".format(pos))
  
  for base, count in dic[pos].iteritems(): 
    #print("{}:{}, ".format(base, count), end='')
    if base=='-': continue
    total+=count
  limit = 0.25*total
  d = dict((k, v) for k, v in dic[pos].items() if v >0 and v >= limit and k in "ACGT") 
  length = len(d)
  if debug:  
    if length>1:
      print("{}{}{}".format(colours['blue'], d, col_end), end='')
    else:
      print("* {}*".format(d), end='')
  if length>3: return 'N' 
  if length==3:
    if not 'A' in d: return 'B' 
    if not 'T' in d: return 'V'
    if not 'G' in d: return 'H' 
    if not 'C' in d: return 'D'  
  if 'A' in d: 
    if 'T' in d: return 'W'
    if 'G' in d: return 'Y'
    if 'C' in d: return 'K'
    return 'A'
  if 'T' in d: 
    if 'G' in d: return 'M'
    if 'C' in d: return 'R'
    return 'T'
  if 'G' in d: 
    if 'C' in d: return 'S'
    return 'G' 
  if 'C' in d: return 'C'
  return '-'
    

####################################
########    PARSING LOOP    ##################################################
####################################

#Get file from standard input OR from -i option.
seq_file = sys.stdin
if args.input: seq_file = args.input

for line in seq_file:
  if len(line)<3:
    continue
  list = line.rstrip().split('\t')
  name = list[0]
  try: 
    if int(list[1]) & 16: # find SAM bitwise flag marking reverse complement 
      name += 'R'
  except IndexError: 
    warnings.warn("Index Error, line ='{}'".format(line))
  if debug: name += " " + list[5]
  start = int(list[3])-1 
  cigar = pattern.split(list[5] )[:-1]
  sequence = list[9]
  length = len(sequence)
  end = start + length -1
  
#  if debug: 
#    print("Perusing sequence "+name)

# if looking for specific match
  if find_name and name.find(find_name)==-1:
    continue

#Error-checking:
  if len(list[5])<2:
    #warnings.warn("invalid cigar string")
    continue

#Remove read if there's too short overlap with reference sequence
  if not full:
    if (end-10)<reference_start: continue
    if (start+10)>reference_end: continue

#Populate lists and update global vars
  names.append(name)
  startpoints.append(start) 
  min_pos=min(min_pos, start)
  max_pos=max(max_pos, len(sequence)+start)
  sequences.append(sequence)

#Parse CIGAR into variables. String of letters?
  cigar_string=""
  for i, c in zip(cigar[::2], cigar[1::2]): 
    cigar_string += int(i)*c 
  cigars.append(cigar_string)
  endpoints.append(start+len(cigar_string.translate(None,'ID')))

#Find sequence number of any I in string (remember to add start number).
  pos=0
  noninsert_length=0

  while pos<len(cigar_string): 
    c=cigar_string[pos]
    if c=='I': 
      nr=0
      insert_pos=noninsert_length+start
      pos-=1
      while cigar_string[pos+1]=='I': 
        nr+=1
        pos+=1        
      if insert_pos in insertions.keys():
        insertions[insert_pos]=max(nr, insertions[insert_pos])
      else: insertions[insert_pos]=nr
      if debug:
        print("Insertion found in {} at c:{} = pos:{}, nr:{}".format(name, insert_pos, pos, nr))
    else: 
      noninsert_length+=1
    pos+=1
if debug: 
  print("done finding insertions")

#If there are no reads, exit

if len(names)==0:
#  warnings.warn("No reads for {}, skipping.".format(reference_name))
  exit()

####################################
#   Prepare consensus dictionary   #
####################################

# Find total length of insertions
total_insertion_count=0
#print(insertions)
for insert in insertions.values():
  total_insertion_count+=insert

# Define aligned_bases_count 
aligned_bases_count={}
length=max_pos-min_pos+total_insertion_count+1
if debug:
  warnings.warn("name={}, reads={}, length={}, insertions={}".format(reference_name, len(names), length, total_insertion_count))

for i in range(0, length):
  aligned_bases_count[i]={'A':0,'T':0,'G':0,'C':0,'N':0,'-':0, '+':0}


#################################################
################  PROCESSING  ###################
#################################################

# PROCESS REFERENCE
print_sequence(reference_name, bases=reference, start=reference_start, end=reference_end, consensus=0)

# PROCESS EXTERNAL CONSENSUS
#if args.consensus_file:
#  print_sequence("External_Consensus", bases=external_consensus_, start=reference_start, end=reference_end)

# PROCESS REST OF SEQUENCES
for i, name in enumerate(names):

  if debug:
    print("{}\n{}".format(sequences[i], cigars[i]))
    print("start:{}, end:{}, length:{}, covers: {} ".format(startpoints[i], endpoints[i], endpoints[i] - startpoints[i], len(sequences[i]) ))
    print("i={}".format(i))  
  print_sequence(name, start=startpoints[i], end=endpoints[i], bases=sequences[i], cigar=cigars[i])

######################################################################
###############    OUTPUT DICTIONARY OF ALIGNMENTS     ###############
######################################################################
 
print()
if debug and positioning: 
  print("Exiting early")
  sys.exit()

##################
### FUNCTION PRINT_MODIFIED_SEQUENCE

def print_modified_sequence(name, seq, color=1):
  print("{}>{}{}".format(colours['blue'], name, col_end))
  if debug: 
    if name in names:
      index = names.index(name)
      print("DNA Sequence: "+sequences[index])
      print("CIGAR string: "+cigars[index])
      print("start:{}, end:{}, covers:{}, length: {} ".format(startpoints[index], endpoints[index], endpoints[index] - startpoints[index], len(sequences[index]) ))
  if len(seq)>len(consensus_sequence):
    warnings.warn("Sequence is longer than consensus sequence ({}>{}): \n{}\n{}\n,len(aligned_bases_count)={}".format(len(seq), len(consensus_sequence), seq, consensus_sequence,len(aligned_bases_count)))
  for pos, base in enumerate(seq):
    if gaps_with_inserts and (pos % 10)==0 and pos!=0: print(' ', end='')
    if pos>=len(consensus_sequence):
      continue
#    if debug:
#      print("{}{}{}".format(colours['blue'], consensus_sequence[pos-1], col_end), end='')
    if color==0:
      print(base, end='')
    elif base=='+':
      print("{}{}{}".format(colours['yellow'], '-', col_end), end='')
    elif consensus_sequence[pos]=='+':
      print("{}{}{}".format(colours['yellow'], base, col_end), end='')
    elif base=='-':    
      print(base, end='')
    elif base==consensus_sequence[pos]:    
      print(base, end='')
    else:
      print("{}{}{}".format(colours['green'], base, col_end), end='')
  print()

#Creating count of bases in given position
#aligned_bases_count={}
#for i, seq in aligned_bases.iteritems():
#  for j, base in enumerate(seq):  
#    if not j in aligned_bases_count.keys():  
#      aligned_bases_count[j]={'A':0,'T':0,'G':0,'C':0,'N':0,'-':0}
#    if not base in aligned_bases_count[j].keys(): 
#      aligned_bases_count[j][base]=0
#    aligned_bases_count[j][base]+=1


#Using count to create consensus sequence
consensus_sequence=""
con_base=""

for pos in range(1, len(aligned_bases_count)):
  con_base=choose_base(pos, aligned_bases_count)
  consensus_sequence+=con_base
  #print(con_base, end='') #REMOVE
if len(consensus_sequence) != (max_pos-min_pos + total_insertion_count): 
  warnings.warn("{} consensus sequence is {} bases long (should be {} bases.), len(aligned_bases_count)={}".format(reference_name, len(consensus_sequence), (max_pos-min_pos + total_insertion_count), len(aligned_bases_count)))

#Printing consensus sequence
print_modified_sequence("Consensus_Sequence", consensus_sequence)

#Printing cut consensus sequence 
cut_consensus=consensus_sequence

print(colours["blue"] + ">Cut_consensus"+ col_end)
i=-1
for j, base in enumerate(cut_consensus):
  if gaps_with_inserts and j>0 and j % 10 == 0: print(' ', end='')
  if j<(reference_start - min_pos):
    print('-', end='')
    continue
  if base=='-' or base =='+': #Catch insertions. 
    print('-', end='')
    continue
  if j>(reference_end - min_pos + total_insertion_count -1):
    print('-', end='')
    continue
  if aligned_bases[reference_name][j]=='-': #This REALLY ought not be needed. And yet it is.
    print('-', end='')
    continue
  try:
    print(cut_consensus[j], end='') 
  except IndexError: 
    warnings.warn("Error in {}, j={}. i={} len(cut_consensus)={}, len(consensus)={}, reference_start={}, reference_end={}, (reference_end - reference_start)={}, min_pos={}, max_pos={}, total_length={}".format(reference_name, j, i, len(cut_consensus), len(consensus_sequence), reference_start, reference_end, (reference_end - reference_start), min_pos, max_pos, (max_pos - min_pos) ))
    continue
print()

#Printing aligned reference_sequence
print(">Aligned reference")
print(aligned_bases[reference_name])

#Printing original reference_sequence
#print(">Original_reference")
#print(reference)

#Printing EXTERNAL consensus sequence 
if args.consensus_file:
  print(colours["blue"] + ">External_Consensus"+ col_end)
  i=-1
  for j, base in enumerate(aligned_bases[reference_name]):
    if gaps_with_inserts and j>0 and j % 10 == 0: print(' ', end='')
    if base=='-' or base =='+': 
      print('-', end='')
      continue
    i+=1
    if i>len(external_consensus_exon)-1: 
      print('-', end='')
      continue
    try:
      print(external_consensus_exon[i], end='') 
    except IndexError: 
      warnings.warn("Error in {}, external_reference print at position {}".format(reference_name, i))
    continue
print()

if args.speedy: 
    sys.exit()

#printing sequences, coloring bases depending on whether they differ from consensus.
sequence_nr=len(names_dict)
#print("sequence_nr is {}".format(sequence_nr))
#print("names_dict is {}".format(len(names_dict)))
#print("aligned_bases is {}".format(len(aligned_bases)))
for i in range(1,sequence_nr):
  #print("i is {}".format(i))
  seq = aligned_bases[names_dict[i]] 
  print_modified_sequence(names_dict[i], seq)

