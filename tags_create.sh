#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate obitools
## USE OLIGOTAG TO OBTAIN A LIST OF TAGS USING THE DIFFERENT SPECIFICATIONS

DIFFERENCES=4
TAG_LENGTH=8
MIN_TAGS=40
HOMOPLOIDY=1
MAX_NUM_GC=5
OUTPUT='tags_list.txt'

oligotag -s ${TAG_LENGTH} -d ${DIFFERENCES} -f ${MIN_TAGS} -p ${HOMOPLOIDY} -g ${MAX_NUM_GC} > "${OUTPUT}"

conda deactivate
## SELECT FROM THE TAGS OBTAINED IN $OUTPUT, THE ONES NOT ENDING AS THE PRIMER_FORWARD OR PRIMER_REVERSE
### PYTHON 

python3 <<END 
from Bio.Seq import Seq

PRIMER_FORWARD_BACT=Seq("AACMGGATTAGATACCCKG")
PRIMER_REVERSE_BACT=Seq("ACGTCATCCCCACCTTCC")
PRIMER_FORWARD_FUNGI=Seq("GCATCGATGAAGAACGCAGC")
PRIMER_REVERSE_FUNGI=Seq("TCCTCCGCTTATTGATATGC")

bact_f=PRIMER_FORWARD_BACT[0]
bact_r=PRIMER_REVERSE_BACT[0]
fungi_f=PRIMER_FORWARD_FUNGI[0]
fungi_r=PRIMER_REVERSE_FUNGI[0]

file = 'tags_list.txt'

with open('tags_list.txt') as f:
    tags = f.readlines()
    tags = [tag.strip().upper() for tag in tags]
    forward_tags = [tag.strip() for tag in tags if tag[-1] != bact_f and tag[-1] != fungi_f]
    reverse_tags = [tag.strip() for tag in tags if tag[-1] != bact_r and tag[-1] != fungi_r]
with open ('tags_forward.txt', 'w') as tg:
    tg.writelines('\n'.join(forward_tags))        
with open ('tags_reverse.txt', 'w') as tr:
    tr.writelines('\n'.join(reverse_tags))


## CREATE SPACERS USING 1 OR 3 NUCLEOTIDES THAT ARE NOT THE BEGINING OF THE TAG

from itertools import product

nucleotides = ['A', 'T', 'G', 'C']
single_nucl_spacer = [str(Seq(nuc)) for nuc in nucleotides]

# Generate three-length nucleotide sequences
two_nucl_spacer = [str(Seq(''.join(seq))) for seq in product(nucleotides, repeat=2)]
three_nucl_spacer = [str(Seq(''.join(seq))) for seq in product(nucleotides, repeat=3)]
# Save spacers in a file

with open ('nucl_spacers.txt','w') as nucl_sp:
    nucl_sp.writelines('\n'.join(single_nucl_spacer))
    nucl_sp.writelines('\n'.join(two_nucl_spacer))
    nucl_sp.writelines('\n'.join(three_nucl_spacer))

END 
