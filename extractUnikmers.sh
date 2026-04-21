#!/bin/bash
# ------- PRE-REQUISITES
# Install the following packages in your conda environment before running this script
# jellyfish (channel: bioconda), xlsxwriter (channel: conda-forge)
# command to install them:
# conda install -c bioconda jellyfish
# conda install -c conda-forge xlsxwriter

# ------- SPECIFICATION
READ_FILE="reads.fasta"  # REPLACE THIS with the path to your input reads file
k=21

# Output files will be generated in the current directory
OUT_COUNT="k${k}_count.jf"
OUT_HISTO="k${k}_counts.histo"
OUT_XCEL="k${k}_distribution.xlsx"
UNIKMERS="k${k}_unikmers.kmers"

# the count and histo commands together generate the frequency distribution of kmers
# in the $OUT_HISTO file 
jellyfish count -m $k -C -o $OUT_COUNT -c 3 -s 10000000 --disk -t 16 $READ_FILE
jellyfish histo -o $OUT_HISTO -v $OUT_COUNT

# the python script generates an excel (.xlsx) file to visualize the distribution
# using this excel file, you can find out the [mean - 3*sd, mean + 3*sd] window
# containing the expected unikmers
# the python script also outputs the calculated mean, sd, lower bound and upper bound
# values from the input distribution file
OUTPUT=$(python3 ./freq_distribution_kmers.py $OUT_HISTO $OUT_XCEL)

echo $OUTPUT

# extract the lower and upper bound values from the output
lower=$(echo "$OUTPUT" | sed -n '2p' | cut -d':' -f2 | xargs)
upper=$(echo "$OUTPUT" | sed -n '3p' | cut -d':' -f2 | xargs)

# extract the unikmers using the lower and upper bound values
jellyfish dump -c -t -L $lower -U $upper $OUT_COUNT | awk '{print $1}' > $UNIKMERS
