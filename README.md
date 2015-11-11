# TelMe
Telomere counting from NGS data

## Explanation of the output file


- SAMPLE	name of the sample
- TARGET	path to the bed file
- CMD		full command-line executed

These are the headers for data based on individual RG-tags, one line per RG-tag								
- RG		rg-tag as in the bam file
- ALL 		# of total reads examined
- ONTARGET	# of reads deemed on target
- TELOMERE 	# of reads deemed telomeric
- CHRM		# of reads mapped to chrM
- TEL-UNMAPPED # of reads deemed telomeric that were also unmapped
- CNT_TEL_LEN 	# of telomere repeats counted
- RTPMO-w/offt 	Relative-telomere (reads) per million off-target reads
- RTPKM-w/chrM	Relative-telomere per thousand chrM reads
- RTPMO-w/offtCNT 	Relative-telomere (length) per million off-target reads
- FRACTION_DUPS 	fraction of duplication
- TELOMERE_SEQ_FOUND ordered list of relative fraction of all the subtelomeric repeats identified 

At the last line, there will be a summary over all the RG-tags in the bam file

- Sample name
- MEAN+-SD represent the values in the adjacent columns for each calculated value for each of RTPMO-w/offt, RTPKM-w/chrM, RTPMO-w/offtCNT
- Mean of all 
- SD of all 
- ONTARGET-RATE 	fraction of all the reads deemed on target
- DUP-RATE 	fraction of all the reads marked as duplicates

And the last value is the ordered-list of  subtelomeric sequences with relative fractions	
