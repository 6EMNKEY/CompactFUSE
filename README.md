# CompactFUSE
Pipeline for integration of diferent gene fusion detection tools, providing a user frendly environement.

### Tools used :

- JAFFA
- Chimerascan
- FusionCatcher
- FusionBloom

### break_point_clean_analisys.py

The script does the following

a) by making calls to the BAM alignment file generated by Arriba for a given sample using samtools 1.13 [96] command
samtools view and passing the fusion breakpoint as an argument; this is done
for the two genes:

### samtools view <PATH to .bam> chrN:base-base+1 > filename.txt

This command returns a file with the alignment information of all fragments that span across the region that is selected (from base to base+1 of chrN). 

b) Since we were only looking for spanning reads and not spanning pairs or their mates, we select only the reads that start before the break point and end
after the break point. 

c) Then the genomic sequence from 180 bp before breakpoint to 180 bp after breakpoint is retrieved for each gene from the assembly of the genome, and separated into post breakpoint sequence and pre breakpoint sequence. 
d) The read alignment information is used to perform the same operation on the reads, separating them into the pre and post breakpoint regions. 

e) These breakpoint regions are aligned using a smith and waterman based algorithm to both their non fusion version (the continuation of the gene)
12 and their fusion version ( the other gene sequence, from the breakpoint to end or start). 

*) Also we need to point out that the breakpoints region reads in the 3’ side of the 5’ fusion gene have to be reversed along with the sequences of
comparison for a proper alignment.

  
