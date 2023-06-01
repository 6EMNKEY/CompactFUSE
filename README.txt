#	     _____ ____  __  __ _____        _____ _______ ______ _    _  _____ ______ 
#	    / ____/ __ \|  \/  |  __ \ /\   / ____|__   __|  ____| |  | |/ ____|  ____|
#	   | |   | |  | | \  / | |__) /  \ | |       | |  | |__  | |  | | (___ | |__   
#	   | |   | |  | | |\/| |  ___/ /\ \| |       | |  |  __| | |  | |\___ \|  __|  
#	   | |___| |__| | |  | | |  / ____ \ |____   | |  | |    | |__| |____) | |____ 
#	    \_____\____/|_|  |_|_| /_/    \_\_____|  |_|  |_|     \____/|_____/|______|
#   	                                                                             
                                                                             
##################################################################################################
				CompactFuse by Nicolas Ares
			   in collaboration with Sequentia Biotech SL
	A pipeline for the detection of fusion genes in cancer, along with the prioritization
	of the given calls using the softwares of FusionCatcher, StarFusion, JAFFA and Arriba
	done by:
		
		- ndaniel(https://github.com/ndaniel/fusioncatcher) 
		- Trinity CTAT Project (https://github.com/STAR-Fusion/STAR-Fusion)
		- Oshlack( https://github.com/Oshlack/JAFFA )
		- Shurig(https://github.com/suhrig/arriba) 


##################################################################################################

#BASIC FUNCTIONING OF THE PIPELINE 
	
	0.5. Trimming: Optional trimming of the reads is possible and sometimes can greatly improve
	the results.
	1. Fusion Calling : The Choosen algorithms from the 4 avaliable calling softwares will 
	produce a set of calls.
	2. Standarization : The Calls are standarized and intersected into a table.
	3. Prioritization : A score is assigned to each of the calls so as to reduce noise and 
	false positives.
		3.1 Mitelman database check : Calls are cross validated using Mitelman Database of
		curated fusion genes for a given tissue and Disease.
		3.2 Fusion gene breakpoints are checked for being in repetitive zones, in pseudo 
		genes or being transposones. Using UCSC, biomart and Dfam.
		3.3 Mapping analisis is preformed: 
			3.3.1 The reads in the bam are checked in the breakpoints for possible 
			false positives in the spanning reads (reads that span in both genes at
			the same time.
			3.3.2 Analisis is done for the spanning pairs ( reads that have one of the
			mates at one gene and the other mate at another gene)
	4. Visualization: A shinny app is provided for a better viewing of the results, including
	a visualization of reads at breakpoint, a plot of transposone prediction, and a VCF search
	feature, that allows for a better search in a structural VCF file.  


#INSTALLATION

	Inside bin folder:
	
		bin/CompactFuse_install.sh
		
	-tools flag can be provided with the values a(arriba),j(jaffa),f(fusion catcher),s(starfusion)
		
		bin/CompactFuse_install.sh -tools aj would install arriba and jaffa.
	
	if -tools is not provided all of the callers will be installed 

#INSTALLING JAFFA

	At the end of the installation Jaffa outputs the following:
	
Checking that all required tools were installed:
bpipe looks like it has been installed
velveth looks like it has been installed
velvetg looks like it has been installed
oases looks like it has been installed
trimmomatic looks like it has been installed
samtools looks like it has been installed
bowtie2 looks like it has been installed
blat looks like it has been installed
dedupe looks like it has been installed
reformat looks like it has been installed
extract_seq_from_fasta looks like it has been installed
make_simple_read_table looks like it has been installed
blastn looks like it has been installed
minimap2 looks like it has been installed
process_transcriptome_align_table looks like it has been installed
make_3_gene_fusion_table looks like it has been installed
**********************************************************
All commands installed successfully!

If one of them (XXX) did not install succesfully 
	- check if in tools/bin/ there is the binary XXX 
	- If it is there check that in the file tools.groovy the path points to the installation
	- If not then install it from the url in install_linux64.sh
	- move it to tools/bin/ and add the path to tools.groovy
	
#INSTALLING FUSIONCATCHER

	Fusion catcher installs from the script tools/install_tools.sh
	A correct installation would have the following inside tools directory
	
	1.2-r101c.tar.gz                  fatotwobit
	2.7.2b.tar.gz                     install_tools.sh
	BBMap_38.44                       liftover
	BBMap_38.44.tar.gz                parallel
	README                            parallel-20201222
	STAR-2.7.2b                       parallel-20201222.tar.bz2
	bbmap                             picard
	blat                              pigz
	bowtie                            pigz-2.4.tar.gz
	bowtie-1.2.3-linux-x86_64         seqtk
	bowtie-1.2.3-linux-x86_64.zip     seqtk-1.2-r101c
	bowtie2                           sratoolkit
	bowtie2-2.3.5.1-linux-x86_64      sratoolkit.2.9.6-centos_linux64
	bowtie2-2.3.5.1-linux-x86_64.zip  sratoolkit.2.9.6-centos_linux64.tar.gz
	fastqtk                           star
	fastqtk-0.27                      v0.27.tar.gz
	
	If anything went wrong during the installation we recomend manual instalation copying the commands in 
	instal_tools.sh.
	
#INSTALLING STARFUSION

	Star uses a make file for installing the binaries, It is important to know that it will call many other
	make files that you might want to check in case of error.
	
	A usual error prompts when R libraries are missing, install said librarie and rerrun the command:	
		bin/CompactFuse_install.sh -tools s

#INSTALLING ARRIBA
	
	Arriba is the easiest of the installations, for checking if it was correctlly installed you can use the
	test reads in the test folder.
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

#Basic usage:
# - Select tools (default all).
# - tools a, j, s, f
# - Select analisis (default complete)
# - analisis complete, prioritization, calling
# If only prioritization is selected: A .txt file has to be generated containing the names of all samples.
# - Select File
# - reads Abs/File (File Contains a subfile for each fasta File/S1/S1_1 File/S1/S1_2)
# - 
#




#Tool installer.
# -tools a,f,s,j default all
