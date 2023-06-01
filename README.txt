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
