import subprocess as sub
import time 
import pandas as pd
import os.path
import numpy as np
from termcolor import colored


#SMITH AND WATERMAN ALORITHM 

# Implementing the Smith Waterman local alignment
def smith_waterman(seq1, seq2):
    # Generating the empty matrices for storing scores and tracing
    row = len(seq1) + 1
    col = len(seq2) + 1
    matrix = np.zeros(shape=(row, col), dtype=int)  
    tracing_matrix = np.zeros(shape=(row, col), dtype=int)  
    
    # Initialising the variables to find the highest scoring cell
    max_score = -1
    max_index = (-1, -1)
    
    # Calculating the scores for all cells in the matrix
    for i in range(1, row):
        for j in range(1, col):
            # Calculating the diagonal score (match score)
            match_value = 1 if seq1[i - 1] == seq2[j - 1] else -1
            diagonal_score = matrix[i - 1, j - 1] + match_value
            
            # Calculating the vertical gap score
            vertical_score = matrix[i - 1, j] -1
            
            # Calculating the horizontal gap score
            horizontal_score = matrix[i, j - 1] -1
            
            # Taking the highest score 
            matrix[i, j] = max(0, diagonal_score, vertical_score, horizontal_score)
            
            # Tracking where the cell's value is coming from    
            if matrix[i, j] == 0: 
                tracing_matrix[i, j] = 0
                
            elif matrix[i, j] == horizontal_score: 
                tracing_matrix[i, j] = 1
                
            elif matrix[i, j] == vertical_score: 
                tracing_matrix[i, j] = 2
                
            elif matrix[i, j] == diagonal_score: 
                tracing_matrix[i, j] = 3 
                
            # Tracking the cell with the maximum score
            if matrix[i, j] >= max_score:
                max_index = (i,j)
                max_score = matrix[i, j]
    
    # Initialising the variables for tracing
    aligned_seq1 = ""
    aligned_seq2 = ""   
    current_aligned_seq1 = ""   
    current_aligned_seq2 = ""  
    (max_i, max_j) = max_index
    
    # Tracing and computing the pathway with the local alignment
    while tracing_matrix[max_i, max_j] != 0:
        if tracing_matrix[max_i, max_j] == 3:
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = seq2[max_j - 1]
            max_i = max_i - 1
            max_j = max_j - 1
            
        elif tracing_matrix[max_i, max_j] == 2:
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = '-'
            max_i = max_i - 1    
            
        elif tracing_matrix[max_i, max_j] == 1:
            current_aligned_seq1 = '-'
            current_aligned_seq2 = seq2[max_j - 1]
            max_j = max_j - 1
            
        aligned_seq1 = aligned_seq1 + current_aligned_seq1
        aligned_seq2 = aligned_seq2 + current_aligned_seq2
    
    # Reversing the order of the sequences
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]
    
    return aligned_seq1, aligned_seq2


#SORTING AND INDEXING BAM IF NOT DONE BEFORE
def sort_indexBamIfNecessary():
    if not os.path.isfile("AlignedS.out.bam"):
        sub.run("samtools sort edgreen/Arriba/MCF_7.SRR064286/Aligned.out.bam -o AlignedS.out.bam" , shell = True)
    if not os.path.isfile("AlignedS.out.bam.bai"):
        sub.run("samtools index AlignedS.out.bam" , shell = True)

#USING SAMTOOLS FOR THE SELECTION OF SPANNING FRAGMENTS AT A GIVEN POSITION
def make_breakpointSpanningFragFile(df, i):
    target = df["Fusion"][i]
    if "," in target:
        target = str(i)
    break1 = df["base1"][i]
    break2 = df["base2"][i]
    ######################
    print("Analising fusion: ",target)
    print("[*] Computing samtools")
    ######################
    pos1 = str(df["chrom1"][i]) +":"+str(break1-100)+"-"+str(break1+1)
    pos2 = str(df["chrom2"][i]) +":"+str(break2-1)+"-"+str(break2+100)
    sub.run("samtools view AlignedS.out.bam "+ pos1 +" -o outpos1_"+str(target)+"_"+str(break1), shell = True)
    sub.run("samtools view AlignedS.out.bam "+ pos2 +" -o outpos2_"+str(target)+"_"+str(break2) , shell = True)

#SUBSELECTION OF ONLY THE SPANNING READS THAT FALL IN THAT POSITION
def subselect_SpanningReads(df,i):
    break1 = df["base1"][i]
    break2 = df["base2"][i]
    target = df["Fusion"][i]
    if "," in target:
        target = str(i)
    outpos1 = pd.read_csv("outpos1_"+str(target)+"_"+str(df["base1"][i]), sep="\t", header=None,on_bad_lines='skip')
    outpos2 = pd.read_csv("outpos2_"+str(target)+"_"+str(df["base2"][i]),sep="\t", header=None,on_bad_lines='skip')
    
    kept1 = []
    for row in outpos1.iterrows():
        if row[1][3] < break1 and row[1][3] + len(row[1][9]) > break1:
            kept1.append(list(row[1]))
    kept1 = pd.DataFrame(kept1)
    kept1.to_csv("outposK1_"+str(target)+"_"+str(df["base1"][i]) ,sep="\t")
    kept2 = []

    for row in outpos2.iterrows():
        comp = (row[1][3] + len(row[1][9])) 
        if row[1][3] < break2 and comp > break2:
            kept2.append(list(row[1]))
    kept2 = pd.DataFrame(kept2)
    kept2.to_csv("outposK2_"+str(target)+"_"+str(df["base2"][i]) ,sep="\t")

    #print("Pos1: "+str(len(kept1))+ " reads kept" +"\n"+"Pos2: "+str(len(kept2)) + " reads kept" )
    return kept1,kept2

#RETURNING RELATIVE POSITIONS FOR ALIGNMENT (VECTOR IS RETRIEVED FROM THE BAM)
def get_relativeCoordinatesVector(vector):
    minum = vector[0]
    for i in vector:
        if i < minum:
            minum = i
    choords = []
    for i in vector:
        choords.append(i-minum)
    return choords, minum

#GETING THE CHROMOSOME START POINTS OF THE PRIMARY ASSEMBLY
def get_ChromosomeStartingPoints():
    startings = sub.run("cat Homo_sapiens.GRCh38.dna.primary_assembly.fa | grep -in '>' ", capture_output=True, shell= True)
    startl = startings.stdout.decode().split()
    dicto = []
    for i in startl:
        if ">" in i:
            b = i.split(">")
            a = b[0]
            b = b[1]
            dicto.append([b,a.strip(":")])
    return dicto

#RETRIEVING 180BP LEFT AND RIGHT OF OUR BREAKPOINT FROM ASSEMBLY
def slice_Assembly_360bp(base,lls):
    c = 181
    t = 0
    selected = ""
    flag = False
    for i in lls[1:]:
        if c + len(i) > base and not flag:
            pos = base -c 
            selected += i[pos:]
            flag=True
        elif c + len(i) == base and not flag:
            selected +=i
            flag=True
            
        elif flag:
            selected+= i
            t+=1
            if t == 5:
                break
        c+= len(i)
    return selected

#RETRIEVING THE SEQUENCES FROM ASSEMBLY
def get_sequencesFromAssembly(minum,chromosome):
    #Technically I want to get from breakpoint :  180bp - breakpoint - 180 bp
    dicto = get_ChromosomeStartingPoints()


    poses = [0,0]
    nextc = False

    for i in dicto:
        if i[0] == chromosome:
            poses[0] = int(i[1])
            nextc = True
        elif nextc:
            poses[1] = int(i[1])
            nextc = False
    lenght = max(poses) - min(poses)

    lines1 = sub.run("tail -n+" +str(poses[0]) + " Homo_sapiens.GRCh38.dna.primary_assembly.fa" + "| head -n"+ str(lenght), capture_output=True, shell=True )
    lls = lines1.stdout.decode().split("\n")
    

    sequence = slice_Assembly_360bp(minum,lls)

    return sequence

#SEPARATING READS INTO THEIR PRE AND POST BREAKPOINT PARTS
def separateReadstoPostPreBreakpoint(readtable, df, i, n):
    breakp = df["base"+str(n)][i]
    aligned_reads = []
    choordvect, minima = get_relativeCoordinatesVector(readtable[3])
    print("Gene"+str(n)+ ": "+ df["Gene"+str(n)][i] +" \n")
    print("Starting position: ", minima)
    print("Max length between read starts: ",max(choordvect)-min(choordvect))
    for i in range(readtable.shape[0]):
        readpos =breakp  - (minima + int(choordvect[i]))
        read = readtable[9][i]
        prebreak = read[:readpos]
        postbreak = read[readpos:]
        aligned_reads.append(" "*int(choordvect[i]) + prebreak+ "***" + postbreak )
        #print(" "*int(choordvect[i]) + prebreak + "     "+ colored(postbreak,"yellow") )

    return aligned_reads, minima


#PRINTING THE ALIGNMENT TO TERMINAL (FOR EASY TOGGLING OFF)
def print_aligned_readsToChr(df,fusion_call,mina,sequenceGene, aligned_reads,n):
    bk1 = df["base"+ n][fusion_call]-mina+180
    print("Break from start:", bk1)
    print(sequenceGene[180:bk1]+ "***" + sequenceGene[bk1:bk1+50])
    for i in aligned_reads:
        print(i)
    return(sequenceGene[:bk1]+ "***" + sequenceGene[bk1:])

#FUNCTION THAT DOES THE MAIN RETRIEVAL OF GENES AND THEIR READS (IN AN ALIGNED WAY)
def RetrieveSpanningReadsAndAssemblyAligned(df,fusion_call):
    mainlist = []

    make_breakpointSpanningFragFile(df,fusion_call)
    isnotMithocondrial = True
    try:
        keptreads1, keptreads2 = subselect_SpanningReads(df,fusion_call)
        if len(keptreads1)> 0:
            
            aligned_reads1, min1 = separateReadstoPostPreBreakpoint(keptreads1,df,fusion_call, 1)
        else:
            min1 = df["base1"][fusion_call]
            print("Gene1: "+ df["Gene1"][fusion_call] +" \n")
            print("Starting position: ", min1)
            print("NO READS MACH THE CRITERIA")
            aligned_reads1 = []
        
        sequenceGene1 = get_sequencesFromAssembly(min1, df["chrom1"][fusion_call])
        #########################################################################3
        sequenceGene1 = print_aligned_readsToChr(df,fusion_call,min1,sequenceGene1, aligned_reads1, "1")
        #########################################################################3

        if len(keptreads2)> 0:
            aligned_reads2, min2 = separateReadstoPostPreBreakpoint(keptreads2,df,fusion_call, 2)
        else:
            min2 = df["base2"][fusion_call]
            print("Gene2: "+ df["Gene2"][fusion_call] +" \n")
            print("Starting position: ", min2)
            print("NO READS MATCH THE CRITERIA")
            aligned_reads2 = []

        sequenceGene2 = get_sequencesFromAssembly(min2,df["chrom2"][fusion_call])
        #########################################################################3
        sequenceGene2 =print_aligned_readsToChr(df,fusion_call,min2,sequenceGene2, aligned_reads2, "2")
        #########################################################################3

    except:
        print("Mithocondrial Chromosome hits can't be analised")
        isnotMithocondrial = False

    if isnotMithocondrial:
        return aligned_reads1, sequenceGene1, aligned_reads2, sequenceGene2
    else:
        return 0,0,0,0

#FUNCTION FOR GENERATING THE ALIGNMENT OF EACH READ TO SEQUENCE
def AlignReads_ToChromosome(Reads, SequenceGene1, SequenceGene2, Orientation):
    Alignment = []
    if Orientation == "R":
        for read in Reads:
            part_to_align ="***"+read.split("***")[1]
            Gene1postBreakpoint = "***"+SequenceGene1.split("***")[1]
            Gene2postBreakpoint = "***"+SequenceGene2.split("***")[1]
            Al1_1,Al1_2 = smith_waterman(part_to_align,Gene1postBreakpoint)
            Al2_1,Al2_2 = smith_waterman(part_to_align,Gene2postBreakpoint)
            Alignment.append([[Al1_1,Al1_2],[Al2_1,Al2_2]])
    if Orientation == "L":
        #IT IS REVERSED TO ACCOUNT FOR THE ALGORITHM ALIGNING FROM THE LEFT
        for read in Reads:
            part_to_align = "***"+read.split("***")[0][::-1]
            Gene2preBreakpoint = "***"+SequenceGene2.split("***")[0][::-1]
            Gene1preBreakpoint = "***"+SequenceGene1.split("***")[0][::-1]
            Al1_1,Al1_2 = smith_waterman(part_to_align,Gene1preBreakpoint)
            Al2_1,Al2_2 = smith_waterman(part_to_align,Gene2preBreakpoint)
            Alignment.append([[Al1_1,Al1_2],[Al2_1,Al2_2]])
    return Alignment

def Print_SW_Alignments(Al):
    c = 1
    print("Alignments to gene")
    for i in Al:
        print(f"Read Number {c}")
        print(i[0][0]) 
        print(i[0][1])
        c+=1
    c = 1
    print("Aligments to Fusion")
    for i in Al:
        print(f"Read Number {c}")
        print(i[1][0]) 
        print(i[1][1])
        c+=1

def main():
    df = pd.read_csv("Standarized_calls.csv")
    sort_indexBamIfNecessary()
    for fusion_call in range(len(df.index)):
        r1,s1,r2,s2 = RetrieveSpanningReadsAndAssemblyAligned(df, fusion_call)
        if r1 != 0:
            Al1 = AlignReads_ToChromosome(r1,s1,s2,"R") 
            print("ALIGNMENTS TO GENE1")
            Print_SW_Alignments(Al1)
            Al2 = AlignReads_ToChromosome(r2,s1,s2,"L") 
            print("ALIGNMENTS TO GENE2")
            Print_SW_Alignments(Al2)
            


if __name__ == "__main__":
    main()
