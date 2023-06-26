import subprocess as sub
import time 
import pandas as pd
import os.path
import numpy as np
from termcolor import colored

#RECONSTRUCT STRINGS
def reconstruct_strings(substrings, G):
    c = 0
    if G == 1:
        for i in substrings:
            substrings[c] = (substrings[c].split("***")[1]).replace("W","")
            c+=1
    else:
        for i in substrings:
            substrings[c] = ((substrings[c].split("***")[0]).replace("W",""))[::-1].strip()
            c+=1
    # Step 1: Identify unique substrings
    unique_substrings = {}
    for substring in substrings:
        if substring not in unique_substrings:
            unique_substrings[substring] = 1
        else:
            unique_substrings[substring] += 1

    # Step 2: Build partial strings
    partial_strings = []
    for substring in unique_substrings:
        partial_strings.append(substring)

    # Step 3: Check for conflicts
    partial_strings_c = partial_strings

    for i in range(len(partial_strings)):
        for j in range(len(partial_strings)):
            if len(partial_strings[i]) < len(partial_strings[j]):
                if partial_strings[i].startswith(partial_strings[j]):
                    #print("Conflicting partial strings:", partial_strings[i], partial_strings[j])
                    partial_strings_c[j] = "n"
                elif partial_strings[j].startswith(partial_strings[i]):
                    #print("Conflicting partial strings:", partial_strings[i], partial_strings[j])
                    partial_strings_c[i] = "n"
    
    #REVERSE READS
    if G == 0:
        for i in range(len(partial_strings_c)):
            partial_strings_c[i] = partial_strings_c[i][::-1]

    # Step 4: Concatenate partial strings
    reconstructed_strings = []
    for s in range(len(partial_strings_c)):
        if partial_strings_c[s] != "n":
            reconstructed_strings.append(partial_strings_c[s])

    return reconstructed_strings


#SMITH AND WATERMAN ALORITHM 
def sequence_similarity(seq1,seq2,l1):
    c = 0
    endpos = 0
    for i in seq1[::-1]:
        if i == "-":
            endpos+=1
        else:
            break

    seq1 = seq1[:len(seq1)-endpos]
    seq2 = seq2[:len(seq2)-endpos]
    la = len(seq1)
    for i in range(la):
        if seq1[i] == seq2[i] and seq1!="-" and seq2!="-":
            c+=1
    if l1 == 0:
        l1 = 50
    return [(c/l1),c]

# Implementing the Smith Waterman local alignment


def Align_two_seqs(readpart, genepart,W):

    readpart = readpart.strip()
    genepart = genepart.strip()
    # Generating the empty matrices for storing scores and tracing
    row = len(readpart) + 1
    col = len(genepart) + 1
    matrix = np.zeros(shape=(row, col), dtype=int)  
    tracing_matrix = np.zeros(shape=(row, col), dtype=int)  
    
    # Initialising the variables to find the highest scoring cell
    max_score = -1
    max_index = (-1, -1)
    # Calculating the scores for all cells in the matrix
    for i in range(1, row):
        for j in range(1, col):
            # Calculating the diagonal score (match score)
            match_value = 1 if readpart[i - 1] == genepart[j - 1] else -1
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
    aligned_readpart = ""
    aligned_genepart = ""   
    current_aligned_readpart = ""   
    current_aligned_genepart = ""  
    (max_i, max_j) = max_index
    score = matrix[max_i][max_j]
    
    # Tracing and computing the pathway with the local alignment
    while tracing_matrix[max_i, max_j] != 0:
        if tracing_matrix[max_i, max_j] == 3:
            current_aligned_readpart = readpart[max_i - 1]
            current_aligned_genepart = genepart[max_j - 1]
            max_i = max_i - 1
            max_j = max_j - 1
            
        elif tracing_matrix[max_i, max_j] == 2:
            current_aligned_readpart = readpart[max_i - 1]
            current_aligned_genepart = '-'
            max_i = max_i - 1    
            
        elif tracing_matrix[max_i, max_j] == 1:
            current_aligned_readpart = '-'
            current_aligned_genepart = genepart[max_j - 1]
            max_j = max_j - 1
            
        aligned_readpart = aligned_readpart + current_aligned_readpart
        aligned_genepart = aligned_genepart + current_aligned_genepart
    
    # Reversing the order of the sequences
    aligned_readpart = aligned_readpart[::-1]
    aligned_genepart = aligned_genepart[::-1]

    score,anchor =  sequence_similarity(aligned_genepart.replace("*",""),aligned_readpart.replace("*",""), len(readpart.replace("*","")))
    return aligned_readpart.replace("*",""), aligned_genepart.replace("*",""), [score,W,anchor]

#SORTING AND INDEXING BAM IF NOT DONE BEFORE
def sort_indexBamIfNecessary(bam):
    sample = "BT_474"
    dataset = "edgreen"
    if bam == 0:
        #ARRIBA
        if not os.path.isfile("AlignedS0.out.bam"):
            if not os.path.isfile("../"+dataset+"/Arriba/"+sample+"/Aligned.out.bam"):
                sub.run("cp ../"+dataset+"/Arriba/"+sample+"/Aligned.sortedByCoord.out.bam ./AlignedS0.out.bam", shell = True)
            else:
                sub.run("samtools sort ../"+dataset+"/Arriba/"+sample+"/Aligned.out.bam -o AlignedS0.out.bam" , shell = True)
        if not os.path.isfile("AlignedSC0.out.bam"):
            sub.run("samtools rmdup -S AlignedS0.out.bam AlignedSC0.out.bam", shell=True)
        if not os.path.isfile("AlignedSC0.out.bam.bai"):
            sub.run("samtools index AlignedSC0.out.bam" , shell = True)
    elif bam == 1:
        #STARFUSION:
        if not os.path.isfile("AlignedS1.out.bam"):
            sub.run("samtools sort ../"+dataset+"/StarFusion/"+sample+"/STAR-Fusion_outdir/Aligned.out.bam -o AlignedS1.out.bam" , shell = True)
        if not os.path.isfile("AlignedSC1.out.bam"):
            sub.run("samtools rmdup -S AlignedS1.out.bam AlignedSC1.out.bam", shell=True)
        if not os.path.isfile("AlignedSC1.out.bam.bai"):
            sub.run("samtools index AlignedSC1.out.bam" , shell = True)

    #!May add JAFFA in the future
    # else:
    #     #JAFFA
    #     if not os.path.isfile("AlignedS1.out.bam"):
    #         sub.run("samtools sort ../"+dataset+"/JAFFA/"+sample+"/SRR064438_discordant_pairs.bam -o AlignedS1.out.bam" , shell = True)
    #     if not os.path.isfile("AlignedSC1.out.bam"):
    #         sub.run("samtools rmdup -S AlignedS1.out.bam AlignedSC1.out.bam", shell=True)
    #     if not os.path.isfile("AlignedSC1.out.bam.bai"):
    #         sub.run("samtools index AlignedSC1.out.bam" , shell = True)
    #     sub.run("samtools view AlignedS1.out.bam | awk -F'\t' '{split($3, a, /[=:]/); $3 = a[2]; print}' > AlignedSC1.out", shell = True)

 
#USING SAMTOOLS FOR THE SELECTION OF SPANNING FRAGMENTS AT A GIVEN POSITION
def make_breakpointSpanningFragFile(df, i, bam):
    target = df["Fusion"][i]

    break1 = df["base1"][i]
    break2 = df["base2"][i]
    ######################
    print("Analising fusion: ",target)
    print("[*] Computing samtools")
    ######################
    pos1 = str(df["chrom1"][i]) +":"+str(break1-150)+"-"+str(break1+150)
    pos2 = str(df["chrom2"][i]) +":"+str(break2-150)+"-"+str(break2+150)
    if bam == 1:
        pos1 = "chr"+ pos1
        pos2 = "chr"+ pos2
    sub.run("samtools view AlignedSC"+str(bam)+".out.bam "+ pos1 +" -o outpos"+str(bam)+"/outpos1_"+str(i)+"_"+str(break1), shell = True)
    sub.run("samtools view AlignedSC"+str(bam)+".out.bam "+ pos2 +" -o outpos"+str(bam)+"/outpos2_"+str(i)+"_"+str(break2) , shell = True)
    # else:
    #     print("doin")
    #     print("cat AlignedSC1.out | grep chr"+str(df["chrom1"][i])+" > outpos"+str(bam)+"/outpos1_"+str(i)+"_"+str(break1))
    #     sub.run("cat AlignedSC1.out | grep chr"+str(df["chrom1"][i])+" > outpos"+str(bam)+"/outpos1_"+str(i)+"_"+str(break1), shell = True)
    #     sub.run("cat AlignedSC1.out | grep chr"+str(df["chrom2"][i])+" > outpos"+str(bam)+"/outpos2_"+str(i)+"_"+str(break1), shell = True)
    #     print("done")

#!Deprecated BY THE MOMENT
def compute_reverse(read):
    nucs = {"A": "T", "G":"C", "C":"G", "T":"A", "N":"N"}
    new = ""
    for i in read:
        new += i
    return new
#!------------------------>
def apply_cigar_insertion_Deletion(read, cig):
    new = ""
    count = 0
    num = ""
    warning = False
    for i in cig:
        if i.isdigit():
            num+=i
        else:
            if i =="M":
                new += read[count:count+int(num)]
                count+= int(num)
            elif i =="D" and int(num)<50:
                new+= int(num)*"D"
            elif i =="I":
                count+=int(num)
            elif i =="N":
                if int(num) > 50:
                    warning = True
            #         new+= int(num)*"D"
            #     else:
            #         return(new)

            num = ""
    if warning:
        new+= "W"
    return new

#Returns an int that alters the starting position of our read
def left_shift_Softclip(cig):
    num = ""
    First = True
    for i in cig:
        if i.isdigit():
            num+=i
        else:
            if i == "S":
                if not First:
                    return 0
                else:
                    return int(num)
    return 0


def correct_necessaryReads(df, gene = 2):
    reversed_df = []
    for row in df.iterrows():
        flag = format(row[1][1],'b')
        read = list(row[1])
        #print(read)
        if len(str(flag)) > 2:
            # if str(flag)[3] == "0":
            #     read[9] = compute_reverse(read[9])
            shift = left_shift_Softclip(read[5])
            if shift <= 50 and gene == 2:
                read[3] = read[3] - shift
            if "I" in read[5] or "D" in read[5] or "N" in read[5]:
                # or "N" in read[5]:
                read[9]=apply_cigar_insertion_Deletion(read[9],read[5])
            
        reversed_df.append(read)
    return pd.DataFrame(reversed_df)

#SUBSELECTION OF ONLY THE SPANNING READS THAT FALL IN THAT POSITION
def subselect_SpanningReads(df,i,bam):
    break1 = df["base1"][i]
    break2 = df["base2"][i]
    target = df["Fusion"][i]

    outpos1 = pd.read_csv("outpos"+str(bam)+"/outpos1_"+str(i)+"_"+str(df["base1"][i]), sep="\t", header=None,on_bad_lines='skip')
    outpos2 = pd.read_csv("outpos"+str(bam)+"/outpos2_"+str(i)+"_"+str(df["base2"][i]),sep="\t", header=None,on_bad_lines='skip')
    
    kept1 = []
    for row in outpos1.iterrows():
        startposi = row[1][3]-left_shift_Softclip(row[1][5])
        if startposi <= break1 and startposi + len(row[1][9]) >= break1:
            kept1.append(list(row[1]))
    kept1 = pd.DataFrame(kept1)
    #kept1 = RemoveDuplicateBAMReads(kept1)
    kept1 = correct_necessaryReads(kept1)
    kept1.to_csv("outpos"+str(bam)+"/outposK1_"+str(i)+"_"+str(df["base1"][i]) ,sep="\t")
    kept2 = []

    for row in outpos2.iterrows():
        startposi = row[1][3]-left_shift_Softclip(row[1][5]) 
        if startposi <= break2 and startposi + len(row[1][9]) >= break2:
            kept2.append(list(row[1]))
    kept2 = pd.DataFrame(kept2)
    #kept2 = RemoveDuplicateBAMReads(kept2)
    kept2 = correct_necessaryReads(kept2)
    kept2.to_csv("outpos"+str(bam)+"/outposK2_"+str(i)+"_"+str(df["base2"][i]) ,sep="\t")

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
    c = 361
    t = 0
    selected = ""
    flag = False
    for i in lls[1:]:
        if c + len(i) > base and not flag:
            pos = base -c 
            selected += i[pos:]
            flag=True
        elif c + len(i) == base and not flag:
            pos = base -c 
            selected +=i[pos:]
            flag=True 
        elif flag:
            selected+= i
            t+=1
            if t == 13:
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
        if i[0] == str(chromosome):
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
def print_aligned_readsToChr(df,fusion_call,mina,sequenceGene, aligned_reads,n, Printer= True):

    bk1 = df["base"+ n][fusion_call]-mina+360
    fasta = []
    if Printer:
        #with open("Fusions/"+df["Fusion"][fusion_call]+"_"+n+".csv","w") as file:
            print("Break from start:", bk1-360)
            #print(sequenceIntr[360:bk1]+ "***" + sequenceIntr[bk1:bk1+50])
            fasta.append(sequenceGene[360:bk1]+ "***" + sequenceGene[bk1:bk1+50]+"\n")
            print(sequenceGene[360:bk1]+ "***" + sequenceGene[bk1:bk1+50])
            for i in aligned_reads:
                fasta.append(i.strip("W"))
                print(i.strip("W"))
    return(sequenceGene[1:bk1]+ "***" + sequenceGene[bk1:], fasta, bk1-360)

#FUNCITON FOR DUPLICATE REMOVAL

def RemoveDuplicateBAMReads(df, light= False):
    if not light:
        df = df.drop_duplicates(subset=[0], keep=False)
        df = df.reset_index(drop=True)
    else:
        df = df.drop_duplicates(subset=[0, 9], keep = False)
        df = df.reset_index(drop=True)
    return df


def get_introns_exons(minum, bkpoint, chrom):
    
    gtf_df = pd.read_csv("mart_exons.txt", sep='\t', skiprows=1,header=None, low_memory=False)
    # filter to only include features for the chromosome of interest
    chrom_df = gtf_df[gtf_df[2] == chrom]
    chrom_df.sort_values(by=0, inplace=True)

    start_pos =  bkpoint-minum
    # iterate over each 300-base pair interval
    end_pos = bkpoint +50
    output_str = ""
    
    for pos in range(start_pos, end_pos):
        overlapping_features = chrom_df[(chrom_df[0] <= pos ) & (chrom_df[1] >= pos)]
        if len(overlapping_features) == 0:
            # no features overlap with this interval, so it's an intron
            output_str += "I"
        else:
            # at least one feature overlaps with this interval, so it's an exon
            output_str += "E"
        if minum == 1:
            output_str += "***"
        minum-=1
    # print the final output string
    return output_str

#FUNCTION THAT DOES THE MAIN RETRIEVAL OF GENES AND THEIR READS (IN AN ALIGNED WAY)
def RetrieveSpanningReadsAndAssemblyAligned(df,fusion_call,bam):
    mainlist = []

    make_breakpointSpanningFragFile(df,fusion_call,bam)
    isnotMithocondrial = True

    try:
        keptreads1, keptreads2 = subselect_SpanningReads(df,fusion_call,bam)
        # # print(keptreads1.list())
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
        sequenceGene1,fasta1,i1 = print_aligned_readsToChr(df,fusion_call,min1,sequenceGene1, aligned_reads1, "1",)
        introns_exons1 = get_introns_exons(i1,df["base1"][fusion_call],df["chrom1"][fusion_call])
        print(introns_exons1)
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
        sequenceGene2,fasta2,i2 = print_aligned_readsToChr(df,fusion_call,min2,sequenceGene2, aligned_reads2, "2")
        introns_exons2 = get_introns_exons(i2, df["base2"][fusion_call],df["chrom2"][fusion_call])
        print(introns_exons2)
        #########################################################################3

        if "*" not in introns_exons1:
            introns_exons1 = "***"+introns_exons1
        if "*" not in introns_exons2:
            introns_exons2 = "***"+introns_exons2
    except:
        print("Mithocondrial Chromosome hits can't be analised")
        isnotMithocondrial = False

    if isnotMithocondrial:
        return aligned_reads1, sequenceGene1, aligned_reads2, sequenceGene2, fasta1,fasta2,introns_exons1,introns_exons2
    else:
        return 0,'AAA',0,'AAA',[],[],"",""

#FUNCTION FOR GENERATING THE ALIGNMENT OF EACH READ TO SEQUENCE
def AlignReads_ToChromosome(Reads, SequenceGene1, SequenceGene2, Orientation):
    Alignment = []
    Norma = []
    Fusa = []
    if Orientation == "R":
        for read in Reads:
            part_to_align = 50*"*"+read.split("***")[1]
            part_to_align = part_to_align.strip("W")
            W = 1 if "W" in read else 0
            Gene1postBreakpoint = 50*"*"+SequenceGene1.split("***")[1][:len(part_to_align)+round(len(part_to_align)*0.25)]
            Gene2postBreakpoint = 50*"*"+SequenceGene2.split("***")[1][:len(part_to_align)+round(len(part_to_align)*0.25)]
            Al1_1,Al1_2, Norm = Align_two_seqs(part_to_align,Gene1postBreakpoint,W)
            Al2_1,Al2_2, Fus = Align_two_seqs(part_to_align,Gene2postBreakpoint,W)
            Norma.append(Norm)
            Fusa.append(Fus)
            Alignment.append([[Al1_1,Al1_2],[Al2_1,Al2_2]])
    if Orientation == "L":
        #IT IS REVERSED TO ACCOUNT FOR THE ALGORITHM ALIGNING FROM THE LEFT
        for read in Reads:
            W = 1 if "W" in read else 0
            part_to_align = 50*"*"+read.split("***")[0][::-1].strip()
            part_to_align = part_to_align.strip("W")
            Gene2preBreakpoint = 50*"*"+SequenceGene2.split("***")[0][::-1][:len(part_to_align)+round(len(part_to_align)*0.25)]
            Gene1preBreakpoint = 50*"*"+SequenceGene1.split("***")[0][::-1][:len(part_to_align)+round(len(part_to_align)*0.25)]
            Al1_1,Al1_2, Fus = Align_two_seqs(part_to_align,Gene1preBreakpoint,W)
            Al2_1,Al2_2, Norm = Align_two_seqs(part_to_align,Gene2preBreakpoint,W)
            Norma.append(Norm)
            Fusa.append(Fus)
            Alignment.append([[Al1_1,Al1_2],[Al2_1,Al2_2]])

    return Alignment, Norma,Fusa

def Print_SW_Alignments(Al, Norm,Fus):
    c = 1
    print("Alignments to gene")
    for i in Al:
        print(f"Read Number {c}")
        print("G: ",i[0][1])
        print("R: ",i[0][0]) 
        print("Score", Norm[c-1])
        c+=1
    c = 1
    print("Aligments to Fusion")
    for i in Al:
        print(f"Read Number {c}")
        print("G: ",i[1][1])
        print("R: ",i[1][0]) 
        print("Score", Fus[c-1])
        c+=1

def basic_score_check(Norm,Fus):
    c = 0
    anchor = 0
    for i in range(len(Norm)):
        if Norm[i][0] < Fus[i][0] and Fus[i][0] > 0.70 and Fus[i][1] == 0:
            c+=1
            if Fus[i][2]> anchor:
                anchor = Fus[i][2]
        elif Norm[i][0] < Fus[i][0] and Fus[i][0] > 0.80 and Fus[i][1] == 1:
            c+=1
            if Fus[i][2]> anchor:
                anchor = Fus[i][2]
        else:
            pass

    return c, anchor

def main():

    df = pd.read_csv("../Standarized_calls.csv")

    Fusions = []
    Norm1l = []
    Fus1l = []
    Norm2l = []
    Fus2l = []
    count1 = []
    count2 = []
    anchor1 = []
    anchor2 = []
    ReadsG1 = []
    ReadsG2 = []
    SeqGene1 = []
    SeqGene2 = []

    for bams in range(1):
        print("COOMPUTING BAM NUMBER: "+str(bams))
        # * CREATES AlignedSC.bam
        sort_indexBamIfNecessary(bams)
        for fusion_call in range(len(df.index)):
            r1,s1,r2,s2,f1,f2,i1,i2 = RetrieveSpanningReadsAndAssemblyAligned(df, fusion_call,bams)
            c1 = 0
            c2 = 0
            if r1 != 0:
                ReadsG1.append(len(r1)) 
            else:
                ReadsG1.append(0)
            #replace('*', '')
            SeqGene1.append(s1)
            SeqGene2.append(s2)

            if r2 != 0: 
                ReadsG2.append(len(r2))
            else:

                ReadsG2.append(0)
            
            Fusions.append(df["Fusion"][fusion_call])

            if r1 != 0:

                Al1, Norm1, Fus1= AlignReads_ToChromosome(r1,s1,s2,"R")
                consensus = reconstruct_strings(r1,1)
                with open("consensus"+str(fusion_call)+"_1.txt", "a") as file:
                    file.write(df["Fusion" ][fusion_call]+"$")
                    stra = ""
                    for cons in consensus:
                        stra+= cons +","
                    file.write(stra.strip(",")+ ";")  
                for i in range(len(f1[1:])):
                    if Norm1[i][0] < Fus1[i][0] and Fus1[i][0] > 0.70 and Fus1[i][1] == 0:
                        f1[i+1]+=" Fusion"+" scoreF: "+ str(round(Fus1[i][0],2))+" scoreN: "+ str(round(Norm1[i][0],2))+"\n"
                    elif Norm1[i][0] < Fus1[i][0] and Fus1[i][0] > 0.80 and Fus1[i][1] == 1:
                        f1[i+1]+=" Fusion"+" scoreF: "+ str(round(Fus1[i][0],2))+" scoreN: "+ str(round(Norm1[i][0],2))+"\n"
                    else:
                        f1[i+1]+=" Normal"+" scoreF: "+ str(round(Fus1[i][0],2))+" scoreN: "+ str(round(Norm1[i][0],2))+"\n"


                #Print_SW_Alignments(Al1,Norm1,Fus1)
                c1, a1 = basic_score_check(Norm1,Fus1)
                count1.append(c1)
                anchor1.append(a1)
                Norm1l.append(sum(Num[0] for Num in Norm1))
                Fus1l.append(sum(Fu[0] for Fu in Fus1))

            else:
                Norm1l.append("noR")
                Fus1l.append("noR")   
                count1.append("noR")
                anchor1.append("noR")

            with open("Fusions"+str(bams)+"/"+str(fusion_call)+"_"+"1"+".csv","w") as file:
                file.write(i1+ "\n")
                for line in f1:
                    file.write(line)
            if r2 != 0:
                Al2, Norm2, Fus2 = AlignReads_ToChromosome(r2,s1,s2,"L")
                consensus = reconstruct_strings(r2,0)
                with open("consensus"+str(fusion_call)+"_2"+".txt", "a") as file:
                    for cons in consensus:
                        file.write(cons+ ",")
                    file.write("\n")
                for i in range(len(f2[1:])):
                    if Norm2[i][0] < Fus2[i][0] and Fus2[i][0]> 0.70 and Fus2[i][1] == 0:
                        f2[i+1]+=" Fusion"+" scoreF: "+ str(round(Fus2[i][0],2))+" scoreN: "+ str(round(Norm2[i][0],2))+"\n"
                    elif Norm2[i][0] < Fus2[i][0] and Fus2[i][0]> 0.85 and Fus2[i][1] == 1:
                        f2[i+1]+=" Fusion"+" scoreF: "+ str(round(Fus2[i][0]))+" scoreN: "+ str(round(Norm2[i][0]))+"\n"
                    else:
                        f2[i+1]+=" Normal"+" scoreF: "+ str(round(Fus2[i][0]))+" scoreN: "+ str(round(Norm2[i][0]))+"\n"



                #Print_SW_Alignments(Al2,Norm2,Fus2)

                c2 , a2 = basic_score_check(Norm2,Fus2)
                count2.append(c2)
                anchor2.append(a2)
                Norm2l.append(sum(Num[0]for Num in Norm2))
                Fus2l.append(sum(Fu[0] for Fu in Fus2))
            else:
                Norm2l.append("noR")
                Fus2l.append("noR")
                count2.append("noR")
                anchor2.append("noR")
            with open("Fusions"+str(bams)+"/"+str(fusion_call)+"_"+"2"+".csv","w") as file:
                file.write(i2+ "\n")
                for line in f2:
                    file.write(line)
                    
            print("COUNT1: ", c1, "COUNT2: ", c2)
        #print(len(count1), len(count2), len(Fusions))
        pd.DataFrame({'Fusion': Fusions,'Gene1total':ReadsG1, 'Gene2total':ReadsG2 ,'Gene1Normal': Norm1l, 'Gene1Fused':Fus1l, 'Gene2Normal':Norm2l,'Gene2Fused': Fus2l, "Fused1" : count1, "Fused2": count2, "Anchor1": anchor1, "Anchor2" : anchor2 }).to_csv('Results_bam_'+str(bams)+'.csv', index=False)
    pd.DataFrame({'Fusion':Fusions,'SeqGene1': SeqGene1 , 'SeqGene2': SeqGene2}).to_csv('GeneSequences.csv', index=False)

if __name__ == "__main__":
    main()
