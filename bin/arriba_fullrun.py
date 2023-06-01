import subprocess
import time
import sys

with open(sys.argv[1]+"/SRA_list.txt", "r") as file:
	lista = file.readlines()
	srrc = 0
	for i in lista:
		if ">" not in i:
			srrc+= 1


pwd = subprocess.run("pwd", shell= True, capture_output= True)
wk = pwd.stdout.decode().strip("\n")

print("Script run in "+ wk)
path = ""
c = 0
for i in lista:
	if "SRR" not in i.strip():
		path = i.strip()
		print("RUNING ARRIBA FOR "+wk+sys.argv[0]+path)
	else:
                print(path+"."+str(i))
                read1 = sys.argv[1]+ "/" + path + "/" + str(i.strip()) + "_1.fastq.gz"
                read2 = sys.argv[1]+ "/"+ path + "/" + str(i.strip()) + "_2.fastq.gz"
                print("./run_arriba.sh STAR_index_GRCh38_ENSEMBL104 Homo_sapiens.GRCh38.108.chr.gtf Homo_sapiens.GRCh38.dna.primary_assembly.fa database/blacklist_hg38_GRCh38_v2.3.0.tsv database/known_fusions_hg38_GRCh38_v2.3.0.tsv database/protein_domains_hg38_GRCh38_v2.3.0.gff3 10 "+read1+ " "+ read2)
                subprocess.run("./run_arriba.sh STAR_index_GRCh38_ENSEMBL104 Homo_sapiens.GRCh38.108.chr.gtf Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa database/blacklist_hg38_GRCh38_v2.3.0.tsv database/known_fusions_hg38_GRCh38_v2.3.0.tsv database/protein_domains_hg38_GRCh38_v2.3.0.gff3 10 "+read1+ " "+ read2, shell= True)
                time.sleep(2)
                subprocess.run("mkdir "+ path+"."+str(i), shell=True)
                time.sleep(2)
                subprocess.run("mv Aligned.sortedByCoord.out.bam Aligned.sortedByCoord.out.bam.bai fusions.tsv Log.out Log.std.out fusions.discarded.tsv Log.final.out Log.progress.out SJ.out.tab "+ path+"."+str(i), shell=True)
                print("mv " + path+"."+str(i) +"/media/sequentia/sdb1/visitor3/"+sys.argv[1]+ "/Arriba/")
                subprocess.run("mv " + path+"."+str(i.strip("\n")) +" /media/sequentia/sdb1/visitor3/"+sys.argv[1]+ "/Arriba/", shell=True)
                c+= 1
                print(str(i.strip())+ " WAS COMPUTED " + str(c) +" OUT OF " + str(srrc) +" Have been completed")
