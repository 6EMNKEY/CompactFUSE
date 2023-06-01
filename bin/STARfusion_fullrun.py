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
	if ">" in i.strip():
		path = i.strip(">")
		print("RUNING STAR-Fusion FOR "+ sys.argv[1] +path)
	else:
		print(path+"."+str(i))
		read1 = sys.argv[1] + "/" + path + "/" + str(i.strip()) + "_1.fastq.gz"
		read2 = sys.argv[1] + "/" + path + "/" + str(i.strip()) + "_2.fastq.gz"
		print("perl STAR-Fusion --left_fq " +read1+" --right_fq " + read2 + " --genome_lib_dir GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ --FusionInspector validate")
		subprocess.run("perl STAR-Fusion --left_fq " +read1+" --right_fq " + read2 + " --genome_lib_dir GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ --FusionInspector validate", shell= True)
		time.sleep(2)
		subprocess.run("mkdir "+ path+"."+str(i), shell=True)
		time.sleep(120)
		subprocess.run("mv STAR-Fusion_outdir "+ path+"."+str(i.strip("\n")), shell=True)
		time.sleep(10)
		subprocess.run("mv " + path+"."+str(i.strip("\n")) +" /media/sequentia/sdb1/visitor3/"+sys.argv[1]+"/StarFusion", shell=True)
		time.sleep(2)
		c+= 1
		print(str(i.strip())+ " WAS COMPUTED " + str(c) +" OUT OF " + str(srrc) +" Have been completed")
