import subprocess
import time
import sys

with open(sys.argv[1]+ "/SRA_list", "r") as file:
	lista = file.readlines()
	srrc = 0
	for i in lista:
		if "SRR" in i:
			srrc+= 1


pwd = subprocess.run("pwd", shell= True, capture_output= True)
wk = pwd.stdout.decode().strip("\n")

print("Script run in "+ wk)
path = ""
c = 0
for i in lista:
	if "SRR" not in i.strip():
		path = i.strip()
		print("RUNING JAFFA FOR "+wk+sys.argv[0]+path)
	else:
                print(path+"."+str(i))
                read1 = sys.argv[1]+ "/" + path + "/" + str(i.strip()) + "_1.fastq.gz"
                read2 = sys.argv[1]+ "/" + path + "/" + str(i.strip()) + "_2.fastq.gz"
                print("tools/bin/bpipe run JAFFA_hybrid.groovy "+read1+ " "+ read2)
                subprocess.run("tools/bin/bpipe run JAFFA_hybrid.groovy "+read1+ " "+ read2, shell= True)
                subprocess.run("tools/bin/bpipe stop", shell=True)
                time.sleep(2)
                subprocess.run("mkdir "+ path+"."+str(i), shell=True)
                time.sleep(2)
                subprocess.run("mv jaffa_results.csv jaffa_results.fasta "+str(i.strip())+"/* "+ path+"."+str(i), shell=True)
                time.sleep(2)
                subprocess.run("mv " + path+"."+str(i.strip("\n")) +" /media/sequentia/sdb1/visitor3/"+sys.argv[1]+"/JAFFA", shell=True)
                c+= 1
                print(str(i.strip())+ " WAS COMPUTED " + str(c) +" OUT OF " + str(srrc) +" Have been completed")
