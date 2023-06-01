import subprocess
import time
import sys

with open(sys.argv[1]+"/SRA_list", "r") as file:
	lista = file.readlines()
	srrc = 0
	for i in lista:
		if ">" not in i:
			srrc+= 1


pwd = subprocess.check_output("pwd", shell= True)
wk = pwd.decode().strip("\n")

print("Script run in "+ wk)
path = ""
c = 0
for i in lista:
	if "SRR" not in i.strip():
		path = i.strip()
		print("RUNING FusionCatcher FOR "+wk+ sys.argv[0] +path)
	else:
                print(path+"."+str(i))
                read1 = sys.argv[1] + "/" + path + "/" + str(i.strip()) + "_1.fastq.gz"
                read2 = sys.argv[1] + "/" + path + "/" + str(i.strip()) + "_2.fastq.gz"
                print("python2 bin/fusioncatcher.py \ -d data/current/ -i "+ sys.argv[0] + "/" +path + "-o " + path+"."+str(i))
                subprocess.call(["python2 bin/fusioncatcher.py", "-d", "data/current/", "-i",  sys.argv[0] + "/" +path , "-o" , path+"."+str(i)], shell= True)
                time.sleep(2)
                subprocess.call(["mv" , path+"."+str(i.strip("\n")) ," /media/sequentia/sdb1/visitor3/"+sys.argv[1]+"/FusionCatcher"], shell=True)
                c+= 1
                print(str(i.strip())+ " WAS COMPUTED " + str(c) +" OUT OF " + str(srrc) +" Have been completed")

