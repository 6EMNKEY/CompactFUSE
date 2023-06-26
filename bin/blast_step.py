import subprocess as sub
import time 
import pandas as pd
import os.path
import numpy as np
from termcolor import colored


def do_blast(lines,c,n):
    with open("fasta"+str(c)+"_"+str(n), "w") as file:
        t = 0
        for i in lines.split(","):
            file.write(">"+ str(t)+ "\n")
            t+=1
            file.write(i+"\n")

    for i in lines.split(","):
        A = sub.run("blastn -db nt -query "+ " fasta"+str(c)+"_"+str(n) +" -subject_besthit -remote", shell= True, capture_output=True)# -out "+ "blasts" + str(c)+"_"+str(n))
        print(A.stdout.decode())
        print(A.stderr.decode())
def iterate_lines(lines):
    fulltext = ""
    c = 0
    for line in lines:
        working = line.split("$")
        fulltext += "Sequences found in: "+ working[0]+"\n"
        if working[1:].count(",;") != 1 and working[1].count(",") != 1:
            wk1 = working[1].split(";")[0]
            do_blast(wk1,c, 1)
            wk2 = working[1].split(";")[1]
            do_blast(wk2,c, 2)
            break

        
def main():
    with open("consensus.txt", "r") as file:
        lines = file.readlines()
    iterate_lines(lines)



if __name__ == "__main__":
    main()

