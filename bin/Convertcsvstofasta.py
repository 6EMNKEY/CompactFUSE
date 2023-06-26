import subprocess as sub
import time 
import pandas as pd
import os.path
import numpy as np
from termcolor import colored

def split_in_80(line):
    n = 80
    return [line[i:i+n] for i in range(0, len(line), n)]

def convert_file_to_fasta(file,o):
    lines = []
    with open("Fusions/"+file,'r') as readF:
        first = True
        unprocessed = readF.readlines()
        c = 1
        for i in unprocessed:
            if first:
                lines.append(">REFERENCE SEQUENCE GENE "+ o + "\n")
                first = False
            else:
                lines.append(">read "+ str(c)+"\n")
            c+=1
            for j in split_in_80(i):
                lines.append(j)      
    with open("Fastas/"+file,'w') as writeF:
        for i in lines:
            writeF.write(i)

def read_files():
    output = sub.run("ls Fusions | xargs -I % echo %", shell=True, capture_output=True)
    files = output.stdout.decode().split("\n")
    return files

def main():
    files = read_files()
    for i in files[:len(files)-1]:
        print(i)
        if "_1.csv" in i:
            convert_file_to_fasta(i, "1")
        else:
            convert_file_to_fasta(i, "2")


if __name__ == "__main__":
    main()