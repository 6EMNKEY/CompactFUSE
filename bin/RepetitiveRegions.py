import subprocess as sub
import time 
import pandas as pd
import os.path
import numpy as np
from termcolor import colored

def read_data():
    calls = pd.read_csv("Standarized_calls.csv")
    return calls

def subselect_nec(chromosome,bk):
    out = sub.run("grep chr" + str(chromosome) + " simpleRepeats.bed", capture_output=True, shell= True)
    result = 0
    for row in out.stdout.decode().split("\n"):
        order = row.split("\t")
        if order[0] != "":
            if bk > int(order[1]) and bk < int(order[2]):
                result = 1
                break
    return result

def main():
    calls = read_data()
    for i in range(calls.shape[0]):
        r1 = subselect_nec(calls["chrom1"][i],calls["base1"][i])
        if r1 == 1:
            print(calls["Gene1"][i]+"\t"+"Repeated region",end="\t")
        else:
            print(calls["Gene1"][i]+"\t"+"Normal region", end="\t")

        r2 = subselect_nec(calls["chrom2"][i],calls["base2"][i])
        if r2 == 1:
            print(calls["Gene2"][i]+"\t"+"Repeated region")
        else:
            print(calls["Gene2"][i]+"\t"+"Normal region")
            
if __name__ == "__main__":
    main()