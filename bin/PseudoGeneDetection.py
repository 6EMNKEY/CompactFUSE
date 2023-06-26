import subprocess as sub
import time 
import pandas as pd
import os.path
import numpy as np
from termcolor import colored


def read_data():
    calls = pd.read_csv("Standarized_calls.csv")
    return calls

def subselect_pseudo(gene):
    out = sub.run("grep " + gene + " Homo_sapiens.GRCh38.109.chr.gtf", capture_output=True, shell= True)
    result = 0
    if "pseudogene" in out.stdout.decode():
        result = 1
    return result

def main():
    calls = read_data()
    for i in range(calls.shape[0]):
        r1 = subselect_pseudo(calls["Gene1"][i])
        if r1 == 1:
            print(calls["Gene1"][i]+"\t"+"Pseudogene",end="\t")
        else:
            print(calls["Gene1"][i]+"\t"+"Gene", end="\t")

        r2 = subselect_pseudo(calls["Gene2"][i])
        if r2 == 1:
            print(calls["Gene2"][i]+"\t"+"Pseudogene")
        else:
            print(calls["Gene2"][i]+"\t"+"Gene")
            
if __name__ == "__main__":
    main()

