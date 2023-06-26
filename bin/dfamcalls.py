import subprocess as sub
import time 
import pandas as pd
import os.path
import numpy as np
from termcolor import colored
import requests

#!REMEMBER TO REMOVE SAMPNAME DEFAULT
def read_data(sampname = "BT_474"):
    df = pd.read_csv("BREAK/GeneSequences.csv")
    return(df)

def make_call(sequence):
    url = "https://dfam.org/api/searches"
    json_send = {  
        "sequence": sequence,
        "organism": "Homo sapiens",
        "cutoff": "curated"}
    
    response = requests.post(url, json=json_send)
    id_search = response.json()["id"]
    print("Requested search with id:", id_search)
    get_id = {"id": id_search}
    response_real = requests.get(url+"/"+id_search,params=get_id)
    finished = False
    while not finished:
        if "Not finished" in response_real.text:
                time.sleep(2)
                response_real = requests.get(url+"/"+id_search,params=get_id)
        else:
            #print(response_real.text)
            finished = True
    return response_real

def retrieve_important_info(response):
    hits = response.json()["results"][0]["hits"]
    tandem_reps = response.json()["results"][0]["tandem_repeats"]
    results = []
    tandems = []
    c = 0
    if hits != []:
        for hit in hits:
            results.append([])
            results[c].append(hit["type"])
            results[c].append(hit["ali_start"])
            results[c].append(hit["ali_end"])
            c+=1
    c = 0
    if tandem_reps != []:
        for tandem in tandem_reps:
            tandems.append([])
            tandems[c].append(tandem["type"])
            tandems[c].append(tandem["start"])
            tandems[c].append(tandem["end"])
            c+=1
    # print("RESULTS:", results)
    # print("TANDEMS", tandems)
    return results, tandems

def reorder_vals(a,b):
    if int(a) > int(b):
        return b,a
    else:
        return a,b
    
#? I SHOULD IMPLEMENT A WAY TO SEE WHICH KIND OF TANDEM IS DETECTED AND THE VIABILITY ETCETCETC....
def main():
    df = read_data()
    t = 1
    rowlist = []
    for i in range(df.shape[0]):
        print("FUSION NUMBER "+ str(t))
        print("Gene1:",)
        response = make_call(str(df["SeqGene1"][i]).replace("*", ""))
        lengther = len(str(df["SeqGene1"][i]).replace("*", ""))
        res,tan = retrieve_important_info(response)
        try:
            bk = str(df["SeqGene1"][i]).index("***")
            if res != []:
                c = 1
                for j in res:
                    a , b = reorder_vals(j[1],j[2])
                    print(j[0],a,b)
                    rowlist.append([df["Fusion"][i],bk,"Transposon","G1",lengther,j[0], a,b,c])
                    c+=1

            if tan != []:
                for j in tan:
                    a , b = reorder_vals(j[1],j[2])
                    print(j[0],a,b)
                    rowlist.append([df["Fusion"][i],bk,"Tandem","G1",lengther,j[0], a,b,c])
                    c+=1
        except:
            print("Mithocondrial Chromosomes not analised.. SKIPPING")
        print("Gene2:")
        response = make_call(str(df["SeqGene2"][i]).replace("*", ""))
        lengther = len(str(df["SeqGene2"][i]).replace("*", ""))
        res,tan = retrieve_important_info(response)
        try:
            bk = str(df["SeqGene2"][i]).index("***")
            if res != []:
                c = 1
                for j in res:
                    a , b = reorder_vals(j[1],j[2])
                    print(j[0],a,b)
                    rowlist.append([df["Fusion"][i],bk,"Transposon","G2",lengther,j[0], a,b,c])
                    c+=1
            if tan != []:
                c = 1
                for j in tan:
                    a , b = reorder_vals(j[1],j[2])
                    print(j[0],a,b)
                    rowlist.append([df["Fusion"][i],bk,"Tandem","G2",lengther,j[0], a,b,c])
                    c+=1
        except:
            print("Mithocondrial Chromosomes not analised.. SKIPPING")
        t += 1

    towrite = pd.DataFrame(columns=["Fusion", "Breakpoint","Gene","Class","Length","Type", "Start", "End", "height"], data=rowlist)
    towrite.to_csv("DFAMcalls.csv", index=False)

if __name__ == "__main__":
    main()